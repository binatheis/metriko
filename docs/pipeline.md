# The pipeline

metriko implements the QGP (quadrilateral grid parametrization) pipeline described by Campen et al., with a QEx
extraction at the end. Knowing the stages helps interpret the timings, the failures and the quality of the
result.

```
input (validate)  ->  decimation  ->  field  ->  parameterization  ->  tmesh
      ->  quantization  ->  collapse  ->  embedding  ->  tutte  ->  slim  ->  qex  ->  output
```

All the stages below run inside `remesh()`. Each is named by a `RemeshStage` / `metriko.Stage`, which is what an
error reports. The "typical cost" ranges are approximations from our reference set (see
[performance.md](performance.md)); run with `--timings` on your own mesh for the real numbers.

## 1. `input` — validation

Checks the contract: closed, manifold, one component, triangles, valid indices. Rejects with
`std::invalid_argument` / `ValueError`, not with a stage error, because it is a caller error rather than a
pipeline failure. Fast (milliseconds).

## 2. `decimation` — reduce the input to a workable size

Runs `pmp::decimate` when the input reaches `faces_threshold` faces (300k) or when forced. Targets
`max_verts` (90k) or `vertices_percentage`. The mesh that is actually remeshed is returned as
`used_pos`/`used_idx`.

Why it exists: the parameterization is a global linear solve on the vertices, and its cost grows superlinearly.
Dense meshes also hit T-mesh degeneracies more often. Decimating first is what makes the 300 s budget achievable
on scan-sized inputs.

Typical cost: 1-10 s for a 250k-face mesh.

## 3. `field` — 4-RoSy cross field

Computes a smooth 4-fold rotationally symmetric direction field on the faces (smoothest, or curvature-aligned).
This field defines the orientation of the quads. A small number of **singularities** (where the field's rotation
index is non-zero) is normal and necessary on a closed surface.

Typical cost: 5-25 s at 90k vertices.

## 4. `parameterization` — seamless uv map

Integrates the field into a seamless parametrization: a least-squares problem whose unknowns are the complex
values of the map plus the per-face rotations. Solved with a Levenberg-Marquardt loop, each iteration assembling a
sparse Jacobian and factorizing the normal equations with CHOLMOD.

**This is often the most expensive stage** (150-400 s at 90k vertices). Two facts matter for tuning:

* the cost is dominated by the numeric factorization of a large sparse system (~65 % of the LM time on our
  measurements); the symbolic analysis runs once;
* the LM's unknown count is ~18 per vertex in the current formulation. Reducing it would require a different
  formulation, which was investigated and is *not* a safe local change (see `HANDOFF.md` §16).

Typical cost: 20 s at 30k vertices, 130-200 s at 60k, 300-400 s at 90k.

## 5. `tmesh` — motorcycle graph and T-mesh

Traces the integer grid lines across the parametrization (the motorcycle graph) and builds the T-mesh: the
arrangement of grid lines where singularities, seams and boundary act as sources. The T-mesh is what the
quantization will make integer.

Typical cost: 1-10 s.

## 6. `quantization` — make the T-mesh integer

Solves for integer coordinates along the T-mesh arcs so that every quad of the arrangement is a unit square in the
integer grid, while keeping the t-edges near their real-valued lengths. It runs in two steps:

1. a rounding pass that drives the fractional parts to 0 (with a wake/sleep heuristic);
2. a *pricing* loop (column generation): while a negative-cost strip loop exists, apply ±1 along it and re-price.
   This is a Bellman-Ford/SPFA `find_negative_cycle` over the strip graph, and it is where the time goes when many
   t-edges are driven to zero (hundreds of them on hard meshes).

Typical cost: 5-40 s on ordinary meshes; 100-460 s on meshes whose quantization collapses many t-edges. This is
the stage that most often exceeds a budget after the parameterization.

## 7. `collapse` — remove zero-arc t-quads

With integer coordinates, many t-quads collapse to zero width. They have to be removed so the T-mesh becomes a
regular quad grid: t-edges are snapped to their neighbours, zero-width thalfs are collapsed, crossing t-edges are
re-routed. This is the most delicate combinatorial stage.

Typical cost: 20-140 s. Failures here are the "degenerate T-mesh" class reported as
`boundary segment ... lies in no face`, `tquad ... is self-intersected`, `left N zero-width tquads`, etc.

## 8. `embedding` — put the T-mesh back on the surface

Maps every T-mesh vertex to the surface through the parametrization, producing the initial quad layout on the
mesh.

Typical cost: ~1-3 s.

## 9. `tutte` — planar embedding of the quad layout

Computes a planar embedding of each patch (the quad layout is a disk per patch, with the singularities and the
seams as constraints) using a Tutte-style solve. This gives the initial uv layout for SLIM.

Typical cost: seconds to tens of seconds; the linear systems grow with the number of T-mesh vertices.

## 10. `slim` — injective relaxation

Runs SLIM (Shtengel et al. / Rabinovich et al.) on the uv layout: it minimizes a symmetric-Dirichlet-like energy
while keeping every triangle positively oriented, so the final parametrization is injective. Two details matter:

* the iteration cap is `slim_max_iter` (30 by default) and landing on the cap is normal;
* the inner linear solve is a large sparse factorization; the assembly is deterministic.

Typical cost: 30-90 s at 90k vertices.

## 11. `qex` — extract the quad mesh

Turns the integer-grid parametrization into quads with a QEx-style algorithm: q-vertices, q-ports, q-edges and
q-faces are traced along the grid, then the quad soup is welded and relaxed on the surface
(`conformalize_on_surface`).

Some rings may not close (folded or collapsed uv triangles break the tracing). Those rings are dropped — the
libQEx behaviour — and counted in `holes`; the loops they leave are closed with fill triangles so the output is
still watertight. If more than 2 % of the rings fail to close, the run fails instead: that means the
parametrization itself is broken rather than a few local cells.

Typical cost: 2-10 s.

## Where the time goes

On our reference set (idle machine, 90k-vertex cap):

| mesh | total | dominant stage |
|---|---|---|
| manotas_30k | 43 s | field + parameterization |
| azure_30k | 19 s | field |
| gat_30k | 91 s | parameterization |
| drag_30k | 188 s | quantization (98 s) |
| GEO-head_stylized | 123 s | parameterization (with one retry) |
| Azure Alien Head Scan (1) | 103 s | parameterization |
| manotas (999k faces) | 132 s | parameterization |
| gatubelo | 261 s | parameterization + collapse (2 retries) |

See [performance.md](performance.md) for how to read the timings and what to tune.
