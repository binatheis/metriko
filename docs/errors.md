# Errors and what to do about them

Every algorithmic failure is a `metriko::RemeshFailure` (C++) / `metriko.RemeshError` (Python, a
`RuntimeError` subclass) with a **stage** and a human-readable detail:

```
parameterization: the initial integration ran past the stage timeout
```

```python
try:
    Q_pos, Q_idx, patch = metriko.remesh(V, F)
except metriko.RemeshError as err:
    print(err.stage)   # metriko.Stage.parameterization
```

In C++, `catch (const metriko::RemeshFailure& e) { e.stage; }`. Invalid *inputs* (open mesh, non-manifold,
multiple components) are not stage failures: they raise `std::invalid_argument` / `ValueError` immediately.

## The stages

| stage | what it does | what a failure means |
|---|---|---|
| `input` | validation | not a stage failure (see above) |
| `decimation` | reduce the input | `pmp::decimate` could not meet the constraints; relax them or disable decimation |
| `field` | 4-RoSy cross field | the eigen/Laplacian solve broke: degenerate or extremely irregular geometry |
| `parameterization` | seamless uv map | the Levenberg-Marquardt solve failed to converge, or ran past the stage timeout |
| `tmesh` | motorcycle graph + T-mesh | the grid lines could not be traced (degenerate parametrization) |
| `quantization` | integer T-mesh | the pricing loop hit the budget, or the rounding produced an inconsistent grid |
| `collapse` | remove zero-arc t-quads | degenerate T-mesh: a snapped chord with no common face, a self-intersected tquad, leftover zero-width tquads |
| `embedding` | map the T-mesh to the surface | a T-mesh edge could not be placed (inconsistent T-mesh) |
| `tutte` | planar embedding | the patch layout is not a disk (topological inconsistency) |
| `slim` | injective uv relaxation | the relaxation could not keep the map injective |
| `qex` | extract quads | more than 2 % of the quad rings failed to close: the uv map is folded |

## Message shapes

* `the time budget was exceeded during stage X` — the total wall-clock budget ran out while `X` was running.
* `X: <detail>` where the detail says what broke (e.g. `t-mesh collapse left 2 zero-width tquads: [113, 260]`).
* `the initial integration ran past the stage timeout` — a single stage exceeded `stage_timeout_seconds`; the
  attempt is abandoned and the retry ladder coarsens the grid.

## Decision table

| situation | do this |
|---|---|
| `parameterization` or `quantization` hits the budget | raise `time_budget`/`stage_timeout`, or lower `decimate.max_verts` (60k-45k), or raise `scale` |
| `collapse` / `embedding` with a "degenerate T-mesh" message | let the retry ladder work (default), or raise `scale` yourself; this class is mesh-specific and resolution-dependent |
| `tutte` / `slim` "not a disk" / "not injective" | raise `scale` (coarser grid), or try `field="smoothest"` |
| `qex` "rings do not close" | the uv map is folded: raise `scale`, raise `slim_max_iter`, or decimate more |
| `decimation` fails | relax the decimation constraints (`hausdorff`, `normal_dev`, `aspect_ratio`) or disable decimation |
| `input` raises `ValueError` | the mesh is not closed/manifold/single-component; repair it outside metriko |

## The retry ladder

By default (`retry=True`, `max_attempts=3`, `coarsen_factor=1.6`), a *stage* failure restarts the pipeline with a
coarser scale (`scale × 1.6` per attempt) and, when decimation targets the automatic cap, a smaller target
(`× 0.6`). A `time_budget`-exceeded failure also triggers the ladder only if the attempt was cheap enough; a run
that already spent the budget raises immediately.

The reason is empirical: nearly all degenerate configurations are tied to the density of the integer grid, and a
coarser grid escapes them. The cost is a coarser mesh than requested — check `used_scale` (or
`stats["used_scale"]`).

If you need predictable density, set `max_attempts=1` (`--no-retry`) and handle the failure upstream.

## What never happens

metriko never returns a mesh that is open or non-manifold to signal a problem: it either returns a watertight
result or raises. The `holes` field is not an error — it reports how many rings had to be dropped and closed with
fill triangles.
