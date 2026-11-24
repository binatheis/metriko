# Parameters

Every option exists in three places with the same name: `metriko::RemeshOptions` (C++), the `metriko-cli` flag,
and a keyword argument of `metriko.remesh()` (Python). Defaults are the same everywhere unless noted.

Quick map:

* **Shape of the result** — `scale`, `field`, `refine`, `slim_max_iter`, `slim_rel_tol`.
* **Input size** — the `decimate*` group.
* **Time and failure behaviour** — `time_budget_seconds`, `stage_timeout_seconds`, `retry_coarser`,
  `max_attempts`, `coarsen_factor`.
* **Reproducibility and parallelism** — `deterministic`, `threads`.
* **Diagnostics** — `report_timings`, `quiet`.

---

## Shape of the result

### `scale` — `--scale` — `scale=`

Target quad size, as a fraction of the input's **bounding-box diagonal**. Default `0.02`.

* 0.02 → roughly 50 quads across the model.
* Halving it roughly quadruples the quad count; doubling it quarters it.
* Cost grows steeply: the t-mesh (and therefore the quantization and the collapse) scales with the number of
  quads, and the parameterization cost grows superlinearly with the vertex count. Halving the scale is much more
  expensive than a factor of two, and it is much more likely to hit a degenerate configuration.
* On a retry, the *effective* scale is multiplied by `coarsen_factor`; the value actually used for the returned
  mesh is reported as `used_scale`.

Practical advice: change it in steps of √2 and re-check `used_scale`.

### `field` — `--field` — `field=`

Which cross field drives the parametrization. Default `smoothest`.

* `smoothest` — the smoothest 4-RoSy field (least singularities). Best general choice, and the fastest.
* `curvature_aligned` — aligned with the principal curvature directions where they are reliable, smoothest
  elsewhere. Produces more natural edge flow on organic shapes (body scans, faces) at a slightly higher cost and
  with more singularities.

### `refine` — `--no-refine` — `refine=`

Run the surface relaxation (shape-up with squares) on the extracted quad mesh. Default on.

Turning it off returns the raw extraction: the quads follow the integer grid exactly but can be visibly irregular.
Keep it on unless you want the raw grid for inspection.

### `slim_max_iter` / `slim_rel_tol` — `--slim-max-iter` / `--slim-rel-tol`

Iteration cap and relative tolerance of the SLIM relaxation that keeps the parametrization injective.
Defaults `30` and `1e-4`.

SLIM hitting the cap (`slim: 30 of 30 iterations`) is normal on hard meshes and is not an error: it means the
energy was still improving but the cost had to be bounded. Raising the cap improves the parametrization marginally
and costs linearly; lowering it saves time and risks a slightly worse map.

---

## Input size: decimation

Decimation is what keeps the pipeline affordable on dense scans. It runs **once**, before the remeshing, and
the mesh that was actually remeshed comes back as `used_pos`/`used_idx` (or the extra arrays with
`return_decimated=True`).

### `decimate.enabled` — `--decimate` / `--no-decimate` — `decimate=`

Master switch. Default: enabled, but only *triggered* automatically when the input reaches
`faces_threshold` faces. `--decimate` forces it on for smaller meshes; `--no-decimate` disables it entirely.

### `decimate.faces_threshold` — `--decimate-faces-threshold`

Input face count at or above which decimation runs automatically. Default `300000`.

### `decimate.vertices_percentage` — `--decimate-verts-pct`

Percentage of the input vertices to keep. `0` (default) means "use the automatic cap" (`max_verts`). Use this when
you want a specific *relative* reduction instead of an absolute target.

### `decimate.max_verts` — `--decimate-max-verts`

Vertex cap used when `vertices_percentage == 0`. Default `90000`.

This is the most consequential knob in the whole API for large inputs. The parameterization and SLIM cost grow
superlinearly with it, and dense decimations hit T-mesh edge cases more often. 90k is a compromise that keeps a
typical 250k-face scan under the budget. Lower it (45k-60k) if you are fighting the time budget or a pathological
mesh; raise it only if your inputs are small and you need the extra fidelity.

### `decimate.normal_deviation` / `aspect_ratio` / `min_edge_length` / `max_valence` / `hausdorff_error`

Passed through to `pmp::decimate`; `0` means unconstrained. Defaults: all `0` except `normal_deviation`, which the
CLI exposes but leaves at 0.

A typical "preserve the shape" combination is `hausdorff_error` at a fraction of the quad size (in mesh units) plus
`normal_deviation` around 30-45°.

### `decimate.seam_threshold` / `seam_deviation`

Texture-seam preservation of `pmp::decimate`. `seam_threshold < 0` (the default) disables it. Note that the OBJ
loader in the Python front-end only reads positions, so seams are not available there; the CLI has the same
limitation.

---

## Time and failure behaviour

The defaults are tuned for interactive/local use: **300 s total, 150 s per stage, 3 attempts**. The C++ API
defaults are `0` (unlimited) because a library should not impose a time policy — set them explicitly.

### `time_budget_seconds` — `--time-budget` — `time_budget=`

Wall-clock budget for the whole call. The budget is checked at stage boundaries and before starting another
attempt; a stage that starts before the deadline is allowed to finish (with a 10 % grace), so the actual wall time
can exceed the budget by up to about one stage. `0` disables it.

Exceeding it raises `RemeshError` with the stage that was running. The message names the stage, e.g.
`quantization: the time budget was exceeded during stage quantization`.

### `stage_timeout_seconds` — `--stage-timeout` — `stage_timeout=`

Maximum duration of a single stage before it is abandoned. A stage timeout is treated as a *failure of that
attempt* and triggers the coarser retry, so a slow but non-hopeless mesh gets a second chance on a coarser grid
instead of consuming the whole budget. `0` disables it.

This is the knob to raise when a mesh legitimately needs a slow parameterization: a 90k-vertex mesh can need
300-400 s there.

### `retry_coarser` — `--no-retry` — `retry=`

When a stage fails, restart with a coarser configuration instead of giving up. Default on.

Each retry:

* multiplies the scale by `coarsen_factor` (default `1.6`);
* multiplies the decimation target by `0.6` (when decimation is enabled and the target is the automatic cap).

The rationale: dense integer-grid maps are the main source of degenerate configurations (folded uv triangles,
zero-width t-mesh quads, self-intersecting tquads), and those get rarer as the grid coarsens. The returned mesh
reports the scale it was actually built with (`used_scale`).

`--no-retry` is the right choice when you prefer a fast, deterministic *failure* to a slower, coarser success —
for example in a batch job that re-queues failures with different settings.

### `max_attempts` — `--max-attempts`

Number of attempts including the first. Default `3`. With the defaults, attempt 3 uses `scale * 1.6² = 2.56 ×` the
requested scale, so the quad count is ~6.5× lower than requested — check `used_scale` before assuming the density.

### `coarsen_factor` — `--coarsen-factor`

Multiplier applied to the scale per retry. Default `1.6`. Larger jumps converge to a usable mesh faster but
deviate more from the requested density.

---

## Reproducibility and parallelism

### `deterministic` — `--fast` (inverse) — (C++ only; the Python front-end does not expose it yet)

Default `true`. With it, the scheduling-dependent floating-point order is removed (single-threaded CHOLMOD
factorization, deterministic SLIM assembly, Eigen's internal threading off), so two runs of the same mesh on the
same machine and thread count give the same quad/vertex/hole counts.

`deterministic = false` (CLI: `--fast`) lets CHOLMOD use all its threads: about 1.3× faster on the initial solve,
but run-to-run variation comes back (different pricing trajectories can change the result).

### `threads` — `--threads`

Number of OpenMP threads. `0` (default) leaves the runtime default (all cores). `--threads 1` is the way to
reproduce a run on a different machine or under a different load.

---

## Diagnostics

### `report_timings` — `--timings` — `timings=`

Print the wall time of every pipeline stage to stderr, plus internal counters:

```
[metriko] stage decimation          6.008 s
[metriko] stage field              19.538 s
[metriko] stage parameterization  412.666 s
[metriko] naive LM: 52 of 100 iterations, 301.324 s
[metriko] cholmod: 1 analyses, copy 7.9 s, analyze 2.7 s, numeric 203.8 s
[metriko] LM system: 1620038 unknowns, 3600082 residuals ...
[metriko] slim: 4 of 30 iterations, 184657 unknowns, energy 2.1e+05
[metriko] stage quantization      231.261 s
```

It costs nothing when off, and it is the first thing to look at when a mesh is slow.

### `quiet` — `--quiet`

Silence the library's informational output (the stage lines above except under `--timings`, the q-port pairing
summary, etc.). Errors are still raised.

---

## `--patches` / `RemeshResult::val`

Writes one integer per quad: the patch label. `-1` means the extractor dropped rings in that run and the labels
could not be transferred reliably (see [output.md](output.md)).

---

## `--decimate-only`

Decimate with the `--decimate-*` options and write the triangle mesh, without remeshing. Useful to prepare test
meshes, to measure how much decimation costs, or to inspect what the pipeline will actually work on. Unlike
`--decimated-out` (which is only written when the remeshing succeeds), `--decimate-only` always writes the output.
