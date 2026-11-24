# Getting started

## 1. Install

### Python (recommended)

```bash
pip install metriko-<version>-cp312-abi3-<platform>.whl
```

The wheel is `abi3` for CPython ≥ 3.12, so one file works for 3.12, 3.13 and later on the same platform. It
bundles every non-system dependency (OpenBLAS, the SuiteSparse libraries), so no compiler or package manager is
needed.

Requirements at run time: a CPU with AVX2/FMA (the bundled OpenBLAS is built for the `HASWELL` target). On an
older CPU the import fails with an illegal-instruction error; build from source with a different OpenBLAS target
in that case (see [building.md](building.md)).

### CLI

Unpack the archive for your platform (`metriko-cli-windows-x86_64.zip`,
`metriko-cli-linux-x86_64.tar.gz`, `metriko-cli-macos-arm64.tar.gz`) and put `metriko-cli` on your `PATH`. It is
a self-contained executable with the same requirements as the wheel.

### From source

See [building.md](building.md).

## 2. Check that your mesh is in scope

metriko only accepts **closed, manifold, single-component triangle meshes**. It validates this on entry and
refuses anything else with a clear message; it never tries to repair the input.

* Closed: every edge has exactly two incident faces. (A mesh with holes needs to be capped first.)
* Manifold: no edge with more than two faces, no pinch points.
* One component.
* Triangles: quads/ngons are rejected by the C++ API; the Python loader triangulates fans automatically.

Nothing else is required: metriko does not need UVs, normals or a specific scale.

## 3. Remesh

### CLI

```bash
metriko-cli -i input.obj -o quads.obj
```

What you get: `quads.obj` contains the quad mesh (`f a b c d` lines). Read it with any DCC tool or with
`metriko.load_obj`-style readers that accept polygons.

Useful first flags:

```bash
metriko-cli -i input.obj -o quads.obj \
  --scale 0.01 \          # smaller quads (finer). Default 0.02
  --timings \             # print the wall time of every stage
  --patches patches.txt   # one patch label per quad
```

Note: `scale` is a fraction of the **bounding-box diagonal**. Halving it roughly quadruples the quad count, and
strongly increases both the run time and the chance of hitting a degenerate configuration. Start from the
default.

### Python

```python
import numpy as np
import metriko

V, F = metriko.load_obj("input.obj")          # (n,3) float64, (m,3) int32
Q_pos, Q_idx, patch = metriko.remesh(V, F, scale=0.02)

print(Q_idx.shape[0], "quads,", Q_pos.shape[0], "vertices")
metriko.save_obj("quads.obj", Q_pos, Q_idx)
```

`remesh` returns `(Q_pos, Q_idx, patch)`:

* `Q_pos` — `float64 (k, 3)`, the quad vertices;
* `Q_idx` — `int32 (q, 4)`, four vertex indices per quad, counter-clockwise seen from outside;
* `patch` — `int32 (q,)`, the patch (chart) of every quad, or `-1` when the extractor dropped rings (see
  [output.md](output.md)).

With `return_decimated=True` two more arrays are appended — the triangle mesh that was actually remeshed (the
decimated one when decimation ran). With `return_stats=True` a final dict is appended with `used_scale`,
`holes` and `fill_idx`.

### Reading the result

* **`holes`** (or `stats["holes"]`) counts the q-face rings the extractor had to drop. It is informative: those
  rings are closed with fill triangles (`stats["fill_idx"]`), so the mesh is watertight regardless. A value of 0
  means the output is *pure* quads; otherwise it is quad-dominant.
* **`used_scale`**: the scale the returned mesh was actually produced with. It differs from the option when a
  stage failed and the retry ladder coarsened the grid.
* The quad mesh is closed and manifold: 0 boundary edges, 0 non-manifold edges. If you want to verify, count edge
  incidences (each undirected edge must appear exactly twice).

## 4. If it fails

Every failure is stage-tagged:

```
metriko.RemeshError: quantization: the time budget was exceeded during stage quantization
```

In Python, `err.stage` is a `metriko.Stage` enum member (`quantization` here). In C++, catch
`metriko::RemeshFailure` and read `stage`. [errors.md](errors.md) explains what each stage means and what to do.

The most common causes, in order:

1. the mesh is not closed/manifold → rejected at `input`;
2. the mesh is too dense for the budget → the parameterization or the quantization runs long;
3. a degenerate configuration in the T-mesh → a `collapse` or `embedding` error (rare, mesh-specific).

For (2), the retry ladder is already on by default (`--no-retry` disables it): on failure the scale is multiplied
by `1.6` and — when decimation is enabled — the decimation target is cut by 0.6, and the pipeline restarts. The
number of attempts is `--max-attempts` (default 3).
