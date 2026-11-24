# metriko documentation

metriko turns a **closed, manifold, single-component triangle mesh** into a **pure-quad mesh** of a requested
size, using the QGP pipeline (cross field → seamless parametrization → motorcycle graph → T-mesh → integer
quantization → quad extraction) followed by a QEx-style extraction and a light surface relaxation.

This documentation is the user guide. If you are looking for the engineering log (what was measured, which fixes
landed, why some avenues were rejected), read [`HANDOFF.md`](HANDOFF.md).

## Contents

| Document | Read it when you want to… |
|---|---|
| [getting-started.md](getting-started.md) | install it and remesh your first mesh (CLI or Python) |
| [parameters.md](parameters.md) | know what every option does, its default and when to change it |
| [pipeline.md](pipeline.md) | understand the stages, what can fail where, and where the time goes |
| [output.md](output.md) | know exactly what comes back (closedness, quad-dominance, fill triangles, labels) |
| [errors.md](errors.md) | interpret a `RemeshError` / `RemeshFailure` and decide what to do next |
| [performance.md](performance.md) | set time budgets, read `--timings`, and understand what dominates |
| [limitations.md](limitations.md) | check whether your mesh is in scope, and what is known not to work |
| [building.md](building.md) | build from source, package wheels/CLI, or reproduce CI locally |

## What metriko guarantees

Given an input that satisfies the contract below:

* **Input contract**: closed (watertight), manifold, one connected component, triangles only, indexed. Anything
  else is rejected up front (`std::invalid_argument` in C++, `ValueError` in Python) — metriko does not silently
  repair meshes.
* **Output**: a quad-dominant mesh with **no boundary edges and no non-manifold edges** (verified on every tested
  mesh), plus a patch label per quad. When the extractor has to drop a few rings (rare, small), the holes they
  leave are closed with fill triangles; the result is still watertight.
* **Determinism**: with the defaults, two runs of the same mesh on the same machine and the same thread count
  give the same quad count, vertex count and hole count. `--fast` trades that for speed.
* **Bounded time**: the CLI and the Python front-end enforce a 300 s budget and a 150 s per-stage timeout by
  default, and retry a failing stage with a coarser grid instead of running forever.

## What it does not guarantee

* A quad for literally every input: meshes with a degenerate T-mesh (see [limitations.md](limitations.md)) fail
  with a stage-tagged error rather than producing a broken mesh.
* A fixed wall-clock time on an arbitrary machine: the 300 s budget is a *cap*, not a promise. Meshes in the
  150-260 s range on an idle machine can exceed the cap on a loaded one. Use `--time-budget`/`time_budget` to
  match your environment, and `--no-retry` when you would rather fail fast than retry coarser.
* Geometry-preserving remeshing of open surfaces, non-manifold meshes or multi-component scenes.

## The 60-second version

```bash
# CLI: a quad mesh of roughly 0.02 × the bounding-box diagonal
metriko-cli -i input.obj -o quads.obj --timings

# add the patch labels (one integer per quad, same order as the faces)
metriko-cli -i input.obj -o quads.obj --patches patches.txt
```

```python
import metriko

V, F = metriko.load_obj("input.obj")
Q_pos, Q_idx, patch = metriko.remesh(V, F, scale=0.02)
metriko.save_obj("quads.obj", Q_pos, Q_idx)
```

`scale` is the only parameter most users ever touch: it is the target quad size as a fraction of the input's
bounding-box diagonal (0.02 → about 50 quads across the model).

## Licence note

The distributed binaries statically link the GPL-licensed CHOLMOD supernodal modules of SuiteSparse, so they are
GPL-affected; see [`../LICENSE`](../LICENSE) and the note in [`../README.md`](../README.md).
