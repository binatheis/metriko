# metriko
![elephant](example/docs/elephant.jpg)

Metriko is a header-only C++ library for quad meshing. It turns a closed triangle mesh into a quad mesh, every quad labelled with the patch of the T-mesh it lies in, in a few seconds for meshes of 10k–30k vertices and well under a minute for 100k+ vertices, without any commercial solver.

## Documentation

The user guide lives in [`docs/`](docs/):

| document | contents |
|---|---|
| [docs/getting-started.md](docs/getting-started.md) | install, first remesh, reading the result |
| [docs/parameters.md](docs/parameters.md) | every option, its default and when to change it |
| [docs/pipeline.md](docs/pipeline.md) | what each stage does and where the time goes |
| [docs/output.md](docs/output.md) | closedness, quad-dominance, fill triangles, patch labels |
| [docs/errors.md](docs/errors.md) | stage-tagged errors and what to do about them |
| [docs/performance.md](docs/performance.md) | measured timings, budgets, tuning |
| [docs/limitations.md](docs/limitations.md) | input domain and known limits |
| [docs/building.md](docs/building.md) | build from source, package, CI |
| [docs/HANDOFF.md](docs/HANDOFF.md) | engineering log (measurements, fixes, rejected avenues) |
## Why Metriko?

Quad meshing is a technique used to convert triangular meshes into quadrilateral meshes. Quadrilateral meshes are highly desirable in various applications due to their topological and geometric advantages. They perform particularly well in tasks such as animation deformation, uv mapping, and CAD operations.

Many existing quad meshing algorithms rely heavily on Mixed Integer Programming (MIP) or Mixed Integer Non-Linear Programming (MINLP). These approaches typically require commercial solvers such as Gurobi, making them inaccessible to many users. While alternative algorithms exist, some of them are numerically unstable or consume several hours of computation, and few have publicly available open-source implementations.

Metriko implements the Quantized Global Parameterization (QGP) algorithm, providing fast and robust quadrilateral parameterization. The algorithm guarantees a valid result and scales linearly with the number of vertices.

## Pipeline
1. Globally optimal 4-rosy tangent field, combed along a seam
2. Integer grid map (IGM) by iterative rounding, locally injective
3. Motorcycle graph of the IGM, giving the T-mesh
4. Quantization of the T-mesh: integer lengths for every arc, no degenerate patch
5. Collapse of the zero-length arcs and snapping of the T-mesh onto mesh vertices
6. Cut along the T-mesh, Tutte embedding of every patch, then SLIM relaxation
7. Quad extraction (QEx) and conformal relaxation of the quads on the input surface
8. Every quad labelled with the T-mesh patch it lies in

## Usage
Meshes with boundaries are currently not supported (planned). Metriko is header-only. It depends on [libigl](https://libigl.github.io/) (core only, which brings Eigen); if [SuiteSparse](https://github.com/DrTimothyAldenDavis/SuiteSparse) and OpenMP are found by CMake they are used automatically and the linear solves get significantly faster than with Eigen alone. The usage is simple: the library offers a single function "compute_quadrangulation".

```cpp
#include "metriko/lib.h"
Eigen::MatrixXd V;   // vertices
Eigen::MatrixXi F;   // triangles
double scale = 0.01; // quad size
auto res = metriko::compute_quadrangulation(V, F, scale);
res->pos; // quad vertices
res->idx; // quad faces
res->val; // patch index of each quad
```

### Demo
The demo in `example/` shows the input and the quads coloured by patch ([polyscope](https://polyscope.run/), as a submodule).
```
git submodule update --init --recursive
cmake -S example -B example/build
cmake --build example/build --target example
./example/build/example example/models/fertility.obj 0.01
```

## Benchmark
Wall time of `compute_quadrangulation`, including the quad relaxation. Apple M5, Apple clang 17 with `-O2`, SuiteSparse and OpenMP. `scale` is the quad edge length as a fraction of the bounding box diagonal.

| mesh | vertices | triangles | scale | quads | time |
|---|---:|---:|---:|---:|---:|
| spot | 2,930 | 5,856 | 0.01 | 7,652 | 1.0 s |
| fandisk | 6,475 | 12,946 | 0.01 | 10,277 | 0.8 s |
| elephant | 12,362 | 24,732 | 0.01 | 5,796 | 3.0 s |
| fertility | 13,971 | 27,954 | 0.01 | 7,885 | 2.1 s |
| bunny | 14,290 | 28,576 | 0.01 | 8,581 | 2.2 s |
| gargoyle | 25,953 | 51,902 | 0.01 | 8,038 | 4.5 s |
| dancer | 12,482 | 24,964 | 0.005 | 6,765 | 2.2 s |
| botijo | 12,441 | 24,898 | 0.005 | 32,691 | 4.0 s |
| rolling stage | 49,988 | 100,000 | 0.005 | 41,298 | 7.4 s |
| igea | 134,345 | 268,686 | 0.005 | 35,914 | 35.1 s |
| bimba | 112,455 | 224,906 | 0.003 | 86,770 | 38.0 s |

## Distributable binaries (CLI + Python wheel)

This fork also builds `metriko-cli` (CLI11) and a stable-ABI Python wheel (`metriko`,
nanobind, `cp312-abi3`, works on Python ≥ 3.12) on top of `metriko_core`, the same
static library both frontends link against. Dependencies: SuiteSparse/CHOLMOD and
OpenMP via `deps/SuiteSparse` + `tools/build_suitesparse.py`, nanobind/CLI11/OpenBLAS
via Conan, libigl/Eigen via FetchContent, and [pmp-library](https://github.com/pmp-library/pmp-library)
via `deps/pmp-library` for decimation.

```bash
conan install . --profile:host conan/profiles/windows-msvc --profile:build conan/profiles/windows-msvc -s build_type=Release --build=missing
python tools/build_suitesparse.py   # once per platform; installs into .deps/
pip install scikit-build-core && python -m build --wheel --no-build-isolation
cmake --install build/cp312-abi3-win_amd64 --component cli --prefix dist/cli
```

### metriko-cli

```
metriko-cli -i input.obj -o output.obj [--scale 0.02] [--field smoothest|curvature_aligned]
            [--no-refine] [--slim-max-iter 50] [--slim-rel-tol 1e-4] [--threads 0]
            [--patches labels.txt] [--quiet] [--timings] [--time-budget 300]
            [--no-retry] [--coarsen-factor 1.6] [--max-attempts 3]
            [--decimate-only] [decimation flags]
```

### Python

```python
import metriko
V, F = metriko.load_obj("input.obj")
Q_pos, Q_idx, patch = metriko.remesh(V, F, scale=0.02)
```

## Production behaviour: time budget, retries, holes, determinism

- **Time budget.** `time_budget_seconds` (CLI `--time-budget`, Python `time_budget`, default
  **300 s** in both front-ends; `0` = unlimited in the C++ API) is checked between pipeline
  stages: a run stops at the next stage boundary and reports the stage it was in. A single
  stage can still overrun it by its own duration, so treat it as a soft cap on the total.
- **Retries with a coarser scale.** Integer-grid maps get fewer degenerate configurations as
  the grid gets coarser, so a failed stage is retried with `scale * coarsen_factor`
  up to `max_attempts` (default 3, factor 1.6). A retry is only started when the previous
  attempt's cost still fits the remaining budget. `--no-retry` / `retry=False` disables it.
  The scale that produced the returned mesh is reported (`result.used_scale`, CLI output).
- **Stage-tagged errors.** A failure raises `RemeshFailure` (C++), `metriko.RemeshError`
  (Python, subclass of `RuntimeError`, with `.stage` a `metriko.Stage`). The message is
  `"<stage>: <detail>"`. Stages: `field, parameterization, tmesh, quantization, collapse,
  embedding, tutte, slim, qex`.
- **Holes instead of failure.** The quad extractor drops a q-face ring that cannot be closed
  (a folded or collapsed uv triangle) instead of aborting — the reference implementation
  (libQEx) does the same. `result.holes` counts them (CLI prints a warning, Python
  `return_stats=True` gives `{"used_scale", "holes"}`); more than 2% of the rings lost is
  reported as a `qex` failure.
- **Timings.** `--timings` / `timings=True` prints the wall time of every stage to stderr
  (with the SLIM iteration count and energy), which is how to find the hot stage of a mesh.
- **Determinism.** The pipeline is **not bit-reproducible run to run** (OpenMP reductions in
  SLIM and in the quantization's basis loops, multithreaded CHOLMOD/BLAS), and a marginal mesh
  can pass one run and fail the next. Use `--threads 1` / `threads=1` when a repro must be
  exact; it is slower in wall-clock time but deterministic in practice.

## Decimation (optional preprocessing, on by default)

Inputs well beyond the benchmarked ~270k faces can take a very long time or hit edge
cases in the T-mesh collapse, so meshes with **≥ 300k faces** are automatically
decimated to **90k vertices** with `pmp::decimate` (greedy halfedge collapse) before
the pipeline runs. The decimated triangle mesh stays closed, manifold and
single-component, so all algorithm guarantees still hold. `result.used_pos` /
`used_idx` (Python: `return_decimated=True`) expose the mesh that was actually
remeshed. Disable with `--no-decimate` / `decimate=False`; force on smaller meshes
with `--decimate` / `decimate=True`.

| CLI flag / Python kwarg (prefix `decimate-` / `decimate_`) | default | meaning (maps to `pmp::decimate`) |
|---|---|---|
| `verts-pct` / `verts_pct` | 0 | % of input vertices to keep; 0 = use `max-verts` cap |
| `max-verts` / `max_verts` | 90000 | automatic target vertex cap |
| `faces-threshold` / `faces_threshold` | 300000 | input face count that triggers decimation; 0 = always |
| `normal-dev` / `normal_dev` | 0 | max normal deviation in degrees |
| `aspect-ratio` / `aspect_ratio` | 0 | max triangle aspect ratio |
| `min-edge-len` / `min_edge_len` | 0 | minimum target edge length (mesh units) |
| `max-valence` / `max_valence` | 0 | maximum vertex valence |
| `hausdorff` / `hausdorff` | 0 | max deviation from the original surface (mesh units) |
| `seam-thresh` / `seam_thresh` | -1 | texture seam threshold (<0: disabled; needs `h:tex`, not loaded from OBJ) |
| `seam-dev` / `seam_dev` | 0 | max texture seam deviation in degrees |

`0` for a constraint means unconstrained. The CLI also accepts `--decimated-out out.obj`
to save the decimated triangle mesh.

## Licensing note

Binaries built with the default configuration link CHOLMOD's **GPL** supernodal module
(required by `Eigen::CholmodSupernodalLLT`) and are therefore covered by the GPL-2.0+.
Rebuilding without SuiteSparse (Eigen solvers only) removes that obligation; the
upstream library itself stays MPL-2.0.

## References
- Bommes et al. [Mixed-integer quadrangulation](https://doi.org/10.1145/1531326.1531383). SIGGRAPH 2009.
- Bommes et al. [Integer-grid maps for reliable quad meshing](https://doi.org/10.1145/2461912.2462014). SIGGRAPH 2013.
- Campen et al. [Quantized global parametrization](https://doi.org/10.1145/2816795.2818140). SIGGRAPH Asia 2015.
- Ebke et al. [QEx: Robust quad mesh extraction](https://doi.org/10.1145/2508363.2508372). SIGGRAPH Asia 2013.
- Eppstein et al. [Motorcycle graphs: canonical quad mesh partitioning](https://doi.org/10.1111/j.1467-8659.2008.01288.x). SGP 2008.
- Knöppel et al. [Globally optimal direction fields](https://doi.org/10.1145/2461912.2462005). SIGGRAPH 2013.
- Lyon et al. [Parametrization quantization with free boundaries for trimmed quad meshing](https://doi.org/10.1145/3306346.3323019). SIGGRAPH 2019.
- Lyon et al. [Quad layouts via constrained T-mesh quantization](https://doi.org/10.1111/cgf.142634). Eurographics 2021.
- Myles et al. [Robust field-aligned global parametrization](https://doi.org/10.1145/2601097.2601154). SIGGRAPH 2014.
- Rabinovich et al. [Scalable locally injective mappings](https://doi.org/10.1145/2983621). TOG 2017.
- Bouaziz et al. [Shape-Up: Shaping discrete geometry with projections](https://doi.org/10.1111/j.1467-8659.2012.03171.x). SGP 2012.

## Output guarantees (production profile)

- **Watertight, quad-dominant output.** The extractor may drop q-face rings where the uv parametrization
  folds; the holes they leave are closed before returning: a four-edge loop becomes a quad again, a longer one
  is filled with a triangle fan (`result.fill_idx`, CLI writes them into the same OBJ, Python
  `stats["fill_idx"]`). `result.holes` / `stats["holes"]` reports how many rings were dropped, it is
  informative only -- the mesh is closed either way. All tested meshes (manotas, gatubelo, Azure, GEO-head,
  drag) come out with zero boundary edges and zero non-manifold edges.
- **Reproducible by default.** The stages that were scheduling-dependent are deterministic (SLIM assembly,
  CHOLMOD single-threaded factorization, Eigen internal threading off): two runs of the same mesh give the same
  quad count, vertex count and hole count. `--fast` / `deterministic=False` gives CHOLMOD its threads back
  (~1.3x on the initial solve) at the cost of run-to-run variation.
- **Bounded time.** `--time-budget` (300 s in the CLI and Python front-ends) and `--stage-timeout` (150 s) stop
  a run that is not going to finish, and a failed stage is retried with a coarser scale and a 0.6x decimation
  (`--no-retry` disables it). `--timings` prints the wall time of every stage and the SLIM/LM/rounding
  internals.
- **Failures say where.** `RemeshFailure` (C++) / `metriko.RemeshError` (Python, `.stage` is a
  `metriko.Stage`) carries the stage: `field, parameterization, tmesh, quantization, collapse, embedding,
  tutte, slim, qex`.
- **Input domain.** Closed (watertight), manifold, single-component triangle meshes. Anything else is rejected
  up front with `ValueError` / `std::invalid_argument`.
- **Known limit.** Meshes whose T-mesh acquires a snapped "chord" between non-adjacent vertices still fail in
  the collapse (`boundary segment ... lies in no face`); one production mesh (dragonita) hits it. Everything
  else in the test set passes in under 300 s.
