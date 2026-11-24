# AGENTS.md — build & test notes for this fork

## Layout
- `include/metriko/` — upstream header-only library (`lib.h` = full pipeline).
- `include/metriko/api.h` — public `RemeshOptions`/`DecimateOptions`/`RemeshResult`, `remesh()`, `decimate_mesh()`.
- `src/metriko_api.cpp` — input validation (closed, manifold, 1 component, indexed), decimation trigger, option plumbing.
- `src/decimate.cpp` — Eigen ↔ pmp `SurfaceMesh` conversion + `pmp::decimate` call.
- `src/cli/main.cpp` — `metriko-cli` (CLI11, OBJ in/out).
- `python/bindings.cpp` + `python/metriko/__init__.py` — nanobind `metriko._core`, stable ABI.
- `deps/SuiteSparse` (submodule @ e202a5b), `deps/pmp-library` (submodule), `example/polyscope` (submodule).
- `tools/build_suitesparse.py` — builds/installs SuiteSparse into `.deps/suitesparse-<platform>-<arch>/`.
- `conanfile.py`, `conan/profiles/`, `conan.lock` — nanobind 2.12, CLI11 2.6, OpenBLAS 0.3.30 (win/linux), llvm-openmp (macOS).
- `build/dev.bat` — dev-env wrapper (vcvars64 + `CMAKE_BUILD_PARALLEL_LEVEL=6`); `build/uvbuild.bat` — `uv build` wrapper.

## Build (Windows)
```bat
conan install . --profile:host conan/profiles/windows-msvc --profile:build conan/profiles/windows-msvc -s build_type=Release --build=missing
python tools/build_suitesparse.py
build\dev.bat build\uvbuild.bat            :: wheel -> wheelhouse/raw/
build\dev.bat cmake --install build\cp312-abi3-win_amd64 --component cli --prefix dist\cli
uv tool run delvewheel repair wheelhouse\raw\*.whl -w wheelhouse
ctest --test-dir build\cp312-abi3-win_amd64 --output-on-failure
```

## Tests
- `ctest` → `cli.fertility` smoke test (fertility.obj @ scale 0.01).
- `tests/python/test_remesh.py` → run under a venv with `numpy,pytest` and the wheel installed (3.12 and 3.13 both verified).

## Gotchas
- **Reproducible by default**: the scheduling-dependent floating point was removed (SLIM assembly buffers merged in thread order, CHOLMOD factorizes single-threaded via `detail::cholmod_deterministic()`, `EIGEN_DONT_PARALLELIZE`), so the same mesh gives the same quad/vertex/hole counts every run. `--fast` (CLI) / `RemeshOptions::deterministic = false` gives CHOLMOD its threads back (~1.3x on the initial solve) and loses that. `--threads 1` still forces one thread everywhere.
- **Watertight output**: `fill_holes` (patching.h) closes every boundary loop of the extracted quad mesh (4 edges -> quad, longer -> centroid fan); `RemeshResult::fill_idx` carries the triangles, the CLI writes them after the quads, Python exposes them in `stats["fill_idx"]`. `holes` is the informative count of dropped q-face rings -- the output is closed regardless.
- **Diagnosing a slow or failing mesh**: `metriko-cli ... --timings` prints per-stage wall time to stderr (plus the SLIM iteration count and energy and the number of unpaired q-ports). `--time-budget 0` removes the 300 s production cap, `--no-retry` disables the coarser-scale retry, `--decimate-only` writes the decimated mesh without remeshing (fast way to build test meshes: `--decimate --decimate-max-verts 60000 --decimated-out` does not work on failure, `--decimate-only` does).
- **Failure contract**: `remesh()` throws `RemeshFailure{stage}` with `what() == "<stage>: <detail>"`; the stage is one of `RemeshStage` (input, decimation, field, parameterization, tmesh, quantization, collapse, embedding, tutte, slim, qex). The Python side is `metriko.RemeshError.stage` (a `metriko.Stage`).
- **Ninja + Spanish-locale MSVC**: `/showIncludes` output ("Nota: inclusión del archivo:") is not parsed, so header edits do NOT trigger recompiles. After editing a header, touch the .cpp that includes it (e.g. `src/metriko_api.cpp`) or rebuild clean. Setting `VSLANG=1033` before configure should fix it.
- **SuiteSparse on MSVC** needs `SUITESPARSE_C_TO_FORTRAN="(name,NAME) name##_"` (already in build_suitesparse.py) and OpenBLAS `dynamic_arch=False` (`target=HASWELL` → the wheel requires AVX2/FMA-capable CPUs!).
- **CHOLMOD supernodal is GPL** — distributed binaries are GPL-affected (see README).
- `pmp` target does not export its `src/` include dir when `PMP_INSTALL=OFF`; CMakeLists adds it manually.
- Decimation runs when `decimate.enabled && (faces_threshold<=0 || nF >= faces_threshold)`; `pmp::decimate` is triangle-only and calls `garbage_collection()` itself, so returned indices are compact.
- `metriko-cli.exe` in `dist/cli/bin/` is locked while a CLI run is in flight — install/rebuild after killing the run.
- `RemeshResult.used_pos/used_idx` = the mesh actually remeshed (post-decimation).
