# Building from source

## Prerequisites

* CMake ≥ 3.27 and a C++23 compiler (MSVC 2022, GCC 13, Clang 16).
* [Conan 2](https://conan.io) — provides nanobind, CLI11, OpenBLAS (Linux/Windows) or llvm-openmp (macOS).
* Python ≥ 3.12 with `scikit-build-core`, `nanobind==2.12.0` and [uv](https://docs.astral.sh/uv/) (or plain
  `pip wheel`).
* Git submodules: `deps/SuiteSparse`, `deps/pmp-library` (and `example/polyscope` only for the example viewer).

```bash
git clone --recursive <repo> metriko && cd metriko
```

## One-time dependency setup

```bash
conan install . --profile:host conan/profiles/<your-profile> \
               --profile:build conan/profiles/<your-profile> \
               -s build_type=Release --build=missing
python tools/build_suitesparse.py --jobs 6
```

Available profiles: `conan/profiles/windows-msvc`, `linux-gcc`, `macos-armv8`.

`build_suitesparse.py` compiles CHOLMOD (with the GPL supernodal modules and the vendored METIS) and installs it
into `.deps/suitesparse-<platform>-<arch>/`. It is cached: re-running it is a no-op unless `--force` or the
submodule SHA changes.

## Build the wheel and the CLI

Both artefacts come from **one** CMake build (the wheel's scikit-build configuration already builds the
`metriko-cli` target):

```bash
# use the conan toolchain, and point CMake at the SuiteSparse install
export CMAKE_ARGS="-DCMAKE_TOOLCHAIN_FILE=$PWD/build/conan/conan_toolchain.cmake -DCMAKE_PREFIX_PATH=$PWD/.deps/suitesparse-<platform>-<arch>"
uv build --wheel --no-build-isolation -o wheelhouse/raw

# the CLI, from the same build tree
cmake --install build/cp312-abi3-<platform-tag> --component cli --prefix dist/cli
```

`<platform-tag>` is what scikit-build-core chose (`build/{wheel_tag}` from `pyproject.toml`), e.g.
`cp312-abi3-win_amd64`, `cp312-abi3-linux_x86_64`, `cp312-abi3-macosx_11_0_arm64`.

Repair the wheel so it is self-contained:

```bash
# Windows
delvewheel repair wheelhouse/raw/metriko-*.whl -w wheelhouse
# Linux
auditwheel repair wheelhouse/raw/metriko-*.whl -w wheelhouse
# macOS
delocate-wheel -w wheelhouse wheelhouse/raw/metriko-*.whl
```

## Tests

```bash
ctest --test-dir build/cp312-abi3-<platform-tag> --output-on-failure   # CLI end-to-end (fertility, scale 0.01)
python -m pytest tests/python -v                                       # wheel installed in the current venv
```

`tests/python/test_remesh.py` needs a venv with `numpy` and `pytest` and the repaired wheel installed. It covers
the API contract (shapes, dtypes, the stage enum on errors, `return_stats`, decimation) and is the gate we use
before publishing.

## Windows specifics

* `build/dev.bat <command>` wraps `vcvars64` and sets a parallel build level; use it for every CMake command.
* Ninja + a non-English MSVC locale does not parse `/showIncludes` ("Nota: inclusión del archivo:"), so **header
  edits do not trigger recompiles**. Touch the including `.cpp` (e.g. `src/metriko_api.cpp`) after editing a
  header, or configure with `VSLANG=1033`.
* SuiteSparse needs `SUITESPARSE_C_TO_FORTRAN="(name,NAME) name##_"` (already handled by
  `tools/build_suitesparse.py`).
* OpenBLAS is built with `dynamic_arch=False, target=HASWELL`, so the wheel requires AVX2/FMA. Change it in
  `conanfile.py` for a portable build.

## Continuous integration

`.github/workflows/build.yml` builds the wheel **and** the CLI on the three supported platforms, from a single
compilation per platform, with caches and per-platform artefacts. See the comments in that file for the exact
steps, and [the workflow reference](building.md#what-ci-does) below.

### What CI does

| platform | runner | wheel tag | repair |
|---|---|---|---|
| Linux x86_64 | `manylinux_2_28_x86_64` container on `ubuntu-latest` | `manylinux_2_28_x86_64` | `auditwheel` |
| Windows x86_64 | `windows-2022` | `win_amd64` | `delvewheel` |
| macOS arm64 | `macos-15` | `macosx_*_arm64` | `delocate` |

Caching, in order of impact:

1. **ccache** (`hendrikmuhs/ccache-action`) — the C++ objects of SuiteSparse, pmp and the pipeline are the bulk of
   the build time.
2. **`.deps/`** — the compiled SuiteSparse install, keyed by the SuiteSparse submodule SHA plus
   `tools/build_suitesparse.py`.
3. **`~/.conan2`** — the conan packages (nanobind, CLI11, OpenBLAS), keyed by `conan.lock`/`conanfile.py`/profiles.
4. **pip** — the wheel build dependencies.

Each job uploads two artefacts (`metriko-wheel-<platform>`, `metriko-cli-<platform>`) and runs the test suite on
the packaged wheel before uploading.
