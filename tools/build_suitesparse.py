#!/usr/bin/env python3
"""Configure, build and install SuiteSparse (CHOLMOD with the GPL supernodal modules) from deps/SuiteSparse.

Usage (from the repository root, after `conan install`):
    python tools/build_suitesparse.py --jobs 6 --conan-dir build/conan

The install prefix defaults to .deps/suitesparse-<os>-<arch>. A stamp (submodule SHA + options) lets repeated
invocations skip the build when nothing changed.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import re
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SOURCE = ROOT / "deps" / "SuiteSparse"
STAMP_NAME = ".metriko-stamp.json"


def host_tag() -> str:
    system = {"Windows": "windows", "Linux": "linux", "Darwin": "macos"}.get(platform.system(), platform.system().lower())
    machine = platform.machine().lower()
    arch = {"amd64": "x86_64", "x86_64": "x86_64", "arm64": "arm64", "aarch64": "arm64"}.get(machine, machine)
    return f"{system}-{arch}"


def submodule_sha() -> str:
    for cmd in (["git", "-C", str(SOURCE), "rev-parse", "HEAD"], ["git", "-C", str(ROOT), "rev-parse", "HEAD:deps/SuiteSparse"]):
        res = subprocess.run(cmd, capture_output=True, text=True)
        if res.returncode == 0 and res.stdout.strip():
            return res.stdout.strip()
    raise SystemExit("cannot determine the SHA of deps/SuiteSparse (is the submodule initialised?)")


def cmake_set(text: str, name: str) -> str | None:
    m = re.search(rf"set\(\s*{re.escape(name)}\s+(.*?)\)", text, re.S)
    return m.group(1).strip() if m else None


def conan_openblas(conan_dir: Path, build_type: str) -> tuple[str, str] | None:
    """(library file, include dir) of the Conan OpenBLAS, read from the CMakeDeps data files."""
    files = sorted(conan_dir.glob(f"*openblas*-{build_type.lower()}-*-data.cmake"))
    if not files:
        return None
    text = files[0].read_text(encoding="utf-8")
    suffix = build_type.upper()
    pkg = re.search(r"set\(\s*\w+_PACKAGE_FOLDER_" + suffix + r'\s+"([^"]+)"', text)
    libs = re.search(r"set\(\s*\w+_LIBS_" + suffix + r"\s+([^)]*)\)", text)
    libdirs = re.search(r"set\(\s*\w+_LIB_DIRS_" + suffix + r'\s+"([^"]*)"', text)
    incdirs = re.search(r"set\(\s*\w+_INCLUDE_DIRS_" + suffix + r'\s+"([^"]*)"', text)
    if not (pkg and libs and libdirs):
        return None
    package = pkg.group(1)

    def expand(s: str) -> str:
        return s.replace("${openblas_PACKAGE_FOLDER_" + suffix + "}", package).replace("${OpenBLAS_PACKAGE_FOLDER_" + suffix + "}", package)

    libdir = Path(expand(libdirs.group(1)).split(";")[0])
    names = [n.strip('"') for n in libs.group(1).split() if n.strip('"')]
    for name in names:
        for candidate in (f"{name}.lib", f"lib{name}.a", f"{name}.a", f"lib{name}.lib"):
            if (libdir / candidate).is_file():
                inc = expand(incdirs.group(1)).split(";")[0] if incdirs else str(Path(package) / "include")
                return str(libdir / candidate), inc
    return None


def run(cmd: list[str], **kw) -> None:
    print("+", " ".join(cmd), flush=True)
    subprocess.run(cmd, check=True, **kw)


def main() -> int:
    tag = host_tag()
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--build-type", default="Release")
    ap.add_argument("--jobs", type=int, default=6)
    ap.add_argument("--conan-dir", type=Path, default=ROOT / "build" / "conan", help="conan install output folder (CMakeDeps data files)")
    ap.add_argument("--prefix", type=Path, default=ROOT / ".deps" / f"suitesparse-{tag}")
    ap.add_argument("--build-dir", type=Path, default=None)
    ap.add_argument("--blas-lib", type=Path, default=None, help="BLAS+LAPACK library (default: Conan OpenBLAS; ignored on macOS)")
    ap.add_argument("--blas-include", type=Path, default=None)
    ap.add_argument("--no-openmp", action="store_true", help="build SuiteSparse without OpenMP")
    ap.add_argument("--force", action="store_true", help="ignore the stamp and rebuild")
    args = ap.parse_args()

    prefix = args.prefix.resolve()
    build_dir = (args.build_dir or ROOT / ".deps" / f"build-suitesparse-{tag}").resolve()
    is_macos = platform.system() == "Darwin"

    defines: dict[str, str] = {
        "CMAKE_BUILD_TYPE": args.build_type,
        "CMAKE_INSTALL_PREFIX": str(prefix),
        "SUITESPARSE_ENABLE_PROJECTS": "cholmod",
        "CHOLMOD_GPL": "ON",
        "SUITESPARSE_USE_CUDA": "OFF",
        "SUITESPARSE_USE_FORTRAN": "OFF",
        "BUILD_SHARED_LIBS": "OFF",
        "BUILD_STATIC_LIBS": "ON",
        "SUITESPARSE_DEMOS": "OFF",
        "BUILD_TESTING": "OFF",
        "SUITESPARSE_USE_OPENMP": "OFF" if args.no_openmp else "ON",
        "SUITESPARSE_USE_STRICT": "ON",
        "SUITESPARSE_USE_64BIT_BLAS": "OFF",
        "CMAKE_POSITION_INDEPENDENT_CODE": "ON",
    }
    if is_macos:
        defines["BLA_VENDOR"] = "Apple"
    else:
        blas_lib, blas_inc = args.blas_lib, args.blas_include
        if blas_lib is None:
            found = conan_openblas(args.conan_dir, args.build_type)
            if not found:
                raise SystemExit(f"cannot locate the Conan OpenBLAS under {args.conan_dir}; run `conan install` first or pass --blas-lib")
            blas_lib, blas_inc = Path(found[0]), Path(found[1])
        blas_lib = Path(blas_lib).resolve()
        if not blas_lib.is_file():
            raise SystemExit(f"BLAS library not found: {blas_lib}")
        cmake_lib = blas_lib.as_posix()
        defines["BLA_VENDOR"] = "OpenBLAS"
        defines["SUITESPARSE_C_TO_FORTRAN"] = "(name,NAME) name##_"
        defines["BLAS_LIBRARIES"] = cmake_lib
        defines["LAPACK_LIBRARIES"] = cmake_lib
        if blas_inc:
            defines["BLAS_INCLUDE_DIRS"] = Path(blas_inc).resolve().as_posix()
            defines["LAPACK_INCLUDE_DIRS"] = Path(blas_inc).resolve().as_posix()

    stamp = {
        "sha": submodule_sha(),
        "defines": defines,
        "blas_size": Path(defines["BLAS_LIBRARIES"]).stat().st_size if "BLAS_LIBRARIES" in defines else 0,
    }
    stamp["hash"] = hashlib.sha256(json.dumps(stamp, sort_keys=True).encode()).hexdigest()
    stamp_file = prefix / STAMP_NAME
    config_exists = any(prefix.glob("lib*/cmake/CHOLMOD/CHOLMODConfig.cmake")) or any(prefix.glob("lib*/cmake/CHOLMOD/cholmod-config.cmake"))
    if not args.force and stamp_file.is_file() and config_exists:
        if json.loads(stamp_file.read_text()).get("hash") == stamp["hash"]:
            print(f"SuiteSparse up to date at {prefix} (stamp {stamp['hash'][:12]}), skipping")
            return 0

    if build_dir.exists():
        shutil.rmtree(build_dir)
    build_dir.mkdir(parents=True)
    cmd = ["cmake", "-S", str(SOURCE), "-B", str(build_dir), "-G", "Ninja"]
    cmd += [f"-D{k}={v}" for k, v in defines.items()]
    run(cmd)
    run(["cmake", "--build", str(build_dir), "--config", args.build_type, "-j", str(args.jobs)])
    run(["cmake", "--install", str(build_dir), "--config", args.build_type])

    prefix.mkdir(parents=True, exist_ok=True)
    stamp_file.write_text(json.dumps(stamp, indent=2))
    print(f"SuiteSparse installed to {prefix}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
