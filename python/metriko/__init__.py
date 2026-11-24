from __future__ import annotations

import re
from enum import IntEnum

import numpy as np

from . import _core

__version__ = "0.1.0"
__all__ = ["remesh", "load_obj", "save_obj", "Stage", "RemeshError", "__version__"]


class Stage(IntEnum):
    """The pipeline stage that raised :class:`RemeshError`."""

    input = 0
    decimation = 1
    field = 2
    parameterization = 3
    tmesh = 4
    quantization = 5
    collapse = 6
    embedding = 7
    tutte = 8
    slim = 9
    qex = 10


RemeshError = _core.RemeshError


def _error_stage(self):
    """The stage named by the message, as a :class:`Stage`, or ``None``."""
    m = re.match(r"^([a-z_]+): ", str(self))
    if m is None:
        return None
    try:
        return Stage[m.group(1)]
    except KeyError:
        return None


RemeshError.stage = property(_error_stage)


def remesh(
    V,
    F,
    scale: float = 0.02,
    field: str = "smoothest",
    refine: bool = True,
    slim_max_iter: int = 30,
    slim_rel_tol: float = 1e-4,
    threads: int = 0,
    quiet: bool = False,
    decimate=None,
    decimate_verts_pct=None,
    decimate_max_verts=None,
    decimate_faces_threshold=None,
    decimate_normal_dev=None,
    decimate_aspect_ratio=None,
    decimate_min_edge_len=None,
    decimate_max_valence=None,
    decimate_hausdorff=None,
    decimate_seam_thresh=None,
    decimate_seam_dev=None,
    return_decimated: bool = False,
    retry: bool = True,
    max_attempts: int = 3,
    coarsen_factor: float = 1.6,
    timings: bool = False,
    time_budget: float = 300.0,
    stage_timeout: float = 150.0,
    return_stats: bool = False,
):
    """Quad-remesh a closed, manifold, single-component triangle mesh.

    Returns ``(Q_pos float64 (k, 3), Q_idx int32 (q, 4), patch int32 (q,))``,
    where ``patch`` is the patch label of every quad. With
    ``return_decimated=True`` two more arrays are appended:
    ``(V_used, F_used)``, the triangle mesh that was actually remeshed
    (the decimated mesh when decimation ran, otherwise the input).

    Decimation (pmp-library) runs automatically when the input has at least
    ``decimate_faces_threshold`` faces (default 300k). ``decimate=True``
    forces it on smaller meshes, ``decimate=False`` disables it entirely.
    ``decimate_verts_pct`` is the percentage of input vertices to keep;
    ``None`` means an automatic cap of ``decimate_max_verts`` vertices (90k by default).
    ``decimate_normal_dev`` (degrees), ``decimate_aspect_ratio``,
    ``decimate_min_edge_len`` (mesh units), ``decimate_max_valence``,
    ``decimate_hausdorff`` (mesh units), ``decimate_seam_thresh`` and
    ``decimate_seam_dev`` (degrees) map to the corresponding
    ``pmp::decimate`` parameters; 0/None leaves them unconstrained.

    With ``retry=True`` (the default) a failed stage is retried with a coarser
    scale (``scale *= coarsen_factor``), up to ``max_attempts`` attempts;
    dense integer-grid maps are the main source of degenerate configurations.
    ``timings=True`` prints the wall time of every stage to stderr. When the
    algorithm fails it raises :class:`RemeshError`, whose ``stage`` attribute
    (a :class:`Stage`) says where.
    """
    V = np.ascontiguousarray(V, dtype=np.float64)
    F = np.ascontiguousarray(F, dtype=np.int32)
    if V.ndim != 2 or V.shape[1] != 3:
        raise ValueError(f"V must have shape (n, 3), got {V.shape}")
    if F.ndim != 2 or F.shape[1] != 3:
        raise ValueError(f"F must have shape (m, 3), got {F.shape}")
    return _core.remesh(
        V, F, float(scale), str(field), bool(refine), int(slim_max_iter),
        float(slim_rel_tol), int(threads), bool(quiet),
        None if decimate is None else bool(decimate),
        None if decimate_verts_pct is None else float(decimate_verts_pct),
        None if decimate_max_verts is None else int(decimate_max_verts),
        None if decimate_faces_threshold is None else int(decimate_faces_threshold),
        None if decimate_normal_dev is None else float(decimate_normal_dev),
        None if decimate_aspect_ratio is None else float(decimate_aspect_ratio),
        None if decimate_min_edge_len is None else float(decimate_min_edge_len),
        None if decimate_max_valence is None else int(decimate_max_valence),
        None if decimate_hausdorff is None else float(decimate_hausdorff),
        None if decimate_seam_thresh is None else float(decimate_seam_thresh),
        None if decimate_seam_dev is None else float(decimate_seam_dev),
        bool(return_decimated),
        bool(retry), int(max_attempts), float(coarsen_factor), bool(timings), float(time_budget),
        float(stage_timeout), bool(return_stats),
    )


def load_obj(path):
    """Read an OBJ file. Polygons are triangulated as fans; vt/vn are ignored.

    Returns (V float64 (n, 3), F int32 (m, 3)).
    """
    verts: list[list[float]] = []
    faces: list[list[int]] = []
    with open(path, "r", encoding="utf-8") as fh:
        for line in fh:
            parts = line.split()
            if not parts:
                continue
            if parts[0] == "v":
                verts.append([float(x) for x in parts[1:4]])
            elif parts[0] == "f":
                idx = []
                for tok in parts[1:]:
                    i = int(tok.split("/")[0])
                    idx.append(i - 1 if i > 0 else len(verts) + i)
                for k in range(1, len(idx) - 1):
                    faces.append([idx[0], idx[k], idx[k + 1]])
    V = np.asarray(verts, dtype=np.float64).reshape(-1, 3)
    F = np.asarray(faces, dtype=np.int32).reshape(-1, 3)
    return V, F


def save_obj(path, V, Q):
    """Write a quad mesh (V (n, 3), Q (q, 4), 0-based) to an OBJ file."""
    V = np.asarray(V, dtype=np.float64)
    Q = np.asarray(Q, dtype=np.int64) + 1
    with open(path, "w", encoding="utf-8") as fh:
        for x, y, z in V:
            fh.write(f"v {x:.17g} {y:.17g} {z:.17g}\n")
        for a, b, c, d in Q:
            fh.write(f"f {a} {b} {c} {d}\n")
