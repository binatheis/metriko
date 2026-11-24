from pathlib import Path

import numpy as np
import pytest

import metriko

FERTILITY = Path(__file__).resolve().parents[2] / "example" / "models" / "fertility.obj"


@pytest.fixture(scope="module")
def fertility():
    return metriko.load_obj(FERTILITY)


def test_version():
    assert metriko.__version__ == "0.1.0"


def test_remesh_fertility(fertility, tmp_path):
    V, F = fertility
    pos, idx, patch = metriko.remesh(V, F, scale=0.01, threads=8, quiet=True)
    assert pos.dtype == np.float64 and pos.ndim == 2 and pos.shape[1] == 3
    assert idx.dtype == np.int32 and idx.ndim == 2 and idx.shape[1] == 4
    assert patch.dtype == np.int32 and patch.shape == (idx.shape[0],)
    assert idx.shape[0] > 0
    assert idx.min() >= 0 and idx.max() < pos.shape[0]
    assert (patch >= 0).all()
    out = tmp_path / "quads.obj"
    metriko.save_obj(out, pos, idx)
    V2, F2 = metriko.load_obj(out)
    assert V2.shape == pos.shape and F2.shape[0] == 2 * idx.shape[0]


def test_decimate_forced(fertility):
    V, F = fertility
    pos, idx, patch, V_used, F_used = metriko.remesh(
        V, F, scale=0.01, threads=8, quiet=True,
        decimate=True, decimate_verts_pct=50, return_decimated=True)
    assert V_used.shape[0] < V.shape[0]
    assert V_used.shape[0] <= 0.55 * V.shape[0]
    assert F_used.ndim == 2 and F_used.shape[1] == 3
    assert idx.shape[0] > 0 and (patch >= 0).all()


def test_decimate_auto_below_threshold(fertility):
    V, F = fertility
    # 28k faces < 300k threshold: decimation must not run
    pos, idx, patch, V_used, F_used = metriko.remesh(
        V, F, scale=0.01, threads=8, quiet=True, return_decimated=True)
    assert V_used.shape == V.shape and F_used.shape == F.shape


def test_decimate_off(fertility):
    V, F = fertility
    pos, idx, patch, V_used, F_used = metriko.remesh(
        V, F, scale=0.01, threads=8, quiet=True,
        decimate=False, decimate_faces_threshold=0, return_decimated=True)
    assert V_used.shape == V.shape and F_used.shape == F.shape


def test_open_mesh_raises_value_error():
    V = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0]], dtype=np.float64)
    F = np.array([[0, 1, 2]], dtype=np.int32)
    with pytest.raises(ValueError, match="boundary"):
        metriko.remesh(V, F)


def test_bad_shapes(fertility):
    V, F = fertility
    with pytest.raises(ValueError):
        metriko.remesh(V[:, :2], F)
    with pytest.raises(ValueError):
        metriko.remesh(V, F[:, :2])
    with pytest.raises(ValueError):
        metriko.remesh(V, F, field="nope")
    with pytest.raises(ValueError, match="out of range"):
        metriko.remesh(V, F + V.shape[0])


def test_stage_enum_and_error_contract():
    assert int(metriko.Stage.qex) == 10
    err = metriko.RemeshError("qex: the ring does not close")
    assert isinstance(err, RuntimeError)
    assert err.stage is metriko.Stage.qex
    assert metriko.RemeshError("no stage here").stage is None


def test_return_stats(fertility):
    V, F = fertility
    out = metriko.remesh(V, F, scale=0.01, threads=2, quiet=True, return_stats=True)
    assert len(out) == 4
    stats = out[3]
    assert stats["used_scale"] == 0.01
    assert stats["holes"] == 0
