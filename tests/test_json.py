# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

from __future__ import annotations

import io
import gzip
from pathlib import Path
from readline import backend

import numpy as np
import pytest

from xcoll.json import load, dump
from zstandard import backend


def _backend_available(backend: str) -> bool:
    if backend == "json":
        return True
    if backend == "orjson":
        try:
            import orjson  # noqa: F401
            return True
        except Exception:
            return False
    if backend == "msgspec":
        try:
            import msgspec  # noqa: F401
            return True
        except Exception:
            return False
    raise ValueError(backend)


@pytest.mark.parametrize("backends", [
    ("json", "json"), ("json", "orjson"), ("json", "msgspec"),
    ("orjson", "json"), ("orjson", "orjson"), ("orjson", "msgspec"),
    ("msgspec", "json"), ("msgspec", "orjson"), ("msgspec", "msgspec")],
    ids=lambda x: f"{x[0]}-{x[1]}"
)
@pytest.mark.parametrize("stream", ["str", "pathlib", "gz", "bytesio",
                                    "stringio"])
def test_roundtrip_path_str(tmp_path: Path, backends: tuple, stream: str):
    if not _backend_available(backends[0]):
        pytest.skip(f"{backends[0]} not installed")
    if not _backend_available(backends[1]):
        pytest.skip(f"{backends[1]} not installed")

    payload = {
        "int_scalar": np.int64(7),
        "float_scalar": np.float64(3.5),
        "bool_scalar": np.bool_(True),
        "float_arr": np.array([1.0, np.nan, 3.0], dtype=np.float64),
        "int_arr": np.array([1, 2, 3], dtype=np.int64),
        "complex_scalar": 1.25 - 3.0j,
        "complex_arr": (np.array([1.0, 2.0], dtype=np.float64)
                        + 1j * np.array([0.5, -0.25])).astype(np.complex128),
        "nested": {"more": [np.float32(1.0), np.int32(2)]},
    }

    if stream == "str":
        p = str(tmp_path / "lm.json")
    elif stream == "pathlib":
        p = tmp_path / "lm.json"
    elif stream == "gz":
        p = tmp_path / "lm.json.gz"
    elif stream == "bytesio":
        p = io.BytesIO()
    elif stream == "stringio":
        p = io.StringIO()

    dump(payload, p, backend=backends[0])
    if stream in ("bytesio", "stringio"):
        p.seek(0)
    out = load(file=p, backend=backends[1])

    # Scalars
    assert int(out["int_scalar"]) == int(payload["int_scalar"])
    assert float(out["float_scalar"]) == float(payload["float_scalar"])
    assert bool(out["bool_scalar"]) == bool(payload["bool_scalar"])

    # float_arr
    fa = np.asarray(payload["float_arr"])
    fb = np.asarray(out["float_arr"])
    assert np.allclose(fb, fa, equal_nan=True)

    # int_arr
    ia = np.asarray(payload["int_arr"])
    ib = np.asarray(out["int_arr"])
    assert np.array_equal(ib, ia)

    # Complex scalar is always reconstructed
    assert out["complex_scalar"] == payload["complex_scalar"]

    # Complex array is reconstructed (your tag scheme)
    assert isinstance(out["complex_arr"], np.ndarray)
    assert np.array_equal(out["complex_arr"], payload["complex_arr"])

    # Nested list with numpy scalars
    assert float(out["nested"]["more"][0]) == float(payload["nested"]["more"][0])
    assert int(out["nested"]["more"][1]) == int(payload["nested"]["more"][1])

    if stream == "gz":
        # sanity: ensure it is actually gzipped (optional)
        with gzip.open(p, "rb") as f:
            raw = f.read(20)
        assert raw  # file not empty
