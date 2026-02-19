# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2026.                 #
# ######################################### #

from __future__ import annotations

import io
import gzip
import numpy as np
from pathlib import Path

# Different json backends
import json as stdjson
try:
    import msgspec
except (ImportError, ModuleNotFoundError):
    msgspec = None
try:
    import orjson
except (ImportError, ModuleNotFoundError):
    orjson = None


# =========================
# === Sanitisation Tags ===
# =========================

_NAN_TAG = "__nan__"
_POSINF_TAG = "__posinf__"
_NEGINF_TAG = "__neginf__"
_COMPLEX_SCALAR_TAG = "__complex__"
_COMPLEX_ARRAY_TAG = "__complex_ndarray__"

_TAG_BYTES = [tag.encode("utf-8") for tag in (
                _NAN_TAG, _POSINF_TAG, _NEGINF_TAG,
                _COMPLEX_SCALAR_TAG, _COMPLEX_ARRAY_TAG)]


# ================
# === Backends ===
# ================

def _validate_backend(backend):
    if backend == "msgspec":
        if msgspec is None:
            raise ImportError("msgspec is not installed. Install it "
                              "with `pip install msgspec`.")
    elif backend == "orjson":
        if orjson is None:
            raise ImportError("orjson is not installed. Install it "
                              "with `pip install orjson`.")
    elif backend != "json":
        raise ValueError(f"Invalid backend {backend!r}. Must be one of "
                         f"'msgspec', 'orjson', or 'json'.")


# ===================
# === I/O helpers ===
# ===================

def _open_stream(file, mode):
    if isinstance(file, io.IOBase):
        return file, False
    path = Path(file)
    if path.suffix == ".gz":
        return gzip.open(path, mode), True
    return open(path, mode), True

def _read(file):
    fp, close = _open_stream(file, "rb")
    try:
        string = fp.read()
    finally:
        if close:
            fp.close()
    return string

def _write(file, bdata):
    fp, close = _open_stream(file, "wb")
    try:
        if isinstance(fp, io.TextIOBase):
            # fp is a text stream, decode bytes to str before writing
            enc = fp.encoding or "utf-8"
            fp.write(bdata.decode(enc))
        else:
            fp.write(bdata)
    finally:
        if close:
            fp.close()


# ========================
# === Encoding Helpers ===
# ========================

def _needs_special_encoding(obj) -> bool:
    stack = [obj]
    while stack:
        x = stack.pop()
        # complex scalar
        if isinstance(x, complex):
            return True
        # numpy arrays
        if isinstance(x, np.ndarray):
            dt = x.dtype
            if np.issubdtype(dt, np.complexfloating):
                # complex arrays always need tagging
                return True
            if np.issubdtype(dt, np.floating):
                if not np.isfinite(x).all():
                    return True
            continue
        # numpy scalar float
        if isinstance(x, np.floating):
            if not np.isfinite(float(x)):
                return True
            continue
        # python float
        if isinstance(x, float):
            if not np.isfinite(x):
                return True
            continue
        # Go down into containers
        if isinstance(x, dict):
            stack.extend(x.values())
        elif isinstance(x, (list, tuple)):
            stack.extend(x)
    return False

def _replace_nonfinite_floats(obj):
    if isinstance(obj, float):
        if np.isnan(obj):
            return {_NAN_TAG: True}
        if np.isposinf(obj):
            return {_POSINF_TAG: True}
        if np.isneginf(obj):
            return {_NEGINF_TAG: True}
    return obj

def _encode_json(obj, *, hook: bool = True):
    # --- Walk dicts/lists ---
    if not hook:
        if isinstance(obj, dict):
            return {k: _encode_json(v, hook=False) for k, v in obj.items()}
        if isinstance(obj, (list, tuple)):
            return [_encode_json(v, hook=False) for v in obj]

    # --- numpy arrays ---
    if isinstance(obj, np.ndarray):
        dt = obj.dtype
        # complex arrays
        if np.issubdtype(dt, np.complexfloating):
            real = obj.real.tolist()
            imag = obj.imag.tolist()
            if not hook:
                real = _encode_json(real, hook=False)
                imag = _encode_json(imag, hook=False)
            return {
                _COMPLEX_ARRAY_TAG: {
                    "dtype": str(dt),
                    "shape": obj.shape,
                    "real": real,
                    "imag": imag,
                }
            }
        # float arrays
        if np.issubdtype(dt, np.floating):
            obj = obj.tolist()
            return obj if hook else _encode_json(obj, hook=False)
        # everything else
        return obj.tolist()

    # --- numpy scalars ---
    if isinstance(obj, np.integer):
        return int(obj)
    if isinstance(obj, np.bool_):
        return bool(obj)
    if isinstance(obj, np.floating):
        obj = float(obj)
        return obj if hook else _replace_nonfinite_floats(obj)

    # --- python floats ---
    if isinstance(obj, float):
        return obj if hook else _replace_nonfinite_floats(obj)

    # --- complex scalar ---
    if isinstance(obj, complex):
        re, im = obj.real, obj.imag
        if not hook:
            re = _replace_nonfinite_floats(re)
            im = _replace_nonfinite_floats(im)
        return {_COMPLEX_SCALAR_TAG: [re, im]}

    if hook:
        raise TypeError

    return obj


# ========================
# === Decoding Helpers ===
# ========================

def _decode_json(obj, *, hook: bool = True):
    # --- deep recursion over containers ---
    if not hook:
        if isinstance(obj, dict):
            # First try to restore this object itself (hook-mode logic)
            restored = _decode_json(obj, hook=True)
            if restored is not obj:
                return restored
            # Otherwise recurse into values
            return {k: _decode_json(v, hook=False) for k, v in obj.items()}
        if isinstance(obj, list):
            return [_decode_json(v, hook=False) for v in obj]
        return obj

    # --- hook-mode: only dicts are expected here ---
    if not isinstance(obj, dict):
        return obj

    # 1) special floats
    if obj.get(_NAN_TAG) is True:
        return float("nan")
    if obj.get(_POSINF_TAG) is True:
        return float("inf")
    if obj.get(_NEGINF_TAG) is True:
        return float("-inf")

    # 2) complex scalar
    if _COMPLEX_SCALAR_TAG in obj:
        real, imag = obj[_COMPLEX_SCALAR_TAG]
        # in hook mode we can restore scalars directly by calling hook=False
        real = _decode_json(real, hook=False)
        imag = _decode_json(imag, hook=False)
        return complex(real, imag)

    # 3) complex ndarray
    if _COMPLEX_ARRAY_TAG in obj:
        md = obj[_COMPLEX_ARRAY_TAG]
        real_list = _decode_json(md["real"], hook=False)
        imag_list = _decode_json(md["imag"], hook=False)
        real = np.asarray(real_list, dtype=float)
        imag = np.asarray(imag_list, dtype=float)
        arr = real + 1j * imag
        dtype = md.get("dtype")
        if dtype is not None:
            arr = arr.astype(np.dtype(dtype), copy=False)
        shape = md.get("shape")
        if shape is not None:
            arr = arr.reshape(tuple(shape))
        return arr

    return obj


# ================
# === JSON API ===
# ================

def dump(data, file, *, indent=None, backend=None):
    if backend is None:
        if orjson is not None:
            backend = "orjson"
        elif msgspec is not None:
            backend = "msgspec"
        else:
            backend = "json"
    _validate_backend(backend)

    # Encode
    if backend == "orjson":
        # no support for arbitrary indent but always use 2 spaces
        opt = orjson.OPT_INDENT_2 if indent is not None else 0
        if not _needs_special_encoding(data):
            # fast path
            bdata = orjson.dumps(data, option=opt | orjson.OPT_SERIALIZE_NUMPY)
        else:
            data = _encode_json(data, hook=False)
            bdata = orjson.dumps(data, option=opt)

    elif backend == "msgspec":
        if not _needs_special_encoding(data):
            # Cannot work via hook because it will write NaNs as null
            bdata = msgspec.json.encode(data)
        else:
            data = _encode_json(data, hook=False)
            bdata = msgspec.json.encode(data)

    else:
        try:
            s = stdjson.dumps(data, allow_nan=False, default=_encode_json,
                              indent=indent)
        except ValueError:
            data = _encode_json(data, hook=False)
            s = stdjson.dumps(data, allow_nan=False, indent=indent)
        bdata = s.encode("utf-8")

    # Write bytes
    _write(file, bdata)


def load(file=None, string=None, backend=None):
    if file is not None and string is not None:
        raise ValueError("Cannot specify both file and string")
    if file is None and string is None:
        raise ValueError("Must specify either file or string")

    if backend is None:
        if orjson is not None:
            backend = "orjson"
        elif msgspec is not None:
            backend = "msgspec"
        else:
            backend = "json"
    _validate_backend(backend)

    if file is not None:
        string = _read(file)
    if isinstance(string, str):
        # Make string into binary
        string = string.encode("utf-8")

    if backend == 'msgspec':
        # dec_hook dsoes not work; use manual post-parsing
        data = msgspec.json.decode(string)
        if any([tag in string for tag in _TAG_BYTES]):
            data = _decode_json(data, hook=False)
        return data

    elif backend == 'orjson':
        # Parse then post-walk to restore complex tags if needed
        data =  orjson.loads(string)
        if any([tag in string for tag in _TAG_BYTES]):
            data = _decode_json(data, hook=False)
        return data

    else:
        # stdlib fallback
        return stdjson.loads(string.decode("utf-8"), object_hook=_decode_json)
