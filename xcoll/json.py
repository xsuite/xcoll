# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2026.                 #
# ######################################### #

import numpy as np

import xobjects as xo


try:
    import orjson
except (ImportError, ModuleNotFoundError):
    orjson = None
import json as stdjson


def json_load(path):
    try:
        if orjson is None:
            with open(path, "r", encoding="utf-8") as fp:
                return stdjson.load(fp)
        else:
            with open(path, "rb") as fp:
                return orjson.loads(fp.read())
    except FileNotFoundError:
        raise FileNotFoundError(f"File {path} not found.") from None


def json_dump(obj, path, *, indent=None):
    if orjson is None:
        with open(path, "w", encoding="utf-8") as fp:
            stdjson.dump(obj, fp, indent=indent, cls=xo.JEncoder)
    else:
        # orjson ignores `indent`; use OPT_INDENT_2 for pretty
        opt = orjson.OPT_INDENT_2 if indent is not None else 0
        try:
            # Fast path: numeric numpy arrays
            bdata = orjson.dumps(obj, option=opt | orjson.OPT_SERIALIZE_NUMPY)
        except TypeError:
            # Fallback: convert ndarrays (incl. string/object) to lists
            bdata = orjson.dumps(obj, option=opt, default=_orjson_default)
        with open(path, "wb") as fp:
            fp.write(bdata)


def _orjson_default(obj):
    # Mirror JEncoder behaviour
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if np.issubdtype(type(obj), np.integer):
        return int(obj)
    if np.issubdtype(type(obj), np.floating):
        return float(obj)
    if np.issubdtype(type(obj), np.bool_):
        return bool(obj)
    raise TypeError

# def atomic_write_bytes(path: Path, bdata: bytes):
#     path = Path(path)
#     tmp = path.with_suffix(path.suffix + ".tmp")
#     with open(tmp, "wb") as f:
#         f.write(bdata)
#         f.flush()
#         os.fsync(f.fileno())
#     os.replace(tmp, path)
