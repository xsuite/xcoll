# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2026.                 #
# ######################################### #

import numpy as np
from pathlib import Path

import xobjects as xo


try:
    import orjson
except (ImportError, ModuleNotFoundError):
    orjson = None
import json as stdjson


def json_load(path: Path):
    data = path.read_bytes()
    if orjson is None:
        return stdjson.loads(data.decode("utf-8"))
    return orjson.loads(data)


def json_dump(obj, path: Path, *, indent=None):
    if orjson is None:
        return path.write_text(stdjson.dumps(obj, indent=indent, cls=xo.JEncoder), encoding="utf-8")
    opt = 0
    if indent is not None:
        # orjson ignores `indent`; use OPT_INDENT_2 for pretty
        opt |= orjson.OPT_INDENT_2
    try:
        # Fast path: numeric numpy arrays
        bdata = orjson.dumps(obj, option=opt | orjson.OPT_SERIALIZE_NUMPY)
    except TypeError as e:
        # Fallback: convert ndarrays (incl. string/object) to lists
        bdata = orjson.dumps(obj, option=opt, default=_orjson_default)
    path.write_bytes(bdata)


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
