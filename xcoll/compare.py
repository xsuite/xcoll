# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

import numpy as np
from collections.abc import Mapping, Iterator
from dataclasses import is_dataclass, asdict


def _arrays_equal(a, b, *, rtol=1e-7, atol=0.0):
    try:
        a = np.asarray(a)
    except:
        return False
    try:
        b = np.asarray(b)
    except:
        return False

    if a.shape != b.shape:
        return False

    a_floaty = np.issubdtype(a.dtype, np.floating) or np.issubdtype(a.dtype, np.complexfloating)
    b_floaty = np.issubdtype(b.dtype, np.floating) or np.issubdtype(b.dtype, np.complexfloating)

    if a_floaty or b_floaty:
        # allclose can't handle object dtype; fall back to exact for object arrays
        if a.dtype == object or b.dtype == object:
            return np.array_equal(a, b)
        return np.allclose(a, b, rtol=rtol, atol=atol, equal_nan=True)

    return np.array_equal(a, b)


def _is_sequence_like(obj):
    # "sequence-like" is safer than "iterable":
    # it suggests stable ordering + finite length + repeatable access
    if isinstance(obj, (str, bytes, bytearray)):
        return False
    if isinstance(obj, (np.ndarray, Mapping)):
        return False
    return hasattr(obj, "__len__") and hasattr(obj, "__getitem__")


def deep_equal(x, y, *, rtol=1e-7, atol=0.0, compare_iterators=False, max_depth=200):
    def _eq(a, b, depth):
        if a is b:
            return True

        if depth > max_depth:
            raise RecursionError(f"Max depth {max_depth} exceeded (possible cycle).")

        # Cheap type mismatch early
        ta, tb = type(a), type(b)
        if ta is not tb:
            # allow numeric-ish comparisons later; don't immediately fail for numpy scalars etc.
            pass

        # NumPy scalars: compare as Python scalars (fast)
        if isinstance(a, np.generic) or isinstance(b, np.generic):
            if not (isinstance(a, np.generic) and isinstance(b, np.generic)):
                return False
            # float-like uses isclose
            if isinstance(a, (np.floating, np.complexfloating)) or isinstance(b, (np.floating, np.complexfloating)):
                return np.isclose(a, b, rtol=rtol, atol=atol, equal_nan=True)
            return a == b

        # NumPy arrays / array-like
        if isinstance(a, np.ndarray) or isinstance(b, np.ndarray):
            if not (isinstance(a, np.ndarray) and isinstance(b, np.ndarray)):
                return False
            return _arrays_equal(a, b, rtol=rtol, atol=atol)

        # to_dict only for non-basic objects (avoid repeated hasattr on primitives)
        if not isinstance(a, (str, bytes, bytearray, int, float, bool, type(None), tuple, list, dict)):
            td = getattr(a, "to_dict", None)
            if td is not None and callable(td):
                a = td()
        if not isinstance(b, (str, bytes, bytearray, int, float, bool, type(None), tuple, list, dict)):
            td = getattr(b, "to_dict", None)
            if td is not None and callable(td):
                b = td()

        # Mappings (dict and friends)
        if isinstance(a, Mapping) or isinstance(b, Mapping):
            if not (isinstance(a, Mapping) and isinstance(b, Mapping)):
                return False
            if a.keys() != b.keys():
                return False
            return all(_eq(a[k], b[k], depth + 1) for k in a.keys())

        # Dataclasses (common “object containers”)
        if is_dataclass(a) or is_dataclass(b):
            if not (is_dataclass(a) and is_dataclass(b)):
                return False
            return _eq(asdict(a), asdict(b), depth + 1)

        # Sequence-like custom objects (supports “objects with __iter__/__getitem__” case)
        if _is_sequence_like(a) or _is_sequence_like(b):
            if not (_is_sequence_like(a) and _is_sequence_like(b)):
                return False
            if len(a) != len(b):
                return False
            return all(_eq(a[i], b[i], depth + 1) for i in range(len(a)))

        # Iterators/generators: opt-in only (dangerous otherwise)
        if compare_iterators and (isinstance(a, Iterator) or isinstance(b, Iterator)):
            if not (isinstance(a, Iterator) and isinstance(b, Iterator)):
                return False
            # WARNING: consumes iterators
            for aa, bb in zip(a, b):
                if not _eq(aa, bb, depth + 1):
                    return False
            # ensure same length
            try:
                next(a)
                return False
            except StopIteration:
                pass
            try:
                next(b)
                return False
            except StopIteration:
                pass
            return True

        # Fallback: atomic comparison
        try:
            return a == b
        except Exception:
            return False

    return _eq(x, y, 0)
