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
        if depth > max_depth:
            raise RecursionError(f"Max depth {max_depth} exceeded (possible cycle).")

        if hasattr(a, 'to_dict') and callable(a.to_dict):
            a = a.to_dict()
        if hasattr(b, 'to_dict') and callable(b.to_dict):
            b = b.to_dict()
        if hasattr(a, 'tolist') and callable(a.tolist):
            a = a.tolist()
        if hasattr(b, 'tolist') and callable(b.tolist):
            b = b.tolist()

        if a is b:
            return True

        # NumPy arrays / array-like
        if isinstance(a, np.ndarray) or isinstance(b, np.ndarray):
            if not (isinstance(a, np.ndarray) and isinstance(b, np.ndarray)):
                return False
            return _arrays_equal(a, b, rtol=rtol, atol=atol)

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
        return a == b

    return _eq(x, y, 0)
