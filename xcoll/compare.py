# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

import numpy as np
from collections.abc import Mapping, Iterator
from dataclasses import is_dataclass, asdict


def _arrays_equal(a, b, *, rtol=1e-7, atol=0.0, verbose=False):
    try:
        a = np.asarray(a)
    except:
        if verbose: print(f"Failed to convert first argument to array:\n{a}")
        return False
    try:
        b = np.asarray(b)
    except:
        if verbose: print(f"Failed to convert second argument to array:\n{b}")
        return False

    if a.shape != b.shape:
        if verbose: print(f"Array shapes differ:\n{a}\nvs\n{b}")
        return False

    a_floaty = np.issubdtype(a.dtype, np.floating) or np.issubdtype(a.dtype, np.complexfloating)
    b_floaty = np.issubdtype(b.dtype, np.floating) or np.issubdtype(b.dtype, np.complexfloating)

    if a_floaty or b_floaty:
        # allclose can't handle object dtype; fall back to exact for object arrays
        if a.dtype == object or b.dtype == object:
            if verbose: print(f"Comparing object arrays exactly:\n{a}\nvs\n{b}")
            return np.array_equal(a, b)
        if verbose: print(f"Comparing float-like arrays with allclose:\n{a}\nvs\n{b}")
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


def deep_equal(x, y, *, rtol=1e-7, atol=0.0, expand_numpy_and_hybridclass=False,
               compare_iterators=False, max_depth=200, verbose=False):
    def _eq(a, b, depth):
        if a is b:
            return True

        if depth > max_depth:
            raise RecursionError(f"Max depth {max_depth} exceeded (possible cycle).")

        if expand_numpy_and_hybridclass:
            if hasattr(a, 'to_dict') and callable(a.to_dict):
                a = a.to_dict()
            if hasattr(b, 'to_dict') and callable(b.to_dict):
                b = b.to_dict()
            if hasattr(a, 'tolist') and callable(a.tolist):
                a = a.tolist()
            if hasattr(b, 'tolist') and callable(b.tolist):
                b = b.tolist()

        # Cheap type mismatch early
        ta, tb = type(a), type(b)
        if ta is not tb:
            # allow numeric-ish comparisons later; don't immediately fail for numpy scalars etc.
            pass

        # NumPy scalars: compare as Python scalars (fast)
        if isinstance(a, np.generic) or isinstance(b, np.generic):
            if not (isinstance(a, np.generic) and isinstance(b, np.generic)):
                if verbose: print(f"Type mismatch at depth {depth} ({type(a)} vs {type(b)}):\n{a}\nvs\n{b}")
                return False
            # float-like uses isclose
            if isinstance(a, (np.floating, np.complexfloating)) or isinstance(b, (np.floating, np.complexfloating)):
                if verbose: print(f"Comparing float-like numpy scalars at depth {depth}:\n{a}\nvs\n{b}")
                return np.isclose(a, b, rtol=rtol, atol=atol, equal_nan=True)
            return a == b

        # NumPy arrays / array-like
        if isinstance(a, np.ndarray) or isinstance(b, np.ndarray):
            if not (isinstance(a, np.ndarray) and isinstance(b, np.ndarray)):
                if verbose: print(f"Type mismatch at depth {depth} ({type(a)} vs {type(b)}):\n{a}\nvs\n{b}")
                return False
            return _arrays_equal(a, b, rtol=rtol, atol=atol, verbose=verbose)

        # to_dict only for non-basic objects (avoid repeated hasattr on primitives)
        if not expand_numpy_and_hybridclass:
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
                if verbose: print(f"Type mismatch at depth {depth} ({type(a)} vs {type(b)}):\n{a}\nvs\n{b}")
                return False
            if a.keys() != b.keys():
                if verbose: print(f"Mapping keys differ at depth {depth}:\n{a.keys()}\nvs\n{b.keys()}")
                return False
            return all(_eq(a[k], b[k], depth + 1) for k in a.keys())

        # Dataclasses (common “object containers”)
        if is_dataclass(a) or is_dataclass(b):
            if not (is_dataclass(a) and is_dataclass(b)):
                if verbose: print(f"Type mismatch at depth {depth} ({type(a)} vs {type(b)}):\n{a}\nvs\n{b}")
                return False
            return _eq(asdict(a), asdict(b), depth + 1)

        # Sequence-like custom objects (supports “objects with __iter__/__getitem__” case)
        if _is_sequence_like(a) or _is_sequence_like(b):
            if not (_is_sequence_like(a) and _is_sequence_like(b)):
                if verbose: print(f"Type mismatch at depth {depth} ({type(a)} vs {type(b)}):\n{a}\nvs\n{b}")
                return False
            if len(a) != len(b):
                if verbose: print(f"Sequence lengths differ at depth {depth} ({len(a)} vs {len(b)}):\n{a} vs {b}")
                return False
            return all(_eq(a[i], b[i], depth + 1) for i in range(len(a)))

        # Iterators/generators: opt-in only (dangerous otherwise)
        if compare_iterators and (isinstance(a, Iterator) or isinstance(b, Iterator)):
            if not (isinstance(a, Iterator) and isinstance(b, Iterator)):
                if verbose: print(f"Type mismatch at depth {depth} ({type(a)} vs {type(b)}):\n{a}\nvs\n{b}")
                return False
            # WARNING: consumes iterators
            for aa, bb in zip(a, b):
                if not _eq(aa, bb, depth + 1):
                    if verbose: print(f"Iterator values differ at depth {depth}:\n{aa}\nvs\n{bb}")
                    return False
            # ensure same length
            try:
                next(a)
                if verbose: print(f"Iterator lengths differ at depth {depth}:\n{a}\nvs\n{b}")
                return False
            except StopIteration:
                pass
            try:
                next(b)
                if verbose: print(f"Iterator lengths differ at depth {depth}:\n{a}\nvs\n{b}")
                return False
            except StopIteration:
                pass
            return True

        # Fallback: atomic comparison
        try:
            if verbose: print(f"Comparing atomic values at depth {depth}:\n{a}\nvs\n{b}")
            return a == b
        except Exception:
            return False

    return _eq(x, y, 0)
