"""Portable initial conditions for the standalone and LHC Everest workloads.

The schema is explicit (rather than a pickle or a dump of the current Xobject
layout). All arrays are stored in primary-ID order and in canonical little-endian
numeric dtypes. Reference and candidate both track a reconstruction of these
inputs, including the original random-generator state and tracking start.
"""
import hashlib

import numpy as np

# Public Particles constructor fields shared by supported Xtrack versions.
FLOAT_FIELDS = (
    'x', 'px', 'y', 'py', 'zeta', 'delta', 's', 'p0c',
    'beta0', 'gamma0', 'pzeta', 'rpp', 'rvv', 'chi', 'charge_ratio',
    'weight', 'ax', 'ay', 'spin_x', 'spin_y', 'spin_z',
    'anomalous_magnetic_moment',
)
INT_FIELDS = ('particle_id', 'parent_particle_id', 'state', 'at_element',
              'at_turn', 'pdg_id')
RNG_FIELDS = ('_rng_s1', '_rng_s2', '_rng_s3', '_rng_s4')
SCALAR_FLOAT_FIELDS = ('mass0', 'q0', 't_sim', 'time_s', 'reference_turns')
SCALAR_INT_FIELDS = ('start_tracking_at_element',)
ARRAY_DTYPES = {**dict.fromkeys(FLOAT_FIELDS, '<f8'),
                **dict.fromkeys(INT_FIELDS, '<i8'),
                **dict.fromkeys(RNG_FIELDS, '<u4')}
SCALAR_DTYPES = {**dict.fromkeys(SCALAR_FLOAT_FIELDS, '<f8'),
                 **dict.fromkeys(SCALAR_INT_FIELDS, '<i8')}
DTYPES = {**ARRAY_DTYPES, **SCALAR_DTYPES}


def initial_snapshot(particles, n):
    """Copy the input before tracking, with no dependency on buffer layout."""
    host = particles._context.nparray_from_context_array
    ids = np.asarray(host(particles.particle_id))
    keep = (ids >= 0) & (ids < n)
    order = np.argsort(ids[keep])
    if not np.array_equal(ids[keep][order], np.arange(n)):
        raise ValueError('Expected exactly one initial particle per primary ID')
    arrays = {f'initial_{key}': np.array(host(getattr(particles, key)),
                                        dtype=dtype)[keep][order]
              for key, dtype in ARRAY_DTYPES.items()}
    arrays.update({f'initial_{key}': np.array(getattr(particles, key), dtype=dtype)
                   for key, dtype in SCALAR_DTYPES.items()})
    return validate_initial(arrays, n)


def validate_initial(snapshot, n):
    """Validate numeric types, shapes and the exact primary-ID alignment."""
    if not isinstance(n, (int, np.integer)) or n < 2:
        raise ValueError('Expected at least two initial primaries')
    result = {}
    for key, dtype in DTYPES.items():
        name = f'initial_{key}'
        arr = np.asarray(snapshot[name])
        expected = (n,) if key in ARRAY_DTYPES else ()
        if arr.shape != expected or arr.dtype != np.dtype(dtype):
            raise ValueError(f'Invalid initial field {name}: {arr.shape}, {arr.dtype}')
        if key == 'reference_turns':
            # NaN is Xtrack's documented uninitialised reference-turn counter.
            if np.isinf(arr):
                raise ValueError(f'Invalid initial field {name}')
            if np.isnan(arr):
                arr = np.array(np.nan, dtype=dtype)  # canonical NaN bit pattern
        elif not np.all(np.isfinite(arr)):
            raise ValueError(f'Non-finite initial field {name}')
        result[name] = arr
    if not np.array_equal(result['initial_particle_id'], np.arange(n)):
        raise ValueError('Initial primary IDs are not in canonical order')
    if np.any(result['initial_state'] <= 0):
        raise ValueError('Expected only active initial primaries')
    if np.any(np.logical_and.reduce([result[f'initial_{k}'] == 0 for k in RNG_FIELDS])):
        raise ValueError('Uninitialised initial RNG state')
    return result


def initial_digest(snapshot, n):
    """Hash all initial fields (including names, dimensions and scalar state)."""
    arrays = validate_initial(snapshot, n)
    h = hashlib.sha256()
    for key in sorted(arrays):
        a = arrays[key]
        h.update(key.encode('ascii') + b'\0')
        h.update(a.dtype.str.encode('ascii') + b'\0')
        h.update(str(a.shape).encode('ascii') + b'\0')
        h.update(a.tobytes(order='C'))
    return h.hexdigest()


def restore_initial(snapshot, n, context):
    """Reconstruct through Xtrack's public API and verify the round trip."""
    import xpart as xp

    arrays = validate_initial(snapshot, n)
    kwargs = {key: (arrays[f'initial_{key}'].item() if key in SCALAR_DTYPES
                    else arrays[f'initial_{key}'].copy()) for key in DTYPES}
    p = xp.Particles.from_dict(kwargs, _context=context, load_rng_state=True)
    if initial_digest(initial_snapshot(p, n), n) != initial_digest(arrays, n):
        raise ValueError('Initial particle reconstruction changed stored inputs')
    return p


def read_initial(snapshot, info):
    """Detect corruption or a missing initial-particle payload before tracking."""
    arrays = validate_initial(snapshot, info['n'])
    if initial_digest(arrays, info['n']) != info['input_sha256']:
        raise ValueError(f'{info["case"]}: initial-particle digest mismatch')
    return arrays
