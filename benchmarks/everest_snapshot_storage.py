"""Lossless shared initial inputs and compact final observables for Everest.

The bank is a normal NPZ archive of content-addressed chunks. Chunks are
shared only when their dtype, shape and bytes agree; no RNG is regenerated.
"""
import hashlib
import io
import json
from pathlib import Path
from zipfile import ZipFile, ZIP_DEFLATED

import numpy as np

from everest_initial_particles import (DTYPES, SCALAR_DTYPES, initial_digest,
                                     validate_initial)

BANK_PATH = Path('inputs') / 'arrays.npz'
CHUNK_SIZE = 20_000
SURVIVOR_FIELDS = ('x', 'px', 'y', 'py', 'delta', 'energy')
FINAL_FIELDS = (*SURVIVOR_FIELDS, 'zeta', 's', 'state', 'at_element',
                'particle_id', 'kick_x')


def chunk_key(array):
    h = hashlib.sha256()
    h.update(array.dtype.str.encode('ascii') + b'\0')
    h.update(str(array.shape).encode('ascii') + b'\0')
    h.update(array.tobytes(order='C'))
    return h.hexdigest()


def store_initial(directory, initial, n):
    """Append only new chunks; return a JSON-serialisable field index."""
    initial = validate_initial(initial, n)
    path = Path(directory) / BANK_PATH
    path.parent.mkdir(parents=True, exist_ok=True)
    refs = {}
    with ZipFile(path, 'a', compression=ZIP_DEFLATED, compresslevel=6) as bank:
        existing = set(bank.namelist())
        for name, array in initial.items():
            chunks = ([array] if array.ndim == 0 else
                      (array[i:i + CHUNK_SIZE] for i in range(0, n, CHUNK_SIZE)))
            keys = []
            for chunk in chunks:
                key = chunk_key(chunk)
                member = key + '.npy'
                if member not in existing:
                    data = io.BytesIO()
                    np.save(data, chunk, allow_pickle=False)
                    bank.writestr(member, data.getvalue())
                    existing.add(member)
                keys.append(key)
            refs[name] = keys
    return np.array(json.dumps(refs, sort_keys=True))


def load_initial(directory, snapshot, info):
    """Reassemble and verify the exact original full initial-particle digest."""
    refs = json.loads(snapshot['initial_refs'].item())
    expected = {f'initial_{key}' for key in DTYPES}
    if set(refs) != expected:
        raise ValueError(f'{info["case"]}: invalid shared initial field index')
    n = info['n']
    arrays = {}
    cache = {}
    with np.load(Path(directory) / BANK_PATH, allow_pickle=False) as bank:
        for name, keys in refs.items():
            scalar = name.removeprefix('initial_') in SCALAR_DTYPES
            count = 1 if scalar else (n + CHUNK_SIZE - 1) // CHUNK_SIZE
            if not isinstance(keys, list) or len(keys) != count:
                raise ValueError(f'{info["case"]}: invalid chunk list for {name}')
            chunks = []
            for i, key in enumerate(keys):
                if (not isinstance(key, str) or len(key) != 64
                        or any(c not in '0123456789abcdef' for c in key)):
                    raise ValueError(f'{info["case"]}: invalid chunk key for {name}')
                if key not in cache:
                    chunk = bank[key]
                    if chunk_key(chunk) != key:
                        raise ValueError(f'{info["case"]}: corrupt initial chunk {name}/{key}')
                    cache[key] = chunk
                chunk = cache[key]
                shape = () if scalar else (min(CHUNK_SIZE, n - i * CHUNK_SIZE),)
                dtype = np.dtype(DTYPES[name.removeprefix('initial_')])
                if chunk.shape != shape or chunk.dtype != dtype:
                    raise ValueError(f'{info["case"]}: corrupt initial chunk {name}/{key}')
                chunks.append(chunk)
            arrays[name] = chunks[0] if scalar else np.concatenate(chunks)
    arrays = validate_initial(arrays, n)
    if initial_digest(arrays, n) != info['input_sha256']:
        raise ValueError(f'{info["case"]}: initial-particle digest mismatch')
    return arrays


def compact_final(full, n):
    """Convert ID-ordered full output to state + survivor-only coordinates."""
    state = np.asarray(full['state'])
    at_element = np.asarray(full['at_element'])
    if state.shape != (n,) or at_element.shape != (n,):
        raise ValueError('Invalid final state or element array shape')
    if not np.array_equal(full['particle_id'], np.arange(n)):
        raise ValueError('Final primary IDs are not in canonical order')
    alive = state > 0
    result = dict(state=state, at_element=at_element)
    for key in SURVIVOR_FIELDS:
        value = np.asarray(full[key])
        if value.shape != (n,):
            raise ValueError(f'Invalid final field {key}')
        result[f'survivor_{key}'] = value[alive]
    return result


def survivor_data(snapshot, initial, n):
    """Return survivor IDs and coordinates, including the exact derived kick."""
    state = np.asarray(snapshot['state'])
    at_element = np.asarray(snapshot['at_element'])
    if state.shape != (n,) or at_element.shape != (n,):
        raise ValueError('Invalid final state or element array shape')
    ids = np.flatnonzero(state > 0)
    result = {}
    for key in SURVIVOR_FIELDS:
        value = np.asarray(snapshot[f'survivor_{key}'])
        if value.shape != (len(ids),):
            raise ValueError(f'Invalid survivor field {key}')
        result[key] = value
    result['kick_x'] = result['px'] - initial['initial_px'][ids]
    return ids, result
