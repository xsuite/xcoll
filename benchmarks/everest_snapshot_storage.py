"""Lossless shared initial inputs and compact final observables for Everest.

The bank is a normal NPZ archive of content-addressed chunks. Hashes always
refer to the decoded dtype, shape and bytes, never to the encoding. Version 2
adds reversible codecs; legacy unencoded members remain readable.
"""
import hashlib
import io
import json
from pathlib import Path
from zipfile import ZipFile, ZIP_DEFLATED
import zlib

import numpy as np

from everest_initial_particles import (DTYPES, SCALAR_DTYPES, RNG_FIELDS,
                                       initial_digest, validate_initial)

BANK_PATH = Path('inputs') / 'arrays.npz'
INPUT_STORAGE_VERSION = 2
CHUNK_SIZE = 20_000
SURVIVOR_FIELDS = ('x', 'px', 'y', 'py', 'delta', 'energy')
FINAL_FIELDS = (*SURVIVOR_FIELDS, 'zeta', 's', 'state', 'at_element',
                'particle_id', 'kick_x')
# Each suffix is a permanent codec identifier. Never change an existing codec.
CODECS = ('', '.shuffle-v1', '.constant-v1', '.range-v1', '.rng-v1')
COMPRESSION_LEVEL = 6


def chunk_key(array):
    h = hashlib.sha256()
    h.update(array.dtype.str.encode('ascii') + b'\0')
    h.update(str(array.shape).encode('ascii') + b'\0')
    h.update(array.tobytes(order='C'))
    return h.hexdigest()


def frozen_rng_v1(seed, n):
    """Xtrack rng_set as in v0.118.7, base_rng.h (Apache-2.0).

    All arithmetic is unsigned modulo 2**32; six warm-up steps. This is a
    snapshot decoder, not the runtime RNG. Do not delegate to installed Xtrack.
    """
    seed = int(seed) & 0xffffffff
    s = np.arange(n, dtype=np.uint32) + np.uint32(seed)
    s1 = s * np.uint32(69069)
    s1[s1 < 2] += np.uint32(2)
    s2 = s1 * np.uint32(69069)
    s2[s2 < 8] += np.uint32(8)
    s3 = s2 * np.uint32(69069)
    s3[s3 < 16] += np.uint32(16)
    s4 = s3 * np.uint32(69069)

    def taus(s, a, b, c, d):
        return ((s & np.uint32(c)) << np.uint32(d)) ^ (((s << np.uint32(a)) ^ s) >> np.uint32(b))

    for _ in range(6):
        s1 = taus(s1, 13, 19, 4294967294, 12)
        s2 = taus(s2, 2, 25, 4294967288, 4)
        s3 = taus(s3, 3, 11, 4294967280, 17)
        s4 = s4 * np.uint32(1664525) + np.uint32(1013904223)
    return tuple(np.asarray(s, dtype='<u4') for s in (s1, s2, s3, s4))


def shuffle_bytes(array):
    return np.ascontiguousarray(array).view(np.uint8).reshape(-1, array.dtype.itemsize).T.copy()


def unshuffle_bytes(data, dtype, shape):
    dtype = np.dtype(dtype)
    if data.dtype != np.uint8 or data.shape != (dtype.itemsize, int(np.prod(shape))):
        raise ValueError('Invalid byte-shuffled chunk')
    return data.T.copy().reshape(-1).view(dtype).reshape(shape)


def npy_bytes(array):
    data = io.BytesIO()
    np.save(data, array, allow_pickle=False)
    return data.getvalue()


def same_bytes(a, b):
    return a.dtype == b.dtype and a.shape == b.shape and a.tobytes() == b.tobytes()


def encode_chunk(chunk, *, rng=None):
    """Return a codec and its payload, with exact raw fallback.

    rng is (first seed, lane, expected decoded chunk) or None. Detection
    compares actual bytes, so changed RNG initialisation never loses data.
    """
    if rng is not None:
        seed, lane, expected = rng
        if same_bytes(chunk, expected):
            return '.rng-v1', np.array([seed & 0xffffffff, lane], dtype='<u8')
    if chunk.ndim == 1 and len(chunk) > 1:
        scalar = chunk[:1].copy()
        if same_bytes(chunk, np.repeat(scalar, len(chunk))):
            return '.constant-v1', scalar
        if chunk.dtype.kind in 'iu':
            sequence = np.asarray(np.arange(len(chunk), dtype=chunk.dtype) + scalar[0], dtype=chunk.dtype)
            if same_bytes(chunk, sequence):
                return '.range-v1', scalar
        if chunk.dtype.itemsize > 1:
            shuffled = shuffle_bytes(chunk)
            # Compare the actual compressed NPY payloads, not raw byte counts.
            raw_size = len(zlib.compress(npy_bytes(chunk), COMPRESSION_LEVEL))
            shuffled_size = len(zlib.compress(npy_bytes(shuffled), COMPRESSION_LEVEL))
            if shuffled_size + len('.shuffle-v1') < raw_size:
                return '.shuffle-v1', shuffled
    return '', chunk


def decode_chunk(codec, payload, dtype, shape):
    dtype = np.dtype(dtype)
    if codec == '':
        return payload
    if codec == '.shuffle-v1':
        return unshuffle_bytes(payload, dtype, shape)
    if codec in ('.constant-v1', '.range-v1'):
        if len(shape) != 1 or payload.shape != (1,) or payload.dtype != dtype:
            raise ValueError('Invalid constant/range chunk')
        if codec == '.constant-v1':
            return np.repeat(payload, shape[0])
        if dtype.kind not in 'iu':
            raise ValueError('Invalid integer range chunk')
        return np.asarray(np.arange(shape[0], dtype=dtype) + payload[0], dtype=dtype)
    if codec == '.rng-v1':
        if (len(shape) != 1 or dtype != np.dtype('<u4')
                or payload.shape != (2,) or payload.dtype != np.dtype('<u8')):
            raise ValueError('Invalid RNG chunk')
        seed, lane = map(int, payload)
        if seed > 0xffffffff or not 0 <= lane < 4:
            raise ValueError('Invalid RNG seed or lane')
        return frozen_rng_v1(seed, shape[0])[lane]
    raise ValueError(f'Unknown input chunk codec {codec}')


def store_initial(directory, initial, n, *, rng_seed=None):
    """Append only new chunks; return the unchanged field-to-hash index."""
    initial = validate_initial(initial, n)
    path = Path(directory) / BANK_PATH
    path.parent.mkdir(parents=True, exist_ok=True)
    rng = frozen_rng_v1(rng_seed, n) if rng_seed is not None else None
    rng_names = {f'initial_{key}': i for i, key in enumerate(RNG_FIELDS)}
    refs = {}
    with ZipFile(path, 'a', compression=ZIP_DEFLATED, compresslevel=COMPRESSION_LEVEL) as bank:
        existing = set(bank.namelist())
        for name, array in initial.items():
            chunks = ([array] if array.ndim == 0 else
                      (array[i:i + CHUNK_SIZE] for i in range(0, n, CHUNK_SIZE)))
            keys = []
            for i, chunk in enumerate(chunks):
                key = chunk_key(chunk)
                members = [key + codec + '.npy' for codec in CODECS]
                if not any(member in existing for member in members):
                    hint = None
                    if rng is not None and name in rng_names:
                        lane = rng_names[name]
                        start = i * CHUNK_SIZE
                        hint = (int(rng_seed) + start, lane, rng[lane][start:start + len(chunk)])
                    codec, payload = encode_chunk(chunk, rng=hint)
                    # Every codec must preserve the original hash before writing.
                    if not same_bytes(chunk, decode_chunk(codec, payload, chunk.dtype, chunk.shape)):
                        raise ValueError(f'Input codec changed chunk {name}/{key}')
                    member = key + codec + '.npy'
                    bank.writestr(member, npy_bytes(payload))
                    existing.add(member)
                keys.append(key)
            refs[name] = keys
    return np.array(json.dumps(refs, sort_keys=True))


def load_initial(directory, snapshot, info):
    """Read legacy/new codecs and verify exact original full input digest."""
    refs = json.loads(snapshot['initial_refs'].item())
    expected = {f'initial_{key}' for key in DTYPES}
    if set(refs) != expected:
        raise ValueError(f'{info["case"]}: invalid shared initial field index')
    n = info['n']
    arrays = {}
    cache = {}
    with np.load(Path(directory) / BANK_PATH, allow_pickle=False) as bank:
        members = set(bank.files)
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
                shape = () if scalar else (min(CHUNK_SIZE, n - i * CHUNK_SIZE),)
                dtype = np.dtype(DTYPES[name.removeprefix('initial_')])
                if key not in cache:
                    codecs = [codec for codec in CODECS if key + codec in members]
                    if len(codecs) != 1:
                        raise ValueError(f'{info["case"]}: missing or ambiguous initial chunk {name}/{key}')
                    codec = codecs[0]
                    chunk = decode_chunk(codec, bank[key + codec], dtype, shape)
                    if chunk_key(chunk) != key:
                        raise ValueError(f'{info["case"]}: corrupt initial chunk {name}/{key}')
                    cache[key] = chunk
                chunk = cache[key]
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
