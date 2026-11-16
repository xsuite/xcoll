"""Unit tests of the portable compact physics reference format."""
import json
from pathlib import Path
import sys
from types import SimpleNamespace
from zipfile import ZipFile

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'benchmarks'))
from everest_initial_particles import ARRAY_DTYPES, SCALAR_DTYPES, initial_digest
import everest_snapshot_storage as storage


def make_initial(n):
    # Explicit synthetic inputs; the real Xtrack round trip is tested separately.
    arrays = {f'initial_{key}': np.zeros(n, dtype=dtype)
              for key, dtype in ARRAY_DTYPES.items()}
    arrays.update({f'initial_{key}': np.zeros((), dtype=dtype)
                   for key, dtype in SCALAR_DTYPES.items()})
    arrays['initial_particle_id'][:] = np.arange(n)
    arrays['initial_parent_particle_id'][:] = np.arange(n)
    arrays['initial_state'][:] = 1
    arrays['initial_p0c'][:] = 26e9
    arrays['initial_px'][:] = np.arange(n) * 1e-6
    for i in range(1, 5):
        arrays[f'initial__rng_s{i}'][:] = np.arange(n) + i * 100
    return arrays


def input_info(initial, n):
    return dict(n=n, case='test', input_sha256=initial_digest(initial, n))


def test_shared_initial_chunks(tmp_path, monkeypatch):
    monkeypatch.setattr(storage, 'CHUNK_SIZE', 4)
    small, large = make_initial(8), make_initial(12)
    large['initial_y'][:] = 1.0  # a genuine difference must not be deduplicated
    a = storage.store_initial(tmp_path, small, 8)
    b = storage.store_initial(tmp_path, large, 12)
    ar, br = json.loads(a.item()), json.loads(b.item())
    assert ar['initial__rng_s1'] == br['initial__rng_s1'][:2]
    assert ar['initial_px'] == br['initial_px'][:2]
    assert ar['initial_y'] != br['initial_y'][:2]
    # Equal arrays are shared even between different fields, e.g. IDs.
    assert ar['initial_particle_id'] == ar['initial_parent_particle_id']
    with ZipFile(tmp_path / storage.BANK_PATH) as bank:
        expected = {key + '.npy' for refs in (ar, br)
                    for keys in refs.values() for key in keys}
        assert set(bank.namelist()) == expected
        assert len(bank.namelist()) == len(expected)
    for data, refs, n in ((small, a, 8), (large, b, 12)):
        restored = storage.load_initial(tmp_path, {'initial_refs': refs}, input_info(data, n))
        for key in data:
            np.testing.assert_array_equal(data[key], restored[key])
    # A valid but wrong chunk reference must fail the full-input digest check.
    ar['initial_y'] = br['initial_y'][:2]
    with pytest.raises(ValueError, match='digest mismatch'):
        storage.load_initial(tmp_path, {'initial_refs': np.array(json.dumps(ar))},
                             input_info(small, 8))


def make_final(initial, state):
    n = len(state)
    result = dict(state=np.array(state), at_element=np.arange(n), particle_id=np.arange(n))
    result.update({key: np.arange(n, dtype=float) + i
                   for i, key in enumerate((*storage.SURVIVOR_FIELDS, 'zeta', 's'))})
    result['kick_x'] = result['px'] - initial['initial_px']
    return result


@pytest.mark.parametrize('state', [[1, -1, 1, -2], [-1, -1, -1, -1], [1, 1, 1, 1]])
def test_compact_survivors(state):
    n = len(state)
    initial = make_initial(n)
    full = make_final(initial, state)
    compact = storage.compact_final(full, n)
    ids, result = storage.survivor_data(compact, initial, n)
    np.testing.assert_array_equal(ids, np.flatnonzero(full['state'] > 0))
    assert set(compact) == {'state', 'at_element',
                            *(f'survivor_{key}' for key in storage.SURVIVOR_FIELDS)}
    for key in result:
        np.testing.assert_array_equal(result[key], full[key][ids])
    np.testing.assert_array_equal(compact['state'], full['state'])
    np.testing.assert_array_equal(compact['at_element'], full['at_element'])


def test_snapshot_particles_reorders_ids():
    from everest_compare_physics import snapshot_particles
    initial = make_initial(4)
    full = make_final(initial, [1, -1, 1, -2])
    order = np.array([3, 2, 0, 1])
    particles = SimpleNamespace(**{key: value[order] for key, value in full.items()},
                                _context=SimpleNamespace(nparray_from_context_array=np.asarray))
    result = snapshot_particles(particles, 4)
    expected = storage.compact_final(full, 4)
    for key in expected:
        np.testing.assert_array_equal(result[key], expected[key])


def test_migrate_v4(tmp_path):
    from everest_compare_physics import compact_reference, compare, load_manifest
    source = tmp_path / 'old'
    source.mkdir()
    case = 'block-cfc-26gev'
    n = 4
    initial = make_initial(n)
    final = make_final(initial, [1, -1, 1, -2])
    provenance = dict(sha='a' * 40, python='3.12.0', deps={}, platform={'system': 'test'})
    manifest = dict(workload_version=4, cases=[case], **provenance)
    info = dict(input_info(initial, n), case=case, family='block-cfc',
                energy_gev=26, p0c=26e9, workload_version=4, **provenance)
    events = dict(event_code=np.array([16]), event_name=np.array(['Single Diffractive']),
                  event_n=np.array([2]), event_count_by_particle=np.array([[1, 0, 1, 0]], dtype=np.uint32),
                  event_length_code=np.array([], dtype=np.int64),
                  event_length_by_particle=np.empty((0, n)), deposited_energy=np.array(0.0))
    path = source / f'{case}.npz'
    np.savez_compressed(path, **initial, **final, **events, info=np.array(json.dumps(info)))
    (source / 'reference.json').write_text(json.dumps(manifest))
    original = path.read_bytes()
    output = tmp_path / 'new'
    compact_reference(source, output)
    assert path.read_bytes() == original
    assert load_manifest(output) == dict(manifest, workload_version=5)
    assert compare(SimpleNamespace(base=output, candidate=output,
                                   fraction_tol=0.01, ks_tol=0.02, map_tol=0.05)) == 0
    with np.load(output / path.name, allow_pickle=False) as new:
        assert json.loads(new['info'].item()) == dict(info, workload_version=5)
        assert set(new.files) == {'info', 'initial_refs', *events,
                                 *storage.compact_final(final, n)}
    # A corrupted original kick must not silently produce a migrated reference.
    final['kick_x'][0] += 1
    np.savez_compressed(path, **initial, **final, **events, info=np.array(json.dumps(info)))
    invalid = tmp_path / 'invalid'
    with pytest.raises(ValueError, match='changed survivor observable kick_x'):
        compact_reference(source, invalid)
    assert not invalid.exists()
