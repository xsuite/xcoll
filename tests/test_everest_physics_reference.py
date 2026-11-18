# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2026.                 #
# ######################################### #
"""CPU physics regression against an accepted (or explicitly selected) reference.

If the codebase has changed and the references will no longer match, one can
regenerate them by running (from the Xcoll root on a clean checkout):

    ts=$(date "+%Y%m%d-%H%M%S")
    python -m benchmarks.perf.prebuild_cpu --skip-tpsa --skip-synrad --n-threads 12 --clear-existing
    python benchmarks/everest_compare_physics.py collect --checkout "$PWD" --use-prebuilt \
        --output /tmp/everest-physics-$ts \
        --particles 20000 --lhc-particles 20000 --turns 10 \
        --case-particles block-tungsten-26gev=500000 \
        --case-particles crystal-amorphous-26gev=100000 \
        --case-particles block-cfc-7000gev=100000

validate it:

    XCOLL_EVEREST_PHYSICS_REFERENCE=/tmp/everest-physics-$ts \
    XCOLL_EVEREST_PHYSICS_FULL_MATRIX=1 \
    XCOLL_EVEREST_PHYSICS_USE_PREBUILT=1 \
    python -m pytest tests/test_everest_physics_reference.py \
        -v --durations=0

and finally store it permanently if it is correct:

    mv /tmp/everest-physics-$ts tests/data/everest_physics_reference
"""


import os
import re
import sys
import pytest
import numpy as np
from pathlib import Path
from types import SimpleNamespace


import xobjects as xo
import xcoll as xc

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'benchmarks'))
from everest_compare_physics import (
    collect, compare, case_specs, load_manifest, snapshot_info,
    DEFAULT_REFERENCE, ELEMENT_ENERGIES, LHC_ENERGIES,
    environment_metadata, environment_differences, assert_compatible,
    lhc_start_names,
)
from everest_interaction_statistics import paired_statistics
from everest_initial_particles import initial_snapshot, initial_digest, restore_initial, read_initial


@pytest.fixture(scope='session')
def physics_context():
    return xo.ContextCpu()


def pytest_generate_tests(metafunc):
    if 'physics_case' not in metafunc.fixturenames:
        return
    configured = os.environ.get('XCOLL_EVEREST_PHYSICS_REFERENCE')
    path = Path(configured).expanduser().resolve() if configured else DEFAULT_REFERENCE
    if not configured and not path.exists():
        # Bootstrap only: once committed, the default is always exercised.
        metafunc.parametrize('physics_case', [pytest.param(None, marks=pytest.mark.skip(
            reason='No committed physics reference yet; generate one or set '
                   'XCOLL_EVEREST_PHYSICS_REFERENCE'))])
        return
    try:
        manifest = load_manifest(path)
        reporter = metafunc.config.pluginmanager.get_plugin('terminalreporter')
        for difference in environment_differences(manifest, environment_metadata()):
            message = f'Physics reference environment difference (informational): {difference}'
            if reporter is not None:
                reporter.write_line(message)
            else:
                print(message)
        sha = manifest['sha']
        expected_sha = os.environ.get('XCOLL_EVEREST_PHYSICS_REF_SHA')
        if not re.fullmatch('[0-9a-f]{40}', sha):
            raise ValueError(f'Invalid reference commit: {sha}')
        if expected_sha is not None and sha != expected_sha:
            raise ValueError(f'Expected reference {expected_sha}, got {sha}')
        # External references may intentionally contain a subset, but the
        # permanent reference and CI must always cover the complete matrix.
        if not configured or os.environ.get('XCOLL_EVEREST_PHYSICS_FULL_MATRIX') == '1':
            full = {case for case, _ in case_specs(SimpleNamespace(
                cases=None, element_energies=ELEMENT_ENERGIES,
                lhc_energies=LHC_ENERGIES, beam=[1, 2], plane=['H', 'V'],
                crystal_modes=[False, True]))}
            if set(manifest['cases']) != full:
                raise ValueError(f'Incomplete full physics matrix: {full ^ set(manifest["cases"])}')
        cases = [pytest.param((path, case, snapshot_info(path, case, manifest)), id=case)
                 for case in manifest['cases']]
    except (OSError, ValueError, KeyError, TypeError) as exc:
        raise pytest.UsageError(f'Invalid physics reference {path}: {exc}') from exc
    metafunc.parametrize('physics_case', cases)


@pytest.mark.everest
@pytest.mark.longtest
def test_everest_physics_against_reference(physics_case, tmp_path, physics_context):
    reference, case, info = physics_case
    n = info['n']
    energy = info['energy_gev']
    candidate = tmp_path / 'candidate'
    # CI runs this trusted harness from main while importing the separately
    # pinned candidate package. Local tests still use their own checkout.
    checkout = Path(os.environ.get('XCOLL_EVEREST_PHYSICS_CHECKOUT',
                                   Path(xc.__file__).resolve().parents[1])).resolve()
    args = SimpleNamespace(
        checkout=checkout, output=candidate, input_reference=reference,
        data=Path(__file__).resolve().parent/'data', cases=[info['family']],
        particles=n, lhc_particles=n, turns=info.get('turns', 10),
        element_energies=[energy], lhc_energies=[energy],
        beam=[info.get('beam', 1)], plane=[info.get('plane', 'H')],
        crystal_modes=[info.get('crystals', False)], interpolation=info.get('interpolation', 0.2),
        events_per_particle=info.get('events_per_particle', 12),
    )
    collect(args, allow_dirty=True, context=physics_context, use_prebuilt=True)
    result = compare(SimpleNamespace(base=reference, candidate=candidate,
                     fraction_tol=0.01, ks_tol=0.02, map_tol=0.05), cases=[case],
                     allow_installed=not (checkout / '.git').exists())
    assert result == 0, f'Physics regression screening failed for {case} against {info["sha"]}'


# Two tests for the everest interaction statistics functions

def test_paired_interaction_statistics():
    ref = np.zeros(100_000)
    ref[:500] = 1
    candidate = ref.copy()
    candidate[500:510] = 1
    result = paired_statistics(ref, candidate)
    assert np.isclose(result['difference'], 10 / 100_000)
    # The variance comes from the ten changed primaries,
    # not all 1,010 recorded events.
    expected = np.sqrt(0.0001 * (1 - 0.0001) / 99_999)
    assert np.isclose(result['stderr'], expected)


def test_channelling_segment_splitting():
    ref_pid = np.array([0, 1])
    ref_ds = np.array([0.0015, 0.0010])
    new_pid = np.array([0, 0, 0, 1, 1])
    new_ds = np.array([0.0005, 0.0005, 0.0005, 0.0005, 0.0005])
    ref = np.bincount(ref_pid, weights=ref_ds, minlength=100)
    new = np.bincount(new_pid, weights=new_ds, minlength=100)
    result = paired_statistics(ref, new)
    assert result['difference'] < 1e-15


def test_initial_particles_roundtrip(tmp_path, physics_context):
    """Public Particles reconstruction preserves inputs and RNG exactly."""
    import xpart as xp
    p = xp.Particles(p0c=450e9, x=[1e-3, 2e-3, 3e-3],
                     px=[1e-6, -2e-6, 3e-6], delta=[0., 1e-4, -2e-4],
                     s=[10., 10., 10.], at_element=2, at_turn=3,
                     start_tracking_at_element=2, _context=physics_context)
    p._init_random_number_generator(seeds=np.array([1, 2, 3], dtype=np.uint32))
    initial = initial_snapshot(p, 3)
    digest = initial_digest(initial, 3)
    path = tmp_path / 'inputs.npz'
    np.savez_compressed(path, **initial)
    info = {'n': 3, 'case': 'test', 'input_sha256': digest}
    with np.load(path, allow_pickle=False) as z:
        loaded = read_initial(z, info)
    restored = restore_initial(loaded, 3, physics_context)
    assert initial_digest(initial_snapshot(restored, 3), 3) == digest
    assert restored.start_tracking_at_element == 2
    loaded['initial_px'][0] += 1e-9
    with pytest.raises(ValueError, match='digest mismatch'):
        read_initial(loaded, info)


def test_physics_environment_is_informational():
    ref = dict(case='test', n=2, input_sha256='abc', sha='a',
               python='3.12', deps={'xtrack': '1'}, platform={'system': 'Darwin'})
    new = dict(ref, sha='b', python='3.14', deps={'xtrack': '2'},
               platform={'system': 'Linux'})
    assert len(environment_differences(ref, new)) == 3
    assert_compatible(ref, new)
    with pytest.raises(ValueError, match='Incompatible physical inputs'):
        assert_compatible(ref, dict(new, input_sha256='changed'))
    with pytest.raises(ValueError, match='Incompatible physical inputs'):
        assert_compatible(ref, dict(new, n=3))


def test_initial_tracking_location_guard():
    initial = {'initial_at_element': np.array([1, 1]),
               'initial_start_tracking_at_element': np.array(1)}
    line = SimpleNamespace(element_names=['a', 'tcp', 'b'])
    assert lhc_start_names(line, initial) == dict(at_element='tcp', start_tracking_at_element='tcp')
    with pytest.raises(ValueError, match='single valid'):
        lhc_start_names(line, dict(initial, initial_at_element=np.array([1, 2])))
