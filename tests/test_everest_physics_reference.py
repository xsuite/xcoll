# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2026.                 #
# ######################################### #
"""CPU physics regression against an accepted (or explicitly selected) reference.

If the codebase has changed and the references will no longer match, one can
regenerate them by running (from the Xcoll root)

    rm -r tests/data/everest_physics_reference
    python benchmarks/compare_everest_physics.py reference --ref main --use-prebuilt \
        --output tests/data/everest_physics_reference \
        --particles 20000 --lhc-particles 20000 --turns 10 \
        --case-particles block-tungsten-26gev=500000 \
        --case-particles crystal-amorphous-26gev=100000 \
        --case-particles block-cfc-7000gev=100000
"""


import os
import re
import sys
import pytest
import numpy as np
from pathlib import Path
from types import SimpleNamespace


import xobjects as xo

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'benchmarks'))
from compare_everest_physics import (
    collect, compare, case_specs, load_manifest, snapshot_info,
    DEFAULT_REFERENCE, ELEMENT_ENERGIES, LHC_ENERGIES,
)
from everest_interaction_statistics import paired_statistics


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
    args = SimpleNamespace(
        checkout=ROOT, output=candidate, data=ROOT/'tests'/'data', cases=[info['family']],
        particles=n, lhc_particles=n, turns=info.get('turns', 10),
        element_energies=[energy], lhc_energies=[energy],
        beam=[info.get('beam', 1)], plane=[info.get('plane', 'H')],
        crystal_modes=[info.get('crystals', False)], interpolation=info.get('interpolation', 0.2),
        events_per_particle=info.get('events_per_particle', 40),
    )
    # A developer must be able to test uncommitted changes against a clean reference.
    use_prebuilt = os.environ.get('XCOLL_EVEREST_PHYSICS_USE_PREBUILT') == '1'
    collect(args, allow_dirty=True, context=physics_context, use_prebuilt=use_prebuilt)
    result = compare(SimpleNamespace(base=reference, candidate=candidate,
        fraction_tol=0.01, ks_tol=0.02, map_tol=0.05, event_tol=0.03), cases=[case])
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
