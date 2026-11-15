# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2026.                 #
# ######################################### #
"""CPU physics regression against an accepted (or explicitly selected) reference.

If the codebase has changed and the references will no longer match, one can
regenerate them by running (from the Xcoll root)

    rm -r tests/data/everest_physics_reference
    python benchmarks/compare_everest_physics.py reference --ref main --output tests/data/everest_physics_reference
"""


import os
from pathlib import Path
import re
import sys
from types import SimpleNamespace

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'benchmarks'))
from compare_everest_physics import (
    collect, compare, case_specs, load_manifest, snapshot_info,
    DEFAULT_REFERENCE, ELEMENT_ENERGIES, LHC_ENERGIES,
)


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
def test_everest_physics_against_reference(physics_case, tmp_path):
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
    collect(args, allow_dirty=True)
    result = compare(SimpleNamespace(base=reference, candidate=candidate,
        fraction_tol=0.01, ks_tol=0.02, map_tol=0.05, event_tol=0.03), cases=[case])
    assert result == 0, f'Physics regression screening failed for {case} against {info["sha"]}'
