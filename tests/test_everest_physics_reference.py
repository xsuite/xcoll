# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2026.                 #
# ######################################### #
"""Opt-in CPU physics regression against an explicitly generated Git reference."""
import json
import os
from pathlib import Path
import re
import sys
from types import SimpleNamespace

import numpy as np
import pytest


ROOT = Path(__file__).resolve().parents[1]
BENCHMARKS = ROOT / 'benchmarks'
sys.path.insert(0, str(BENCHMARKS))
from compare_everest_physics import collect, compare
from everest_workloads import ELEMENT_CASES


def pytest_generate_tests(metafunc):
    if 'physics_case' not in metafunc.fixturenames:
        return
    configured = os.environ.get('XCOLL_EVEREST_PHYSICS_REFERENCE')
    if not configured:
        metafunc.parametrize('physics_case', [pytest.param(None, marks=pytest.mark.skip(
            reason='Set XCOLL_EVEREST_PHYSICS_REFERENCE to a generated reference directory'))])
        return
    path = Path(configured).expanduser().resolve()
    files = sorted(path.glob('*.npz'))
    if not files:
        raise pytest.UsageError(f'No physics reference snapshots in {path}')
    expected_sha = os.environ.get('XCOLL_EVEREST_PHYSICS_REF_SHA')
    cases = []
    shas = set()
    for file in files:
        case = file.stem
        if case not in (*ELEMENT_CASES, 'lhc'):
            raise pytest.UsageError(f'Unexpected physics reference case: {file}')
        with np.load(file, allow_pickle=False) as data:
            info = json.loads(data['info'].item())
        sha = info.get('sha', '')
        if not re.fullmatch('[0-9a-f]{40}', sha):
            raise pytest.UsageError(f'Invalid reference commit in {file}')
        if info.get('case') != case or info.get('workload_version') != 2:
            raise pytest.UsageError(f'Incompatible reference schema/case in {file}; regenerate it')
        shas.add(sha)
        cases.append(pytest.param((path, case, info), id=case))
    if len(shas) != 1 or (expected_sha is not None and shas != {expected_sha}):
        raise pytest.UsageError(f'Unexpected or mixed reference commit(s): {shas}; expected {expected_sha}')
    metafunc.parametrize('physics_case', cases)


@pytest.mark.everest
@pytest.mark.longtest
def test_everest_physics_against_reference(physics_case, tmp_path):
    reference, case, info = physics_case
    # The accepted reference fixes the inputs; don't silently change them in a test.
    n = info['n']
    candidate = tmp_path / 'candidate'
    args = SimpleNamespace(
        checkout=ROOT, output=candidate, data=ROOT/'tests'/'data', cases=[case],
        particles=n, lhc_particles=n, turns=info.get('turns', 10),
        beam=info.get('beam', 1), plane=info.get('plane', 'H'),
        crystals=info.get('crystals', False), interpolation=info.get('interpolation', 0.2),
        events_per_particle=info.get('events_per_particle', 40),
    )
    # A developer must be able to test uncommitted changes against a clean reference.
    collect(args, allow_dirty=True)
    result = compare(SimpleNamespace(base=reference, candidate=candidate,
        fraction_tol=0.01, ks_tol=0.02, map_tol=0.05, event_tol=0.03), cases=[case])
    assert result == 0, f'Physics regression screening failed for {case} against {info["sha"]}'
