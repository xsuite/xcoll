"""Regression tests for benchmark loss counting, version metadata, and RF ordering.

No FLUKA, Geant4, GPU, or Xsuite installation required for these checks.
"""
from __future__ import annotations

import ast
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

BENCH = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(BENCH))

from perf.particle_metrics import count_lost_primaries


def particles(state, particle_id):
    return SimpleNamespace(state=np.array(state), particle_id=np.array(particle_id))


def test_unused_particle_capacity_does_not_count_as_loss():
    example = particles(
        [1] * 500 + [-999999999] * 1500,
        list(range(500)) + [-999999999] * 1500,
    )
    assert count_lost_primaries(example, 500) == 0


def test_only_primary_losses_counted_not_secondaries_or_spare_capacity():
    example = particles([1, -335, 1, -333, -999999999], [0, 1, 2, 100, -999999999])
    assert count_lost_primaries(example, 3) == 1


def test_duplicate_or_missing_primary_ids_are_not_silently_ignored():
    with pytest.raises(ValueError, match='distinct original primaries'):
        count_lost_primaries(particles([1, 1, 1, -999999999],
                                       [0, 0, 2, -999999999]), 3)


def test_rf_build_particles_precedes_prepare():
    source = ast.parse((BENCH / 'bench_rf_sweep.py').read_text())
    fn = next(n for n in source.body if isinstance(n, ast.FunctionDef)
              and n.name == 'bench_rf_sweep_tracking_lhc')
    lines = {}
    for node in ast.walk(fn):
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute):
            if node.func.attr in ('build_particles', 'prepare'):
                lines[node.func.attr] = node.lineno
    assert lines['build_particles'] < lines['prepare']


def test_version_is_imported_not_distribution_metadata(monkeypatch):
    import conftest

    original_version = conftest.metadata.version
    monkeypatch.setattr(conftest.metadata, 'version',
                        lambda name: 'stale-installed-version' if name == 'xcoll'
                        else original_version(name))
    monkeypatch.setattr(conftest, 'import_module',
                        lambda name: SimpleNamespace(__version__='0.13.0+working-tree'))
    monkeypatch.setattr(conftest, 'external_versions', lambda: ({}, {}))
    monkeypatch.setattr(conftest, 'machine_name', lambda: 'BlackMetal3')
    data = {}
    conftest.pytest_benchmark_update_json(None, None, data)
    assert data['xsuite_versions']['xcoll'] == '0.13.0+working-tree'
    assert data['machine_name'] == 'BlackMetal3'


def test_corrected_physics_protocol_separates_old_baselines():
    # Old results counted empty slots as losses. Do not compare their
    # physics observables against the corrected primary-only metrics.
    for name in ('bench_geant4.py', 'bench_fluka.py'):
        source = (BENCH / name).read_text()
        assert 'protocol=4, particles=N, engine=' in source
