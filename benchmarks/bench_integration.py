"""Installation and line-tool benchmarks on the real LHC Run 3 sequence.

The machine, collimator database, and benchmarked operation are visible in each
benchmark. No external scenario module is needed. Loading input data is excluded
from the timed installation measurements, and material tracking is not involved.
"""

import hashlib
from pathlib import Path

import numpy as np
import pytest
import xtrack as xt
import xcoll as xc

DATA = Path(__file__).resolve().parents[1] / 'tests' / 'data'
ROUNDS = 5


def _require_inputs(*paths):
    """Small local utility: missing optional test data should skip cleanly."""
    for path in paths:
        if not path.is_file():
            pytest.skip(f'LHC input dataset is unavailable: {path}')


@pytest.mark.quick
@pytest.mark.benchmark(group='LHC Run 3 | installation')
def bench_lhc_install_everest_b1(benchmark):
    """Install the B1 IR7 Everest collimation system, including the crystals."""
    sequence = DATA / 'sequence_lhc_run3_b1.json'
    database = DATA / 'colldb_lhc_run3_ir7.yaml'
    _require_inputs(sequence, database)
    sequence_hash = hashlib.sha256(sequence.read_bytes()).hexdigest()
    database_hash = hashlib.sha256(database.read_bytes()).hexdigest()

    def setup():
        # Fresh, uninstalled line for every round; none of this is timed.
        env = xt.load(sequence)
        line = env['lhcb1']
        colldb = xc.CollimatorDatabase.from_yaml(database, beam=1, ignore_crystals=False)
        original_length = line.get_length()
        return (line, colldb, original_length), {}

    def install(line, colldb, original_length):
        colldb.install_everest_collimators(line=line)
        return line, colldb, original_length

    line, colldb, original_length = benchmark.pedantic(
        install, setup=setup, rounds=ROUNDS, warmup_rounds=1,
    )
    installed = len(line.xcoll.collimators)
    assert installed == len(colldb.collimator_names)
    assert installed >= 10
    assert np.isclose(line.get_length(), original_length)
    assert any(isinstance(line[name], xc.EverestCrystal)
               for name in colldb.collimator_names)
    benchmark.extra_info.update({
        'protocol': 2, 'phase': 'installation', 'beam': 1, 'crystals': True,
        'physics_case': 'install-everest-ir7', 'installed_collimators': installed,
        'sequence_sha256': sequence_hash,
        'collimator_database_sha256': database_hash,
        'timed': 'colldb.install_everest_collimators(line=line)',
    })


@pytest.mark.full
@pytest.mark.benchmark(group='LHC Run 3 | installation')
def bench_lhc_install_everest_b2(benchmark):
    """Install the B2 IR7 Everest collimation system without crystals."""
    sequence = DATA / 'sequence_lhc_run3_b2.json'
    database = DATA / 'colldb_lhc_run3_ir7.yaml'
    _require_inputs(sequence, database)
    sequence_hash = hashlib.sha256(sequence.read_bytes()).hexdigest()
    database_hash = hashlib.sha256(database.read_bytes()).hexdigest()

    def setup():
        env = xt.load(sequence)
        line = env['lhcb2']
        colldb = xc.CollimatorDatabase.from_yaml(database, beam=2, ignore_crystals=True)
        original_length = line.get_length()
        return (line, colldb, original_length), {}

    def install(line, colldb, original_length):
        colldb.install_everest_collimators(line=line)
        return line, colldb, original_length

    line, colldb, original_length = benchmark.pedantic(
        install, setup=setup, rounds=ROUNDS, warmup_rounds=1,
    )
    installed = len(line.xcoll.collimators)
    assert installed == len(colldb.collimator_names)
    assert installed >= 10
    assert np.isclose(line.get_length(), original_length)
    assert all(not isinstance(line[name], xc.EverestCrystal)
               for name in colldb.collimator_names)
    benchmark.extra_info.update({
        'protocol': 2, 'phase': 'installation', 'beam': 2, 'crystals': False,
        'physics_case': 'install-everest-ir7', 'installed_collimators': installed,
        'sequence_sha256': sequence_hash,
        'collimator_database_sha256': database_hash,
        'timed': 'colldb.install_everest_collimators(line=line)',
    })


@pytest.mark.quick
@pytest.mark.benchmark(group='LHC Run 3 | line tools')
def bench_lhc_collimator_family_and_gap_lookup(benchmark):
    """Read collimator family membership and configured gaps on a real line."""
    sequence = DATA / 'sequence_lhc_run3_b1.json'
    database = DATA / 'colldb_lhc_run3_ir7.yaml'
    _require_inputs(sequence, database)
    env = xt.load(sequence)
    line = env['lhcb1']
    colldb = xc.CollimatorDatabase.from_yaml(database, beam=1, ignore_crystals=True)
    colldb.install_everest_collimators(line=line)
    expected_names = set(colldb.collimator_names)

    def read_configuration():
        access = line.xcoll.collimators
        return access.families, access.gap

    families, gaps = benchmark.pedantic(
        read_configuration, rounds=7, iterations=100, warmup_rounds=1,
    )
    assert isinstance(families, dict) and families
    assert isinstance(gaps, dict)
    assert set(name for names in families.values() for name in names) == expected_names
    assert set(gaps) == expected_names
    benchmark.extra_info.update({
        'protocol': 2, 'phase': 'line-tools', 'beam': 1, 'crystals': False,
        'physics_case': 'family-gap-snapshot',
        'num_collimators': len(expected_names),
        'num_families': len(families),
        'sequence_sha256': hashlib.sha256(sequence.read_bytes()).hexdigest(),
        'collimator_database_sha256': hashlib.sha256(database.read_bytes()).hexdigest(),
        'timed': 'line.xcoll.collimators.families and .gap (per iteration)',
    })
