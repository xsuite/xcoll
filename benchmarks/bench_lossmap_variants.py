"""LHC loss-map analysis without aperture interpolation.

This complements bench_lhc_lossmap.py, where interpolation=0.2 m. Tracking,
model preparation, and particle copies are deliberately excluded from timing.
The physical input population is identical in every measured round.
"""

import hashlib
from pathlib import Path

import numpy as np
import pytest
import xobjects as xo
import xtrack as xt
import xcoll as xc

DATA = Path(__file__).resolve().parents[1] / 'tests' / 'data'


@pytest.mark.full
@pytest.mark.benchmark(group='LHC Run 3 | analysis')
def bench_lhc_lossmap_without_interpolation(benchmark):
    """LossMap construction for 3k B1 particles/3 turns, interpolation disabled."""
    sequence = DATA / 'sequence_lhc_run3_b1.json'
    database = DATA / 'colldb_lhc_run3_ir7.yaml'
    if not sequence.is_file() or not database.is_file():
        pytest.skip('LHC Run 3 datasets not available')

    env = xt.load(sequence)
    line = env['lhcb1']
    colldb = xc.CollimatorDatabase.from_yaml(database, beam=1, ignore_crystals=True)
    colldb.install_everest_collimators(line=line)
    assert not np.any(line.check_aperture().has_aperture_problem)
    line.build_tracker(_context=xo.ContextCpu())
    line.xcoll.collimators.assign_optics()

    # generate_pencil uses NumPy's global random generator.
    old_state = np.random.get_state()
    try:
        np.random.seed(20261109)
        particles = line['tcp.c6l7.b1'].generate_pencil(3_000)
    finally:
        np.random.set_state(old_state)
    particles._init_random_number_generator(
        seeds=24681357 + np.arange(3_000, dtype=np.uint32))
    line.xcoll.scattering.enable()
    line.track(particles, num_turns=3)
    line.xcoll.scattering.disable()
    assert np.count_nonzero(np.asarray(particles.state) <= 0) > 0

    # LossMap consumes/reset energy-deposition counters and may change loss
    # locations. Preserve both the counters and the untouched particle state.
    reference_counters = {
        name: (float(element._acc_ionisation_loss),
               float(element._acc_ionisation_loss_sec))
        for name, element in zip(line.element_names, line.elements)
        if hasattr(element, '_acc_ionisation_loss')
    }

    def restore_counters():
        for name, (primary, secondary) in reference_counters.items():
            element = line[name]
            element._acc_ionisation_loss = primary
            element._acc_ionisation_loss_sec = secondary

    def make_lossmap(p):
        return xc.LossMap(line=line, part=p, line_is_reversed=False,
                          interpolation=False, verbose=False)

    restore_counters()
    reference = make_lossmap(particles.copy())  # Warm-up outside timing
    restore_counters()

    def setup():
        restore_counters()
        return (particles.copy(),), {}

    result = benchmark.pedantic(
        make_lossmap, setup=setup, rounds=5, warmup_rounds=1,
    )
    assert result.num_initial == 3_000
    assert result.num_collimator_losses == reference.num_collimator_losses
    assert result.num_aperture_losses == reference.num_aperture_losses
    benchmark.extra_info.update({
        'protocol': 2, 'phase': 'lossmap-only', 'beam': 1, 'plane': 'H',
        'particles': 3000, 'turns': 3, 'crystals': False,
        'interpolation_m': False,
        'coordinate_seed': 20261109, 'particle_seed': 24681357,
        'physics_case': 'lossmap-without-interpolation',
        'reference_collimator_losses': int(reference.num_collimator_losses),
        'reference_aperture_losses': int(reference.num_aperture_losses),
        'sequence_sha256': hashlib.sha256(sequence.read_bytes()).hexdigest(),
        'collimator_database_sha256': hashlib.sha256(database.read_bytes()).hexdigest(),
        'timed': 'xc.LossMap(..., interpolation=False)',
    })
