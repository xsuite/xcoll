"""Real LHC Run 3, IR7 Everest benchmarks using Xcoll's test datasets.

No scenario module is imported. All physics setup is explicit below.
Only tracking and/or LossMap are timed; sequence loading, installation,
Twiss, pencil generation, compilation, and per-round particle copying are not.
"""

import hashlib
from pathlib import Path

import numpy as np
import pytest
import xobjects as xo
import xtrack as xt
import xcoll as xc

DATA = Path(__file__).resolve().parents[1] / "tests" / "data"
SEED_COORDS = 20261109
SEED_PARTICLES = 24681357


def _sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _build_lhc_case(n_particles, n_turns, include_crystals):
    """Prepare an actual LHC simulation and pretracked particles once.

    Reuses release/v0.13.0 tests/data inputs; NOT an artificial toy ring.
    This helper stays in this benchmark file for readability.
    """
    beam = 1
    sequence = DATA / "sequence_lhc_run3_b1.json"
    database = DATA / "colldb_lhc_run3_ir7.yaml"
    if not sequence.exists() or not database.exists():
        pytest.skip(f"LHC benchmark data not found in {DATA}")

    env = xt.load(sequence)
    line = env["lhcb1"]
    colldb = xc.CollimatorDatabase.from_yaml(
        database, beam=beam, ignore_crystals=not include_crystals)
    colldb.install_everest_collimators(line=line)
    aperture = line.check_aperture()
    assert not np.any(aperture.has_aperture_problem)

    line.build_tracker(_context=xo.ContextCpu())
    line.xcoll.collimators.assign_optics()
    if include_crystals:
        line.xcoll.collimators.align_to_beam_divergence()

    tcp = "tcp.c6l7.b1"
    # generate_pencil itself uses NumPy's global RNG.
    previous_state = np.random.get_state()
    try:
        np.random.seed(SEED_COORDS)
        p0 = line[tcp].generate_pencil(n_particles)
    finally:
        np.random.set_state(previous_state)
    p0._init_random_number_generator(
        seeds=SEED_PARTICLES + np.arange(n_particles, dtype=np.uint32))

    line.xcoll.scattering.enable()
    pretracked = p0.copy()
    line.track(pretracked, num_turns=n_turns)  # warm-up / reference only
    n_lost = int(np.count_nonzero(np.asarray(pretracked.state) <= 0))
    assert n_lost > 0, "No beam losses: this is not a useful loss-map workload"

    # LossMap reads/reset element-specific energy-deposition counters.
    # Preserve values measured during the reference tracking pass.
    counters = {
        name: (float(el._acc_ionisation_loss),
               float(el._acc_ionisation_loss_sec))
        for name, el in zip(line.element_names, line.elements)
        if hasattr(el, "_acc_ionisation_loss")
    }

    def restore(reference):
        for name, (primary, secondary) in reference.items():
            line[name]._acc_ionisation_loss = primary
            line[name]._acc_ionisation_loss_sec = secondary

    # Prewarm the LossLocationRefinement code, which can be costly. It may
    # change loss locations in the particles, hence the copy. It also resets
    # counters, hence the restore immediately afterwards.
    reference = xc.LossMap(line=line, part=pretracked.copy(),
                           line_is_reversed=False, interpolation=0.2,
                           verbose=False)
    restore(counters)
    assert reference.num_losses > 0

    hashes = {"sequence_sha256": _sha256(sequence),
              "collimator_database_sha256": _sha256(database)}
    return line, p0, pretracked, counters, restore, reference, hashes


@pytest.fixture(scope="module")
def lhc_quick():
    return _build_lhc_case(n_particles=3_000, n_turns=3,
                           include_crystals=False)


@pytest.fixture(scope="module")
def lhc_full():
    return _build_lhc_case(n_particles=12_000, n_turns=10,
                           include_crystals=True)


def _report(benchmark, *, phase, n_particles, n_turns, crystals, ref, hashes):
    benchmark.extra_info.update({
        "protocol": 2,
        "phase": phase,
        "particles": n_particles,
        "turns": n_turns,
        "crystals": crystals,
        "beam": 1,
        "plane": "H",
        "interpolation_m": 0.2,
        "coordinate_seed": SEED_COORDS,
        "particle_seed": SEED_PARTICLES,
        "reference_collimator_losses": int(ref.num_collimator_losses),
        "reference_aperture_losses": int(ref.num_aperture_losses),
        **hashes,
    })


@pytest.mark.quick
@pytest.mark.benchmark(group="LHC Run 3 | tracking")
def bench_lhc_tracking_quick(benchmark, lhc_quick):
    """Real LHC IR7, 3k pencil protons, 3 turns; no loss-map processing."""
    line, p0, _, counters, restore, ref, hashes = lhc_quick
    zeros = {name: (0., 0.) for name in counters}

    def setup():
        restore(zeros)
        return (p0.copy(),), {}

    def track(p):
        line.track(p, num_turns=3)
        return p

    p = benchmark.pedantic(track, setup=setup, rounds=5, warmup_rounds=1)
    assert np.count_nonzero(np.asarray(p.state) <= 0) > 0
    _report(benchmark, phase="tracking", n_particles=3_000, n_turns=3,
            crystals=False, ref=ref, hashes=hashes)


@pytest.mark.full
@pytest.mark.benchmark(group="LHC Run 3 | tracking")
def bench_lhc_tracking_full(benchmark, lhc_full):
    """Real LHC IR7, 12k pencil protons, 10 turns, including bent crystals."""
    line, p0, _, counters, restore, ref, hashes = lhc_full
    zeros = {name: (0., 0.) for name in counters}

    def setup():
        restore(zeros)
        return (p0.copy(),), {}

    def track(p):
        line.track(p, num_turns=10)
        return p

    p = benchmark.pedantic(track, setup=setup, rounds=3, warmup_rounds=1)
    assert np.count_nonzero(np.asarray(p.state) <= 0) > 0
    _report(benchmark, phase="tracking", n_particles=12_000, n_turns=10,
            crystals=True, ref=ref, hashes=hashes)


@pytest.mark.full
@pytest.mark.benchmark(group="LHC Run 3 | analysis")
def bench_lhc_lossmap_full(benchmark, lhc_full):
    """Only LossMap creation, including aperture loss refinement (12k/10t)."""
    line, _, pretracked, counters, restore, ref, hashes = lhc_full

    def setup():
        restore(counters)
        return (pretracked.copy(),), {}

    def make_lossmap(p):
        return xc.LossMap(line=line, part=p, line_is_reversed=False,
                          interpolation=0.2, verbose=False)

    lm = benchmark.pedantic(make_lossmap, setup=setup,
                            rounds=5, warmup_rounds=1)
    assert lm.num_initial == 12_000
    assert lm.num_collimator_losses == ref.num_collimator_losses
    assert lm.num_aperture_losses == ref.num_aperture_losses
    _report(benchmark, phase="lossmap-only", n_particles=12_000, n_turns=10,
            crystals=True, ref=ref, hashes=hashes)


@pytest.mark.full
@pytest.mark.benchmark(group="LHC Run 3 | combined")
def bench_lhc_tracking_and_lossmap_full(benchmark, lhc_full):
    """Tracking + complete LossMap creation, excluding machine preparation."""
    line, p0, _, counters, restore, ref, hashes = lhc_full
    zeros = {name: (0., 0.) for name in counters}

    def setup():
        restore(zeros)
        return (p0.copy(),), {}

    def track_and_analyse(p):
        line.track(p, num_turns=10)
        return xc.LossMap(line=line, part=p, line_is_reversed=False,
                          interpolation=0.2, verbose=False)

    lm = benchmark.pedantic(track_and_analyse, setup=setup,
                            rounds=3, warmup_rounds=1)
    assert lm.num_initial == 12_000
    assert lm.num_losses > 0
    _report(benchmark, phase="tracking+lossmap", n_particles=12_000,
            n_turns=10, crystals=True, ref=ref, hashes=hashes)
