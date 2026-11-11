"""CPU blow-up tracking: common kicks versus independent particle kicks."""

import numpy as np
import pytest
import xobjects as xo
import xpart as xp
import xcoll as xc

N = 100_000
SEED = 20261110


def _particles(context):
    rng = np.random.default_rng(SEED)
    particles = xp.Particles(
        p0c=450e9,
        x=rng.normal(0, 1e-3, N),
        px=rng.normal(0, 3e-6, N),
        y=rng.normal(0, 1e-3, N),
        py=rng.normal(0, 3e-6, N),
        _context=context,
    )
    particles._init_random_number_generator(
        seeds=SEED + np.arange(N, dtype=np.uint32)
    )
    return particles


@pytest.mark.quick
@pytest.mark.benchmark(group="BlowUp | CPU tracking")
@pytest.mark.parametrize("individual", [False, True], ids=["common-kick", "individual-kicks"])
def bench_blowup_cpu(benchmark, individual):
    context = xo.ContextCpu()
    particles0 = _particles(context)
    old_random_state = np.random.get_state()
    np.random.seed(SEED)
    try:
        element = xc.BlowUp(
            plane="H", amplitude=1e-12,
            start_at_turn=0, stop_at_turn=1,
            use_individual_kicks=individual, _context=context,
        )
    finally:
        np.random.set_state(old_random_state)
    element.activate()
    element.track(particles0.copy())  # compile outside timing

    def setup():
        return (particles0.copy(),), {}

    def track(particles):
        element.track(particles)
        return particles

    result = benchmark.pedantic(track, setup=setup, rounds=5, warmup_rounds=1)
    kick = np.asarray(result.px) - np.asarray(particles0.px)
    assert np.any(kick != 0)
    benchmark.extra_info.update(
        protocol=2, particles=N, particle_seed=SEED,
        physics_case="blowup-individual" if individual else "blowup-common",
        individual_kicks=individual,
        rms_kick=float(np.sqrt(np.mean(kick**2))),
    )
