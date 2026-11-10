"""CPU emittance monitor accumulation and subsequent derived statistics."""

import numpy as np
import pytest
import xobjects as xo
import xpart as xp
import xcoll as xc

N = 100_000
SEED = 20261110


def _particles(context):
    rng = np.random.default_rng(SEED)
    return xp.Particles(
        p0c=450e9, x=rng.normal(0, 1e-3, N),
        px=rng.normal(0, 3e-6, N),
        y=rng.normal(0, 8e-4, N), py=rng.normal(0, 4e-6, N),
        zeta=rng.normal(0, .1, N), delta=rng.normal(0, 1e-4, N),
        _context=context,
    )


def _monitor(context, planes):
    return xc.EmittanceMonitor(
        start_at_turn=0, stop_at_turn=1,
        horizontal=True,
        vertical=(planes == "6D"),
        longitudinal=(planes == "6D"),
        _context=context,
    )


@pytest.mark.quick
@pytest.mark.benchmark(group="EmittanceMonitor | CPU accumulation")
@pytest.mark.parametrize("planes", ["H", "6D"])
def bench_emittance_monitor_accumulate(benchmark, planes):
    context = xo.ContextCpu()
    initial = _particles(context)
    monitor = _monitor(context, planes)
    monitor.track(initial.copy())  # compile outside timing
    monitor.reset()

    def setup():
        monitor.reset()
        return (initial.copy(),), {}

    def track(particles):
        monitor.track(particles)
        return particles

    benchmark.pedantic(track, setup=setup, rounds=5, warmup_rounds=1)
    assert int(np.asarray(monitor.data.count)[0]) == N
    benchmark.extra_info.update(
        protocol=2, particles=N, physics_case="emittance-monitor-accumulate",
        planes=planes, recorded_particles=int(np.asarray(monitor.data.count)[0]),
    )


@pytest.mark.full
@pytest.mark.benchmark(group="EmittanceMonitor | analysis")
def bench_emittance_monitor_calculate_horizontal(benchmark):
    context = xo.ContextCpu()
    particles = _particles(context)
    monitor = _monitor(context, "H")
    monitor.track(particles)  # data accumulation outside measurement

    def setup():
        # Reset cached flags, not the accumulated moments. Emittance calculation
        # should perform a fresh reduction each round, not return cached results.
        monitor.data.cached[0] = 0
        return (), {}

    result = benchmark.pedantic(lambda: monitor.gemitt_x, setup=setup,
                                rounds=10, warmup_rounds=1)
    assert np.isfinite(np.asarray(result)).all()
    benchmark.extra_info.update(
        protocol=2, particles=N, physics_case="emittance-monitor-reduce",
        planes="H", geometric_emittance=float(np.asarray(result).flat[0]),
    )
