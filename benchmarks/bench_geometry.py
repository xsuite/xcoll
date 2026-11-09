"""CPU geometry baselines.

These are material-free BlackAbsorber/BlackCrystal workloads. They benchmark
current jaw/tilt/crystal geometry + Xtrack dispatch, NOT NewSegments internals.
NewSegments must get its own benchmarks after its API stabilises.
"""

import numpy as np
import pytest
import xobjects as xo
import xpart as xp
import xcoll as xc

N = 1_000_000
ROUNDS = 5
SEED = 20261109


@pytest.mark.quick
@pytest.mark.benchmark(group="geometry CPU | absorber")
def bench_geometry_absorber_parallel_hit_1M(benchmark):
    """Straight two-jaw absorber, 1 million particles in left material at entrance."""
    context = xo.ContextCpu()
    absorber = xc.BlackAbsorber(length=0.6, jaw=[0.001, -0.001],
                                _context=context)
    rng = np.random.default_rng(SEED)
    p0 = xp.Particles(p0c=450e9, x=rng.normal(1.1e-3, 4e-6, N),
                      px=rng.normal(0, 3e-6, N), y=np.zeros(N),
                      py=np.zeros(N), _context=context)
    assert np.all(np.asarray(p0.x) > absorber.jaw_LU)
    absorber.track(p0.copy())  # kernel compilation outside timing

    def setup():
        return (p0.copy(),), {}

    def track(p):
        absorber.track(p)
        return p

    p = benchmark.pedantic(track, setup=setup, rounds=ROUNDS, warmup_rounds=1)
    lost = int(np.count_nonzero(np.asarray(p.state) <= 0))
    assert lost == N
    benchmark.extra_info.update(protocol=2, geometry="parallel-hit",
                                particles=N, lost_fraction=lost / N)


@pytest.mark.full
@pytest.mark.benchmark(group="geometry CPU | absorber")
def bench_geometry_absorber_parallel_miss_1M(benchmark):
    """Straight absorber, 1 million particles well inside the unobstructed opening."""
    context = xo.ContextCpu()
    absorber = xc.BlackAbsorber(length=0.6, jaw=[0.001, -0.001],
                                _context=context)
    rng = np.random.default_rng(SEED)
    p0 = xp.Particles(p0c=450e9, x=rng.normal(0, 0.1e-3, N),
                      px=rng.normal(0, 3e-6, N), y=np.zeros(N),
                      py=np.zeros(N), _context=context)
    assert np.max(np.abs(np.asarray(p0.x))) < 0.001
    absorber.track(p0.copy())

    def setup():
        return (p0.copy(),), {}

    def track(p):
        absorber.track(p)
        return p

    p = benchmark.pedantic(track, setup=setup, rounds=ROUNDS, warmup_rounds=1)
    lost = int(np.count_nonzero(np.asarray(p.state) <= 0))
    assert lost == 0
    benchmark.extra_info.update(protocol=2, geometry="parallel-miss",
                                particles=N, lost_fraction=0.0)


@pytest.mark.full
@pytest.mark.benchmark(group="geometry CPU | absorber")
def bench_geometry_absorber_tilted_1M(benchmark):
    """A tilted jaw, particles near the x=1 mm boundary with finite px."""
    context = xo.ContextCpu()
    absorber = xc.BlackAbsorber(length=0.6, jaw=[0.001, -0.001],
                                tilt=2e-5, _context=context)
    rng = np.random.default_rng(SEED)
    p0 = xp.Particles(p0c=450e9, x=rng.normal(0.001, 10e-6, N),
                      px=rng.normal(0, 20e-6, N), y=np.zeros(N),
                      py=np.zeros(N), _context=context)
    absorber.track(p0.copy())

    def setup():
        return (p0.copy(),), {}

    def track(p):
        absorber.track(p)
        return p

    p = benchmark.pedantic(track, setup=setup, rounds=ROUNDS, warmup_rounds=1)
    lost = int(np.count_nonzero(np.asarray(p.state) <= 0))
    assert 0 < lost < N
    benchmark.extra_info.update(protocol=2, geometry="tilted-two-jaw",
                                particles=N, lost_fraction=lost / N,
                                tilt_radians=0.002)


@pytest.mark.full
@pytest.mark.benchmark(group="geometry CPU | crystal")
def bench_geometry_black_crystal_curved_1M(benchmark):
    """Curved BlackCrystal geometry with a +1 mm Si-like crystal outline.

    BlackCrystal has no material scattering: this measures current crystal
    geometry and boundary handling independently of the Everest engine.
    """
    context = xo.ContextCpu()
    crystal = xc.BlackCrystal(length=0.002, bending_angle=149e-6,
                              width=0.002, height=0.05, side="+",
                              jaw=0.001, _context=context)
    rng = np.random.default_rng(SEED)
    p0 = xp.Particles(p0c=400e9, x=rng.normal(0.001, 0.4e-3, N),
                      px=rng.normal(0, 30e-6, N),
                      y=rng.normal(0, 1e-3, N),
                      py=rng.normal(0, 5e-6, N), _context=context)
    crystal.track(p0.copy())

    def setup():
        return (p0.copy(),), {}

    def track(p):
        crystal.track(p)
        return p

    p = benchmark.pedantic(track, setup=setup, rounds=ROUNDS, warmup_rounds=1)
    lost = int(np.count_nonzero(np.asarray(p.state) <= 0))
    assert 0 < lost < N
    benchmark.extra_info.update(protocol=2, geometry="curved-crystal",
                                particles=N, lost_fraction=lost / N)
