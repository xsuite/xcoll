"""FLUKA coupling tracking. Requires the separately prepared FLUKA environment."""

import numpy as np
import pytest
import xpart as xp
import xtrack as xt
import xcoll as xc

N = 500
SEED = 20261110


@pytest.fixture(scope="module")
def fluka_case():
    if not xc.fluka.interface.compiled:
        pytest.skip("FLUKA interface not compiled; run the FLUKA setup workflow")
    engine = xc.fluka.engine
    if engine.is_running():
        engine.stop(clean=True)
    coll = xc.FlukaCollimator(length=0.6, jaw=0.001, assembly="hilumi_tcppm")
    engine.particle_ref = xt.Particles.reference_from_pdg_id(
        pdg_id="proton", p0c=6.8e12,
    )
    rng = np.random.default_rng(SEED)
    particles = xp.build_particles(
        x=rng.normal(0.002, 1e-3, N),
        px=rng.normal(0, 5e-6, N),
        y=rng.normal(0, 1e-3, N),
        py=rng.normal(0, 5e-6, N),
        particle_ref=engine.particle_ref, _capacity=N * 4,
    )
    try:
        engine.start(elements=coll, seed=1993, verbose=False)
        yield coll, particles
    finally:
        engine.stop(clean=True)


@pytest.mark.fluka
@pytest.mark.benchmark(group="FLUKA | tracking")
def bench_fluka_collimator_500(benchmark, fluka_case):
    coll, initial = fluka_case
    # Engine start and input generation are outside the timed region.
    coll.track(initial.copy())  # warm up Python/C interfaces

    def setup():
        return (initial.copy(),), {}

    def track(particles):
        coll.track(particles)
        return particles

    result = benchmark.pedantic(track, setup=setup, rounds=3, warmup_rounds=0)
    lost = int(np.count_nonzero(np.asarray(result.state) <= 0))
    assert lost > 0, "FLUKA produced no material losses"
    benchmark.extra_info.update(
        protocol=2, physics_case="fluka-hilumi-tcppm",
        engine="fluka", particles=N, particle_seed=SEED,
        lost_fraction=lost / N,
    )
