"""BDSIM/Geant4 coupling tracking on CPU, excluding engine start-up."""

import numpy as np
import pytest
import xpart as xp
import xtrack as xt
import xcoll as xc

N = 500
SEED = 20261110


@pytest.fixture(scope="module")
def geant4_case():
    if not xc.geant4.interface.compiled:
        pytest.skip("BDSIM/Geant4 interface not compiled; run setup workflow")
    engine = xc.geant4.engine
    if engine.is_running():
        engine.stop(clean=True)
    coll = xc.Geant4Collimator(length=0.6, jaw=0.001, material="Ti")
    engine.particle_ref = xt.Particles("proton", p0c=6.8e12)
    rng = np.random.default_rng(SEED)
    particles = xp.build_particles(
        x=rng.normal(0.0011, 1e-4, N),
        px=rng.normal(0, 5e-6, N),
        y=rng.normal(0, 1e-3, N),
        py=rng.normal(0, 5e-6, N),
        particle_ref=engine.particle_ref, _capacity=N * 4,
    )
    try:
        engine.start(elements=coll, seed=1993)
        yield coll, particles
    finally:
        engine.stop(clean=True)


@pytest.mark.geant4
@pytest.mark.benchmark(group="Geant4 | tracking")
def bench_geant4_collimator_500(benchmark, geant4_case):
    coll, initial = geant4_case
    coll.track(initial.copy())  # warm up; engine initialisation excluded

    def setup():
        return (initial.copy(),), {}

    def track(particles):
        coll.track(particles)
        return particles

    result = benchmark.pedantic(track, setup=setup, rounds=3, warmup_rounds=0)
    lost = int(np.count_nonzero(np.asarray(result.state) <= 0))
    assert lost > 0, "Geant4 produced no material losses"
    benchmark.extra_info.update(
        protocol=2, physics_case="geant4-ti-jaw",
        engine="geant4", particles=N, particle_seed=SEED,
        lost_fraction=lost / N,
    )
