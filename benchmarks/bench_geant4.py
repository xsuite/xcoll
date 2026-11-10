"""BDSIM/Geant4 benchmark matrix; Geant4Crystal is not implemented in Xcoll.

The supported Geant4CollimatorTip case exercises layered material geometry.
Engine start-up is amortised across the complete module (not timed).
"""
import numpy as np
import pytest
import xpart as xp
import xtrack as xt
import xcoll as xc

N = 500
SEED = 20261110
ROUNDS = 3

CASES = ['titanium-hit', 'titanium-miss', 'titanium-grazing',
         'tungsten-hit', 'beryllium-hit', 'copper-hit', 'tungsten-tip-hit']


def _make_case(name):
    material = ({'titanium': xc.materials.Titanium, 'tungsten': xc.materials.Tungsten,
                 'beryllium': xc.materials.Beryllium, 'copper': xc.materials.Copper}[name.split('-')[0]])
    if name == 'tungsten-tip-hit':
        return xc.Geant4CollimatorTip(length=0.6, jaw=0.001,
                                     material=xc.materials.Copper,
                                     tip_material=xc.materials.Tungsten,
                                     tip_thickness=0.005)
    return xc.Geant4Collimator(length=0.6, jaw=0.001, material=material)


def _particles(name, reference):
    rng = np.random.default_rng(SEED)
    if name.endswith('miss'):
        x = rng.normal(0.0, 5e-5, N)
    elif name.endswith('grazing'):
        x = rng.normal(0.001, 2e-6, N)
    else:
        x = rng.normal(0.0011, 1e-5, N)
    return xp.build_particles(x=x, px=rng.normal(0, 3e-6, N),
                              y=rng.normal(0, 1e-4, N), py=rng.normal(0, 3e-6, N),
                              particle_ref=reference, _capacity=N*4)


@pytest.fixture(scope='module')
def geant4_cases():
    if not xc.geant4.interface.ready:
        pytest.fail('Geant4/BDSIM interface is not ready; run its setup workflow')
    engine = xc.geant4.engine
    if engine.is_running():
        engine.stop(clean=True)
    elements = {name: _make_case(name) for name in CASES}
    engine.particle_ref = xt.Particles('proton', p0c=6.8e12)
    initial = {name: _particles(name, engine.particle_ref) for name in CASES}
    try:
        engine.start(elements=list(elements.values()), seed=1993)
        yield elements, initial
    finally:
        if engine.is_running():
            engine.stop(clean=True)


@pytest.mark.geant4
@pytest.mark.benchmark(group='Geant4 | tracking')
@pytest.mark.parametrize('case', CASES, ids=CASES)
def bench_geant4_tracking(benchmark, geant4_cases, case):
    elements, initial = geant4_cases
    coll = elements[case]
    p0 = initial[case]
    coll.track(p0.copy())

    def setup():
        return (p0.copy(),), {}

    def track(p):
        coll.track(p)
        return p

    result = benchmark.pedantic(track, setup=setup, rounds=ROUNDS, warmup_rounds=0)
    lost = int(np.count_nonzero(np.asarray(result.state) <= 0))
    if case.endswith('miss'):
        assert lost == 0
    benchmark.extra_info.update(
        protocol=3, particles=N, engine='geant4', physics_case=case,
        material=case.split('-')[0], geometry=('tip' if 'tip' in case else 'collimator'),
        coordinate_seed=SEED, particle_seed=1993, lost_fraction=lost/N,
    )
