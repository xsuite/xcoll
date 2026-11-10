"""FLUKA benchmark matrix: hit/miss/grazing, materials, and a Si crystal.

All input geometries and FLUKA engine initialisation are outside timed code.
The optional named FEDB assembly is omitted when unavailable; generic cases
remain mandatory once the FLUKA environment is configured.
"""
import numpy as np
import pytest
import xpart as xp
import xtrack as xt
import xcoll as xc

N = 500
SEED = 20261110
ROUNDS = 3

CASES = [
    'cfc-hit', 'cfc-miss', 'cfc-grazing',
    'tungsten-hit', 'copper-hit', 'fedb-hilumi-hit',
    'silicon-crystal-hit', 'silicon-crystal-miss',
]


def _make_case(name):
    if name.startswith('silicon-crystal'):
        return xc.FlukaCrystal(length=0.004, material=xc.materials.Silicon,
                               side='+', jaw=0.001, bending_angle=149e-6)
    if name == 'fedb-hilumi-hit':
        return xc.FlukaCollimator(length=0.6, jaw=0.001, assembly='hilumi_tcppm')
    mat = (xc.materials.Tungsten if name.startswith('tungsten') else
           xc.materials.Copper if name.startswith('copper') else
           xc.materials.CarbonFibreComposite)
    return xc.FlukaCollimator(length=0.6, jaw=0.001, material=mat)


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
def fluka_cases():
    if not xc.fluka.interface.ready:
        pytest.fail('FLUKA interface is not ready; run the FLUKA setup workflow')
    engine = xc.fluka.engine
    if engine.is_running():
        engine.stop(clean=True)
    names = [s for s in CASES if not s.startswith('fedb-')]
    available = set(xc.fluka.assemblies)
    if 'hilumi_tcppm' in available:
        names.append('fedb-hilumi-hit')
    elements = {name: _make_case(name) for name in names}
    engine.particle_ref = xt.Particles.reference_from_pdg_id('proton', p0c=6.8e12)
    initial = {name: _particles(name, engine.particle_ref) for name in names}
    try:
        engine.start(elements=list(elements.values()), seed=1993, verbose=False)
        yield elements, initial
    finally:
        if engine.is_running():
            engine.stop(clean=True)


@pytest.mark.fluka
@pytest.mark.benchmark(group='FLUKA | tracking')
@pytest.mark.parametrize('case', CASES, ids=CASES)
def bench_fluka_tracking(benchmark, fluka_cases, case):
    elements, initial = fluka_cases
    if case not in elements:
        pytest.skip(f'Optional FLUKA FEDB assembly unavailable: {case}')
    coll = elements[case]
    p0 = initial[case]
    coll.track(p0.copy())  # warm-up outside timing

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
        protocol=3, particles=N, engine='fluka', physics_case=case,
        material=('silicon' if 'crystal' in case else case.split('-')[0]),
        geometry=('crystal' if 'crystal' in case else 'collimator'),
        coordinate_seed=SEED, particle_seed=1993,
        lost_fraction=lost/N,
    )
