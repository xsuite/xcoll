"""Opt-in CUDA benchmarks; run locally with ``pytest -m gpu``.

Timing includes a CUDA device synchronisation so asynchronous launches do
not produce misleadingly small runtime measurements. Input preparation,
particle copies, compilation, and initialisation are outside timing.
"""
import numpy as np
import pytest
import xobjects as xo
import xpart as xp
import xcoll as xc

SEED = 20261110
ROUNDS = 5


@pytest.fixture(scope='module')
def gpu():
    cp = pytest.importorskip('cupy', reason='GPU benchmarks require CuPy')
    try:
        if cp.cuda.runtime.getDeviceCount() < 1:
            pytest.skip('No CUDA device available')
        context = xo.ContextCupy()
    except Exception as error:
        pytest.skip(f'CUDA context is unavailable: {error}')
    return cp, context


def _particles(cp, context, case, npart):
    rng = np.random.default_rng(SEED)
    if case.endswith('miss'):
        x = rng.normal(0, 5e-5, npart)
    elif case.endswith('grazing'):
        x = rng.normal(0.001, 2e-6, npart)
    else:
        x = rng.normal(0.0011, 1e-5, npart)
    particles = xp.Particles(p0c=450e9, x=x, px=rng.normal(0, 3e-6, npart),
                             y=rng.normal(0, 1e-4, npart), py=rng.normal(0, 3e-6, npart),
                             _context=context)
    particles._init_random_number_generator(seeds=13579 + np.arange(npart, dtype=np.uint32))
    return particles


@pytest.mark.gpu
@pytest.mark.benchmark(group='Everest GPU | tracking')
@pytest.mark.parametrize('case', ['cfc-hit', 'cfc-miss', 'cfc-grazing', 'tungsten-hit'])
@pytest.mark.parametrize('npart', [20, 100, 500, 2000, 5000, 20_000, 50_000, 200_000, 1_000_000])
def bench_gpu_everest(benchmark, gpu, case, npart):
    cp, context = gpu
    material = (xc.materials.Tungsten if case.startswith('tungsten')
                else xc.materials.CarbonFibreComposite)
    element = xc.EverestCollimator(length=0.6, jaw=[0.001, -0.001],
                                  material=material, _context=context)
    initial = _particles(cp, context, case, npart)

    def track(p):
        element.track(p)
        cp.cuda.runtime.deviceSynchronize()
        return p

    track(initial.copy())  # compile and warm up outside timing

    def setup():
        return (initial.copy(),), {}

    result = benchmark.pedantic(track, setup=setup, rounds=ROUNDS, warmup_rounds=1)
    lost = int(cp.count_nonzero(result.state <= 0).get())
    if case.endswith('miss'):
        assert lost == 0
    properties = cp.cuda.runtime.getDeviceProperties(cp.cuda.Device().id)
    name = properties.get('name', b'unknown GPU')
    if isinstance(name, bytes):
        name = name.decode(errors='replace')
    benchmark.extra_info.update(protocol=3, particles=npart, physics_case='gpu-'+case,
                                context='ContextCupy', gpu_name=name,
                                coordinate_seed=SEED, particle_seed=13579,
                                lost_fraction=lost/npart,
                                timed='element.track(particles) + CUDA device synchronisation')
