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


# Shared workloads also used by the standalone crossover and physics scripts.
try:
    from .everest_workloads import (CRYSTAL_CASES, make_element,
        make_particles, reset_counters, host, COORD_SEED, PARTICLE_SEED)
except ImportError:
    from everest_workloads import (CRYSTAL_CASES, make_element,
        make_particles, reset_counters, host, COORD_SEED, PARTICLE_SEED)


@pytest.mark.gpu
@pytest.mark.benchmark(group='Everest GPU | block and crystal')
@pytest.mark.parametrize('case', ['block-cfc', *('crystal-'+c for c in CRYSTAL_CASES)])
@pytest.mark.parametrize('npart', [1_000, 20_000, 200_000])
def bench_gpu_block_crystal(benchmark, gpu, case, npart):
    cp, context = gpu
    element = make_element(case, context)
    initial = make_particles(case, npart, context)
    initial_px = host(initial, 'px').copy()

    def setup():
        reset_counters([element])
        p = initial.copy()
        cp.cuda.runtime.deviceSynchronize()
        return (p,), {}

    def track(p):
        element.track(p)
        cp.cuda.runtime.deviceSynchronize()
        return p

    track(*setup()[0])
    result = benchmark.pedantic(track, setup=setup, rounds=ROUNDS, warmup_rounds=1)
    state = host(result, 'state')
    ids = host(result, 'particle_id')
    alive = (state > 0) & (ids >= 0) & (ids < npart)
    kick = host(result, 'px')[alive] - initial_px[ids[alive]]
    lost = int(np.count_nonzero(state <= 0))
    if case.endswith('miss'):
        assert lost == 0
    gpu_name = cp.cuda.runtime.getDeviceProperties(cp.cuda.Device().id)['name']
    if isinstance(gpu_name, bytes):
        gpu_name = gpu_name.decode(errors='replace')
    benchmark.extra_info.update(protocol=3, particles=npart, physics_case='gpu-'+case,
        context='ContextCupy', gpu_name=gpu_name, coordinate_seed=COORD_SEED,
        particle_seed=PARTICLE_SEED, lost_fraction=lost/npart,
        timed='element.track(particles) + CUDA device synchronisation')
    if case=='crystal-channeling':
        benchmark.extra_info['channeling_proxy_fraction_of_survivors'] = (
            float(np.mean(kick>74.5e-6)) if kick.size else None)
    # Match the existing configured observable key exactly.
    benchmark.extra_info['survivors_with_kick_gt_0.1urad'] = (
        float(np.mean(np.abs(kick)>0.1e-6)) if kick.size else None)


@pytest.mark.gpu
@pytest.mark.benchmark(group='InteractionRecord GPU | data extraction')
@pytest.mark.parametrize('format', ['full-jaw', 'full-lattice', 'selected-jaw'])
def bench_gpu_interaction_record_to_pandas(benchmark, gpu, format):
    cp, context = gpu
    n = 50_000
    coll = xc.TransparentCollimator(name='TCP', length=0.6,
                                    jaw=[0.001,-0.001], _context=context)
    rng = np.random.default_rng(SEED)
    p = xp.Particles(p0c=450e9, x=rng.normal(1.025e-3,2e-6,n),
        px=rng.normal(0,2e-6,n), y=np.zeros(n), py=np.zeros(n), _context=context)
    selected = format=='selected-jaw'
    frame = 'lattice' if format=='full-lattice' else 'jaw'
    rec = xc.InteractionRecord(elements=[coll], names=['TCP'], num_rows=3*n+100,
        columns=['particle_id_before','x_before'] if selected else None,
        record_impacts=True, record_exits=True, record_scatterings=False)
    try:
        coll.track(p)
        cp.cuda.runtime.deviceSynchronize()
        recorded = int(rec._index.num_recorded)
        assert recorded == 2*n < rec.num_rows
        result = benchmark.pedantic(rec.to_pandas, kwargs={'frame':frame},
                                    rounds=ROUNDS, warmup_rounds=1)
        assert len(result)==recorded
        assert 'interaction_type' in result.columns
        name = cp.cuda.runtime.getDeviceProperties(cp.cuda.Device().id)['name']
        if isinstance(name, bytes):
            name = name.decode(errors='replace')
        benchmark.extra_info.update(protocol=3, particles=n,
            physics_case='gpu-transparent-jaw-record-export', context='ContextCupy',
            gpu_name=name, record_mode='impacts-exits',
            record_columns='selected' if selected else 'all', record_frame=frame,
            record_capacity=rec.num_rows, recorded_interactions=recorded,
            output_columns=len(result.columns), coordinate_seed=SEED,
            timed='record.to_pandas(frame=...) including device-to-host transfers')
    finally:
        rec.stop(elements=[coll])
