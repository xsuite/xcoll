"""FLUKA benchmark matrix: hit/miss/grazing, materials, and a Si crystal.

All input geometries and FLUKA engine initialisation are outside timed code.
The optional named FEDB assembly is omitted when unavailable; generic cases
remain mandatory once the FLUKA environment is configured.
"""
from pathlib import Path
import numpy as np
import pytest
import xpart as xp
import xtrack as xt
import xcoll as xc

try:
    from .perf.particle_metrics import count_lost_primaries
except ImportError:
    from perf.particle_metrics import count_lost_primaries

N = 500
SEED = 20261110
ROUNDS = 3

CASES = [
    'cfc-hit', 'cfc-miss', 'cfc-grazing',
    'tungsten-hit', 'copper-hit', 'fedb-hilumi-hit',
    'silicon-crystal-hit', 'silicon-crystal-miss',
]
CRYSTAL_CASES = tuple(case for case in CASES if case.startswith('silicon-crystal'))
COLLIMATOR_CASES = tuple(case for case in CASES if case not in CRYSTAL_CASES)


def _make_case(name):
    if name in CRYSTAL_CASES:
        return xc.FlukaCrystal(length=0.004, material=xc.materials.Silicon,
                               side='+', jaw=0.001, bending_angle=149e-6)
    if name == 'fedb-hilumi-hit':
        return xc.FlukaCollimator(length=0.6, jaw=0.001, assembly='hilumi_tcppm')
    material = (xc.materials.Tungsten if name.startswith('tungsten') else
                xc.materials.Copper if name.startswith('copper') else
                xc.materials.CarbonFibreComposite)
    return xc.FlukaCollimator(length=0.6, jaw=0.001, material=material)


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


def _crystal_missing_in_log(folder):
    """Recognise the specific FLUKA CRYINI configuration error, not every crash.

    Engine.start() may already have called stop() after failure, so inspect
    the explicitly supplied working directory rather than engine.cwd.
    """
    for log in Path(folder).rglob('*'):
        if not log.is_file() or log.suffix.lower() not in {'.log', '.out', '.err'}:
            continue
        if log.stat().st_size > 4_000_000:
            continue
        content = log.read_text(encoding='utf-8', errors='replace').lower()
        if 'crystal not found' in content and 'cryini' in content:
            return True
    return False


@pytest.fixture(scope='module')
def fluka_cases(tmp_path_factory):
    """Lazily start FLUKA once per geometry class, never mix crystals with jaws."""
    try:
        ready = xc.fluka.interface.ready
    except (ImportError, OSError, RuntimeError, ValueError) as exc:
        pytest.skip(f'FLUKA environment unavailable: {exc}')
    if not ready:
        pytest.skip('FLUKA environment not ready (FLUKA, FlukaIO, LineBuilder, or compiled interface missing)')
    engine = xc.fluka.engine
    if engine.is_running():
        engine.stop(clean=True)

    state = {'group': None, 'elements': {}, 'initial': {}, 'crystal_skip': None}
    available_assemblies = set(xc.fluka.assemblies)

    def get_case(case):
        if case == 'fedb-hilumi-hit' and 'hilumi_tcppm' not in available_assemblies:
            pytest.skip('Optional FLUKA FEDB assembly hilumi_tcppm is unavailable')

        group = 'crystal' if case in CRYSTAL_CASES else 'collimator'
        if group == 'crystal' and state['crystal_skip']:
            pytest.skip(state['crystal_skip'])

        if state['group'] != group:
            if engine.is_running():
                engine.stop(clean=True)
            state.update(group=None, elements={}, initial={})
            names = (list(CRYSTAL_CASES) if group == 'crystal' else
                     [name for name in COLLIMATOR_CASES
                      if name != 'fedb-hilumi-hit' or 'hilumi_tcppm' in available_assemblies])
            elements = {name: _make_case(name) for name in names}
            reference = xt.Particles.reference_from_pdg_id('proton', p0c=6.8e12)
            engine.particle_ref = reference
            initial = {name: _particles(name, reference) for name in names}
            # Xcoll appends _0000 to an existing cwd, so pass a new child
            # directory and keep the actual server logs at this exact path.
            run_dir = tmp_path_factory.mktemp(f'fluka-bench-{group}') / 'run'
            try:
                # Keep input/output on failure so we can diagnose CRYINI.
                engine.start(elements=list(elements.values()), seed=1993,
                             verbose=False, cwd=run_dir, clean=False)
            except Exception as exc:
                if group == 'crystal' and _crystal_missing_in_log(run_dir):
                    reason = (f'FLUKA crystal initialisation unavailable: '
                              f'CRYINI / Crystal not found. See {run_dir}')
                    state['crystal_skip'] = reason
                    pytest.skip(reason)
                raise RuntimeError(f'FLUKA {group} engine failed during startup. '
                                   f'Preserved output in {run_dir}') from exc
            state.update(group=group, elements=elements, initial=initial)

        return state['elements'][case], state['initial'][case]
    try:
        yield get_case
    finally:
        if engine.is_running():
            engine.stop(clean=True)


@pytest.mark.fluka
@pytest.mark.benchmark(group='FLUKA | tracking')
@pytest.mark.parametrize('case', CASES, ids=CASES)
def bench_fluka_tracking(benchmark, fluka_cases, case):
    coll, p0 = fluka_cases(case)
    coll.track(p0.copy())  # warm-up outside timing

    def setup():
        return (p0.copy(),), {}

    def track(p):
        coll.track(p)
        return p

    result = benchmark.pedantic(track, setup=setup, rounds=ROUNDS, warmup_rounds=0)
    lost = count_lost_primaries(result, N)
    if case.endswith('miss'):
        assert lost == 0
    benchmark.extra_info.update(
        protocol=4, particles=N, engine='fluka', physics_case=case,
        material=('silicon' if 'crystal' in case else case.split('-')[0]),
        geometry=('crystal' if 'crystal' in case else 'collimator'),
        coordinate_seed=SEED, particle_seed=1993,
        lost_fraction=lost/N,
    )
