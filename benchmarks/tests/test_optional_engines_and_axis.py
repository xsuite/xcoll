"""Structural and isolated engine fixture tests; no Xsuite/FLUKA/BDSIM needed."""
from __future__ import annotations

import importlib.util
import sys
import types
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
BENCH = ROOT / 'benchmarks'


def _import_benchmark(path, monkeypatch, *, ready=True, crash_crystal=False, crash_collimator=False):
    class FakeElement:
        def __init__(self, **kwargs):
            self.crystal = 'bending_angle' in kwargs
        def track(self, particles):
            return None

    class FakeEngine:
        def __init__(self):
            self.running = False
            self.groups = []
            self.particle_ref = None
        def is_running(self):
            return self.running
        def stop(self, clean=True):
            self.running = False
        def start(self, *, elements, seed, cwd=None, **kwargs):
            crystal = any(e.crystal for e in elements)
            self.groups.append('crystal' if crystal else 'collimator')
            if (crystal and crash_crystal) or (not crystal and crash_collimator):
                Path(cwd).mkdir(parents=True, exist_ok=True)
                Path(cwd, 'rfluka.log').write_text(
                    'Abort called from CRYINI reason Crystal not found Run stopped!'
                    if crystal else 'FLUKA failed for an unrelated reason')
                raise RuntimeError('The flukaserver died')
            self.running = True

    engine = FakeEngine()
    fake_xc = types.ModuleType('xcoll')
    fake_xc.fluka = types.SimpleNamespace(
        interface=types.SimpleNamespace(ready=ready),
        engine=engine, assemblies={'hilumi_tcppm'},
    )
    fake_xc.geant4 = types.SimpleNamespace(interface=types.SimpleNamespace(ready=ready), engine=engine)
    fake_xc.materials = types.SimpleNamespace(
        Silicon=object(), Tungsten=object(), Copper=object(),
        CarbonFibreComposite=object(), Titanium=object(), Beryllium=object(),
    )
    fake_xc.FlukaCrystal = FakeElement
    fake_xc.FlukaCollimator = FakeElement
    fake_xc.Geant4Collimator = FakeElement
    fake_xc.Geant4CollimatorTip = FakeElement
    fake_xt = types.ModuleType('xtrack')
    fake_xt.Particles = types.SimpleNamespace(
        reference_from_pdg_id=lambda *args, **kwargs: object()
    )
    fake_xp = types.ModuleType('xpart')
    fake_xp.build_particles = lambda **kwargs: object()
    for key, obj in [('xcoll', fake_xc), ('xtrack', fake_xt), ('xpart', fake_xp)]:
        monkeypatch.setitem(sys.modules, key, obj)
    spec = importlib.util.spec_from_file_location('external_bench_mock', path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod, engine


def test_fluka_crystal_failure_does_not_poison_collimator_cases(monkeypatch, tmp_path):
    mod, engine = _import_benchmark(BENCH / 'bench_fluka.py', monkeypatch,
                                    crash_crystal=True)
    counter = iter(range(8))
    class PathFactory:
        def mktemp(self, name):
            folder = tmp_path / f'{name}-{next(counter)}'
            folder.mkdir()
            return folder
    fixture = mod.fluka_cases.__wrapped__(PathFactory())
    get_case = next(fixture)
    assert get_case('cfc-hit')[0].crystal is False
    assert get_case('tungsten-hit')[0].crystal is False
    assert engine.groups == ['collimator']
    with pytest.raises(pytest.skip.Exception, match='CRYINI'):
        get_case('silicon-crystal-hit')
    with pytest.raises(pytest.skip.Exception, match='CRYINI'):
        get_case('silicon-crystal-miss')
    assert engine.groups == ['collimator', 'crystal']
    # Collimator group can be restarted after crystal failure.
    assert get_case('cfc-miss')[0].crystal is False
    assert engine.groups == ['collimator', 'crystal', 'collimator']
    with pytest.raises(StopIteration):
        next(fixture)


@pytest.mark.parametrize('name,fixture', [('fluka', 'fluka_cases'),
                                          ('geant4', 'geant4_cases')])
def test_missing_engine_skips(monkeypatch, tmp_path, name, fixture):
    mod, engine = _import_benchmark(BENCH / f'bench_{name}.py', monkeypatch, ready=False)
    if name == 'fluka':
        args = (tmp_path,)
    else:
        args = ()
    g = getattr(mod, fixture).__wrapped__(*args)
    with pytest.raises(pytest.skip.Exception, match='not ready'):
        next(g)
    assert not engine.groups


def test_unrelated_fluka_crash_is_not_silently_skipped(monkeypatch, tmp_path):
    mod, engine = _import_benchmark(BENCH / 'bench_fluka.py', monkeypatch,
                                    crash_collimator=True)
    class PathFactory:
        def mktemp(self, name):
            folder = tmp_path / name
            folder.mkdir()
            return folder
    fixture = mod.fluka_cases.__wrapped__(PathFactory())
    get_case = next(fixture)
    with pytest.raises(RuntimeError, match='collimator engine failed'):
        get_case('cfc-hit')
    assert engine.groups == ['collimator']
    with pytest.raises(StopIteration):
        next(fixture)


def test_crystal_log_detector_requires_specific_error(monkeypatch, tmp_path):
    mod, _ = _import_benchmark(BENCH / 'bench_fluka.py', monkeypatch)
    log = tmp_path / 'rfluka.log'
    log.write_text('unrelated server error')
    assert not mod._crystal_missing_in_log(tmp_path)
    log.write_text('Abort called from CRYINI reason Crystal not found Run stopped!')
    assert mod._crystal_missing_in_log(tmp_path)


def test_optional_markers_collected_by_default():
    ini = (BENCH / 'pytest.ini').read_text()
    assert 'not fluka and not geant4 and not gpu' not in ini
    assert all(f'    {x}:' in ini for x in ('fluka', 'geant4', 'gpu'))


def test_runtime_axis_is_sticky_and_separate():
    html = (BENCH / 'perf/templates/dashboard.html').read_text()
    assert '.trend-y-axis{flex:' in html
    assert 'position:sticky;left:0' in html
    assert 'const yAxis=svg(' in html
    assert 'headers(pic,rr,L,0)' in html


def test_external_workflows_do_not_publish_when_every_case_skips():
    for engine in ('fluka', 'geant4'):
        workflow = (ROOT / '.github/workflows' / f'xcoll-performance-{engine}.yml').read_text()
        assert 'id: check_measurements' in workflow
        assert 'steps.check_measurements.outputs.measured' in workflow
        assert f'needs.benchmark_{engine}.outputs.measured' in workflow