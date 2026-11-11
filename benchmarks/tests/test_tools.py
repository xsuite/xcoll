"""Infrastructure tests: no Xsuite or pytest-benchmark installation required."""

from __future__ import annotations

import ast
import json
import sys
from pathlib import Path

import pytest

BENCH = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(BENCH))

from perf import metrics, compare_physics, report, results as results_cli, history
import conftest
from perf.prebuild_cpu import select_kernel_names


def entry(*, name='bench_a', lost=0.3, particles=100000, **changes):
    extra = {'protocol': 2, 'particles': particles, 'lost_fraction': lost,
             'survivors_with_kick_gt_0.1urad': 0.2}
    extra.update(changes)
    return {'name': name, 'fullname': f'test.py::{name}',
            'stats': {'median': 0.01, 'stddev': 0.0001, 'rounds': 5},
            'extra_info': extra}


def doc(entries=None, *, sha='f'*40, branch='main', date='2026-11-09T12:00:00'):
    return {'datetime': date,
            'machine_info': {'node': 'runner-123', 'system': 'Linux', 'processor': 'Xeon'},
            'testbed': 'github-ubuntu-24.04-cpu-informational',
            'commit_info': {'id': sha, 'branch': branch},
            'benchmarks': entries if entries is not None else [entry()]}


def test_generic_config_not_code_specific(tmp_path):
    config = tmp_path / 'metric.toml'
    config.write_text('''[comparison]\nworkload_fields = ["order", "algorithm"]\n'''
                      '''[observables.symplecticity_error]\nsource = "symp_err"\nabsolute_tolerance = 0.0001\n''')
    specs = metrics.load_config(config)
    before = entry(order=5, algorithm='A', symp_err=.0002)
    after = entry(order=5, algorithm='A', symp_err=.00025)
    assert compare_physics.compare_entry(before, after, specs).status == 'OK'
    after['extra_info']['symp_err'] = .1
    assert compare_physics.compare_entry(before, after, specs).status == 'CHANGED'


def test_physics_ok_and_changed():
    assert compare_physics.compare_documents(doc(), doc([entry(lost=.305)]))[0].status == 'OK'
    assert compare_physics.compare_documents(doc(), doc([entry(lost=.4)]))[0].status == 'CHANGED'


def test_missing_and_workload_mismatch():
    assert compare_physics.compare_documents(doc(), doc([]))[0].status == 'MISSING'
    assert compare_physics.compare_documents(doc(), doc([entry(particles=300)]))[0].status == 'INCOMPATIBLE'


def test_unverified_physics():
    a=entry()
    b=entry()
    del a['extra_info']['lost_fraction']
    assert compare_physics.compare_documents(doc([a]), doc([b]))[0].status == 'UNVERIFIED'


def test_normalised_counts():
    a=entry(reference_collimator_losses=100)
    b=entry(reference_collimator_losses=150)
    assert compare_physics.compare_documents(doc([a]),doc([b]))[0].status == 'OK'
    b['extra_info']['reference_collimator_losses'] = 10000
    assert compare_physics.compare_documents(doc([a]),doc([b]))[0].status == 'CHANGED'


def test_dashboard_uses_schema_and_sanitises_script(tmp_path):
    files=[]
    for i, host in enumerate(['runner-123','runner-999']):
        d=doc(date=f'2026-11-0{i+8}T12:00:00', sha=str(i+1)*40)
        d['machine_info']['node']=host
        d['benchmarks'][0]['extra_info']['unexpected'] = '</script><script>alert(1)</script>'
        path=tmp_path/f'run-{i}.json'
        path.write_text(json.dumps(d))
        files.append(path)
    runs=report.read_runs(files)
    assert runs[0]['machine'] == runs[1]['machine']
    assert runs[0]['entries'][0]['observables']['lost_fraction'] == .3
    assert 'collimator_losses_per_particle' not in runs[0]['entries'][0]['observables']
    dest=report.create_dashboard(files, tmp_path/'site'/'index.html')
    html=dest.read_text()
    assert '</script><script>alert(1)' not in html
    assert 'function physics(e)' in html
    assert 'All benchmarks · chronological overview' in html
    assert 'Compare two sets' in html
    assert 'runner-999' not in html


def test_cli_list_tag_compare_resolution(tmp_path, capsys):
    directory=tmp_path/'results'
    directory.mkdir()
    for i in (1,2):
        (directory/f'{i:04}_sample.json').write_text(json.dumps(doc(
            date=f'2026-11-0{i+7}T12:00:00', branch=f'feature-{i}')))
    assert results_cli.main(['--results',str(directory),'list']) == 0
    out=capsys.readouterr().out
    assert 'feature-1' in out and 'feature-2' in out
    assert results_cli.main(['--results',str(directory),'tag','1','release-baseline']) == 0
    rows=results_cli.entries(directory)
    labels=results_cli._load_labels(directory)
    assert results_cli.resolve('release-baseline', directory, rows, labels).name == '0002_sample.json'
    assert results_cli.main(['--results',str(directory),'untag','release-baseline']) == 0


def test_archive_preserves_tags_and_files(tmp_path):
    directory=tmp_path/'results'; storage=directory/'example-machine'; storage.mkdir(parents=True)
    for i in range(1,6):
        (storage/f'{i:04d}_sample.json').write_text(json.dumps(doc(date=f'2026-11-{i:02}T12:00:00')))
    results_cli.main(['--results',str(directory),'tag','5','baseline'])
    moved=results_cli.archive_results(directory, threshold=4)
    assert len(moved) == 5
    assert len(list(storage.glob('*.json'))) == 0
    assert len(results_cli.entries(directory)) == 5
    labels=results_cli._load_labels(directory)
    assert results_cli.resolve('baseline',directory,results_cli.entries(directory),labels).exists()


def test_xdist_guard():
    class Fake:
        def __init__(self, workers=None, tx=None):
            self.workers,self.tx=workers,tx
        def getoption(self, name, default=None):
            return {'numprocesses': self.workers, 'tx': self.tx}.get(name,default)
    for option in (2,'auto',1,'logical'):
        with pytest.raises(pytest.UsageError):
            conftest.pytest_configure(Fake(option))
    with pytest.raises(pytest.UsageError):
        conftest.pytest_configure(Fake(None,['ssh=host']))
    for option in (None,0,'0'):
        conftest.pytest_configure(Fake(option))


def test_history_pack(tmp_path):
    original=tmp_path/'raw'/'machine';original.mkdir(parents=True)
    (original/'0001_file.json').write_text(json.dumps(doc()))
    file=history.pack(tmp_path/'raw',tmp_path/'history','123','1','c'*40,'main','full')
    assert json.loads(file.read_text())['ci']['source_branch'] == 'main'
    with pytest.raises(FileExistsError):
        history.pack(tmp_path/'raw',tmp_path/'history','123','1','c'*40,'main','full')


def test_paths_are_relative_to_package_not_cwd(tmp_path, monkeypatch):
    from perf.paths import BENCHMARK_DIR, results_dir, config_file
    monkeypatch.chdir(tmp_path)
    assert config_file() == BENCHMARK_DIR / 'metrics.toml'
    assert results_dir() == BENCHMARK_DIR / 'results'
    assert metrics.load_config()['fields']


def test_dashboard_accepts_result_directory(tmp_path):
    p=tmp_path/'history'/'runs'
    p.mkdir(parents=True)
    (p/'one.json').write_text(json.dumps(doc()))
    (p/'ignore.json').write_text('{"not-a-benchmark": true}')
    dest=report.create_dashboard([p],tmp_path/'site'/'index.html')
    html=dest.read_text()
    assert 'test.py::bench_a' in html
    assert 'Performance history' in html


def test_revision_section_labels_present():
    template = (BENCH / 'perf' / 'templates' / 'dashboard.html').read_text()
    assert 'function revisionSections(rr)' in template
    assert 'function headers(pic,rr,L,x0=0)' in template
    assert 'sec.branch' in template and 'short(sec.sha)' in template
    assert 'function labelMachine(r)' in template


def test_new_benchmarks_use_serial_pedantic_setup():
    for name in ['bench_integration.py', 'bench_lossmap_variants.py']:
        content = (BENCH / name).read_text()
        assert 'benchmark.pedantic(' in content
        assert 'pytest.mark.' in content
        assert 'def setup():' in content


def test_new_files_are_valid_python():
    for name in ("import", "materials", "rf_sweep", "blowup", "emittance_monitor", "fluka", "geant4", "gpu"):
        path = BENCH / f"bench_{name}.py"
        tree = ast.parse(path.read_text())
        assert any(isinstance(n, ast.FunctionDef) and n.name.startswith("bench_") for n in tree.body)


def test_external_suites_have_separate_marks():
    for engine in ("fluka", "geant4", "gpu"):
        tree = ast.parse((BENCH / f"bench_{engine}.py").read_text())
        assert f"pytest.mark.{engine}" in ast.unparse(tree)


def test_setup_requires_configured_existing_paths(monkeypatch):
    from benchmarks.perf.setup_external import _required

    monkeypatch.delenv("XCOLL_MISSING_TEST_PATH", raising=False)
    with pytest.raises(RuntimeError, match="Missing"):
        _required("XCOLL_MISSING_TEST_PATH")
    monkeypatch.setenv("XCOLL_MISSING_TEST_PATH", "/a/nonexistent/xcoll/benchmark/path")
    with pytest.raises(FileNotFoundError):
        _required("XCOLL_MISSING_TEST_PATH")


def test_fluka_geant4_gpu_skip_when_unavailable():
    config = (BENCH / "pytest.ini").read_text()
    assert '-m "not fluka and not geant4 and not gpu"' not in config
    for mark in ('fluka:', 'geant4:', 'gpu:'):
        assert mark in config


def test_history_accepts_external_profiles():
    src = (BENCH / "perf/history.py").read_text()
    assert '"fluka", "geant4", "gpu"' in src
    ast.parse(src)


def test_excludes_tpsa_by_name_and_config():
    defs = [
        ("default_no_config", {"config": {}}),
        ("tpsa_base_config", {"config": {"XTRACK_TPSA_TRACK": True}}),
        ("renamed_tpsa_kernel", {"config": {"XTRACK_TPSA_TRACK": True}}),
        ("all_with_synrad", {"config": {"XTRACK_MULTIPOLE_NO_SYNRAD": False}}),
    ]
    assert select_kernel_names(defs, skip_tpsa=True, skip_synrad=True) == ["default_no_config"]
    assert len(select_kernel_names(defs)) == 4

