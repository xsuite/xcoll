"""Integration tests for storage, daily dashboard, and CI metadata (no Xsuite needed)."""
from __future__ import annotations

import copy
import gzip
import json
import sys
from pathlib import Path

import pytest
import yaml

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from benchmarks.perf import report, rolling_physics, results


def sample(*, profile='full', trigger='schedule', branch='main', sha='a'*40,
           time='2026-11-10T10:00:00', name='bench_example', lost=0.3,
           testbed='cern-openstack-xcoll-cpu-v1'):
    return {
        'datetime': time, 'testbed': testbed,
        'machine_info': {'node': 'changed-fqdn', 'system': 'Linux', 'processor': 'Xeon'},
        'ci': {'run_id': '42800000000', 'run_attempt': 1,
               'source_sha': sha, 'source_branch': branch,
               'profile': profile, 'trigger': trigger,
               'physical_host': 'cern-openstack-xcoll-vm-v1'},
        'benchmarks': [{
            'name': name, 'fullname': f'test.py::{name}', 'group': 'test',
            'stats': {'median': .01, 'stddev': .0001, 'rounds': 4, 'data': [.01, .02]},
            'extra_info': {'protocol': 3, 'particles': 500,
                           'physics_case': name, 'lost_fraction': lost,
                           'raw_physics': [1, 2, 3]},
        }],
    }


def write(tmp_path, name, data):
    p = tmp_path / name
    p.write_text(json.dumps(data))
    return p


def test_online_groups_only_scheduled_profiles(tmp_path):
    p = []
    for profile,sha,time in [('full','a'*40,'2026-11-10T04:30:00'),
                             ('fluka','b'*40,'2026-11-10T06:30:00'),
                             ('geant4','c'*40,'2026-11-10T09:30:00')]:
        d = sample(profile=profile, sha=sha, time=time, name=f'bench_{profile}',
                   testbed=f'cern-openstack-xcoll-{profile}-v1')
        p.append(write(tmp_path,profile+'.json', d))
    manual = sample(profile='quick', trigger='workflow_dispatch',
                    name='manual_case', time='2026-11-10T10:30:00')
    p.append(write(tmp_path,'manual.json',manual))
    assert len(report.read_runs(p)) == 4
    merged = report.read_runs(p, online=True)
    assert len(merged) == 2
    daily = next(r for r in merged if r['profile'] == 'daily')
    assert {e['name'] for e in daily['entries']} == {f'test.py::bench_{x}' for x in ('full','fluka','geant4')}
    assert len(daily['components']) == 3
    assert daily['sha'] == 'mixed'
    assert {e['profile'] for e in daily['entries']} == {'full','fluka','geant4'}


def test_local_hostname_not_part_of_machine_key(tmp_path):
    a = sample(trigger='workflow_dispatch', testbed='unlabelled-local')
    b = copy.deepcopy(a)
    b['machine_info']['node'] = 'different-network-name'
    x = report.read_runs([write(tmp_path,'a.json',a),write(tmp_path,'b.json',b)])
    assert x[0]['machine'] == x[1]['machine']
    assert len(report.merge_daily(x)) == 2


def test_rolling_baseline_changes_and_compacts(tmp_path):
    first=sample()
    f=write(tmp_path,'original.json',first)
    base=tmp_path/'baseline.json.gz'
    out1=tmp_path/'compact-1.json'
    statuses=rolling_physics.process(f,baseline_out=base,history_out=out1)
    assert next(iter(statuses.values()))['status']=='BASELINE'
    assert base.is_file()
    old=gzip.decompress(base.read_bytes())
    assert b'raw_physics' not in old
    old_result=json.loads(out1.read_text())
    bench=old_result['benchmarks'][0]
    assert 'raw_physics' not in bench['extra_info']
    assert 'data' not in bench['stats']
    # Original artefact remains unchanged; only the published compact copy changes.
    assert 'raw_physics' in json.loads(f.read_text())['benchmarks'][0]['extra_info']
    other=sample(lost=.7)
    other['ci']['run_id']='42800000001'
    f2=write(tmp_path,'next.json',other)
    out2=tmp_path/'compact-2.json'
    statuses=rolling_physics.process(f2,previous=base,baseline_out=tmp_path/'next.gz',history_out=out2)
    assert next(iter(statuses.values()))['status']=='CHANGED'
    assert len(next(iter(statuses.values()))['details']) > 0
    assert json.loads(out2.read_text())['benchmarks'][0]['extra_info']['physics_check']['status']=='CHANGED'


def test_manual_run_does_not_change_baseline(tmp_path):
    f=write(tmp_path,'manual.json',sample(trigger='workflow_dispatch',profile='quick'))
    base=tmp_path/'baseline.json.gz'
    st=rolling_physics.process(f,baseline_out=base,history_out=tmp_path/'published.json')
    assert not base.exists()
    assert next(iter(st.values()))['status']=='NOT_CHECKED'


def test_normalised_baseline_not_double_normalised(tmp_path):
    first=sample(lost=.3)
    first['benchmarks'][0]['extra_info']['reference_collimator_losses']=50
    initial=write(tmp_path,'initial.json',first)
    baseline=tmp_path/'baseline.gz'
    rolling_physics.process(initial,baseline_out=baseline)
    nextdata=copy.deepcopy(first)
    nextdata['benchmarks'][0]['extra_info']['reference_collimator_losses']=50
    f=write(tmp_path,'next.json',nextdata)
    statuses=rolling_physics.process(f,previous=baseline)
    assert next(iter(statuses.values()))['status']=='OK'


def test_dashboard_valid_html_and_sanitises(tmp_path):
    d=sample(trigger='workflow_dispatch')
    d['benchmarks'][0]['extra_info']['raw_physics']='</script><script>alert(1)</script>'
    result=write(tmp_path,'x.json',d)
    html=report.create_dashboard([result],tmp_path/'index.html').read_text()
    assert 'Compare two sets' in html
    assert 'compareTable' in html
    assert '</script><script>alert(1)' not in html
    assert 'id="machine"' not in html  # no machine selection


def test_results_list_prints_filename(tmp_path,capsys):
    d=sample()
    root=tmp_path/'results';root.mkdir()
    write(root,'0001_bench.json',d)
    assert results.main(['--results',str(root),'list'])==0
    out=capsys.readouterr().out
    assert 'Filename' in out and '0001_bench.json' in out


def test_workflow_triggers_paths_and_retrieval():
    workflows=ROOT/'.github/workflows'
    for name in ['xcoll-performance-fluka.yml','xcoll-performance-geant4.yml']:
        text=(workflows/name).read_text()
        d=yaml.load(text,Loader=yaml.BaseLoader)
        assert 'schedule' in d['on'] and 'workflow_dispatch' in d['on']
        assert 'push' not in d['on']
        assert 'runner.temp' not in text
        assert '--skip-tpsa --skip-synrad --n-threads 2' in text
    for name in ['xcoll-performance.yml','xcoll-performance-publish.yml']:
        text=(workflows/name).read_text()
        assert 'github-token: ${{ github.token }}' in text
        assert 'actions: read' in text
        assert 'rolling_physics' in text and '--online' in text


def test_gpu_and_engine_case_matrix():
    text=(ROOT/'benchmarks/pytest.ini').read_text()
    assert 'not gpu' not in text and 'gpu:' in text
    for name in ('gpu','fluka','geant4'):
        src=(ROOT/f'benchmarks/bench_{name}.py').read_text()
        assert f'@pytest.mark.{name}' in src
    assert 'FlukaCrystal' in (ROOT/'benchmarks/bench_fluka.py').read_text()
    assert 'Geant4CollimatorTip' in (ROOT/'benchmarks/bench_geant4.py').read_text()
