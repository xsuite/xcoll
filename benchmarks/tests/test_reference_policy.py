"""Regression tests for the heatmap reference, manual physics, and paired PRs."""
import copy
import json
import shutil
import subprocess
from pathlib import Path

import pytest

from benchmarks.perf import pr_compare, rolling_physics


def sample(*, trigger='schedule', profile='full', median=.1, lost=.2, testbed='cern-openstack-xcoll-cpu-v1'):
    return {'datetime': '2026-11-12T03:17:00', 'testbed': testbed,
            'ci': {'trigger': trigger, 'profile': profile, 'source_branch': 'main', 'source_sha': 'a'*40},
            'benchmarks': [{'name': 'case', 'fullname': 'bench_example.py::bench_case',
                            'stats': {'median': median},
                            'extra_info': {'particles': 500, 'protocol': 4, 'lost_fraction': lost}}]}


def record(tmp_path, name, data):
    path = tmp_path / name
    path.write_text(json.dumps(data))
    return path


def test_manual_quick_reads_full_reference_without_writing(tmp_path):
    baseline = tmp_path / 'full.json.gz'
    rolling_physics.process(record(tmp_path, 'first.json', sample()), baseline_out=baseline)
    original = baseline.read_bytes()
    manual = sample(trigger='workflow_dispatch', profile='quick', lost=.8, testbed='github-hosted')
    result = rolling_physics.process(record(tmp_path, 'manual.json', manual),
                                     previous=baseline, baseline_out=baseline)
    assert next(iter(result.values()))['status'] == 'CHANGED'
    assert baseline.read_bytes() == original


def test_manual_missing_or_wrong_profile(tmp_path):
    baseline = tmp_path / 'full.json.gz'
    rolling_physics.process(record(tmp_path, 'first.json', sample()), baseline_out=baseline)
    manual = sample(trigger='workflow_dispatch', profile='fluka')
    path = record(tmp_path, 'manual.json', manual)
    assert next(iter(rolling_physics.process(path, previous=baseline).values()))['status'] == 'NOT_CHECKED'
    assert next(iter(rolling_physics.process(path).values()))['status'] == 'NOT_CHECKED'


def test_pr_pair_and_scheduled_physics(tmp_path):
    baseline_path = tmp_path / 'full.json.gz'
    rolling_physics.process(record(tmp_path, 'daily.json', sample(lost=.2)), baseline_out=baseline_path)
    baseline = rolling_physics._read_baseline(baseline_path)
    base = sample(median=.1, lost=.8)
    candidate = sample(median=.11, lost=.8)
    row, = pr_compare.compare(base, candidate, baseline)
    assert row['delta'] == pytest.approx(10)
    assert row['pair'].status == 'OK'
    assert row['scheduled'].status == 'CHANGED'
    incompatible = copy.deepcopy(candidate)
    incompatible['benchmarks'][0]['extra_info']['particles'] = 100
    row, = pr_compare.compare(base, incompatible, baseline)
    assert row['delta'] is None
    assert row['pair'].status == 'INCOMPATIBLE'
    assert row['scheduled'].status == 'INCOMPATIBLE'
    assert 'Physics vs scheduled' in pr_compare.markdown([row], 'a'*40, 'b'*40, scheduled_available=True)


def test_pr_new_missing_cases():
    base, candidate = sample(), sample()
    candidate['benchmarks'][0]['fullname'] = 'new'
    rows = pr_compare.compare(base, candidate)
    assert {r['timing_status'] for r in rows} == {'New in PR', 'Missing in PR'}
    assert all(r['delta'] is None for r in rows)


@pytest.mark.skipif(shutil.which('node') is None, reason='Node.js not available for JS policy test')
def test_dashboard_reference_policy():
    template = (Path(__file__).resolve().parents[1] / 'perf/templates/dashboard.html').read_text()
    # Exercise the actual browser reference function without a DOM dependency.
    script = template[template.index('function pointTime('):template.index('function selectedReference(')]
    checks = """
const assert = require('node:assert/strict');
function compatible(a,b){return a.run.machine===b.run.machine&&a.e.signature===b.e.signature;}
function p(id,day,time,median,trigger='schedule',machine='vm',signature='s',branch='main'){
 return {run:{id,datetime:day+'T23:00:00',trigger,machine,branch},e:{median,signature,datetime:day+'T'+time+':00'}};
}
const a=p('a','2026-11-10','04:00',100), b=p('b','2026-11-11','04:00',110);
const m=p('m','2026-11-11','12:00',150,'workflow_dispatch');
const c=p('c','2026-11-12','04:00',110);
const x=p('x','2026-11-11','05:00',1000,'schedule','other');
const y=p('y','2026-11-11','05:00',1000,'schedule','vm','different');
const all=[a,b,m,c,x,y];
assert.equal(heatReference(c,[c],all,'scheduled').value,110);
assert.equal(heatReference(b,[b,c],all,'scheduled').value,100);
assert.equal(heatReference(m,[m],all,'scheduled').value,110);
assert.equal(heatReference(a,[a],all,'scheduled'),null);
assert.equal(heatReference(c,[b,c],all,'first').value,110);
assert.equal(heatReference(c,[a,b,m,c,x,y],all,'mean').value,117.5);
assert.equal(heatReference(c,[a,b,m,c,x,y],all,'median').value,110);
assert.equal(heatReference(c,[a,c],all,'mean').value,105);
const local=p('local','2026-11-10','06:00',80,'workflow_dispatch','local');
const later=p('later','2026-11-12','06:00',90,'workflow_dispatch','local');
assert.equal(heatReference(later,[later],[local,later],'scheduled').value,80);
"""
    subprocess.run(['node','-e',script+checks], check=True, capture_output=True, text=True)
