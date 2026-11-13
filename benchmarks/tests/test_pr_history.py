"""Trusted PR history storage and standalone dashboard tests."""
import json
from pathlib import Path

import pytest

from benchmarks.perf import pr_store, report

SHA = 'a' * 40


def sample():
    return {'schema': 1, 'pr': 241, 'run_id': '12345', 'attempt': 1,
            'branch': 'feature', 'base_sha': SHA, 'candidate_sha': SHA,
            'head_sha': SHA, 'datetime': '2026-11-12T11:00:00',
            'machine': 'GitHub Ubuntu', 'scheduled_reference': 'b' * 64,
            'rows': [{'name': 'bench.py::benchmark_one', 'base': 0.1,
                      'candidate': 0.11, 'delta': 10.0,
                      'timing_status': 'Comparable', 'pair': {'status': 'OK',
                      'summary': 'Within tolerance', 'details': []}}]}


def test_pr_history_validation():
    valid = pr_store.validate(sample(), run_id=12345, attempt=1)
    assert valid['rows'][0]['delta'] == 10.0
    with pytest.raises(ValueError, match='match'):
        pr_store.validate(sample(), run_id=54321, attempt=1)
    bad = sample()
    bad['rows'][0]['delta'] = float('nan')
    with pytest.raises(ValueError, match='Non-finite'):
        pr_store.validate(bad, run_id=12345, attempt=1)
    bad = sample()
    bad['rows'] *= 2
    with pytest.raises(ValueError, match='Duplicate'):
        pr_store.validate(bad, run_id=12345, attempt=1)


def test_read_empty_and_saved_pr_history(tmp_path):
    assert report.read_pr_records(tmp_path/'nonexistent') == []
    folder = tmp_path/'pr-runs'/'241'
    folder.mkdir(parents=True)
    (folder/'12345-1.json').write_text(json.dumps(sample()), encoding='utf-8')
    records = report.read_pr_records(tmp_path/'pr-runs')
    assert len(records) == 1
    assert records[0]['pr'] == 241


def test_dashboard_contains_pr_tab_and_sanitises_records(tmp_path):
    saved = sample()
    saved['rows'][0]['name'] = '</script><script>alert(1)</script>'
    folder = tmp_path/'pr-runs'/'241'
    folder.mkdir(parents=True)
    (folder/'12345-1.json').write_text(json.dumps(saved), encoding='utf-8')
    timing = {'datetime': '2026-11-12T11:00:00',
              'machine_info': {'system': 'Linux', 'processor': 'Xeon'},
              'benchmarks': [{'name': 'ordinary', 'stats': {'median': 0.01}}]}
    path = tmp_path/'run.json'
    path.write_text(json.dumps(timing), encoding='utf-8')
    html = report.create_dashboard([path], tmp_path/'index.html',
                                   pr_runs=tmp_path/'pr-runs').read_text()
    assert 'id="prMode"' in html
    assert 'id="prPanel"' in html
    assert 'id="prRun"' in html
    assert '</script><script>alert(1)' not in html
    assert '\\u003c/script>' in html


def test_compiler_and_version_provenance_kept():
    data = sample()
    data['base_compiler'] = {'cc_version': 'GCC 16'}
    data['candidate_versions'] = {'xcoll': '0.14.0'}
    saved = pr_store.validate(data, run_id=12345, attempt=1)
    assert saved['base_compiler']['cc_version'] == 'GCC 16'
    assert saved['candidate_versions']['xcoll'] == '0.14.0'
