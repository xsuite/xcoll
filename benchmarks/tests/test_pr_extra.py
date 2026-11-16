"""Unit tests for the trusted PR resolver and dashboard publication records."""
import json

import pytest

from benchmarks.perf import pr_extra as p

BASE, HEAD, MERGE = 'a' * 40, 'b' * 40, 'c' * 40


def sample_pr():
    return dict(number=12, state='open', draft=False, mergeable=True,
                base=dict(ref='main', sha=BASE, repo={'full_name': p.REPO}),
                head=dict(ref='contributor', sha=HEAD), merge_commit_sha=MERGE)


def test_resolve_requires_exact_reviewed_head(monkeypatch):
    pr = sample_pr()
    monkeypatch.setattr(p, 'api', lambda _: pr)
    monkeypatch.setattr(p, 'parents', lambda _: [BASE, HEAD])
    assert p.current_merge(12, HEAD)[1:] == (BASE, HEAD, MERGE)
    with pytest.raises(ValueError, match='head changed'):
        p.current_merge(12, 'd' * 40)
    pr['draft'] = True
    with pytest.raises(ValueError, match='non-draft'):
        p.current_merge(12, HEAD)


def test_admin_gate_fails_closed(monkeypatch):
    monkeypatch.setattr(p, 'api', lambda _: {'permission': 'write'})
    with pytest.raises(ValueError, match='not a repository administrator'):
        p.assert_admin('developer')
    monkeypatch.setattr(p, 'api', lambda _: {'permission': 'admin'})
    p.assert_admin('administrator')
    with pytest.raises(ValueError, match='Invalid GitHub actor'):
        p.assert_admin('bad/name')


def test_matrix_partition():
    groups = [p.matrix_cases(k) for k in p.MATRIX]
    assert [len(g) for g in groups] == [24, 8, 8]
    assert set.union(*groups) == p.expected_cases()
    assert sum(map(len, groups)) == len(p.expected_cases()) == 40


def test_junit_parses_checks_and_rejects_duplicates(tmp_path):
    name = 'block-cfc-26gev'
    path = tmp_path / 'physics-pr.xml'
    case = f'''<testcase name="test_everest_physics_against_reference[{name}]" time="0.3">
    <failure message="screening failed"/><system-out>  survival fraction difference  0.02 limit=0.01 REVIEW
    state -1 fraction difference 0 limit=0.01 OK
    </system-out></testcase>'''
    path.write_text('<testsuites><testsuite>' + case + '</testsuite></testsuites>')
    row = p.xml_rows(path, {name})[name]
    assert row['status'] == 'FAILED'
    assert row['seconds'] == 0.3
    assert row['checks'][0] == dict(label='survival fraction difference', value=.02, limit=.01, status='REVIEW')
    assert row['checks'][1]['status'] == 'OK'
    path.write_text('<testsuite>' + case + case + '</testsuite>')
    with pytest.raises(ValueError, match='duplicate'):
        p.xml_rows(path, {name})
    path.write_text('<!DOCTYPE test [<!ENTITY x "data">]><testsuite/>')
    with pytest.raises(ValueError, match='DTD'):
        p.xml_rows(path, {name})


def test_missing_matrix_report_is_not_success(tmp_path, monkeypatch):
    monkeypatch.setattr(p, 'download', lambda *args, **kwargs: False)
    monkeypatch.setattr(p, 'reference_sha', lambda _: BASE)
    ctx = dict(base_sha=BASE, pr=12, kind='physics')
    run = dict(id=100, run_attempt=1, conclusion='success', updated_at='2026-11-16T08:00:00Z')
    record = p.physics_record(ctx, run, tmp_path)
    assert record['conclusion'] == 'failure'
    assert len(record['rows']) == 40
    assert all(r['status'] == 'NOT_RUN' for r in record['rows'])


def test_context_parent_and_run_validation(tmp_path, monkeypatch):
    ctx = dict(schema=1, pr=12, branch='contributor', kind='physics', base_sha=BASE,
               head_sha=HEAD, candidate_sha=MERGE, run_id='100', attempt=1)
    path = tmp_path / 'pr-context.json'
    path.write_text(json.dumps(ctx))
    run = dict(id=100, run_attempt=1, name='Xcoll physics (PR)', event='pull_request',
               head_sha=HEAD, pull_requests=[{'number': 12}])
    monkeypatch.setattr(p, 'api', lambda _: sample_pr())
    monkeypatch.setattr(p, 'parents', lambda _: [BASE, HEAD])
    assert p.load_context(path, run) == ctx
    with pytest.raises(ValueError, match='another workflow attempt'):
        p.load_context(path, dict(run, run_attempt=2))
    monkeypatch.setattr(p, 'parents', lambda _: [BASE, 'd' * 40])
    with pytest.raises(ValueError, match='merge parents'):
        p.load_context(path, run)


def test_external_record_cannot_change_approved_revision(tmp_path, monkeypatch):
    context = dict(schema=1, pr=12, branch='contributor', kind='fluka',
                   base_sha=BASE, head_sha=HEAD, candidate_sha=MERGE,
                   run_id='100', attempt=1)
    record = dict(context, rows=[dict(name='bench_fluka', base=.1, candidate=.11,
                                     delta=10., timing_status='Comparable')])
    artifact = tmp_path / 'pr-comparison.json'
    artifact.write_text(json.dumps(record))
    monkeypatch.setattr(p, 'download', lambda *args, **kwargs: True)
    run = dict(id=100, run_attempt=1)
    clean = p.benchmark_record(context, run, tmp_path)
    assert clean['profile'] == 'fluka' and clean['kind'] == 'benchmark'
    assert clean['rows'][0]['delta'] == 10.
    artifact.write_text(json.dumps(dict(record, head_sha='d' * 40)))
    with pytest.raises(ValueError, match='head_sha differs'):
        p.benchmark_record(context, run, tmp_path)
