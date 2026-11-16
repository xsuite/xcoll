"""Unit tests for the trusted PR resolver and dashboard publication records."""
import base64
import json
import subprocess

import pytest

from benchmarks.perf import pr_extra as p

BASE, HEAD, MERGE = 'a' * 40, 'b' * 40, 'c' * 40


def sample_pr():
    return dict(number=12, state='open', draft=False, mergeable=True,
                base=dict(ref='main', sha=BASE, repo={'full_name': p.REPO}),
                head=dict(ref='contributor', sha=HEAD), merge_commit_sha=MERGE)


def test_resolve_requires_exact_reviewed_head(monkeypatch):
    pr = sample_pr()
    monkeypatch.setattr(p, 'api', lambda path: (
        {'object': {'sha': BASE}} if path == 'git/ref/heads/main' else pr))
    monkeypatch.setattr(p, 'parents', lambda _: [BASE, HEAD])
    assert p.current_merge(12, HEAD)[1:] == (BASE, HEAD, MERGE)
    with pytest.raises(ValueError, match='head changed'):
        p.current_merge(12, 'd' * 40)
    pr['draft'] = True
    with pytest.raises(ValueError, match='non-draft'):
        p.current_merge(12, HEAD)



def test_current_merge_rejects_stale_pr_base(monkeypatch):
    pr = sample_pr()
    new_base = 'd' * 40
    monkeypatch.setattr(p, 'api', lambda path: (
        {'object': {'sha': new_base}} if path == 'git/ref/heads/main' else
        {'object': {'sha': MERGE}} if path == 'git/ref/pull/12/merge' else pr))
    monkeypatch.setattr(p, 'parents', lambda _: [BASE, HEAD])
    monkeypatch.setattr(p.time, 'sleep', lambda _: None)
    with pytest.raises(ValueError, match='Refusing to test an older base'):
        p.current_merge(12)


def test_current_merge_waits_for_fresh_merge(monkeypatch):
    pr = sample_pr()
    new_base, new_merge = 'd' * 40, 'e' * 40
    # The PR's base payload and mergeability initially refer to old main.
    stale = dict(pr, mergeable=False)
    fresh = dict(pr, merge_commit_sha=new_merge)
    responses = iter([stale, pr, fresh])
    monkeypatch.setattr(p, 'api', lambda path: (
        {'object': {'sha': new_base}} if path == 'git/ref/heads/main' else
        {'object': {'sha': MERGE}} if path == 'git/ref/pull/12/merge' else next(responses)))
    monkeypatch.setattr(p, 'parents', lambda merge: (
        [BASE, HEAD] if merge == MERGE else [new_base, HEAD]))
    monkeypatch.setattr(p.time, 'sleep', lambda _: None)
    # Even if base.sha is stale, the verified parents and branch tip are enough.
    assert p.current_merge(12, HEAD)[1:] == (new_base, HEAD, new_merge)



def test_current_merge_uses_fresh_merge_ref_despite_stale_metadata(monkeypatch):
    new_base, new_merge = 'd' * 40, 'e' * 40
    pr = sample_pr()
    monkeypatch.setattr(p, 'api', lambda path: (
        {'object': {'sha': new_base}} if path == 'git/ref/heads/main' else
        {'object': {'sha': new_merge}} if path == 'git/ref/pull/12/merge' else pr))
    monkeypatch.setattr(p, 'parents', lambda merge: (
        [BASE, HEAD] if merge == MERGE else [new_base, HEAD]))
    assert p.current_merge(12, HEAD)[1:] == (new_base, HEAD, new_merge)


def test_current_merge_retries_if_main_moves_during_resolution(monkeypatch):
    new_base, new_merge = 'd' * 40, 'e' * 40
    main_tips = iter([BASE, new_base, new_base, new_base])
    prs = iter([sample_pr(), dict(sample_pr(), merge_commit_sha=new_merge)])
    monkeypatch.setattr(p, 'api', lambda path: (
        {'object': {'sha': next(main_tips)}} if path == 'git/ref/heads/main' else
        {'object': {'sha': MERGE}} if path == 'git/ref/pull/12/merge' else next(prs)))
    monkeypatch.setattr(p, 'parents', lambda merge: (
        [BASE, HEAD] if merge == MERGE else [new_base, HEAD]))
    monkeypatch.setattr(p.time, 'sleep', lambda _: None)
    assert p.current_merge(12)[1:] == (new_base, HEAD, new_merge)


def test_reference_manifest_required(monkeypatch):
    source = 'd' * 40
    encoded = base64.b64encode(json.dumps({'sha': source}).encode()).decode()
    monkeypatch.setattr(p.subprocess, 'run', lambda *a, **kw:
        subprocess.CompletedProcess(a, 0, stdout=json.dumps({'content': encoded})))
    assert p.reference_sha(BASE) == source
    monkeypatch.setattr(p.subprocess, 'run', lambda *a, **kw:
        subprocess.CompletedProcess(a, 1, stderr='gh: Not Found (HTTP 404)'))
    with pytest.raises(ValueError, match='No committed physics reference'):
        p.reference_sha(BASE)
    monkeypatch.setattr(p.subprocess, 'run', lambda *a, **kw:
        subprocess.CompletedProcess(a, 1, stderr='gh: Forbidden (HTTP 403)'))
    with pytest.raises(ValueError, match='Cannot read baseline manifest'):
        p.reference_sha(BASE)


@pytest.mark.parametrize('event_name', ['workflow_dispatch', 'pull_request'])
def test_resolve_physics_stops_before_context_when_reference_missing(tmp_path, monkeypatch, event_name):
    env = dict(GITHUB_REPOSITORY=p.REPO, KIND='physics', PR_NUMBER='12',
               GITHUB_EVENT_NAME=event_name, GITHUB_REF='refs/heads/main',
               GITHUB_SHA=('d' * 40 if event_name == 'workflow_dispatch' else MERGE),
               GITHUB_RUN_ID='100', GITHUB_RUN_ATTEMPT='1',
               RUNNER_TEMP=str(tmp_path), GITHUB_OUTPUT=str(tmp_path / 'output'))
    for key, value in env.items():
        monkeypatch.setenv(key, value)
    pr = sample_pr()
    monkeypatch.setattr(p, 'existing_merge', lambda *a: (pr, BASE, HEAD, MERGE))
    monkeypatch.setattr(p, 'tools_checkout_sha', lambda: 'd' * 40)
    monkeypatch.setattr(p, 'main_contains', lambda _: True)
    monkeypatch.setattr(p, 'event', lambda: {'pull_request': pr})
    monkeypatch.setattr(p, 'parents', lambda _: [BASE, HEAD])

    def missing(base):
        assert base == 'd' * 40
        raise ValueError('No committed physics reference')

    monkeypatch.setattr(p, 'reference_sha', missing)
    with pytest.raises(ValueError, match='No committed physics reference'):
        p.resolve()
    assert not (tmp_path / 'pr-context.json').exists()
    assert not (tmp_path / 'output').exists()


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


def test_existing_merge_does_not_wait_for_current_main(monkeypatch):
    pr = sample_pr()
    calls = []

    def api(path):
        calls.append(path)
        return {'object': {'sha': MERGE}} if path == 'git/ref/pull/12/merge' else pr

    monkeypatch.setattr(p, 'api', api)
    monkeypatch.setattr(p, 'parents', lambda _: [BASE, HEAD])
    monkeypatch.setattr(p.time, 'sleep', lambda _: pytest.fail('Must not poll'))
    assert p.existing_merge(12)[1:] == (BASE, HEAD, MERGE)
    assert calls == ['pulls/12', 'git/ref/pull/12/merge']
    # The actual merge first parent wins over stale PR API base metadata.
    pr['base']['sha'] = 'e' * 40
    assert p.existing_merge(12)[1:] == (BASE, HEAD, MERGE)


def test_existing_merge_rejects_old_head_and_conflict(monkeypatch):
    pr = sample_pr()
    monkeypatch.setattr(p, 'api', lambda path: (
        {'object': {'sha': MERGE}} if path == 'git/ref/pull/12/merge' else pr))
    monkeypatch.setattr(p, 'parents', lambda _: [BASE, 'd' * 40])
    with pytest.raises(ValueError, match='does not match current PR head'):
        p.existing_merge(12)
    monkeypatch.setattr(p, 'parents', lambda _: [BASE, HEAD])
    pr['mergeable'] = False
    with pytest.raises(ValueError, match='merge conflicts'):
        p.existing_merge(12)


@pytest.mark.parametrize('event_name', ['workflow_dispatch', 'pull_request'])
def test_resolve_separate_tools_and_reference(tmp_path, monkeypatch, event_name):
    tools, source = 'd' * 40, 'e' * 40
    env = dict(GITHUB_REPOSITORY=p.REPO, KIND='physics', PR_NUMBER='12',
               GITHUB_EVENT_NAME=event_name, GITHUB_REF='refs/heads/main',
               GITHUB_SHA=(tools if event_name == 'workflow_dispatch' else MERGE),
               GITHUB_RUN_ID='100', GITHUB_RUN_ATTEMPT='1',
               RUNNER_TEMP=str(tmp_path), GITHUB_OUTPUT=str(tmp_path / 'output'))
    for key, value in env.items():
        monkeypatch.setenv(key, value)
    pr = sample_pr()
    monkeypatch.setattr(p, 'existing_merge', lambda *a: (pr, BASE, HEAD, MERGE))
    monkeypatch.setattr(p, 'event', lambda: {'pull_request': pr})
    monkeypatch.setattr(p, 'parents', lambda _: [BASE, HEAD])
    monkeypatch.setattr(p, 'tools_checkout_sha', lambda: tools)
    monkeypatch.setattr(p, 'main_contains', lambda commit: commit in (tools, BASE))

    def reference_sha(commit):
        assert commit == tools
        return source

    monkeypatch.setattr(p, 'reference_sha', reference_sha)
    p.resolve()
    ctx = json.loads((tmp_path / 'pr-context.json').read_text())
    assert (ctx['base_sha'], ctx['head_sha'], ctx['candidate_sha']) == (BASE, HEAD, MERGE)
    assert (ctx['tools_sha'], ctx['reference_sha']) == (tools, source)
    assert f'tools={tools}\nreference={source}' in (tmp_path / 'output').read_text()
    # A separate reference/tooling revision survives publisher validation.
    monkeypatch.setattr(p, 'api', lambda _: pr)
    run = dict(id=100, run_attempt=1, name='Xcoll physics (PR)', event=event_name,
               head_sha=env['GITHUB_SHA'], head_branch='main', pull_requests=[{'number': 12}])
    assert p.load_context(tmp_path / 'pr-context.json', run) == ctx
    monkeypatch.setattr(p, 'download', lambda *args, **kwargs: False)
    record = p.physics_record(ctx, dict(run, conclusion='failure',
                             updated_at='2026-11-16T10:00:00Z'), tmp_path)
    assert (record['base_sha'], record['tools_sha'], record['reference_sha']) == (BASE, tools, source)
    # The artifact cannot substitute its own reference SHA.
    (tmp_path / 'pr-context.json').write_text(json.dumps(dict(ctx, reference_sha=BASE)))
    with pytest.raises(ValueError, match='Reference source differs'):
        p.load_context(tmp_path / 'pr-context.json', run)


def test_main_contains_checks_ancestry(monkeypatch):
    tip = 'd' * 40
    monkeypatch.setattr(p, 'main_sha', lambda: tip)
    monkeypatch.setattr(p, 'api', lambda path: {'merge_base_commit': {'sha': BASE}})
    assert p.main_contains(tip)
    assert p.main_contains(BASE)
    assert not p.main_contains(HEAD)
