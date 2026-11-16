"""Trusted PR revision resolution and publication of supplementary PR results.

Run resolve/publish on GitHub-hosted runners using main's code only. Candidate
artifacts are data; never import or execute candidate code in the publisher.
"""
from __future__ import annotations

import argparse
import base64
import json
import math
import os
from pathlib import Path
import re
import subprocess
import tempfile
import time
import xml.etree.ElementTree as ET

from .pr_store import validate as validate_benchmarks

REPO = 'xsuite/xcoll'
SHA = re.compile(r'[0-9a-f]{40}')
ENGINES = ('fluka', 'geant4')
WORKFLOWS = {'Xcoll physics (PR)': 'physics',
             'Xcoll external performance (PR)': 'external'}
MATRIX = ('elements', 'lhc-b1', 'lhc-b2')
MAX_XML = 20_000_000


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha(value):
    require(isinstance(value, str) and SHA.fullmatch(value), 'Invalid SHA')
    return value


def number(value):
    value = str(value)
    require(re.fullmatch(r'[1-9][0-9]{0,11}', value), 'Invalid positive integer')
    return int(value)


def gh(*args):
    return subprocess.check_output(['gh', *args], text=True)


def api(path):
    return json.loads(gh('api', f'repos/{REPO}/{path}'))


def event():
    return json.loads(Path(os.environ['GITHUB_EVENT_PATH']).read_text())


def parents(commit):
    return [p['sha'] for p in api(f'git/commits/{sha(commit)}')['parents']]


def eligible(pr):
    return (pr['state'] == 'open' and not pr['draft']
            and pr['base']['repo']['full_name'] == REPO
            and pr['base']['ref'] == 'main'
            and not pr['head']['ref'].startswith('release/'))


def assert_admin(actor):
    require(re.fullmatch(r'[A-Za-z0-9-]{1,39}', actor or ''), 'Invalid GitHub actor')
    require(api(f'collaborators/{actor}/permission')['permission'] == 'admin',
            f'{actor} is not a repository administrator')


def main_sha():
    """Resolve the actual branch tip independently of the PR metadata."""
    return sha(api('git/ref/heads/main')['object']['sha'])


def current_merge(pr_number, approved_head=None):
    """Wait for GitHub's test merge to have the current main and PR head parents."""
    # A pull request's base.sha and merge_commit_sha can both be stale after
    # main advances. Their agreement alone does not establish current main.
    for attempt in range(20):
        base = main_sha()
        pr = api(f'pulls/{pr_number}')
        require(eligible(pr), 'PR must be open, non-draft, and target main (not release/)')
        head = sha(pr['head']['sha'])
        reported_base = sha(pr['base']['sha'])
        if approved_head is not None:
            require(head == approved_head, 'PR head changed since approval; review and dispatch again')
        merge = pr.get('merge_commit_sha')
        if merge and pr.get('mergeable') is True:
            merge = sha(merge)
            if parents(merge) == [base, head] and main_sha() == base:
                return pr, base, head, merge
        # A conflict computed against a stale base need not be a conflict
        # against current main. Wait for GitHub to refresh mergeability too.
        if pr.get('mergeable') is False and reported_base == base and main_sha() == base:
            raise ValueError('PR has merge conflicts against current main')
        # The Git ref can already be refreshed while the PR payload is stale.
        # Never use it without checking both exact parents and the branch tip.
        try:
            merge_ref = sha(api(f'git/ref/pull/{pr_number}/merge')['object']['sha'])
        except subprocess.CalledProcessError:
            # The test merge ref may not exist while GitHub is computing it.
            merge_ref = None
        if (merge_ref and (merge_ref != merge or pr.get('mergeable') is not True)
                and parents(merge_ref) == [base, head] and main_sha() == base):
            return pr, base, head, merge_ref
        print(f'Waiting for PR #{pr_number} merge against main {base}; '
              f'PR reports base {reported_base}, merge {merge}, merge ref {merge_ref} '
              f'(attempt {attempt + 1}/20)', flush=True)
        if attempt < 19:
            time.sleep(3)
    raise ValueError(
        f'GitHub did not produce a merge commit with current main {base} '
        f'and PR head {head} as parents. The PR merge ref may be stale; '
        'check mergeability and dispatch again. Refusing to test an older base.')


def resolve():
    require(os.environ['GITHUB_REPOSITORY'] == REPO, 'Unexpected repository')
    kind = os.environ['KIND']
    require(kind in ('physics', *ENGINES), 'Invalid PR run kind')
    pr_number = number(os.environ['PR_NUMBER'])
    manual = os.environ['GITHUB_EVENT_NAME'] == 'workflow_dispatch'
    if manual:
        require(os.environ['GITHUB_REF'] == 'refs/heads/main', 'Dispatch from main only')
        approved = None
        if kind in ENGINES:
            # Workflow dispatch alone requires write, not admin, access.
            for actor in {os.environ['GITHUB_ACTOR'], os.environ['TRIGGERING_ACTOR']}:
                assert_admin(actor)
            approved = sha(os.environ['APPROVED_HEAD'])
        pr, base, head, merge = current_merge(pr_number, approved)
    else:
        require(kind == 'physics' and os.environ['GITHUB_EVENT_NAME'] == 'pull_request',
                'External PR benchmarks are manual only')
        pr = event()['pull_request']
        require(pr['number'] == pr_number and eligible(pr), 'Ineligible PR')
        base, head, merge = sha(pr['base']['sha']), sha(pr['head']['sha']), sha(os.environ['GITHUB_SHA'])
        require(parents(merge) == [base, head], 'PR merge parents differ from event')
    if kind == 'physics':
        # Fail on the hosted resolver before scheduling any matrix jobs.
        # The manifest's SHA is the accepted physics source, not the base SHA.
        print(f'Using committed physics reference {reference_sha(base)}', flush=True)
    context = dict(schema=1, kind=kind, pr=pr_number, branch=pr['head']['ref'],
                   base_sha=base, head_sha=head, candidate_sha=merge,
                   run_id=str(number(os.environ['GITHUB_RUN_ID'])),
                   attempt=number(os.environ['GITHUB_RUN_ATTEMPT']))
    path = Path(os.environ['RUNNER_TEMP']) / 'pr-context.json'
    path.write_text(json.dumps(context, indent=2) + '\n')
    with open(os.environ['GITHUB_OUTPUT'], 'a') as f:
        print(f'base={base}\ncandidate={merge}\nhead={head}', file=f)
        print(f'pr={pr_number}', file=f)
        print(f"authorised_attempt={context['attempt']}", file=f)
    print(f'PR #{pr_number}: base {base}, head {head}, merge {merge}, kind {kind}')


def load_context(path, run):
    require(path.stat().st_size <= 16_000, 'Oversized context artifact')
    obj = json.loads(path.read_text())
    require(obj.get('schema') == 1, 'Invalid context schema')
    pr = number(obj['pr'])
    require(type(obj['pr']) is int, 'PR number must be an integer')
    require(obj['run_id'] == str(run['id']) and obj['attempt'] == run['run_attempt'],
            'Context belongs to another workflow attempt')
    kind = obj['kind']
    require((run['name'] == 'Xcoll physics (PR)' and kind == 'physics') or
            (run['name'] == 'Xcoll external performance (PR)' and kind in ENGINES),
            'Unexpected workflow/kind')
    base, head, merge = (sha(obj[k]) for k in ('base_sha', 'head_sha', 'candidate_sha'))
    require(parents(merge) == [base, head], 'Recorded merge parents do not match')
    actual = api(f'pulls/{pr}')
    require(actual['base']['ref'] == 'main' and actual['base']['repo']['full_name'] == REPO,
            'PR does not target main')
    require(obj['branch'] == actual['head']['ref'], 'PR branch differs')
    if actual['head']['sha'] != head:
        # Historical records remain valid after a contributor pushes again.
        associated = api(f'commits/{head}/pulls')
        require(any(p['number'] == pr for p in associated), 'Head not associated with PR')
    if run['event'] == 'pull_request':
        require(kind == 'physics' and run['head_sha'] in (head, merge),
                'PR workflow revision differs from artifact')
        prs = run.get('pull_requests') or []
        require(not prs or any(p['number'] == pr for p in prs), 'Workflow PR differs')
    else:
        require(run['event'] == 'workflow_dispatch' and run['head_branch'] == 'main',
                'Not a manual main workflow')
    # Canonical record: discard all other artifact fields.
    return dict(schema=1, kind=kind, pr=pr, branch=actual['head']['ref'],
                base_sha=base, head_sha=head, candidate_sha=merge,
                run_id=str(run['id']), attempt=run['run_attempt'])


def download(run, name, directory, *, optional=False):
    """Exact name and separate destination; only the expected data file is read."""
    directory.mkdir(parents=True, exist_ok=True)
    result = subprocess.run(['gh', 'run', 'download', str(run['id']), '--repo', REPO,
                             '--name', name, '--dir', str(directory)],
                            capture_output=True, text=True)
    if result.returncode and not optional:
        raise RuntimeError(result.stderr)
    return result.returncode == 0


def expected_cases():
    families = ('block-cfc', 'block-tungsten', 'collimator-cfc', 'collimator-tungsten',
                'crystal-channeling', 'crystal-volume-reflection', 'crystal-amorphous', 'crystal-miss')
    return {f'{f}-{e}gev' for f in families for e in (26, 450, 7000)} | {
        f'lhc-b{b}-{p}-{c}-{e}gev' for b in (1, 2) for p in ('h', 'v')
        for c in ('crystals', 'no-crystals') for e in (450, 7000)}


def matrix_cases(matrix):
    return {c for c in expected_cases() if (not c.startswith('lhc-') if matrix == 'elements'
                                           else c.startswith(matrix + '-'))}


def xml_rows(path, expected):
    require(path.stat().st_size <= MAX_XML, 'Oversized JUnit report')
    raw = path.read_text(encoding='utf-8')
    require('<!DOCTYPE' not in raw.upper() and '<!ENTITY' not in raw.upper(), 'DTD/entities forbidden')
    rows = {}
    for node in ET.fromstring(raw).iter('testcase'):
        match = re.fullmatch(r'test_everest_physics_against_reference\[([^\]]+)\]', node.get('name', ''))
        if not match:
            continue  # helper tests still affect the workflow conclusion
        name = match[1]
        require(name in expected and name not in rows, 'Unexpected or duplicate physics case')
        status, details = 'PASSED', ''
        for tag, label in (('failure', 'FAILED'), ('error', 'ERROR'), ('skipped', 'SKIPPED')):
            entry = node.find(tag)
            if entry is not None:
                status = label
                # Do not publish arbitrary full logs, only a bounded diagnostic.
                details = (entry.get('message') or entry.text or '').strip()[:500]
                out = node.findtext('system-out', '')
                flags = [line.strip() for line in out.splitlines() if line.rstrip().endswith('REVIEW')]
                if flags:
                    details = ('; '.join(flags))[:1000]
                break
        # Keep the numerical screening results, not arbitrary console output.
        checks = []
        pattern = re.compile(r'\s*(.{1,200}?)\s+([+\-0-9.eEinfnaINFNA]+)\s+limit=([+\-0-9.eE]+)\s+(OK|REVIEW)\s*')
        for line in node.findtext('system-out', '').splitlines():
            found = pattern.fullmatch(line)
            if found:
                label, value, limit, outcome = found.groups()
                try:
                    v, lim = float(value), float(limit)
                except ValueError:
                    continue
                if math.isfinite(lim) and lim >= 0:
                    checks.append(dict(label=label, value=v if math.isfinite(v) else None,
                                       limit=lim, status=outcome))
            if len(checks) >= 150:
                break
        seconds = float(node.get('time', '0'))
        require(math.isfinite(seconds) and seconds >= 0, 'Invalid JUnit duration')
        rows[name] = dict(name=name, status=status, seconds=seconds, details=details, checks=checks)
    return rows


def reference_sha(base):
    # The trusted base manifest is mandatory. CI must never generate a reference.
    path = f'contents/tests/data/everest_physics_reference/reference.json?ref={sha(base)}'
    proc = subprocess.run(['gh', 'api', f'repos/{REPO}/{path}'], capture_output=True, text=True)
    if proc.returncode:
        if 'HTTP 404' in proc.stderr:
            raise ValueError(
                f'No committed physics reference at PR base {base}. '
                'Update the PR against main with the accepted reference, '
                'or manually dispatch against current main. CI will not generate one.')
        raise ValueError('Cannot read baseline manifest: ' + proc.stderr)
    return sha(json.loads(base64.b64decode(json.loads(proc.stdout)['content']))['sha'])


def physics_record(context, run, directory):
    rows = {}
    suffix = f"{run['id']}-{run['run_attempt']}"
    for matrix in MATRIX:
        folder = directory / matrix
        if download(run, f'xcoll-pr-physics-{suffix}-{matrix}', folder, optional=True):
            path = folder / 'physics-pr.xml'
            if path.is_file():
                rows.update(xml_rows(path, matrix_cases(matrix)))
    for name in expected_cases() - rows.keys():
        rows[name] = dict(name=name, status='NOT_RUN', seconds=None, details='No completed test report', checks=[])
    conclusion = run['conclusion']
    if conclusion == 'success' and any(r['status'] != 'PASSED' for r in rows.values()):
        conclusion = 'failure'
    return dict(context, datetime=run['updated_at'].removesuffix('Z'), conclusion=conclusion,
                reference_sha=reference_sha(context['base_sha']), rows=[rows[k] for k in sorted(rows)])


def benchmark_record(context, run, directory):
    suffix = f"{run['id']}-{run['run_attempt']}"
    download(run, f'xcoll-external-pr-comparison-{suffix}', directory)
    path = directory / 'pr-comparison.json'
    require(path.stat().st_size <= 3_000_000, 'Oversized benchmark comparison')
    clean = validate_benchmarks(json.loads(path.read_text()), run_id=run['id'], attempt=run['run_attempt'])
    for key in ('pr', 'branch', 'base_sha', 'head_sha', 'candidate_sha'):
        require(clean[key] == context[key], f'External comparison {key} differs from resolved context')
    return dict(clean, profile=context['kind'], kind='benchmark')


def publish():
    run = event()['workflow_run']
    require(run['name'] in WORKFLOWS and run['repository']['full_name'] == REPO,
            'Unexpected source workflow')
    expected_path = {
        'physics': '.github/workflows/xcoll-physics-pr.yml',
        'external': '.github/workflows/xcoll-performance-external-pr.yml',
    }[WORKFLOWS[run['name']]]
    require(run['path'].split('@', 1)[0] == expected_path, 'Unexpected source workflow path')
    suffix = f"{number(run['id'])}-{number(run['run_attempt'])}"
    with tempfile.TemporaryDirectory(dir=os.environ['RUNNER_TEMP']) as tmp:
        directory = Path(tmp)
        # Draft/skipped runs have no resolver artifact, and must not be published.
        if not download(run, f'xcoll-pr-context-{suffix}', directory / 'context', optional=True):
            print('No resolved PR context; nothing to publish')
            return
        ctx = load_context(directory / 'context' / 'pr-context.json', run)
        if ctx['kind'] == 'physics':
            record = physics_record(ctx, run, directory)
        elif run['conclusion'] == 'success':
            record = benchmark_record(ctx, run, directory / 'benchmark')
        else:
            print('External benchmark did not complete successfully; see source workflow')
            return
        dest = Path(os.environ['RUNNER_TEMP']) / 'pr-extra-record.json'
        dest.write_text(json.dumps(record, indent=2, allow_nan=False) + '\n')
        with open(os.environ['GITHUB_OUTPUT'], 'a') as f:
            print(f"pr={ctx['pr']}\nrecord={dest}\nkind={ctx['kind']}", file=f)
        # A manual dispatch has its own check suite on main. Attach an explicit
        # status to the tested merge only if it is still the current PR merge.
        if ctx['kind'] == 'physics' and run['event'] == 'workflow_dispatch':
            try:
                pr = api(f"pulls/{ctx['pr']}")
                merge_ref = sha(api(f"git/ref/pull/{ctx['pr']}/merge")['object']['sha'])
                if (pr['state'] == 'open' and pr['head']['sha'] == ctx['head_sha']
                        and main_sha() == ctx['base_sha']
                        and merge_ref == ctx['candidate_sha']):
                    state = 'success' if record['conclusion'] == 'success' else 'failure'
                    gh('api', '--method', 'POST', f"repos/{REPO}/statuses/{ctx['candidate_sha']}",
                       '-f', f'state={state}', '-f', 'context=Xcoll physics (manual)',
                       '-f', f"target_url=https://github.com/{REPO}/actions/runs/{run['id']}",
                       '-f', 'description=Everest reference validation (manual)')
            except subprocess.CalledProcessError:
                print('::warning::Could not verify/publish manual commit status; retaining dashboard record')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command', choices=('resolve', 'publish'))
    args = parser.parse_args()
    {'resolve': resolve, 'publish': publish}[args.command]()


if __name__ == '__main__':
    main()
