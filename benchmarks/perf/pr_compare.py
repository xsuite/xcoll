"""Read-only, paired PR comparison (base and candidate on one runner).

Physics is screened against both the main measurement and the latest daily
reference. Timing is compared only between the paired local measurements.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from .compare_physics import compare_entry
from .metrics import load_config, workload
from .rolling_physics import _read_baseline


def _load(directory):
    paths = sorted(Path(directory).rglob('*.json'))
    if len(paths) != 1:
        raise ValueError(f'Expected exactly one benchmark JSON in {directory}, got {len(paths)}')
    data = json.loads(paths[0].read_text(encoding='utf-8'))
    if not isinstance(data.get('benchmarks'), list) or not data['benchmarks']:
        raise ValueError(f'No measurements in {paths[0]}')
    return data


def _entries(data):
    return {e.get('fullname') or e['name']: e for e in data['benchmarks']}


def _time(entry):
    try:
        value = float(entry['stats']['median'])
    except (KeyError, TypeError, ValueError):
        return None
    return value if math.isfinite(value) and value > 0 else None


def _escape(text):
    return str(text).replace('&', '&amp;').replace('<', '&lt;').replace('>', '&gt;').replace('|', '&#124;').replace('\n', ' ')


def compare(base, candidate, scheduled=None, *, config=None):
    """Return rows with timing changes and two independent physics checks."""
    config = config or load_config()
    left, right = _entries(base), _entries(candidate)
    previous = _entries(scheduled) if scheduled else {}
    rows = []
    for name in sorted(left.keys() | right.keys()):
        a, b = left.get(name), right.get(name)
        at, bt = _time(a), _time(b)
        comparable = a is not None and b is not None and workload(a, config) == workload(b, config)
        delta = 100 * (bt / at - 1) if comparable and at is not None and bt is not None else None
        pair = compare_entry(a, b, config) if a is not None and b is not None else None
        daily = compare_entry(previous[name], b, config) if b is not None and name in previous else None
        rows.append({'name': name, 'base': at, 'candidate': bt, 'delta': delta,
                     'pair': pair, 'scheduled': daily,
                     'timing_status': 'Comparable' if comparable else ('New in PR' if a is None else 'Missing in PR' if b is None else 'Incompatible workload')})
    return rows


def markdown(rows, base_sha, candidate_sha, *, scheduled_available):
    def ms(value):
        return f'{value * 1000:.4g}' if value is not None else '—'
    def status(result):
        return result.status if result is not None else 'NOT_CHECKED'
    lines = ['## Paired Xcoll CPU benchmarks', '',
             f'Main: `{base_sha}` · PR merge: `{candidate_sha}`', '',
             'Both revisions ran sequentially on the same GitHub-hosted runner, with separate environments and freshly compiled kernels. Physics changes are informational; missing or incompatible measurements are not agreement.', '',
             'Scheduled CPU physics reference: ' + ('available' if scheduled_available else 'not yet available'), '',
             '| Benchmark | Main (ms) | PR (ms) | Runtime change | Physics vs main | Physics vs scheduled |',
             '|---|---:|---:|---:|---|---|']
    for row in rows:
        delta = row['delta']
        change = f'{delta:+.2f}%' if delta is not None else row['timing_status']
        lines.append(f"| {_escape(row['name'])} | {ms(row['base'])} | {ms(row['candidate'])} | {change} | {status(row['pair'])} | {status(row['scheduled'])} |")
    lines += ['', '### Physics details', '']
    for row in rows:
        for key, label in [('pair', 'vs main'), ('scheduled', 'vs scheduled')]:
            result = row[key]
            if result is not None and result.status != 'OK':
                lines.append(f"- {_escape(row['name'])} ({label}): **{result.status}** — {_escape(result.summary)}")
                lines.extend('  - ' + _escape(d) for d in result.differences)
    return '\n'.join(lines) + '\n'


def comparison_record(rows, base, candidate, *, pr, branch, head_sha,
                      base_sha, candidate_sha, run_id, attempt,
                      scheduled_reference=''):
    """Retain exactly the paired comparison performed in this workflow.

    The publisher only validates and archives this document; it must not
    recompute physics against what might be a newer daily baseline.
    """
    from .metrics import load_config, observables, workload

    cfg = load_config()
    a = _entries(base)
    b = _entries(candidate)

    def physics(result):
        if result is None:
            return None
        return {'status': result.status, 'summary': result.summary,
                'details': list(result.differences)}

    def deviation(entry):
        return entry.get('stats', {}).get('stddev') if entry else None

    output = []
    for row in rows:
        name = row['name']
        left, right = a.get(name), b.get(name)
        output.append({
            'name': name, 'base': row['base'], 'candidate': row['candidate'],
            'delta': row['delta'], 'timing_status': row['timing_status'],
            'base_stddev': deviation(left), 'candidate_stddev': deviation(right),
            'base_workload': workload(left, cfg) if left else None,
            'candidate_workload': workload(right, cfg) if right else None,
            'base_observables': observables(left, cfg) if left else None,
            'candidate_observables': observables(right, cfg) if right else None,
            'pair': physics(row['pair']), 'scheduled': physics(row['scheduled']),
        })
    return {
        'schema': 1, 'pr': pr, 'branch': branch, 'run_id': str(run_id),
        'attempt': attempt, 'datetime': candidate.get('datetime', ''),
        'base_sha': base_sha, 'candidate_sha': candidate_sha,
        'head_sha': head_sha, 'machine': candidate.get('machine_name') or
        (candidate.get('machine_info') or {}).get('processor') or 'GitHub hosted',
        'scheduled_reference': scheduled_reference,
        'base_versions': base.get('xsuite_versions') or {},
        'candidate_versions': candidate.get('xsuite_versions') or {},
        'base_compiler': base.get('compiler') or {},
        'candidate_compiler': candidate.get('compiler') or {},
        'rows': output,
    }


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--base', required=True, type=Path)
    p.add_argument('--candidate', required=True, type=Path)
    p.add_argument('--scheduled', type=Path)
    p.add_argument('--base-sha', required=True)
    p.add_argument('--candidate-sha', required=True)
    p.add_argument('--output', required=True, type=Path)
    p.add_argument('--json-output', type=Path)
    p.add_argument('--pr-number', type=int)
    p.add_argument('--branch')
    p.add_argument('--head-sha')
    p.add_argument('--run-id')
    p.add_argument('--attempt', type=int)
    a = p.parse_args(argv)
    baseline = _read_baseline(a.scheduled)
    if baseline is not None and baseline.get('profile') != 'full':
        raise ValueError('Expected scheduled full CPU baseline')
    left, right = _load(a.base), _load(a.candidate)
    rows = compare(left, right, baseline)
    a.output.write_text(markdown(rows, a.base_sha, a.candidate_sha,
                                 scheduled_available=baseline is not None), encoding='utf-8')
    if a.json_output is not None:
        from hashlib import sha256
        if not all((a.pr_number, a.branch, a.head_sha, a.run_id, a.attempt)):
            p.error('JSON publication needs the PR and run metadata')
        digest = (sha256(a.scheduled.read_bytes()).hexdigest() if a.scheduled
                  and a.scheduled.is_file() else '')
        document = comparison_record(
            rows, left, right, pr=a.pr_number, branch=a.branch,
            head_sha=a.head_sha, base_sha=a.base_sha,
            candidate_sha=a.candidate_sha, run_id=a.run_id,
            attempt=a.attempt, scheduled_reference=digest,
        )
        a.json_output.write_text(json.dumps(document, indent=2, allow_nan=False)
                                 + '\n', encoding='utf-8')


if __name__ == '__main__':
    main()
