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


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--base', required=True, type=Path)
    p.add_argument('--candidate', required=True, type=Path)
    p.add_argument('--scheduled', type=Path)
    p.add_argument('--base-sha', required=True)
    p.add_argument('--candidate-sha', required=True)
    p.add_argument('--output', required=True, type=Path)
    a = p.parse_args(argv)
    baseline = _read_baseline(a.scheduled)
    if baseline is not None and baseline.get('profile') != 'full':
        raise ValueError('Expected scheduled full CPU baseline')
    rows = compare(_load(a.base), _load(a.candidate), baseline)
    a.output.write_text(markdown(rows, a.base_sha, a.candidate_sha,
                                 scheduled_available=baseline is not None), encoding='utf-8')


if __name__ == '__main__':
    main()
