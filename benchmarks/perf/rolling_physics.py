"""Compare benchmarks to the last daily physics baseline; only daily runs update it.

Only compact physical observables currently exist in the benchmark JSON.
Checks therefore screen scalar observables using metrics.toml, not particle
arrays, phase-space distributions, or full scientific reproducibility.

The previous baseline is a replaceable, compressed snapshot in the history
branch. Permanent run records retain *only* comparison outcomes and configured
scalar observables. The upload-artifact archive still contains raw results.
"""
from __future__ import annotations

import argparse
import gzip
import json
from pathlib import Path

from .compare_physics import compare_entry
from .metrics import load_config, observables, workload

MAX_BASELINE_BYTES = 1_000_000


def _name(entry):
    return entry.get('fullname') or entry.get('name') or 'unnamed'


def _read_baseline(path):
    if path is None or not Path(path).is_file():
        return None
    with gzip.open(path, 'rt', encoding='utf-8') as stream:
        obj = json.load(stream)
    if not isinstance(obj, dict) or not isinstance(obj.get('benchmarks'), list):
        raise ValueError('Invalid baseline file')
    return obj


def _snapshot(current, config):
    """Save only validated scalar physics and workload fields for next run."""
    names = config['fields']
    specs = config['observables']
    compact = []
    for item in current['benchmarks']:
        extra = item.get('extra_info') or {}
        fields = dict(extra.get('workload') or {key: extra.get(key) for key in names})
        obs = observables(item, config)
        compact.append({
            'name': item.get('name'), 'fullname': item.get('fullname'),
            'extra_info': {'workload': fields, 'observables': obs, 'observables_pre_normalized': True},
        })
    return {'testbed': current.get('testbed'),
            'profile': (current.get('ci') or {}).get('profile'),
            'benchmarks': compact, 'schema': 1}


def annotate(current, previous=None, config=None, *, daily=False):
    config = config or load_config()
    old = {_name(e): e for e in previous['benchmarks']} if previous else {}
    statuses = {}
    for item in current.get('benchmarks') or []:
        name = _name(item)
        if previous is None:
            status = {'status': 'BASELINE' if daily else 'NOT_CHECKED',
                      'summary': ('First scheduled measurement; no previous baseline' if daily
                                  else 'No scheduled physics baseline available'), 'details': []}
        elif name not in old:
            status = {'status': 'BASELINE' if daily else 'NOT_CHECKED',
                      'summary': 'No scheduled reference for this benchmark', 'details': []}
        else:
            comparison = compare_entry(old[name], item, config)
            status = {'status': comparison.status, 'summary': comparison.summary,
                      'details': list(comparison.differences)}
        item.setdefault('extra_info', {})['physics_check'] = status
        statuses[name] = status
    return statuses


def compact_history(current, config):
    """Remove non-physics bulk and retain the numbers needed for future plots."""
    for item in current.get('benchmarks') or []:
        extra = item.get('extra_info') or {}
        supported = {key for key in config['fields']}
        supported.update(spec.get('source', name)
                         for name, spec in config['observables'].items())
        supported.update(spec['denominator'] for spec in config['observables'].values()
                         if 'denominator' in spec)
        supported.update({'workload', 'observables', 'physics_check', 'protocol'})
        item['extra_info'] = {k: v for k, v in extra.items() if k in supported}
        item.pop('data', None)
        if isinstance(item.get('stats'), dict):
            item['stats'].pop('data', None)
    return current


def process(record, *, previous=None, baseline_out=None, history_out=None, config=None):
    config = config or load_config()
    current = json.loads(Path(record).read_text(encoding='utf-8'))
    ci = current.get('ci') or {}
    daily = (ci.get('trigger') == 'schedule' and ci.get('source_branch') == 'main'
             and ci.get('profile') in ('full', 'fluka', 'geant4'))
    reference_profile = 'full' if ci.get('profile') == 'quick' else ci.get('profile')
    old = _read_baseline(previous)
    # Physics can be compared across machines; only a scheduled baseline writer
    # must have the same testbed. Workload compatibility is checked per entry.
    if old and (old.get('profile') != reference_profile or
                (daily and old.get('testbed') != current.get('testbed'))):
        old = None
    statuses = annotate(current, old, config, daily=daily)
    if daily and baseline_out is not None:
        new = _snapshot(current, config)
        payload = json.dumps(new, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()
        compressed = gzip.compress(payload, compresslevel=9, mtime=0)
        if len(compressed) > MAX_BASELINE_BYTES:
            raise ValueError(f'Physics baseline is {len(compressed)} bytes (limit {MAX_BASELINE_BYTES})')
        path = Path(baseline_out)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(compressed)
    if history_out is not None:
        path = Path(history_out)
        path.parent.mkdir(parents=True, exist_ok=True)
        compact_history(current, config)
        path.write_text(json.dumps(current, indent=2, ensure_ascii=True) + '\n', encoding='utf-8')
    return statuses


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('record', type=Path)
    p.add_argument('--previous', type=Path)
    p.add_argument('--baseline-out', type=Path)
    p.add_argument('--history-out', type=Path)
    p.add_argument('--config', type=Path)
    a = p.parse_args(argv)
    statuses = process(a.record, previous=a.previous, baseline_out=a.baseline_out,
                       history_out=a.history_out, config=load_config(a.config))
    for name, status in statuses.items():
        print(f'{status["status"]:12} {name}: {status["summary"]}')
        for diff in status['details']:
            print('              ' + diff)
    # Physics regressions should be red-flagged, not suppress measurements.
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
