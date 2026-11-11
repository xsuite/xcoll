"""Build a self-contained, chronological benchmark dashboard from JSON results.

From the benchmarks directory:
    python -m perf.report results -o results/dashboard.html

The input may be a JSON file, a directory, or a list of paths. Paths are
resolved from the invoking directory; project defaults live in perf.paths.

Normal/local reports retain one column per JSON file. ``--online`` merges
scheduled full/FLUKA/Geant4 runs from the same physical CERN VM and UTC day,
keeping each contributing revision/timestamp in its entry metadata.
"""
from __future__ import annotations

import argparse
import json
import math
from collections import defaultdict
from pathlib import Path

from .metrics import load_config, limits, observables, workload
from .paths import results_dir


ONLINE_PROFILES = frozenset({'full', 'fluka', 'geant4'})


def _files(paths):
    files = []
    for item in paths:
        item = Path(item)
        if item.is_dir():
            files.extend(p for p in item.rglob('*.json') if p.is_file())
        elif item.is_file():
            files.append(item)
        else:
            raise FileNotFoundError(f'No results found at {item}')
    return sorted(set(files))


def _hardware(machine, testbed, physical_host, explicit_name=None):
    cpu = machine.get('cpu') or {}
    processor = (machine.get('processor') or cpu.get('brand_raw') or
                 cpu.get('brand') or 'unknown CPU')
    # Prefer the configured OS computer name captured when the benchmark ran.
    if isinstance(explicit_name, str) and explicit_name.strip():
        return explicit_name.strip()
    if physical_host:
        return physical_host
    if testbed.startswith('cern-openstack-xcoll-'):
        return 'cern-openstack-xcoll-vm-v1'
    # Exclude the volatile network host name: a laptop may move networks.
    return f"{machine.get('system', '?')} | {processor}"


def read_runs(paths, config=None, *, online=False):
    config = config or load_config()
    runs = []
    for path in _files(paths):
        data = json.loads(path.read_text(encoding='utf-8'))
        if not isinstance(data, dict) or not isinstance(data.get('benchmarks'), list):
            # An arbitrary JSON (e.g., another tool's manifest) is not a run.
            continue
        machine = data.get('machine_info') or {}
        ci = data.get('ci') or {}
        commit = data.get('commit_info') or {}
        testbed = data.get('testbed') or 'unlabelled-local'
        physical_host = (ci.get('physical_host') or
                         ('cern-openstack-xcoll-vm-v1' if testbed.startswith('cern-openstack-xcoll-') else ''))
        machine_key = _hardware(machine, testbed, physical_host, data.get('machine_name'))
        sha = ci.get('source_sha') or commit.get('id') or 'unknown'
        branch = ci.get('source_branch') or commit.get('branch') or 'unknown'
        stamp = str(data.get('datetime') or '')
        profile = ci.get('profile') or ''
        # Xcoll release/version is stored by pytest-benchmark's metadata hook.
        release = str(ci.get('release') or (data.get('xsuite_versions') or {}).get('xcoll') or 'unknown')
        entries = []
        for b in data['benchmarks']:
            stats = b.get('stats') or {}
            try:
                median = float(stats['median'])
            except (KeyError, TypeError, ValueError):
                continue
            if not math.isfinite(median) or median <= 0:
                continue
            extra = b.get('extra_info') or {}
            obs = observables(b, config)
            entries.append({
                'name': b.get('fullname') or b.get('name') or 'unnamed',
                'group': b.get('group') or 'Ungrouped',
                'median': median, 'stddev': stats.get('stddev'),
                'rounds': stats.get('rounds'),
                'signature': workload(b, config),
                'observables': obs, 'tolerances': limits(config, obs),
                # Only compact physics status, never detailed particle arrays.
                'physics': extra.get('physics_check'),
                'extra': {k: v for k, v in extra.items() if k not in {'physics', 'raw_physics', 'particle_data'}},
                'sha': sha, 'profile': profile, 'datetime': stamp,
            })
        if entries:
            runs.append({
                'id': f"{ci['run_id']}-{ci.get('run_attempt', '')}" if ci.get('run_id') else path.name,
                'datetime': stamp, 'sha': sha, 'branch': branch,
                'release': release, 'machine': machine_key, 'testbed': testbed,
                'host': '', 'physical_host': physical_host,
                'trigger': ci.get('trigger') or '',
                'versions': data.get('xsuite_versions') or {},
                'profile': profile, 'entries': entries,
                'components': [{'profile': profile, 'sha': sha, 'datetime': stamp,
                                'release': release}],
            })
    runs.sort(key=lambda r: (r['datetime'], r['id']))
    return merge_daily(runs) if online else runs


def merge_daily(runs):
    """Group only scheduled CERN profiles, never ordinary local/manual results."""
    buckets = defaultdict(list)
    singles = []
    for run in runs:
        if (run['trigger'] == 'schedule' and run['profile'] in ONLINE_PROFILES and
                run['physical_host'] and run['branch'] == 'main' and
                len(run['datetime']) >= 10):
            key = (run['datetime'][:10], run['physical_host'], run['branch'])
            buckets[key].append(run)
        else:
            singles.append(run)
    for (date, host, branch), matches in buckets.items():
        # One result per profile per day; on retries keep the latest successful.
        most_recent = {r['profile']: r for r in matches}
        chosen = sorted(most_recent.values(), key=lambda r: (r['datetime'], r['profile']))
        entries = []
        components = []
        for r in chosen:
            entries.extend(r['entries'])
            components.extend(r['components'])
        singles.append({
            'id': f'daily:{date}:{host}:{branch}', 'datetime': max(r['datetime'] for r in chosen),
            'sha': chosen[-1]['sha'] if len({r['sha'] for r in chosen}) == 1 else 'mixed',
            'branch': branch, 'release': (chosen[-1]['release']
                if len({r['release'] for r in chosen}) == 1 else 'mixed'),
            'machine': chosen[-1]['machine'], 'testbed': 'cern-daily',
            'host': '', 'physical_host': host, 'trigger': 'schedule',
            'versions': chosen[-1]['versions'], 'profile': 'daily',
            'components': components, 'entries': entries,
        })
    return sorted(singles, key=lambda r: (r['datetime'], r['id']))


def create_dashboard(paths, output, config=None, *, online=False):
    runs = read_runs(paths, config=config, online=online)
    if not runs:
        raise ValueError('No valid benchmark results found')
    template = (Path(__file__).resolve().parent / 'templates' / 'dashboard.html').read_text(encoding='utf-8')
    # Only insert JSON in JS source, never raw HTML. Escaping '<' prevents
    # values such as '</script>' in benchmark extra_info closing the script.
    serialized = json.dumps(runs, ensure_ascii=True, separators=(',', ':')).replace('<', '\\u003c')
    html = template.replace('/*__BENCHMARK_DATA__*/null', serialized)
    if html == template:
        raise RuntimeError('Dashboard template missing its data placeholder')
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(html, encoding='utf-8')
    return output


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('paths', nargs='*', type=Path)
    p.add_argument('-o', '--output', type=Path, default=results_dir() / 'dashboard.html')
    p.add_argument('--config', type=Path)
    p.add_argument('--online', action='store_true', help='Group scheduled CERN CPU/FLUKA/Geant4 per UTC day')
    a = p.parse_args(argv)
    try:
        output = create_dashboard(a.paths or [results_dir()], a.output, load_config(a.config), online=a.online)
    except (FileNotFoundError, ValueError) as err:
        p.error(str(err))
    print(f'Dashboard: {output}')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
