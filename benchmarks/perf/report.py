"""Build a self-contained, chronological benchmark dashboard from JSON results.

From the benchmarks directory:
    python -m perf.report results -o results/dashboard.html

The input may be a JSON file, a directory, or a list of paths. Paths are
resolved from the invoking directory; project defaults live in perf.paths.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from .metrics import load_config, limits, observables, workload
from .paths import results_dir


def _files(paths):
    files = []
    for item in paths:
        item = Path(item)
        if item.is_dir():
            files.extend(p for p in item.rglob('*.json') if p.is_file())
        elif item.is_file():
            files.append(item)
        else:
            raise FileNotFoundError(f'No benchmark results found at {item}')
    return sorted(set(files))


def read_runs(paths, config=None):
    config = config or load_config()
    runs = []
    for path in _files(paths):
        data = json.loads(path.read_text(encoding='utf-8'))
        if not isinstance(data, dict) or not isinstance(data.get('benchmarks'), list):
            # An arbitrary JSON (e.g., another tool's manifest) is not a run.
            continue
        machine = data.get('machine_info') or {}
        cpu = machine.get('cpu') or {}
        processor = machine.get('processor') or cpu.get('brand_raw') or cpu.get('brand') or 'unknown CPU'
        host = machine.get('node') or 'unknown host'
        testbed = data.get('testbed') or 'unlabelled-local'
        # A named benchmark testbed must identify a stable machine/flavour.
        # Do not use the transient GitHub-hosted instance hostname.
        machine_key = f"{testbed} | {machine.get('system', '?')} | {processor}"
        if testbed == 'unlabelled-local':
            machine_key += f' | {host}'
        ci = data.get('ci') or {}
        commit = data.get('commit_info') or {}
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
                'median': median,
                'stddev': stats.get('stddev'),
                'rounds': stats.get('rounds'),
                'signature': workload(b, config),
                'observables': obs,
                'tolerances': limits(config, obs),
                'extra': extra,
            })
        if not entries:
            continue
        runs.append({
            'id': ci.get('run_id', '') + '-' + str(ci.get('run_attempt', '')) if ci.get('run_id') else path.name,
            'datetime': str(data.get('datetime') or ''),
            'sha': ci.get('source_sha') or commit.get('id') or 'unknown',
            'branch': ci.get('source_branch') or commit.get('branch') or 'unknown',
            'machine': machine_key,
            'testbed': testbed,
            'versions': data.get('xsuite_versions') or {},
            'profile': ci.get('profile') or '',
            'entries': entries,
        })
    return sorted(runs, key=lambda r: (r['datetime'], r['id']))


def create_dashboard(paths, output, config=None):
    runs = read_runs(paths, config=config)
    if not runs:
        raise ValueError('No valid benchmark result files found')
    template = (Path(__file__).resolve().parent / 'templates' / 'dashboard.html').read_text(encoding='utf-8')
    # Only insert JSON in JS source, never raw HTML. Escaping '<' prevents
    # values such as '</script>' in benchmark extra_info closing the script.
    serialized = json.dumps(runs, ensure_ascii=True, separators=(',', ':')).replace('<', '\\u003c')
    html = template.replace('/*__BENCHMARK_DATA__*/null', serialized)
    if html == template:
        raise RuntimeError('Dashboard template is missing its data placeholder')
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(html, encoding='utf-8')
    return output


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('paths', nargs='*', type=Path, help='JSON files or directories (default: results/)')
    p.add_argument('-o', '--output', type=Path, default=results_dir() / 'dashboard.html')
    p.add_argument('--config', type=Path, help='Project metrics.toml (default: benchmarks/metrics.toml)')
    args = p.parse_args(argv)
    try:
        output = create_dashboard(args.paths or [results_dir()], args.output, load_config(args.config))
    except (FileNotFoundError, ValueError) as e:
        p.error(str(e))
    print(f'Dashboard: {output}')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
