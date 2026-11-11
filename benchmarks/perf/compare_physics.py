"""Compare project-defined benchmark observables (not a distributional physics test).

    python -m perf.compare_physics OLD.json NEW.json
    python -m perf.compare_physics OLD.json NEW.json --config benchmarks/metrics.toml
"""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass
from pathlib import Path
from .metrics import load_config, workload, observables, too_different


@dataclass(frozen=True)
class Result:
    name: str
    status: str
    summary: str
    differences: tuple[str, ...] = ()


def read_json(path):
    data = json.loads(Path(path).read_text(encoding='utf-8'))
    if not isinstance(data, dict) or not isinstance(data.get('benchmarks'), list):
        raise ValueError(f'Not a pytest-benchmark result: {path}')
    return data


def _name(entry):
    return entry.get('fullname') or entry.get('name') or 'unnamed'


def compare_entry(old, new, config):
    name = _name(new)
    if workload(old, config) != workload(new, config):
        return Result(name, 'INCOMPATIBLE', 'Benchmark workload parameters differ')
    a, b = observables(old, config), observables(new, config)
    if not a and not b:
        return Result(name, 'UNVERIFIED', 'No configured physical observables saved')
    missing = sorted(set(a) ^ set(b))
    if missing:
        return Result(name, 'UNVERIFIED', 'Observables missing in one run: ' + ', '.join(missing))
    bad = []
    for key in sorted(a):
        spec = config['observables'][key]
        allowed = spec['absolute_tolerance'] + spec.get('relative_tolerance', 0) * abs(a[key])
        if spec.get('binomial_sigmas') is not None:
            # Account for shot noise in small stochastic particle samples.
            n = workload(new, config).get(spec.get('trials_field', 'particles'))
            if n is not None and float(n) > 0 and 0 <= a[key] <= 1:
                sigma = math.sqrt(a[key] * (1 - a[key]) / float(n))
                allowed = max(allowed, float(spec['binomial_sigmas']) * sigma)
        if abs(b[key] - a[key]) > allowed:
            bad.append(f'{key}: {a[key]:.6g} -> {b[key]:.6g}; '
                       f'delta={b[key] - a[key]:+.6g}, allowed={allowed:.6g}')
    if bad:
        return Result(name, 'CHANGED', 'One or more physical metrics outside tolerance', tuple(bad))
    return Result(name, 'OK', f'{len(a)} observables within tolerance')


def compare_documents(old, new, config=None):
    config = config or load_config()
    first = {_name(entry): entry for entry in old['benchmarks']}
    second = {_name(entry): entry for entry in new['benchmarks']}
    result = []
    for name in sorted(set(first) | set(second)):
        if name not in first:
            result.append(Result(name, 'MISSING', 'Not in baseline'))
        elif name not in second:
            result.append(Result(name, 'MISSING', 'Not in candidate'))
        else:
            result.append(compare_entry(first[name], second[name], config))
    return result


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('baseline', type=Path)
    p.add_argument('candidate', type=Path)
    p.add_argument('--config', type=Path)
    p.add_argument('--allow-unverified', action='store_true')
    a = p.parse_args(argv)
    cfg = load_config(a.config)
    results = compare_documents(read_json(a.baseline), read_json(a.candidate), cfg)
    if not results:
        print('No benchmarks found')
        return 2
    for result in results:
        print(f'{result.status:12} {result.name}: {result.summary}')
        for d in result.differences:
            print(' '*14 + d)
    failed = {r.status for r in results}
    print('Checks are a screening tool; full scientific validation belongs in tests.')
    if {'CHANGED', 'INCOMPATIBLE'} & failed:
        return 1
    if {'MISSING', 'UNVERIFIED'} & failed and not a.allow_unverified:
        return 2
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
