"""Project-neutral benchmark metadata and observable specifications.

Each Xsuite package can supply its own metrics.toml. Neither the comparator
nor the HTML dashboard contains project-specific physics field names.
"""

from __future__ import annotations

import math
import tomllib
from pathlib import Path

from .paths import config_file

DEFAULT_CONFIG = config_file()


def load_config(path=None):
    path = Path(path) if path is not None else DEFAULT_CONFIG
    with path.open('rb') as f:
        data = tomllib.load(f)
    fields = data.get('comparison', {}).get('workload_fields', [])
    if not isinstance(fields, list) or not all(isinstance(v, str) for v in fields):
        raise ValueError('comparison.workload_fields must be a list of strings')
    specifications = data.get('observables', {})
    for name, spec in specifications.items():
        if 'absolute_tolerance' not in spec:
            raise ValueError(f'Missing absolute_tolerance for {name}')
        if spec['absolute_tolerance'] < 0 or spec.get('relative_tolerance', 0) < 0 or spec.get('binomial_sigmas', 0) < 0:
            raise ValueError(f'Negative tolerance for {name}')
    return {'fields': tuple(fields), 'observables': specifications}


def workload(entry, config):
    """Explicit nested workload wins; otherwise select the configured flat keys."""
    extra = entry.get('extra_info') or {}
    explicit = extra.get('workload')
    if explicit is not None:
        if not isinstance(explicit, dict):
            raise ValueError('extra_info.workload must be a mapping')
        return explicit
    return {name: extra.get(name) for name in config['fields']}


def observables(entry, config):
    """Return measured values and tolerances, independent of domain names.

    A metric is included only when present in the benchmark output. Missing
    metrics in one of two runs are detected by the comparator.
    """
    extra = entry.get('extra_info') or {}
    embedded = extra.get('observables') or {}
    values = {}
    for name, spec in config['observables'].items():
        source = spec.get('source', name)
        value = embedded.get(name, extra.get(source))
        if value is None:
            continue
        try:
            value = float(value)
            denominator = spec.get('denominator')
            if denominator is not None and not (name in embedded and extra.get("observables_pre_normalized")):
                d = float(extra[denominator])
                if d <= 0 or not math.isfinite(d):
                    raise ValueError(f'Invalid denominator {denominator} for {name}')
                value /= d
        except (ValueError, TypeError, KeyError) as exc:
            raise ValueError(f'Invalid observable {name}: {exc}') from exc
        if not math.isfinite(value):
            raise ValueError(f'Non-finite observable: {name}')
        values[name] = value
    # Explicit observables are also supported for projects without a registered
    # flat-field source. Their tolerances still belong to the project config.
    return values


def limits(config, observed):
    return {key: {'absolute': float(config['observables'][key]['absolute_tolerance']),
                  'relative': float(config['observables'][key].get('relative_tolerance', 0))}
            for key in observed}


def too_different(a, b, spec):
    """Absolute tolerance plus relative tolerance times baseline magnitude."""
    return abs(b - a) > (float(spec['absolute_tolerance']) +
                         float(spec.get('relative_tolerance', 0)) * abs(a))
