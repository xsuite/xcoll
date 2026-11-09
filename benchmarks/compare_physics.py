"""Compare Xsuite benchmark physical observables in pytest-benchmark JSON files.

Usage: python benchmarks/compare_physics.py OLD.json NEW.json

This is a *companion* to `pytest-benchmark compare`, not a monkey-patch.
The upstream CLI only compares timings. We intentionally leave it unchanged.

A statistical distribution test still belongs in the functional tests; these
thresholds are screening tolerances rather than proof of physics equivalence.
"""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass
from pathlib import Path

# Maximum allowed ABSOLUTE change in a fraction, not a percentage increase.
# A 1M-particle measurement can later use much tighter tolerances, but the
# LHC quick suite has only 3,000 particles, and its physics is stochastic.
FRACTION_TOLERANCES = {
    "lost_fraction": 0.01,
    "survivors_with_kick_gt_0.1urad": 0.02,
    "channeling_proxy_fraction_of_survivors": 0.02,
    "reference_collimator_losses_per_particle": 0.02,
    "reference_aperture_losses_per_particle": 0.02,
}

# Fields defining the workload. These must match for comparisons to be valid.
WORKLOAD_FIELDS = (
    "protocol", "particles", "turns", "crystals", "beam", "plane",
    "interpolation_m", "physics_case", "geometry", "tilt_radians",
    "sequence_sha256", "collimator_database_sha256",
)


@dataclass(frozen=True)
class Result:
    name: str
    status: str
    summary: str
    differences: tuple[str, ...] = ()


def _bench_name(entry):
    # pytest-benchmark uses both forms depending on invocation and version.
    # The short name is better for stable cross-checkout comparisons.
    return entry.get("name") or entry.get("fullname") or "unnamed"


def read_json(path: Path):
    data = json.loads(Path(path).read_text(encoding="utf-8"))
    if not isinstance(data, dict) or not isinstance(data.get("benchmarks"), list):
        raise ValueError(f"Not a pytest-benchmark result: {path}")
    return data


def _observables(extra):
    """Translate v1's stored `extra_info` into comparable fractions."""
    observed = {}
    for name in ("lost_fraction", "survivors_with_kick_gt_0.1urad",
                 "channeling_proxy_fraction_of_survivors"):
        if extra.get(name) is not None:
            observed[name] = float(extra[name])
    n = extra.get("particles")
    if n is not None and float(n) > 0:
        for name in ("reference_collimator_losses", "reference_aperture_losses"):
            if extra.get(name) is not None:
                observed[name + "_per_particle"] = float(extra[name]) / float(n)
    # Future benchmarks can provide their own explicit numerical observables,
    # as long as a tolerance is registered above for each quantity.
    for name, value in (extra.get("physics_observables") or {}).items():
        if name in FRACTION_TOLERANCES:
            observed[name] = float(value)
    return observed


def compare_entry(old, new) -> Result:
    name = _bench_name(new)
    before = old.get("extra_info") or {}
    after = new.get("extra_info") or {}
    different_workloads = [key for key in WORKLOAD_FIELDS
                           if before.get(key) != after.get(key)]
    if different_workloads:
        return Result(name, "INCOMPATIBLE", "Workload changed: " + ", ".join(different_workloads))
    a = _observables(before)
    b = _observables(after)
    if not a and not b:
        return Result(name, "UNVERIFIED", "Neither run records recognised physics observables")
    missing = sorted(set(a) ^ set(b))
    if missing:
        return Result(name, "UNVERIFIED", "Observables missing in one run: " + ", ".join(missing))
    bad = []
    for metric in sorted(a):
        tolerance = FRACTION_TOLERANCES[metric]
        x, y = a[metric], b[metric]
        if not math.isfinite(x) or not math.isfinite(y):
            bad.append(f"{metric}: non-finite measurement")
        elif abs(y - x) > tolerance:
            bad.append(f"{metric}: {x:.5g} -> {y:.5g} (Δ={y-x:+.5g}; limit ±{tolerance:g})")
    if bad:
        return Result(name, "CHANGED", "One or more observables exceeded tolerances", tuple(bad))
    return Result(name, "OK", f"{len(a)} physical observable(s) within tolerances")


def compare_documents(old, new):
    first = {_bench_name(b): b for b in old["benchmarks"]}
    second = {_bench_name(b): b for b in new["benchmarks"]}
    results = []
    for name in sorted(set(first) | set(second)):
        if name not in first:
            results.append(Result(name, "MISSING", "Benchmark not present in baseline"))
        elif name not in second:
            results.append(Result(name, "MISSING", "Benchmark not present in candidate"))
        else:
            results.append(compare_entry(first[name], second[name]))
    return results


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("baseline", type=Path)
    parser.add_argument("candidate", type=Path)
    parser.add_argument("--allow-unverified", action="store_true",
                        help="Do not fail on missing benchmarks/observables (exploratory use only)")
    args = parser.parse_args(argv)
    old, new = read_json(args.baseline), read_json(args.candidate)
    results = compare_documents(old, new)
    if not results:
        print("No benchmarks found!")
        return 2
    print("Physics comparison (stored benchmark.extra_info):")
    for r in results:
        print(f"  {r.status:12} {r.name}: {r.summary}")
        for detail in r.differences:
            print(f"               {detail}")
    counts = {status: sum(r.status == status for r in results)
              for status in ("OK", "CHANGED", "MISSING", "INCOMPATIBLE", "UNVERIFIED")}
    print("Summary:", ", ".join(f"{key}={v}" for key, v in counts.items()))
    print("These are screening tolerances, not a replacement for scientific validation tests.")
    if counts["CHANGED"] or counts["INCOMPATIBLE"]:
        return 1
    if not args.allow_unverified and (counts["MISSING"] or counts["UNVERIFIED"]):
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
