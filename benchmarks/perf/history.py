"""Make one immutable history record from a pytest-benchmark autosave run.

This script NEVER pushes to GitHub. CI performs Git operations separately,
with a write token limited to the publication job. CI invokes this trusted
script from main, not from the benchmarked feature branch.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def pack(input_dir, output_dir, run_id, attempt, sha, branch, profile, dependencies=None, trigger="", physical_host=""):

    sources = sorted(Path(input_dir).rglob("*.json"))
    if len(sources) != 1:
        raise ValueError(f"Expected exactly one pytest-benchmark autosave JSON in {input_dir}, found {len(sources)}")
    data = json.loads(sources[0].read_text(encoding="utf-8"))
    if not isinstance(data.get("benchmarks"), list) or not data["benchmarks"]:
        raise ValueError("Benchmark result is empty or invalid; refusing publication")
    sha = str(sha)
    if len(sha) != 40 or not all(c in "0123456789abcdefABCDEF" for c in sha):
        raise ValueError("Expected a full 40-character Git commit SHA")
    if not (str(run_id).isdigit() and str(attempt).isdigit()):
        raise ValueError("Run ID/attempt must be numeric")
    data["ci"] = {
        "run_id": str(run_id), "run_attempt": int(attempt),
        "source_sha": sha, "source_branch": str(branch), "profile": str(profile),
        "trigger": str(trigger), "physical_host": str(physical_host),
    }
    if dependencies is not None and Path(dependencies).exists():
        data["ci"]["dependencies"] = Path(dependencies).read_text(encoding="utf-8").splitlines()
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    out = out_dir / f"run-{run_id}-{attempt}.json"
    if out.exists():
        raise FileExistsError(f"Refusing to replace an existing history record: {out}")
    out.write_text(json.dumps(data, ensure_ascii=True, indent=2) + "\n", encoding="utf-8")
    return out


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input", required=True, type=Path)
    p.add_argument("--output", required=True, type=Path)
    p.add_argument("--run-id", required=True)
    p.add_argument("--attempt", required=True)
    p.add_argument("--sha", required=True)
    p.add_argument("--branch", required=True)
    p.add_argument("--profile", choices=["quick", "full", "fluka", "geant4", "gpu"], required=True)
    p.add_argument("--dependencies", type=Path)
    p.add_argument("--trigger", default="")
    p.add_argument("--physical-host", default="")
    a = p.parse_args(argv)
    print(pack(a.input, a.output, a.run_id, a.attempt,
               a.sha, a.branch, a.profile, a.dependencies,
               a.trigger, a.physical_host))


if __name__ == "__main__":
    main()
