"""Convenience command: native pytest-benchmark timing table + physics check.

    python -m perf.compare OLD.json NEW.json

The original `pytest-benchmark compare` executable is not modified. Use
`--no-timings` to run physics comparison without pytest-benchmark installed.
"""

import argparse
from pathlib import Path
import subprocess
import sys

from .compare_physics import main as physics_main


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("baseline", type=Path)
    parser.add_argument("candidate", type=Path)
    parser.add_argument("--no-timings", action="store_true")
    parser.add_argument("--allow-unverified", action="store_true")
    args = parser.parse_args(argv)
    if not args.no_timings:
        try:
            completed = subprocess.run(
                ["pytest-benchmark", "compare", str(args.baseline), str(args.candidate)],
                check=False,
            )
        except FileNotFoundError:
            print("pytest-benchmark not installed; use --no-timings for physics-only comparison", file=sys.stderr)
            return 2
        if completed.returncode:
            return completed.returncode
    print("\n" + "=" * 70, flush=True)
    params = [str(args.baseline), str(args.candidate)]
    if args.allow_unverified:
        params.append("--allow-unverified")
    return physics_main(params)


if __name__ == "__main__":
    raise SystemExit(main())
