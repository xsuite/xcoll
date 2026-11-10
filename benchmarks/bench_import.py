"""Cold-process import latency (the operating system's file cache may be warm)."""

import os
import subprocess
from statistics import median
import sys

import pytest

# Import the package in a fresh interpreter for every measurement. It does not
# matter whether pytest itself has imported xcoll during collection.
_COMMAND = (
    "import time; t = time.perf_counter_ns(); "
    "import xcoll; print((time.perf_counter_ns() - t) / 1e9)"
)


@pytest.mark.quick
@pytest.mark.benchmark(group="Xcoll | startup")
def bench_import_xcoll_fresh_process(benchmark):
    """Python start-up + xcoll import; import-only duration is recorded separately."""
    child_timings = []

    def import_once():
        env = os.environ.copy()
        env["PYTHONDONTWRITEBYTECODE"] = "1"
        result = subprocess.run(
            [sys.executable, "-B", "-c", _COMMAND],
            check=True, capture_output=True, text=True, env=env,
        )
        seconds = float(result.stdout.strip().splitlines()[-1])
        child_timings.append(seconds)
        return seconds

    import_once()  # warm file cache (but still a new interpreter)
    child_timings.clear()
    seconds = benchmark.pedantic(import_once, rounds=7, iterations=1)
    assert seconds > 0
    benchmark.extra_info.update(
        protocol=2, physics_case="fresh-process-import", fresh_interpreter=True,
        child_import_median_s=median(child_timings),
        timed="spawn Python interpreter, then import xcoll",
    )
