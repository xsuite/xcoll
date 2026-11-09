"""Optional metadata only: pytest-benchmark remains the entire benchmark runner."""

import os
from importlib import metadata


def pytest_benchmark_update_json(config, benchmarks, output_json):
    """Version provenance without replacing pytest-benchmark's own JSON."""
    versions = {}
    for package in ("xcoll", "xtrack", "xobjects", "xpart", "xdeps", "numpy"):
        try:
            versions[package] = metadata.version(package)
        except metadata.PackageNotFoundError:
            versions[package] = None
    output_json["xsuite_versions"] = versions
    output_json["testbed"] = os.environ.get("BENCHMARK_TESTBED", "unlabelled-local")
    output_json["benchmark_protocol"] = 2
