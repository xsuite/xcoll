"""Optional metadata only: pytest-benchmark remains the entire benchmark runner."""

import os
from importlib import metadata
import pytest

try:
    # When pytest imports conftest as benchmarks.conftest (repository root).
    from .perf.results import archive_results
    from .perf.paths import results_dir
    from .perf.dependencies import external_versions
    from .perf.machine import machine_name
except ImportError:
    # Also support executing `pytest` from inside benchmarks/.
    from perf.results import archive_results
    from perf.paths import results_dir
    from perf.dependencies import external_versions
    from perf.machine import machine_name


def pytest_configure(config):
    # pytest.ini cannot itself reject -n. Check controller configuration
    # before collection and before any expensive kernels are compiled.
    workers = config.getoption('numprocesses', default=None)
    tx = config.getoption('tx', default=None)
    if workers not in (None, 0, '0') or tx:
        raise pytest.UsageError(
            'Performance benchmarks must run serially; do not use pytest-xdist '
            '(-n/--numprocesses or --tx). Use -n 0 to explicitly disable xdist.'
        )


def pytest_sessionstart(session):
    # pytest-benchmark's file-store lookup only matches four-digit counters.
    # Rotate long-lived local result directories before hitting 10000.
    archive_results(results_dir(), threshold=9000)


@pytest.hookimpl(optionalhook=True)
def pytest_benchmark_update_json(config, benchmarks, output_json):
    """Version provenance without replacing pytest-benchmark's own JSON."""
    versions = {}
    for package in ('xcoll', 'xtrack', 'xobjects', 'xpart', 'xdeps', 'numpy'):
        try:
            versions[package] = metadata.version(package)
        except metadata.PackageNotFoundError:
            versions[package] = None
    external, sources = external_versions()
    # Keep this inside xsuite_versions so the existing dashboard displays it.
    versions.update(external)
    output_json['xsuite_versions'] = versions
    output_json['external_version_sources'] = sources
    output_json['machine_name'] = machine_name()
    output_json['testbed'] = os.environ.get('BENCHMARK_TESTBED', 'unlabelled-local')
    output_json['benchmark_protocol'] = 2
