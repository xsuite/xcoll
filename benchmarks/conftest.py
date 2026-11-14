"""Optional metadata only: pytest-benchmark remains the entire benchmark runner."""

import os
import pytest
import subprocess
from pathlib import Path
from importlib import import_module, metadata

try:
    # When pytest imports conftest as benchmarks.conftest (repository root).
    from .perf.results import archive_results
    from .perf.paths import results_dir
    from .perf.dependencies import compiler_metadata, external_versions
    from .perf.machine import machine_name
    from .perf.provenance import git_state, check_state
except ImportError:
    # Also support executing `pytest` from inside benchmarks/.
    from perf.results import archive_results
    from perf.paths import results_dir
    from perf.dependencies import compiler_metadata, external_versions
    from perf.machine import machine_name
    from perf.provenance import git_state, check_state


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
    config = session.config

    # Only enforce provenance for sessions that will
    # save benchmark results.
    autosave = config.getoption(
        'benchmark_autosave', default=None
    )

    if (autosave
            or config.getoption('benchmark_save', default=None)
            or config.getoption('benchmark_json', default=None)):

        root = Path(__file__).resolve().parents[1]

        try:
            start = git_state(root)
        except (
            OSError,
            subprocess.CalledProcessError,
            RuntimeError,
        ) as exc:
            raise pytest.UsageError(
                f'Cannot capture benchmark Git provenance: {exc}'
            ) from exc

        # pytest-benchmark computes this tag even earlier,
        # during option setup.
        if autosave and str(autosave).split('_', 1)[0] != start.commit:
            raise pytest.UsageError(
                'Git HEAD changed before benchmark session start; rerun.'
            )

        config._xcoll_benchmark_git_start = (root, start)

    # pytest-benchmark's file-store lookup only matches four-digit counters.
    # Rotate long-lived local result directories before hitting 10000.
    archive_results(results_dir(), threshold=9000)


@pytest.hookimpl(optionalhook=True, trylast=True)
def pytest_benchmark_update_commit_info(config, commit_info):
    """Abort before pytest-benchmark writes any results."""

    snapshot = getattr(
        config, '_xcoll_benchmark_git_start', None
    )
    if snapshot is None:
        return

    root, start = snapshot

    try:
        check_state(
            start,
            git_state(root),
            commit_info.get('id'),
        )
    except (
        ValueError,
        OSError,
        subprocess.CalledProcessError,
        RuntimeError,
    ) as exc:
        pytest.exit(str(exc), returncode=1)


@pytest.hookimpl(optionalhook=True)
def pytest_benchmark_update_json(config, benchmarks, output_json):
    """Version provenance without replacing pytest-benchmark's own JSON."""
    versions = {}
    for package in ('xcoll', 'xtrack', 'xobjects', 'xpart', 'xdeps', 'xfields'):
        try:
            versions[package] = str(import_module(package).__version__)
        except:
            try:
                versions[package] = metadata.version(package)
            except metadata.PackageNotFoundError:
                versions[package] = None
    for package in ('numpy', 'scipy', 'pandas'):
        try:
            versions[package] = metadata.version(package)
        except metadata.PackageNotFoundError:
            versions[package] = None
    external, sources = external_versions()
    # Keep this inside xsuite_versions so the existing dashboard displays it.
    versions.update(external)
    output_json['xsuite_versions'] = versions
    output_json['external_version_sources'] = sources
    output_json['compiler'] = compiler_metadata()
    output_json['machine_name'] = machine_name()
    output_json['testbed'] = os.environ.get('BENCHMARK_TESTBED', 'unlabelled-local')
    output_json['benchmark_protocol'] = 2
