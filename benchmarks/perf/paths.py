"""Resolve project-owned files relative to the perf package, not the shell cwd.

Put this package in ``<project>/benchmarks/perf`` alongside metrics.toml and
pytest.ini. No project name or absolute path is embedded in the tools.
"""
from __future__ import annotations
from pathlib import Path

BENCHMARK_DIR = Path(__file__).resolve().parent.parent


def config_file(path: str | Path | None = None) -> Path:
    return Path(path).expanduser().resolve() if path is not None else BENCHMARK_DIR / "metrics.toml"


def results_dir(path: str | Path | None = None) -> Path:
    return Path(path).expanduser().resolve() if path is not None else BENCHMARK_DIR / "results"
