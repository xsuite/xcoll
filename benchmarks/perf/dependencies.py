"""Optional external-software provenance for Xsuite benchmarks.

All probes are best-effort, bounded, and side-effect-free.  Unavailable tools
produce ``None``; no physics software is imported or initialised here.

Example (from the repository root)::

    python -m benchmarks.perf.dependencies

Version overrides are useful for copied binary installations without a version
command or a Git checkout: BENCHMARK_FLUKA_VERSION,
BENCHMARK_LINEBUILDER_VERSION, BENCHMARK_FLUKAIO_VERSION,
BENCHMARK_GEANT4_VERSION, and BENCHMARK_BDSIM_VERSION.
"""

from __future__ import annotations

import json
import os
from pathlib import Path
import re
import shutil
import subprocess


EXTERNAL_NAMES = ("fluka", "linebuilder", "flukaio", "geant4", "bdsim")


def _run(command: list[str]) -> str | None:
    """Read a short version string, never failing the benchmark session."""
    try:
        proc = subprocess.run(
            command, capture_output=True, text=True, check=False, timeout=3,
        )
    except (OSError, ValueError, subprocess.TimeoutExpired):
        return None
    if proc.returncode != 0:
        return None
    for line in proc.stdout.splitlines():
        line = line.strip()
        if line:
            return line[:256]
    return None


def _git_revision(path: str | None, *, name: str) -> str | None:
    """Use a Git revision only if the path belongs to the intended checkout."""
    if not path:
        return None
    location = Path(path).expanduser()
    if location.is_file():
        location = location.parent
    if not location.is_dir():
        return None
    root = _run(["git", "-C", str(location), "rev-parse", "--show-toplevel"])
    if not root:
        return None
    root_path = Path(root).resolve()
    location = location.resolve()
    # An arbitrary parent Git repository (e.g. Xcoll) must not be attributed
    # to FlukaIO or LineBuilder merely because a binary lives inside it.
    names = {
        "linebuilder": {"linebuilder", "fluka-linebuilder"},
        "flukaio": {"flukaio", "fluka-io"},
    }
    if root_path != location and root_path.name.lower() not in names[name]:
        return None
    version = _run(["git", "-C", str(root_path), "describe", "--tags", "--always", "--dirty"])
    return f"git:{version}" if version else None


def _fluka_version_from_path() -> str | None:
    # FLUPRO is the usual FLUKA installation variable.  The Xcoll interface
    # can also be configured with an rfluka executable on PATH.
    paths = (
        os.environ.get("FLUPRO"),
        os.environ.get("FLUKA_ROOT"),
        os.environ.get("FLUKA_PATH"),
        str(Path(shutil.which("rfluka")).resolve()) if shutil.which("rfluka") else None,
    )
    for raw in paths:
        if not raw:
            continue
        for component in reversed(Path(raw).parts):
            match = re.fullmatch(r"fluka[-_]?([0-9]+)[._-]([0-9]+)(?:[._-]([0-9]+))?", component, re.I)
            if match:
                return ".".join(group for group in match.groups() if group is not None)
    return None


def external_versions() -> tuple[dict[str, str | None], dict[str, str | None]]:
    """Return (versions, provenance), with no required external dependencies."""
    versions: dict[str, str | None] = {}
    sources: dict[str, str | None] = {}
    for name in EXTERNAL_NAMES:
        override = os.environ.get(f"BENCHMARK_{name.upper()}_VERSION", "").strip()
        if override:
            versions[name], sources[name] = override, "environment override"
            continue
        if name == "geant4":
            value, source = _run(["geant4-config", "--version"]), "geant4-config --version"
        elif name == "bdsim":
            value, source = _run(["bdsim", "--version"]), "bdsim --version"
        elif name == "fluka":
            value, source = _fluka_version_from_path(), "FLUKA installation directory"
        elif name == "linebuilder":
            value = _git_revision(
                os.environ.get("LINEBUILDER_PATH") or os.environ.get("LB_PATH"),
                name="linebuilder",
            )
            source = "LineBuilder Git checkout"
        else:  # flukaio
            value = _git_revision(
                os.environ.get("FLUKAIO_PATH") or os.environ.get("FLUKAIO_ROOT"),
                name="flukaio",
            )
            source = "FlukaIO Git checkout"
        versions[name] = value
        sources[name] = source if value is not None else None
    return versions, sources


def main() -> None:
    versions, sources = external_versions()
    print(json.dumps({"external_versions": versions, "sources": sources}, indent=2))


if __name__ == "__main__":
    main()
