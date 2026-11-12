"""Compile an external Xcoll interface for the *current* Python environment.

Machine-level installations are configured through environment variables.
The build occurs inside the disposable benchmark Conda environment.

Run from the repository root:
    python -m benchmarks.perf.setup_external fluka
    python -m benchmarks.perf.setup_external geant4
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path


def _required(name: str) -> Path:
    raw = os.getenv(name)
    if not raw:
        raise RuntimeError(f"Missing {name}; configure the runner environment or GitHub Actions variable")
    path = Path(raw).expanduser().resolve()
    if not path.exists():
        raise FileNotFoundError(f"{name}: {path} does not exist")
    return path


def setup_fluka():
    import xcoll as xc

    root = _required("XCOLL_FLUKA_ROOT")
    interface = xc.fluka.interface
    interface.fluka = root / "bin" / "rfluka"
    interface.flukaserver = _required("XCOLL_FLUKA_SERVER")
    interface.linebuilder = _required("XCOLL_LINEBUILDER_ROOT")
    lib = _required("XCOLL_FLUKAIO_LIB")
    for executable in (interface.fluka, interface.flukaserver):
        if not Path(executable).exists():
            raise FileNotFoundError(f"FLUKA executable not found: {executable}")
    interface.compile(flukaio_lib=lib, verbose=True)
    fedb = os.getenv("XCOLL_FLUKA_FEDB")
    if fedb:
        interface.import_fedb(_required("XCOLL_FLUKA_FEDB"),
                              verbose=True, overwrite=True)
    print("FLUKA interface compiled for the active Python environment")


def setup_geant4():
    import xcoll as xc

    # A conda environment activated with bdsim-g4 provides the executable
    # paths and Geant4 configuration; no external setup scripts are needed.
    import shutil

    for executable in ("geant4-config", "bdsim"):
        if shutil.which(executable) is None:
            raise RuntimeError(f"{executable} not on PATH in activated Conda environment")
    xc.geant4.interface.compile(verbose=True, verbose_compile_output=False)
    print("Geant4/BDSIM interface compiled for the active Python environment")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("engine", choices=["fluka", "geant4"])
    args = parser.parse_args(argv)
    {"fluka": setup_fluka, "geant4": setup_geant4}[args.engine]()


if __name__ == "__main__":
    main()
