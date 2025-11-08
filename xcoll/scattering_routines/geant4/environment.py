# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

import os
import requests
from subprocess import run

from ..environment import BaseEnvironment
from ...general import _pkg_root
try:
    from xaux import FsPath  # TODO: once xaux is in Xsuite keep only this
except (ImportError, ModuleNotFoundError):
    from ...xaux import FsPath


class Geant4Environment(BaseEnvironment):
    _read_only_paths = {'bdsim': 0, 'geant4': 0}

    def __init__(self):
        super().__init__()
        self._in_constructor = True
        self._geant4 = None
        self._bdsim = None
        self._in_constructor = False
        cmd = run(['which', 'geant4-config'], capture_output=True)
        if cmd.returncode == 0:
            path = FsPath(cmd.stdout.decode().strip())
            if path.exists():
                self._geant4 = path
        cmd = run(['which', 'bdsim'], capture_output=True)
        if cmd.returncode == 0:
            path = FsPath(cmd.stdout.decode().strip())
            if path.exists():
                self._bdsim = path

    @property
    def compiled(self):
        if self.geant4 is None or self.bdsim is None:
            return False
        so = list((self.data_dir).glob('g4interface.*so'))
        if len(so) == 1 and so[0].exists():
            return True
        return False

    def compile(self, verbose=True):
        # Check all dependencies
        self.assert_make_installed(verbose=verbose)
        self.assert_gcc_installed(verbose=verbose)
        if self.geant4 is None:
            raise RuntimeError("Could not find Geant4 installation! Please install Geant4.")
        if self.bdsim is None:
            raise RuntimeError("Could not find BDSIM installation! Please install BDSIM.")
        try:
            import pybind11 # noqa F401
        except (ModuleNotFoundError, ImportError):
            pybind_repo = 'https://github.com/pybind/pybind11.git'
            try:
                not_found = requests.get(pybind_repo).status_code == 404
            except (requests.exceptions.HTTPError, requests.exceptions.ConnectionError):
                not_found = True
            if not_found:
                raise RuntimeError("Could not find pybind11! Please pip install "
                                  f"pybind11, or make sure {pybind_repo} is accessible.")

        # Copy the C source files to a temporary directory and compile it
        dest = (self.temp_dir / 'xcoll_bdsim').resolve()
        dest.mkdir(parents=True, exist_ok=True)
        for path in (_pkg_root / 'scattering_routines' / 'geant4' / 'scattering_src').glob('*'):
            if path.name == '.git' or path.name == 'docs' or path.name == 'tests' \
            or path.name == 'samples':
                continue
            else:
                FsPath(path).copy_to(dest, method='mount')
        cwd = FsPath.cwd()
        os.chdir(dest)

        # Compile
        cmd = run(['cmake', '-S', '.', '-B', 'build'], capture_output=True)
        if cmd.returncode != 0:
            stderr = cmd.stderr.decode('UTF-8').strip().split('\n')
            os.chdir(cwd)
            raise RuntimeError(f"Failed to build Xcoll-BDSIM interface!\nError given is:\n{stderr}")
        cmd = run(['cmake', '--build', 'build'], capture_output=True)
        if cmd.returncode == 0:
            if verbose:
                print("Compiled Xcoll-BDSIM interface successfully.")
        else:
            stderr = cmd.stderr.decode('UTF-8').strip().split('\n')
            os.chdir(cwd)
            raise RuntimeError(f"Failed to compile Xcoll-BDSIM interface!\nError given is:\n{stderr}")
        # Collect the compiled shared library
        so = list((dest / 'build').glob('g4interface.*so'))
        if len(so) > 1:
            raise RuntimeError(f"Compiled into multiple g4interface shared libraries!")
        if len(so) == 0:
            raise RuntimeError(f"Failed Xcoll-BDSIM compilation! No shared library found in "
                             + f"{self.data_dir.as_posix()}!")
        so = FsPath(so[0])
        so.move_to(self.data_dir / so.name)
        if verbose:
            print(f"Created Xcoll-BDSIM shared library in {so}.")
        # Clean up the temporary directory
        self.temp_dir = None
        os.chdir(cwd)
