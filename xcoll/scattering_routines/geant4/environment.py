# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

import os
from subprocess import run

try:
    from xaux import FsPath  # TODO: once xaux is in Xsuite keep only this
except (ImportError, ModuleNotFoundError):
    from ...xaux import FsPath

from ..environment import BaseEnvironment
from ...general import _pkg_root


class Geant4Environment(BaseEnvironment):
    _read_only_paths = {'bdsim': 0, 'geant4': 0, 'collimasim': 0}

    def __init__(self):
        cmd = run(['which', 'geant4-config'], capture_output=True)
        self._geant4 = FsPath(cmd.stdout.decode().strip()) if cmd.returncode == 0 else None
        cmd = run(['which', 'bdsim'], capture_output=True)
        self._bdsim = FsPath(cmd.stdout.decode().strip()) if cmd.returncode == 0 else None
        super().__init__()

    @property
    def geant4(self):
        return self._geant4

    @property
    def bdsim(self):
        return self._bdsim

    @property
    def collimasim(self):
        return self._collimasim

    @property
    def compiled(self):
        return self.geant4 is not None and self.bdsim is not None and self.collimasim is not None

    def compile(self, verbose=False):
        # Check all dependencies
        self.assert_make_installed(verbose=verbose)
        self.assert_gcc_installed(verbose=verbose)
        if self.geant4 is None:
            raise RuntimeError("Could not find Geant4 installation! Please install Geant4.")
        if self.bdsim is None:
            raise RuntimeError("Could not find BDSIM installation! Please install BDSIM.")

        # Copy the C source files to a temporary directory and compile it
        dest = (self.temp_dir / 'xcoll_bdsim').resolve()
        dest.mkdir(parents=True, exist_ok=True)
        for path in (_pkg_root / 'scattering_routines' / 'geant4' / 'scattering_src').glob('*'):
            if path.name == '.git' or path.name == 'docs' or path.name == 'tests' \
            or path.name == 'samples':
                continue
            else:
                path.copy_to(dest, method='mount')
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
        so = list(dest.glob('g4interface.*so'))
        if len(so) > 1:
            raise RuntimeError(f"Compiled into multiple g4interface shared libraries!")
        if len(so) == 0:
            raise RuntimeError(f"Failed Xcoll-BDSIM compilation! No shared library found in "
                             + f"{self.data_dir.as_posix()}!")
        so = so[0]
        so.move_to(self.data_dir / so.name)
        if verbose:
            print(f"Created Xcoll-BDSIM shared library in {so}.")
        # Clean up the temporary directory
        self.temp_dir = None
