# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

import os
import sys
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
        self._geant4_sourced = False
        self._bdsim_sourced = False
        try:
            cmd = run(['which', 'geant4-config'], capture_output=True)
        except FileNotFoundError:
            pass
        else:
            if cmd.returncode == 0:
                path = FsPath(cmd.stdout.decode().strip())
                if path.exists():
                    self._geant4 = path
                    self._geant4_sourced = True
        try:
            cmd = run(['which', 'bdsim'], capture_output=True)
        except FileNotFoundError:
            pass
        else:
            if cmd.returncode == 0:
                path = FsPath(cmd.stdout.decode().strip())
                if path.exists():
                    self._bdsim = path
                    self._bdsim_sourced = True

    @property
    def compiled(self):
        if self.geant4 is None or self.bdsim is None:
            return False
        try:
            from g4interface import XtrackInterface
            return True
        except (ModuleNotFoundError, ImportError) as error:
            return False

    @property
    def ready(self):
        return super().ready and self._geant4_sourced and self._bdsim_sourced

    def compile(self, verbose=True, verbose_compile_output=False):
        # Check all dependencies
        self.assert_installed('make', verbose=verbose)
        self.assert_installed('cmake', verbose=verbose)
        self.assert_gcc_installed(verbose=verbose)
        self.assert_gxx_installed(verbose=verbose)
        self.assert_geant4_installed()
        self.assert_bdsim_installed()
        bdsim_version = self.get_bdsim_version()
        if verbose:
            print(f"BDSIM version: {bdsim_version}")

        # Check pybind11 is installed
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
        self._adapt_source_to_bdsim_version(bdsim_version, verbose)

        # Configure
        ctab = '    '
        cmd = run(['cmake', '-S', '.', '-B', 'build'], capture_output=True)
        if cmd.returncode != 0:
            stderr = cmd.stderr.decode('UTF-8').strip()
            os.chdir(cwd)
            raise RuntimeError(f"Failed to build Xcoll-BDSIM interface!\nError given is:\n{stderr}")
        if verbose_compile_output:
            print()
            print("CMake: Configuring")
            print(ctab + cmd.stdout.decode('UTF-8').strip().replace('\n', f'\n{ctab}'))
            if cmd.stderr:
                print(ctab + cmd.stderr.decode('UTF-8').strip().replace('\n', f'\n{ctab}'))
            print()

        # Build
        cmd = run(['cmake', '--build', 'build'], capture_output=True)
        if cmd.returncode != 0:
            stderr = cmd.stderr.decode('UTF-8').strip()
            os.chdir(cwd)
            raise RuntimeError(f"Failed to compile Xcoll-BDSIM interface!\nError given is:\n{stderr}")
        if verbose_compile_output:
            print("CMake: Building")
            print(ctab + cmd.stdout.decode('UTF-8').strip().replace('\n', f'\n{ctab}'))
            if cmd.stderr:
                print(ctab + cmd.stderr.decode('UTF-8').strip().replace('\n', f'\n{ctab}'))
            print()
        if verbose:
            print("Compiled Xcoll-BDSIM interface successfully.")

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
            print(f"Created Xcoll-BDSIM shared library in {self.data_dir / so.name}.")
        # Clean up the temporary directory
        self.temp_dir = None
        os.chdir(cwd)

    def assert_geant4_installed(self):
        if self.geant4 is None:
            raise RuntimeError("Could not find Geant4 installation! Please install Geant4.")
        if not self._geant4_sourced:
            raise RuntimeError(f"Geant4 installation found in {self.geant4} "
                               f"but not active! Please source environment.")

    def assert_bdsim_installed(self):
        if self.bdsim is None:
            raise RuntimeError("Could not find BDSIM installation! Please install BDSIM.")
        if not self._bdsim_sourced:
            raise RuntimeError(f"BDSIM installation found in {self.bdsim} "
                               f"but not active! Please source environment.")

    def assert_environment_ready(self):
        self.assert_geant4_installed()
        self.assert_bdsim_installed()
        super().assert_environment_ready()


    def get_bdsim_version(self):
        cmd = run(['bdsim', '--version'], capture_output=True)
        if cmd.returncode != 0:
            stderr = cmd.stderr.decode('UTF-8').strip()
            raise RuntimeError(f"Failed to run 'bdsim --version'!\nError given is:\n{stderr}")
        return cmd.stdout.decode('UTF-8').strip()

    def bdsim_older_than(self, bdsim_version=None, compare_version='1.7.7.develop'):
        if bdsim_version is None:
            bdsim_version = self.get_bdsim_version()
        n_ver = sum([10**(3*(2-i))*int(j) for i, j in enumerate(bdsim_version.strip().split('.')[:3])])
        if 'develop' in bdsim_version:
            n_ver += 0.5
        n_ver_cmp = sum([10**(3*(2-i))*int(j) for i, j in enumerate(compare_version.strip().split('.')[:3])])
        if 'develop' in compare_version:
            n_ver_cmp += 0.5
        return n_ver < n_ver_cmp


    def _adapt_source_to_bdsim_version(self, bdsim_version, verbose):
        version_mark = '// BDSIM >= '
        adapted_any = False

        # Adapt BDSXtrackInterface.hh
        adapted = False
        with open('BDSXtrackInterface.hh', 'r') as file:
            filedata = file.read()
        if self.bdsim_older_than(bdsim_version, '1.7.7.develop'):
            filedata = filedata.replace('#include "BDSLinkBunch.hh"', '#include "BDSBunchSixTrackLink.hh"')
            filedata = filedata.replace('BDSLinkBunch* stp = nullptr;', 'BDSBunchSixTrackLink* stp = nullptr;')
            adapted = True
        new_filedata = []
        for line in filedata.split('\n'):
            if version_mark in line:
                if self.bdsim_older_than(bdsim_version, line.split(version_mark)[1]):
                    adapted = True
                    continue
            new_filedata.append(line)
        if adapted:
            adapted_any = True
            with open('BDSXtrackInterface.hh', 'w') as file:
                file.write('\n'.join(new_filedata))

        # Adapt BDSXtrackInterface.cpp
        adapted = False
        with open('BDSXtrackInterface.cpp', 'r') as file:
            filedata = file.read()
        if self.bdsim_older_than(bdsim_version, '1.7.7.develop'):
            filedata = filedata.replace('stp = new BDSLinkBunch();', 'stp = new BDSBunchSixTrackLink();')
            adapted = True
        new_filedata = []
        for line in filedata.split('\n'):
            if version_mark in line:
                if self.bdsim_older_than(bdsim_version, line.split(version_mark)[1]):
                    adapted = True
                    continue
            new_filedata.append(line)
        if adapted:
            adapted_any = True
            with open('BDSXtrackInterface.cpp', 'w') as file:
                file.write('\n'.join(new_filedata))

        if verbose and adapted_any:
                print("Adapted source code to BDSIM version.")
