# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

import importlib
import os
import shutil
import subprocess
import sys
from importlib import machinery
from pathlib import Path
from subprocess import run

try:
    from xaux import FsPath  # TODO: once xaux is in Xsuite keep only this
except (ImportError, ModuleNotFoundError):
    from ...xaux import FsPath

from ..environment import BaseEnvironment


class Geant4Environment(BaseEnvironment):
    _read_only_paths = {'bdsim': 0, 'geant4': 0}
    _package_name = 'collimasim'
    _module_name = 'g4interface'

    def __init__(self):
        self._collimasim = None
        cmd = run(['which', 'geant4-config'], capture_output=True)
        self._geant4 = FsPath(cmd.stdout.decode().strip()) if cmd.returncode == 0 else None
        cmd = run(['which', 'bdsim'], capture_output=True)
        self._bdsim = FsPath(cmd.stdout.decode().strip()) if cmd.returncode == 0 else None
        super().__init__()
        self._detect_collimasim_package()

    @property
    def geant4(self):
        return self._geant4

    @property
    def bdsim(self):
        return self._bdsim

    @property
    def collimasim(self):
        return self._collimasim

    @collimasim.setter
    def collimasim(self, value):
        self._collimasim = FsPath(value) if value else None

    @property
    def compiled(self):
        return (self.geant4 is not None and
                self.bdsim is not None and
                self.collimasim is not None and
                self._extension_exists(self.collimasim))

    @property
    def collimasim_python_path(self):
        if not self.collimasim:
            return None
        return self.collimasim.parent

    def _detect_collimasim_package(self):
        package_dir = self._default_package_dir
        if self._extension_exists(package_dir):
            self._collimasim = package_dir
        elif package_dir.exists():
            # Package without extension. Treat as not compiled yet.
            self._collimasim = package_dir
            if not self._extension_exists(package_dir):
                self._collimasim = None

    @property
    def _source_dir(self):
        return FsPath(Path(__file__).parent / 'collimasim_src').resolve()

    @property
    def _default_package_dir(self):
        return FsPath(self.data_dir / 'geant4' / self._package_name).resolve()

    @property
    def _build_dir(self):
        return FsPath(self.data_dir / 'geant4' / 'build').resolve()

    def _extension_exists(self, directory):
        if not directory:
            return False
        directory = FsPath(directory)
        if not directory.exists():
            return False
        for suffix in machinery.EXTENSION_SUFFIXES:
            candidate = directory / f'{self._module_name}{suffix}'
            if candidate.exists():
                return True
        # Fallback glob in case extensions follow custom naming
        pattern = f'{self._module_name}*'
        return any((directory / name).is_file() and name.startswith(self._module_name)
                   and any(name.endswith(sfx) for sfx in machinery.EXTENSION_SUFFIXES)
                   for name in os.listdir(directory))

    def _ensure_python_package(self, package_dir):
        package_dir = FsPath(package_dir)
        package_dir.mkdir(parents=True, exist_ok=True)
        init_file = package_dir / '__init__.py'
        if not init_file.exists():
            init_file.write_text("from .g4interface import XtrackInterface\n\n__all__ = ['XtrackInterface']\n")

    def ensure_python_path(self):
        package_path = self.collimasim_python_path
        if not package_path:
            return None
        path_str = package_path.as_posix()
        if path_str not in sys.path:
            sys.path.insert(0, path_str)
        return path_str

    def compile(self, *, force=False, build_type='Release', cmake_args=None, build_args=None):
        if cmake_args is None:
            cmake_args = []
        if build_args is None:
            build_args = []

        if self.geant4 is None:
            raise RuntimeError('Geant4 not found in PATH. Unable to compile interface.')
        if self.bdsim is None:
            raise RuntimeError('BDSIM not found in PATH. Unable to compile interface.')

        package_dir = self._default_package_dir
        build_dir = self._build_dir

        if force and build_dir.exists():
            shutil.rmtree(build_dir)
        if force and package_dir.exists():
            for item in package_dir.iterdir():
                if item.is_file() and any(item.name.endswith(suffix) for suffix in machinery.EXTENSION_SUFFIXES):
                    item.unlink()

        if not force and self._extension_exists(package_dir):
            self._collimasim = package_dir
            self._ensure_python_package(package_dir)
            return package_dir

        if shutil.which('cmake') is None:
            raise RuntimeError('CMake executable not found in PATH. Cannot compile Geant4 interface.')

        try:
            import pybind11
        except ImportError as exc:
            raise RuntimeError('pybind11 is required to compile the Geant4 interface.') from exc

        self._ensure_python_package(package_dir)
        build_dir.mkdir(parents=True, exist_ok=True)

        source_dir = self._source_dir
        if not source_dir.exists():
            raise FileNotFoundError(f'Geant4 interface sources not found at {source_dir}')

        cmake_cache_args = [
            f'-DCMAKE_LIBRARY_OUTPUT_DIRECTORY={package_dir.as_posix()}',
            f'-DPYTHON_EXECUTABLE={sys.executable}',
            f'-Dpybind11_DIR={pybind11.get_cmake_dir()}',
            f'-DCMAKE_BUILD_TYPE={build_type}',
        ]
        cmake_cache_args.extend(cmake_args)

        configure_cmd = ['cmake', source_dir.as_posix()] + cmake_cache_args
        subprocess.run(configure_cmd, cwd=build_dir, check=True)

        parallel_arg = ['--parallel'] if shutil.which('cmake') else []
        build_cmd = ['cmake', '--build', build_dir.as_posix(), '--config', build_type]
        if parallel_arg:
            build_cmd += parallel_arg
        build_cmd += build_args
        subprocess.run(build_cmd, cwd=build_dir, check=True)

        if not self._extension_exists(package_dir):
            raise RuntimeError('Compilation finished but the Geant4 interface module was not produced.')

        self._collimasim = package_dir
        self.save()
        return package_dir

    def load_collimasim(self):
        try:
            if not self._extension_exists(self.collimasim):
                self.compile()
        except Exception as exc:
            raise ImportError('Failed to prepare the Geant4 interface module.') from exc
        self.ensure_python_path()
        return importlib.import_module(self._package_name)
