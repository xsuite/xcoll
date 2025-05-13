# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

import os
import sys
import json
from math import floor, log10
from subprocess import run, PIPE

try:
    from xaux import FsPath  # TODO: once xaux is in Xsuite keep only this
except (ImportError, ModuleNotFoundError):
    from ...xaux import FsPath

from ...general import _pkg_root


class FlukaEnvironment:
    _default_fluka_eos_path = FsPath('/eos/project/f/flukafiles/fluka-coupling').resolve()

    def __init__(self, *args, **kwargs):
        self._make_installed = False
        self._cmake_installed = False
        self._gcc_installed = False
        self._gfortran_installed = False
        self._old_sys_path = None
        self._old_os_env = None
        self._overwritten_paths = {}
        self.fluka = kwargs.pop('fluka', None)
        self.flukaserver = kwargs.pop('flukaserver', None)
        self.flair = kwargs.pop('flair', None)
        self.linebuilder = kwargs.pop('linebuilder', None)
        if not (self.fedb / 'index.json').exists():
            self._generate_xcoll_index()

    def __del__(self):
        if self._old_sys_path:
            sys.path = self._old_sys_path
        if self._old_os_env:
            os.environ = self._old_os_env


    # =============
    # === Paths ===
    # =============

    @property
    def fluka(self):
        return self._fluka

    @fluka.setter
    def fluka(self, val):
        if val is None:
            self._fluka = (self._default_fluka_eos_path / 'fluka4-5.0' / 'bin' / 'rfluka').resolve()
        else:
            val = FsPath(val)
            self._brute_force_path(val.parents[1])
            self._fluka = val

    @fluka.deleter
    def fluka(self):
        self.fluka = None

    @property
    def flukaserver(self):
        return self._flukaserver

    @flukaserver.setter
    def flukaserver(self, val):
        if val is None:
            self._flukaserver = (self._default_fluka_eos_path / 'fluka_coupling' / 'fluka' / 'flukaserver').resolve()
        else:
            val = FsPath(val)
            self._brute_force_path(val.parent)
            self._flukaserver = val

    @flukaserver.deleter
    def flukaserver(self):
        self.flukaserver = None

    @property
    def flair(self):
        return self._flair

    @flair.setter
    def flair(self, val):
        if val is None:
            self._flair = (self._default_fluka_eos_path / 'flair-3.4' /  'flair').resolve()
        else:
            val = FsPath(val)
            self._brute_force_path(val.parent)
            self._flair = val

    @flair.deleter
    def flair(self):
        self.flair = None

    @property
    def linebuilder(self):
        return self._linebuilder

    @linebuilder.setter
    def linebuilder(self, val):
        if val is None:
            # linebuilder = (_linebuilder_coupling / 'linebuilder').resolve()
            # TODO if MR on gitlab accepted, update path
            self._linebuilder = FsPath("/eos/project/c/collimation-team/software/fluka_coupling_tmp_patch_xsuite").resolve()
        else:
            val = FsPath(val)
            self._brute_force_path(val)
            self._linebuilder = val

    @linebuilder.deleter
    def linebuilder(self):
        self.linebuilder = None

    @property
    def fedb(self):
        return (_pkg_root / 'scattering_routines' / 'fluka' / 'fedb').resolve()

    @property
    def index(self):
        with open(self.fedb / 'index.json', 'r') as fid:
            data = json.load(fid)
        return data

    # ======================
    # === Public Methods ===
    # ======================

    def compile(self, flukaio_path=None, verbose=False):
        self.test_make(verbose=verbose)
        self.test_gcc(verbose=verbose)
        self.test_gfortran(verbose=verbose)
        if flukaio_path is None:
            flukaio_path = FsPath('/eos/project-c/collimation-team/software/flukaio').resolve()
        else:
            flukaio_path = FsPath(flukaio_path).resolve()
        if not flukaio_path.exists():
            raise FileNotFoundError(f"Could not find FlukaIO path {flukaio_path}!")
        self._brute_force_path(flukaio_path)
        # TODO: xaux fails if copying from EOS
        # flukaio_path.copy_to(_pkg_root / 'scattering_routines' / 'fluka')
        import shutil
        shutil.copytree(flukaio_path, _pkg_root / 'scattering_routines' / 'fluka' / 'flukaio',
                        dirs_exist_ok=True, ignore=shutil.ignore_patterns('.git'))
        cwd = FsPath.cwd()
        so = (_pkg_root / 'scattering_routines' / 'fluka').glob('pyflukaf.*so')
        for file in so:
            file.unlink()
        os.chdir(_pkg_root / 'scattering_routines' / 'fluka' / 'flukaio')
        cmd = run(['make', 'libs', 'BUILD64=Y'], stdout=PIPE, stderr=PIPE)
        if cmd.returncode == 0:
            print("Compiled FlukaIO successfully.")
        else:
            stderr = cmds.stderr.decode('UTF-8').strip().split('\n')
            os.chdir(cwd)
            raise RuntimeError(f"Failed to compile FlukaIO!\nError given is:\n{stderr}")
        os.chdir(_pkg_root / 'scattering_routines' / 'fluka' / 'FORTRAN_src')
        mod = (_pkg_root / 'scattering_routines' / 'fluka' / 'FORTRAN_src').glob('*.{mod,o}')
        for file in mod:
            file.unlink()
        cmd = run(['gfortran', '-fpic', '-c', 'core_tools.f90', 'constants.f90', 'strings.f90',
                   'mod_alloc.f90', 'common_modules.f90', 'string_tools.f90', 'mod_units.f90',
                   'pdgid.f90', 'mod_fluka.f90'], stdout=PIPE, stderr=PIPE)
        if cmd.returncode == 0:
            print("Compiled FORTRAN source successfully.")
        else:
            stderr = cmds.stderr.decode('UTF-8').strip().split('\n')
            os.chdir(cwd)
            raise RuntimeError(f"Failed to compile FORTRAN source!\nError given is:\n{stderr}")
        cmd = run(['f2py', '-m', 'pyflukaf', '-c', 'pyfluka.f90',
                   'core_tools.o', 'constants.o', 'strings.o', 'mod_alloc.o',
                   'common_modules.o', 'string_tools.o', 'mod_units.o',
                   'pdgid.o', 'mod_fluka.o', '../flukaio/lib/libFlukaIO64.a',
                   '--backend', 'distutils', '--fcompiler=gfortran'], stdout=PIPE, stderr=PIPE)
        if cmd.returncode == 0:
            print("Linked pyflukaf successfully.")
        else:
            stderr = cmd.stderr.decode('UTF-8').strip().split('\n')
            os.chdir(cwd)
            raise RuntimeError(f"Failed to link pyflukaf!\nError given is:\n{stderr}")
        os.chdir(cwd)
        # Clean up
        mod = (_pkg_root / 'scattering_routines' / 'fluka' / 'FORTRAN_src').glob('*.{mod,o}')
        for file in mod:
            file.unlink()
        shutil.rmtree(_pkg_root / 'scattering_routines' / 'fluka' / 'flukaio')
        so = list((_pkg_root / 'scattering_routines' / 'fluka' / 'FORTRAN_src').glob('pyflukaf.*so'))
        if len(so) > 1:
            raise RuntimeError(f"Found multiple pyflukaf shared libraries!")
        if len(so) == 0:
            raise RuntimeError(f"Failed pyFLUKA compilation! No shared library found in "
                             + f"{(_pkg_root / 'scattering_routines' / 'fluka').as_posix()}!")
        so = so[0]
        so = so.rename(_pkg_root / 'scattering_routines' / 'fluka' / so.name)
        print(f"Created pyFLUKA shared library in {so}.")

    def test_make(self, verbose=False):
        try:
            cmd = run(["make", "--version"], stdout=PIPE, stderr=PIPE)
        except FileNotFoundError:
            self._make_installed = False
            raise RuntimeError("Could not find make installation!")
        if cmd.returncode == 0:
            self._make_installed = True
        else:
            stderr = cmd.stderr.decode('UTF-8').strip().split('\n')
            self._make_installed = False
            raise RuntimeError(f"Could not run make! Verify its installation.\nError given is:\n{stderr}")
        try:
            cmd = run(["cmake", "--version"], stdout=PIPE, stderr=PIPE)
        except FileNotFoundError:
            self._cmake_installed = False
            raise RuntimeError("Could not find cmake installation!")
        if cmd.returncode == 0:
            self._cmake_installed = True
        else:
            stderr = cmd.stderr.decode('UTF-8').strip().split('\n')
            self._cmake_installed = False
            raise RuntimeError(f"Could not run cmake! Verify its installation.\nError given is:\n{stderr}")

    def test_gcc(self, verbose=False):
        try:
            cmd = run(["gcc", "-dumpversion"], stdout=PIPE, stderr=PIPE)
        except FileNotFoundError:
            self._gcc_installed = False
            raise RuntimeError("Could not find gcc installation! Need gcc 9 or higher.")
        if cmd.returncode == 0:
            version = cmd.stdout.decode('UTF-8').strip().split('\n')[0]
            if int(version.split('.')[0]) < 9:
                self._gcc_installed = False
                raise RuntimeError(f"Need gcc 9 or higher, but found gcc {version}!")
            self._gcc_installed = True
            if verbose:
                cmd2 = run(["which", "gcc"], stdout=PIPE, stderr=PIPE)
                if cmd2.returncode == 0:
                    file = cmd2.stdout.decode('UTF-8').strip().split('\n')[0]
                    print(f"Found gcc {version} in {file}", flush=True)
                else:
                    stderr = cmd2.stderr.decode('UTF-8').strip().split('\n')
                    raise RuntimeError(f"Could not run `which gcc`!\nError given is:\n{stderr}")
        else:
            stderr = cmd.stderr.decode('UTF-8').strip().split('\n')
            self._gcc_installed = False
            raise RuntimeError(f"Could not run gcc! Verify its installation.\nError given is:\n{stderr}")

    def test_gfortran(self, verbose=False):
        try:
            cmd = run(["gfortran", "-dumpversion"], stdout=PIPE, stderr=PIPE)
        except FileNotFoundError:
            self._gfortran_installed = False
            raise RuntimeError("Could not find gfortran installation! Need gfortran 9 or higher.")
        if cmd.returncode == 0:
            version = cmd.stdout.decode('UTF-8').strip().split('\n')[0]
            if int(version.split('.')[0]) < 9:
                self._gfortran_installed = False
                raise RuntimeError(f"Need gfortran 9 or higher, but found gfortran {version}!")
            self._gfortran_installed = True
            if verbose:
                cmd2 = run(["which", "gfortran"], stdout=PIPE, stderr=PIPE)
                if cmd2.returncode == 0:
                    file = cmd2.stdout.decode('UTF-8').strip().split('\n')[0]
                    print(f"Found gfortran {version} in {file}", flush=True)
                else:
                    stderr = cmd2.stderr.decode('UTF-8').strip().split('\n')
                    raise RuntimeError(f"Could not run `which gfortran`!\nError given is:\n{stderr}")
        else:
            stderr = cmd.stderr.decode('UTF-8').strip().split('\n')
            self._gfortran_installed = False
            raise RuntimeError(f"Could not run gfortran! Verify its installation.\nError given is:\n{stderr}")


    def import_fedb(self, fedb_path, overwrite=True):
        fedb_path = FsPath(fedb_path)
        if not fedb_path.exists():
            raise FileNotFoundError(f"Could not find FEDB path {fedb_path}!")
        with open(self.fedb / 'index.json', 'r') as fid:
            data = json.load(fid)
        for directory in ['assemblies', 'bodies', 'materials', 'regions', 'stepsizes']:
            for file in fedb_path.glob(f'{directory}/*'):
                if file.name == 'materials.inp':
                    # TODO: later we want this to be a materials DB in Xcoll and import potential new
                    # materials from this file into the Xcoll materials DB, and regenerate the file
                    continue
                if f'{directory}/{file.name}' in data['xcoll_files']:
                    raise FileExistsError(f"Trying to import file {file.name} from {fedb_path}, "
                                         + "but it already exists in the Xcoll FEDB!")
                if f'{directory}/{file.name}' in data['external_files']:
                    if not overwrite:
                        print(f"Trying to import file {file.name} from {fedb_path}, "
                              + "but it already exists in the FEDB! Skipped.")
                        continue
                    else:
                        # Remove the old file
                        path = self.fedb / directory / file.name
                        path.unlink()
                        file.copy_to(self.fedb / directory / file.name)
                        data['external_files'].append(f'{directory}/{file.name}')
        with open(self.fedb / 'index.json', 'w') as fid:
            json.dump(index, fid, indent=4)


    def set_fluka_environment(self):
        self._brute_force_path(self.fluka.parents[1])
        self._brute_force_path(self.flukaserver.parent)

    def set_flair_environment(self):
        self._brute_force_path(self.flair.parent)

    def set_fedb_environment(self):
        self._old_sys_path = sys.path.copy()
        self._old_os_env = os.environ.copy()
        self._brute_force_path(self.fedb_coupling)
        self._brute_force_path(self.linebuilder)
        os.environ['FEDB_PATH'] = self.fedb.as_posix()
        os.environ['LB_PATH'] = self.linebuilder.as_posix()
        # Brute-force the system paths
        sys.path.append(self.fedb.as_posix())
        sys.path.append((self.fedb / "tools").as_posix())
        sys.path.append((self.fedb / "tools" / "materials" / "cables").as_posix())
        sys.path.append((self.linebuilder / "src").as_posix())
        sys.path.append((self.linebuilder / "lib").as_posix())

    def unset_fedb_environment(self):
        sys.path = self._old_sys_path
        self._old_sys_path = None
        os.environ = self._old_os_env
        self._old_os_env = None


    def run_flair(self, input_file=None):
        if input_file is None:
            return
        input_file = FsPath(input_file)
        try:
            self.set_flair_environment()
        except FileNotFoundError:
            print("Flair not found. Cannot view input file.")
            return
        cmd = run([self.flair.as_posix(), input_file.as_posix()],
                  stdout=PIPE, stderr=PIPE)
        if cmd.returncode != 0:
            stderr = cmd.stderr.decode('UTF-8').strip().split('\n')
            raise RuntimeError(f"Failed to run flair on input file {input_file}!\n"
                             + f"Error given is:\n{stderr}")


    def test_assembly(self, fedb_series, fedb_tag, *, show=True, keep_files=False):
        if show:
            try:
                self.set_flair_environment()
            except FileNotFoundError:
                print("Flair not found. Cannot view assembly.")
                return
        else:
            keep_files=True
        self.set_fedb_environment()
        cmd = run(['python', self.fedb / 'tools' / 'test_assembly.py', fedb_series,
                   fedb_tag], stdout=PIPE, stderr=PIPE)
        self.unset_fedb_environment()
        if cmd.returncode != 0:
            stderr = cmd.stderr.decode('UTF-8').strip().split('\n')
            raise RuntimeError(f"Failed to run flair on input file {input_file}!\n"
                             + f"Error given is:\n{stderr}")
        file = FsPath.cwd() / f"{fedb_series}_{fedb_tag}.inp"
        if not file.exists():
            raise FileNotFoundError(f"Temporary input file {file} not generated!")
        logfile = file.with_suffix('.log')
        if not logfile.exists():
            raise FileNotFoundError(f"Temporary log file {logfile} not generated!")
        self.run_flair(file)
        if not keep_files:
            file.unlink()
            logfile.unlink()


    # =======================
    # === Private Methods ===
    # =======================

    def _add_to_index(self, directory, file):
        file = FsPath(file)
        with open(self.fedb / 'index.json', 'r') as fid:
            data = json.load(fid)
        data['xcoll_files'].append(f'{directory}/{file.name}')
        with open(self.fedb / 'index.json', 'w') as fid:
            json.dump(data, fid, indent=4)

    def _generate_xcoll_index(self):
        index = {'xcoll_files': [], 'external_files': []}
        if not (self.fedb / 'stepsizes').exists():
            (self.fedb / 'stepsizes').mkdir()
        for directory in ['assemblies', 'bodies', 'materials', 'regions', 'stepsizes']:
             for file in self.fedb.glob(f'{directory}/*'):
                index['xcoll_files'].append(f'{directory}/{file.name}')
        with open(self.fedb / 'index.json', 'w') as fid:
            json.dump(index, fid, indent=4)

    def _brute_force_path(self, path):
        if not path.exists():
            raise FileNotFoundError(f"Could not find path {path}!")
        try:
            cmd = run(['tree', path.as_posix()], stdout=PIPE, stderr=PIPE)
        except FileNotFoundError:
            # No tree executable. Return as no more can be done
            return
        if cmd.returncode != 0:
            stderr = cmd.stderr.decode('UTF-8').strip().split('\n')
            raise RuntimeError(f"Could not resolve {path} tree!\nError given is:\n{stderr}")


def format_fluka_float(value):
    if abs(value) > 9.99e99:
        raise ValueError(f"Value {value} is too large for FLUKA.")
    elif abs(value) < 1e-99:
        return f"{'0.0':>10}"
    if 1.e-4 <= value <= 999999999. or -1.e-4 >= value >= -99999999.:
        max_digits = 8 if value > 0 else 7
        n_decimal_digits = int(max_digits - max(floor(log10(abs(value))), 0))
        value_string = f"{value:10.{n_decimal_digits}f}"
        value_splitted = value_string.split('.')
        if len(value_splitted) == 1:
            return f"{value_string[1:]}."
        else:
            decimals = value_splitted[1]
            while decimals[-1] == '0':
                decimals = decimals[:-1]
                if decimals == '':
                    decimals = '0'
                    break
            value_string = f"{value_splitted[0]}.{decimals}"
            return f"{value_string:>10}"
    else:
        max_digits = 4 if value > 0 else 3
        value_string = f"{value:.{max_digits}E}"
        value_splitted = value_string.split('E')
        value_splitted_mant = value_splitted[0].split('.')
        if len(value_splitted_mant) == 1:
            return value_string
        else:
            decimals = value_splitted_mant[1]
            dot = '.'
            while decimals[-1] == '0':
                decimals = decimals[:-1]
                if decimals == '':
                    dot = ''
                    break
            value_string = f"{value_splitted_mant[0]}{dot}{decimals}E{value_splitted[1]}"
            return f"{value_string:>10}"
