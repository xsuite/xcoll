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

from ..environment import BaseEnvironment
from ...general import _pkg_root


_FORTRAN_SRC   = FsPath(_pkg_root / 'scattering_routines' / 'fluka' / 'FORTRAN_src').resolve()
_FEDB_TEMPLATE = FsPath(_pkg_root / 'scattering_routines' / 'fluka' / 'fedb').resolve()


class FlukaEnvironment(BaseEnvironment):
    # The paths to be set. The value is the parent depth that needs to be brute-forced (0 = file itself, None = no brute-force)
    _paths = {'fluka': 2, 'flukaserver': 1, 'flair': 1, 'linebuilder': 0}
    _read_only_paths = {'fedb': 0}

    def __init__(self):
        super().__init__()
        # Initialize the FEDB
        self._init_fedb()

    @property
    def fedb(self):
        return self._data_dir / 'fedb'

    @property
    def compiled(self):
        so = list(self.data_dir.glob('pyflukaf.*so'))
        return len(so) > 0

    def compile(self, flukaio_path=None, verbose=False):
        # Check all dependencies
        self.assert_installed('make', verbose=verbose)
        self.assert_installed('cmake', verbose=verbose)
        self.assert_gcc_installed(verbose=verbose)
        self.assert_gfortran_installed(verbose=verbose)

        # Check the provided FlukaIO path
        if flukaio_path is None:
            raise ValueError("FlukaIO path must be provided!")
        flukaio_path = FsPath(flukaio_path).resolve()
        if not flukaio_path.exists():
            raise FileNotFoundError(f"Could not find FlukaIO path {flukaio_path}!")
        self.brute_force_path(flukaio_path)
        # Copy the FlukaIO files to a temporary directory and compile it
        dest = (self.temp_dir / 'flukaio').resolve()
        dest.mkdir(parents=True, exist_ok=True)
        for path in flukaio_path.glob('*'):
            if path.name == '.git' or path.name == 'docs' or path.name == 'tests' \
            or path.name == 'samples':
                continue
            else:
                path.copy_to(dest, method='mount')
        cwd = FsPath.cwd()
        os.chdir(dest)
        cmd = run(['make', 'clean'], stdout=PIPE, stderr=PIPE)
        cmd = run(['make', 'libs', 'BUILD64=Y'], stdout=PIPE, stderr=PIPE)
        if cmd.returncode == 0:
            if verbose:
                print("Compiled FlukaIO successfully.")
        else:
            stderr = cmd.stderr.decode('UTF-8').strip().split('\n')
            os.chdir(cwd)
            raise RuntimeError(f"Failed to compile FlukaIO!\nError given is:\n{stderr}")
        # Copy the FORTRAN source files to the temporary directory and compile it
        dest = (self.temp_dir / 'FORTRAN_src').resolve()
        dest.mkdir(parents=True, exist_ok=True)
        for path in _FORTRAN_SRC.glob('*'):
            path.copy_to(dest, method='mount')
        os.chdir(dest)
        cmd = run(['gfortran', '-fpic', '-c', 'core_tools.f90', 'constants.f90', 'strings.f90',
                   'mod_alloc.f90', 'common_modules.f90', 'string_tools.f90', 'mod_units.f90',
                   'pdgid.f90', 'mod_fluka.f90'], stdout=PIPE, stderr=PIPE)
        if cmd.returncode == 0:
            if verbose:
                print("Compiled FORTRAN source successfully.")
        else:
            stderr = cmd.stderr.decode('UTF-8').strip().split('\n')
            os.chdir(cwd)
            raise RuntimeError(f"Failed to compile FORTRAN source!\nError given is:\n{stderr}")
        # Link the FORTRAN source files and the FlukaIO library to create the pyflukaf shared library
        cmd = run(['f2py', '-m', 'pyflukaf', '-c', 'pyfluka.f90',
                   'core_tools.o', 'constants.o', 'strings.o', 'mod_alloc.o',
                   'common_modules.o', 'string_tools.o', 'mod_units.o',
                   'pdgid.o', 'mod_fluka.o', '../flukaio/lib/libFlukaIO64.a',
                   '--backend', 'distutils', '--fcompiler=gfortran'], stdout=PIPE, stderr=PIPE)
        if cmd.returncode == 0:
            if verbose:
                print("Linked pyflukaf successfully.")
        else:
            stderr = cmd.stderr.decode('UTF-8').strip().split('\n')
            os.chdir(cwd)
            raise RuntimeError(f"Failed to link pyflukaf!\nError given is:\n{stderr}")
        os.chdir(cwd)
        # Collect the compiled shared library
        so = list((self.temp_dir / 'FORTRAN_src').glob('pyflukaf.*so'))
        if len(so) > 1:
            raise RuntimeError(f"Compiled into multiple pyflukaf shared libraries!")
        if len(so) == 0:
            raise RuntimeError(f"Failed pyFLUKA compilation! No shared library found in "
                             + f"{self.data_dir.as_posix()}!")
        so = so[0]
        so.move_to(self.data_dir / so.name)
        if verbose:
            print(f"Created pyFLUKA shared library in {so}.")
        # Clean up the temporary directory
        self.temp_dir = None


    def import_fedb(self, fedb_path, verbose=False, overwrite=False):
        import xcoll as xc
        fedb_path = FsPath(fedb_path)
        if not fedb_path.exists():
            raise FileNotFoundError(f"Could not find FEDB path {fedb_path}!")
        if not self.fedb.exists():
            self.fedb.mkdir(parents=True, exist_ok=True)
        for file in fedb_path.glob(f'assemblies/*.lbp'):
            if (self.fedb / 'assemblies' / file.name).exists() and not overwrite:
                if verbose:
                    print(f"Warning: Assembly {file.name} already exists in Xcoll. "
                         + "Use overwrite=True to overwrite it.")
                continue
            meta = _FEDB_TEMPLATE / 'metadata' / 'coupling' / f'{file.name}.json'
            if meta.exists():
                pro = xc.FlukaAssembly.from_json(meta)
                if verbose:
                    print(f"Imported assembly {pro.fedb_series}_{pro.fedb_tag}")
            else:
                if verbose:
                    print(f"Warning: No metadata found for assembly {file.name} in {fedb_path}!")
                part = file.stem.split('_')
                fedb_series = part[0]
                fedb_tag = '_'.join(part[1:])
                pro = xc.FlukaAssembly(fedb_series, fedb_tag)
            pro.assembly_file = file
        for file in fedb_path.glob(f'bodies/*.bodies'):
            if (self.fedb / 'bodies' / file.name).exists() and not overwrite:
                if verbose:
                    print(f"Warning: Prototype {file.name} already exists in Xcoll. "
                         + "Use overwrite=True to overwrite it.")
                continue
            meta = _FEDB_TEMPLATE / 'metadata' / 'coupling' / f'{file.name}.json'
            if meta.exists():
                pro = xc.FlukaPrototype.from_json(meta)
                if verbose:
                    print(f"Imported prototype {pro.fedb_series}_{pro.fedb_tag}")
            else:
                if verbose:
                    print(f"Warning: No metadata found for prototype {file.name} in {fedb_path}!")
                part = file.stem.split('_')
                fedb_series = part[0]
                fedb_tag = '_'.join(part[1:])
                pro = xc.FlukaPrototype(fedb_series, fedb_tag)
            pro.body_file = file
            pro.material_file = fedb_path / 'materials' / f'{file.stem}.assignmat'
            pro.region_file = fedb_path / 'regions' / f'{file.stem}.regions'
        for file in fedb_path.glob(f'stepsizes/*'):
            file.copy_to(self.fedb / 'stepsizes' / file.name, method='mount')
        tools = self.fedb / 'tools'
        if tools.exists():
            tools.rmtree()
        (_FEDB_TEMPLATE / 'tools').copy_to(tools, method='mount')
        structure = self.fedb / 'structure.py'
        structure.unlink(missing_ok=True)
        (_FEDB_TEMPLATE / 'structure.py').copy_to(structure, method='mount')


    def set_fluka_environment(self):
        self.brute_force_path('fluka')
        self.brute_force_path('flukaserver')

    def set_flair_environment(self):
        self.brute_force_path('flair')

    def set_fedb_environment(self):
        self.store_environment()
        self.brute_force_path(self.fedb)
        self.brute_force_path('linebuilder')
        os.environ['FEDB_PATH'] = self.fedb.as_posix()
        os.environ['LB_PATH'] = self.linebuilder.as_posix()
        # Brute-force the system paths
        sys.path.insert(0, (self.linebuilder / "src").as_posix())
        sys.path.insert(0, (self.linebuilder / "lib").as_posix())
        sys.path.insert(0, self.fedb.as_posix())


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
        self.restore_environment()
        if cmd.returncode != 0:
            stderr = cmd.stderr.decode('UTF-8').strip().split('\n')
            raise RuntimeError(f"Failed to run flair!\nError given is:\n{stderr}")
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

    def _init_fedb(self):
        from xcoll import FlukaPrototype
        if not self.fedb.exists():
            self.fedb.mkdir(parents=True, exist_ok=True)
        for directory in ['assemblies', 'bodies', 'regions', 'materials', 'stepsizes', 'metadata']:
            (self.fedb / directory).mkdir(parents=True, exist_ok=True)
            for f in (_FEDB_TEMPLATE / directory).glob('*.*'):
                link = self.fedb / directory / f.name
                if not link.exists():
                    link.symlink_to(f)
        tools = self.fedb / 'tools'
        if not tools.exists():
            (_FEDB_TEMPLATE / 'tools').copy_to(tools, method='mount')
        structure = self.fedb / 'structure.py'
        if not structure.exists():
            (_FEDB_TEMPLATE / 'structure.py').copy_to(structure, method='mount')
        prototypes = (self.fedb / 'metadata').glob('*_*.json')
        for file in prototypes:
            FlukaPrototype.from_json(file)


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
