# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

import os
import sys
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
    _paths = {'fluka': 2, 'flukaserver': 1, 'linebuilder': 0}
    _optional_paths = {'flair': 1}
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
        try:
            from pyflukaf import track_fluka
            return True
        except (ModuleNotFoundError, ImportError) as error:
            return False

    def compile(self, flukaio_lib=None, flukaio_path=None, verbose=False):
        # Check all dependencies
        self.assert_installed('make', verbose=verbose)
        self.assert_installed('cmake', verbose=verbose)
        self.assert_gcc_installed(verbose=verbose)
        self.assert_gfortran_installed(verbose=verbose)
        cwd = FsPath.cwd()
        self.store_environment()

        # Get FlukaIO
        if flukaio_lib is None:
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
            os.chdir(dest)
            cmd = run(['make', 'clean'], stdout=PIPE, stderr=PIPE)
            cmd = run(['make', 'libs', 'BUILD64=Y'], stdout=PIPE, stderr=PIPE)
            if cmd.returncode == 0:
                if verbose:
                    print("Compiled FlukaIO successfully.")
            else:
                stderr = cmd.stderr.decode('UTF-8').strip()
                os.chdir(cwd)
                raise RuntimeError(f"Failed to compile FlukaIO!\nError given is:\n{stderr}")
            flukaio_lib = dest  / 'lib' / 'libFlukaIO64.a'
        flukaio_lib = FsPath(flukaio_lib).resolve()
        if not flukaio_lib.exists():
            raise FileNotFoundError(f"Could not find FlukaIO library {flukaio_lib}!")

        # Copy the FORTRAN source files to the temporary directory and compile it
        dest = (self.temp_dir / 'FORTRAN_src').resolve()
        dest.mkdir(parents=True, exist_ok=True)
        flukaio_lib.copy_to(dest / 'libFlukaIO.a', method='mount')
        for path in _FORTRAN_SRC.glob('*'):
            path.copy_to(dest, method='mount')
        os.chdir(dest)
        cmd = run(['python', '-m', 'numpy.f2py', '-m', 'pyflukaf', 'pyfluka.f90',
                   '--quiet'], stdout=PIPE, stderr=PIPE)
        if cmd.returncode == 0:
            if verbose:
                print("Created numpy FORTRAN headers.")
        else:
            stderr = cmd.stderr.decode('UTF-8').strip()
            os.chdir(cwd)
            raise RuntimeError(f"Failed to create numpy FORTRAN headers!\nError given is:\n{stderr}")
        cmd = run(['meson', 'setup', 'build'], stdout=PIPE, stderr=PIPE)
        if cmd.returncode == 0:
            if verbose:
                print("Setup meson build successfully.")
        else:
            stdout = cmd.stdout.decode('UTF-8').strip()
            stderr = cmd.stderr.decode('UTF-8').strip()
            os.chdir(cwd)
            raise RuntimeError(f"Failed to setup meson build!\n"
                               f"Output given is:\n{stdout}"
                               f"Error given is:\n{stderr}")
        cmd = run(["meson", "compile", "-C", "build"], stdout=PIPE, stderr=PIPE)
        if cmd.returncode == 0:
            if verbose:
                print(cmd.stdout.decode('UTF-8').strip())
                print("Compiled pyflukaf successfully.")
        else:
            stdout = cmd.stdout.decode('UTF-8').strip()
            stderr = cmd.stderr.decode('UTF-8').strip()
            os.chdir(cwd)
            raise RuntimeError(f"Failed to compile pyflukaf!\n"
                               f"Output given is:\n{stdout}"
                               f"Error given is:\n{stderr}")
        os.chdir(cwd)
        # Collect the compiled shared library
        so = list((self.temp_dir / 'FORTRAN_src' / 'build').glob('pyflukaf.*so'))
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
        self.restore_environment()


    def import_fedb(self, fedb_path, verbose=False, overwrite=False):
        import xcoll as xc
        fedb_path = FsPath(fedb_path)
        if not fedb_path.exists():
            raise FileNotFoundError(f"Could not find FEDB {fedb_path}!")
        self._init_fedb(overwrite=True)  # Always overwrite the build-in FEDB (for in-between modifications in Xcoll)
        # Get assemblies
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
        # Get prototypes
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
        # Get stepsizes (currently not used)
        for file in fedb_path.glob(f'stepsizes/*'):
            file.copy_to(self.fedb / 'stepsizes' / file.name, method='mount')
        # Link the new files into the registry
        self._init_fedb()


    def set_fluka_environment(self):
        self.brute_force_path('fluka')
        self.brute_force_path('flukaserver')

    def set_flair_environment(self):
        self.brute_force_path('flair')

    def set_fedb_environment(self, fedb=None):
        self.store_environment()
        if fedb is None:
            fedb = self.fedb
        self.brute_force_path(fedb)
        self.brute_force_path('linebuilder')
        os.environ['FEDB_PATH'] = fedb.as_posix()
        os.environ['LB_PATH'] = self.linebuilder.as_posix()
        # Brute-force the system paths
        sys.path.insert(0, (self.linebuilder / "src").as_posix())
        sys.path.insert(0, (self.linebuilder / "lib").as_posix())
        sys.path.insert(0, fedb.as_posix())


    def create_temp_fedb(self, assemblies):
        fedb = FsPath('temp_fedb').resolve()
        fedb.mkdir(parents=True)
        (fedb / 'assemblies').mkdir(parents=True)
        (fedb / 'bodies').mkdir(parents=True)
        (fedb / 'regions').mkdir(parents=True)
        (fedb / 'materials').mkdir(parents=True)
        (fedb / 'stepsizes').mkdir(parents=True)
        (fedb / 'tools').symlink_to(self.fedb / 'tools')
        (self.fedb / 'structure.py').copy_to(fedb / 'structure.py')
        (self.fedb / 'materials' / 'materials.inp').copy_to(fedb / 'materials' / 'materials.inp')
        for assm in assemblies:
            assm.populate_into_temp_fedb(fedb)
        return fedb


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
            stderr = cmd.stderr.decode('UTF-8').strip()
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
            stderr = cmd.stderr.decode('UTF-8').strip()
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

    def _init_fedb(self, overwrite=False):
        from xcoll import FlukaPrototype
        self.fedb.mkdir(parents=True, exist_ok=True)
        for directory in ['assemblies', 'bodies', 'regions', 'materials', 'stepsizes', 'metadata']:
            (self.fedb / directory).mkdir(parents=True, exist_ok=True)
            for f in (_FEDB_TEMPLATE / directory).glob('*.*'):
                new_file = self.fedb / directory / f.name
                if new_file.exists() and overwrite:
                    new_file.unlink()
                if not new_file.exists():
                    f.copy_to(new_file, method='mount')
        tools = self.fedb / 'tools'
        if tools.exists() and overwrite:
            tools.rmtree()
        if not tools.exists():
            (_FEDB_TEMPLATE / 'tools').copy_to(tools.parent, method='mount')
        structure = self.fedb / 'structure.py'
        if structure.exists() and overwrite:
            structure.unlink()
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
