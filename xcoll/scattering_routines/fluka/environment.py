# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

import os
import sys
from math import floor, log10
from subprocess import run, PIPE

try:
    from xaux import ClassProperty, ClassPropertyMeta, singleton, FsPath  # TODO: once xaux is in Xsuite keep only this
except (ImportError, ModuleNotFoundError):
    from ...xaux import ClassProperty, ClassPropertyMeta, singleton, FsPath

from ...general import _pkg_root


@singleton(allow_underscore_vars_in_init=False)
class FlukaEnvironment(metaclass=ClassPropertyMeta):
    _default_fluka_eos_path = FsPath('/eos/project/f/flukafiles/fluka-coupling').resolve()

    def __init__(self, *args, **kwargs):
        self._gfortran_installed = False
        self._old_sys_path = None
        self._old_os_env = None
        self._overwritten_paths = {}
        self.fluka = kwargs.pop('fluka', None)
        self.flukaserver = kwargs.pop('flukaserver', None)
        self.flair = kwargs.pop('flair', None)
        self.fedb = kwargs.pop('fedb', None)
        self.linebuilder = kwargs.pop('linebuilder', None)

    def __del__(self):
        self._restore_fedb_base()
        if self._old_sys_path:
            sys.path = self._old_sys_path
        if self._old_os_env:
            os.environ = self._old_os_env


    # =============
    # === Paths ===
    # =============

    @ClassProperty
    def fluka(cls):
        return cls.get_self()._fluka

    @fluka.setter
    def fluka(cls, val):
        self = cls.get_self()
        if val is None:
            self._fluka = (self._default_fluka_eos_path / 'fluka4-4.1' / 'bin' / 'rfluka').resolve()
        else:
            val = FsPath(val)
            self._brute_force_path(val.parents[1])
            self._fluka = val

    @fluka.deleter
    def fluka(cls):
        cls.fluka = None

    @ClassProperty
    def flukaserver(cls):
        return cls.get_self()._flukaserver

    @flukaserver.setter
    def flukaserver(cls, val):
        self = cls.get_self()
        if val is None:
            self._flukaserver = (self._default_fluka_eos_path / 'fluka_coupling' / 'fluka' / 'flukaserver').resolve()
        else:
            val = FsPath(val)
            self._brute_force_path(val.parent)
            self._flukaserver = val

    @flukaserver.deleter
    def flukaserver(cls):
        cls.flukaserver = None

    @ClassProperty
    def flair(cls):
        return cls.get_self()._flair

    @flair.setter
    def flair(cls, val):
        self = cls.get_self()
        if val is None:
            self._flair = (self._default_fluka_eos_path / 'flair-3.3' /  'flair').resolve()
        else:
            val = FsPath(val)
            self._brute_force_path(val.parent)
            self._flair = val

    @flair.deleter
    def flair(cls):
        cls.flair = None

    @ClassProperty
    def linebuilder(cls):
        return cls.get_self()._linebuilder

    @linebuilder.setter
    def linebuilder(cls, val):
        self = cls.get_self()
        if val is None:
            # linebuilder = (_linebuilder_coupling / 'linebuilder').resolve()
            # TODO if MR on gitlab accepted, update path
            self._linebuilder = FsPath("/eos/project/c/collimation-team/software/fluka_coupling_tmp_patch_xsuite").resolve()
        else:
            val = FsPath(val)
            self._brute_force_path(val)
            self._linebuilder = val

    @linebuilder.deleter
    def linebuilder(cls):
        cls.linebuilder = None

    @ClassProperty
    def fedb(cls):
        # This is the user fedb path, only for user-defined assemblies
        return cls.get_self()._fedb

    @fedb.setter
    def fedb(cls, val):
        self = cls.get_self()
        if val is None:
            self._restore_fedb_base()
            self._fedb = None
        else:
            val = FsPath(val)
            self._brute_force_path(val)
            self._fedb = val
            self._sync_user_fedb()

    @fedb.deleter
    def fedb(cls):
        cls.fedb = None

    @ClassProperty
    def fedb_base(cls):
        return (_pkg_root / 'scattering_routines' / 'fluka' / 'fedb').resolve()

    @ClassProperty
    def fedb_xcoll(cls):
        return (cls.fedb_base / 'fedb_xcoll').resolve()

    @ClassProperty
    def fedb_coupling(cls):
        return (cls.fedb_base / 'fedb_coupling').resolve()

    # ======================
    # === Public Methods ===
    # ======================

    @classmethod
    def test_gfortran(cls, verbose=False):
        self = cls.get_self()
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

    @classmethod
    def set_fluka_environment(cls):
        self = cls.get_self()
        self._brute_force_path(self.fluka.parents[1])
        self._brute_force_path(self.flukaserver.parent)

    @classmethod
    def set_flair_environment(cls):
        self = cls.get_self()
        self._brute_force_path(self.flair.parent)

    @classmethod
    def set_fedb_environment(cls):
        self = cls.get_self()
        self._old_sys_path = sys.path.copy()
        self._old_os_env = os.environ.copy()
        self._brute_force_path(self.fedb_coupling)
        self._brute_force_path(self.linebuilder)
        os.environ['FEDB_PATH'] = self.fedb_base.as_posix()
        os.environ['LB_PATH'] = self.linebuilder.as_posix()
        # Brute-force the system paths
        sys.path.append(self.fedb_base.as_posix())
        sys.path.append((self.fedb_base / "tools").as_posix())
        sys.path.append((self.fedb_base / "tools" / "materials" / "cables").as_posix())
        sys.path.append((self.linebuilder / "src").as_posix())
        sys.path.append((self.linebuilder / "lib").as_posix())

    @classmethod
    def unset_fedb_environment(cls):
        self = cls.get_self()
        sys.path = self._old_sys_path
        self._old_sys_path = None
        os.environ = self._old_os_env
        self._old_os_env = None

    @classmethod
    def run_flair(cls, input_file=None):
        self = cls.get_self()
        if input_file is None:
            return
        input_file = FsPath(input_file)
        try:
            self.set_flair_environment()
        except FileNotFoundError:
            print("Flair not found. Cannot view input file.")
            return
        cmd = run([FlukaEnvironment().flair.as_posix(), input_file.as_posix()],
                  stdout=PIPE, stderr=PIPE)
        if cmd.returncode != 0:
            stderr = cmd.stderr.decode('UTF-8').strip().split('\n')
            raise RuntimeError(f"Failed to run flair on input file {input_file}!\n"
                             + f"Error given is:\n{stderr}")

    @classmethod
    def test_assembly(cls, fedb_series, fedb_tag, *, show=True, keep_files=False):
        if show:
            try:
                cls.set_flair_environment()
            except FileNotFoundError:
                print("Flair not found. Cannot view assembly.")
                return
        else:
            keep_files=True
        self = cls.get_self()
        self.set_fedb_environment()
        cmd = run(['python', self.fedb_base / 'tools' / 'test_assembly.py', fedb_series,
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

    def _brute_force_path(self, path):
        if not path.exists():
            raise FileNotFoundError(f"Could not find path {path}!")
        try:
            cmd = run(['tree', path.as_posix()], stdout=PIPE, stderr=PIPE)
        except FileNotFoundError:
            raise RuntimeError(f"Could not find 'tree' exectuable!")
        if cmd.returncode != 0:
            stderr = cmd.stderr.decode('UTF-8').strip().split('\n')
            self._gfortran_installed = False
            raise RuntimeError(f"Could not resolve {path} tree!\nError given is:\n{stderr}")

    def _sync_user_fedb(self):
        if self._fedb:
            for directory in ['assemblies', 'bodies', 'materials', 'regions']:
                for file in self._fedb.glob(f'{directory}/*.lbp'):
                    path = self.fedb_base / directory / file.name
                    if path.exists():
                        if not path.is_symlink():
                            raise FileExistsError(f"File {path} exists but is not a symlink! "
                                                + "Only symlinks are allowed in the FEDB inside "
                                                + "Xcoll. Something is wrong.")
                        self._overwritten_paths[path] = path.readlink().as_posix()
                        path.unlink()
                    path.symlink_to(file)

    def _restore_fedb_base(self):
        for path, target in self._overwritten_paths.items():
            path.unlink()
            path.symlink_to(target)
        self._overwritten_paths = {}


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
