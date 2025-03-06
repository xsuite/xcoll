# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

import sys
import os
from subprocess import run, PIPE

try:
    from xaux import singleton, FsPath  # TODO: once xaux is in Xsuite keep only this
except (ImportError, ModuleNotFoundError):
    from ...xaux import singleton, FsPath

from ...general import _pkg_root


@singleton
class FlukaEnvironment:
    _default_fluka_eos_path = FsPath('/eos/project/f/flukafiles/fluka-coupling').resolve()

    def __init__(self):
        if not self._initialised:
            self._gfortran_installed = False
            self.fluka = None
            self.flukaserver = None
            self.flair = None
            self.fedb = None
            self.linebuilder = None
            self._old_sys_path = None
            self._old_os_env = None
            self._overwritten_paths = {}

    def __del__(self):
        self._restore_base_fedb()
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
            self._fluka = (self._default_fluka_eos_path / 'fluka4-4.1' / 'bin' / 'rfluka').resolve()
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
            self._flair = (self._default_fluka_eos_path / 'flair-3.3' /  'flair').resolve()
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
        # This is the user fedb path, only for user-defined assemblies
        return self._fedb

    @fedb.setter
    def fedb(self, val):
        if val is None:
            self._restore_base_fedb()
            self._fedb = None
        else:
            val = FsPath(val)
            self._brute_force_path(val)
            self._fedb = val
            self._sync_user_fedb()

    @fedb.deleter
    def fedb(self):
        self.fedb = None

    @property
    def fedb_base(self):
        return (_pkg_root / 'scattering_routines' / 'fluka' / 'fedb').resolve()

    @property
    def fedb_xcoll(self):
        return (self.fedb_base / 'fedb_xcoll').resolve()

    @property
    def fedb_coupling(self):
        return (self.fedb_base / 'fedb_coupling').resolve()

    # ======================
    # === Public Methods ===
    # ======================

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
        os.environ['FEDB_PATH'] = self.fedb_base.as_posix()
        os.environ['LB_PATH'] = self.linebuilder.as_posix()
        # Brute-force the system paths
        sys.path.append(self.fedb_base.as_posix())
        sys.path.append((self.fedb_base / "tools").as_posix())
        sys.path.append((self.fedb_base / "tools" / "materials" / "cables").as_posix())
        sys.path.append((self.linebuilder / "src").as_posix())
        sys.path.append((self.linebuilder / "lib").as_posix())

    def unset_fedb_environment(self):
        sys.path = self._old_sys_path
        self._old_sys_path = None
        os.environ = self._old_os_env
        self._old_os_env = None

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

    def _restore_base_fedb(self):
        for path, target in self._overwritten_paths.items():
            path.unlink()
            path.symlink_to(target)
        self._overwritten_paths = {}
