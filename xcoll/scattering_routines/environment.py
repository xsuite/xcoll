# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

import os
import sys
import json
import tempfile
from subprocess import run, PIPE
# from platformdirs import user_config_dir, user_data_dir

from ..general import _pkg_root
try:
    from xaux import FsPath  # TODO: once xaux is in Xsuite keep only this
except (ImportError, ModuleNotFoundError):
    from ..xaux import FsPath


class BaseEnvironment:
    _config_dir = FsPath(_pkg_root / 'config').resolve()
    _data_dir   = FsPath(_pkg_root / 'lib').resolve()
    # _config_dir = FsPath(user_config_dir('xcoll')).resolve()
    # _data_dir   = FsPath(user_data_dir('xcoll')).resolve()
    _paths = {} # The value is the parent depth that needs to be brute-forced (0 = file itself, None = no brute-force)
    _read_only_paths = {}

    def __init__(self, *args, **kwargs):
        self._old_sys_path = None
        self._old_os_env = None
        self._temp_dir = None
        self._in_constructor = True
        for path in self._paths.keys():
            setattr(self, f'_{path}', None)
        self._config_dir.mkdir(parents=True, exist_ok=True)
        self._data_dir.mkdir(parents=True, exist_ok=True)
        self._config_file = self._config_dir / f'{self.__class__.__name__[:-11].lower()}.config.json'
        sys.path.append(self._data_dir.as_posix())
        self.load()
        self._in_constructor = False

    def __del__(self):
        self.restore_environment()
        if self._temp_dir:
            self._temp_dir.cleanup()

    def __repr__(self):
        return f"<{self.__class__.__name__} at {hex(id(self))} (use .show() to see the paths)>"

    def __str__(self):
        res = ["XcollEnvironment"]
        res.append(f"    Configuration file:  {self._config_file.as_posix()}")
        res.append(f"    Configuration dir:   {self._config_dir.as_posix()}")
        res.append(f"    Data dir:            {self._data_dir.as_posix()}")
        if self._temp_dir:
            res.append(f"    Temporary dir:       {self._temp_dir.name}")
        if self.__class__ is not BaseEnvironment:
            res.append("")
            res.append(f"{self.__class__.__name__}")
            for path in self._paths.keys():
                value = getattr(self, f'_{path}', None)
                path = f'{path}:'
                if value is None:
                    res.append(f"    {path:<20} None")
                else:
                    res.append(f"    {path:<20} {value.as_posix()}")
            for path in self._read_only_paths.keys():
                value = getattr(self, path)
                path = f'{path}:'
                if value is None:
                    res.append(f"    {path:<20} None (read-only)")
                else:
                    res.append(f"    {path:<20} {value.as_posix()} (read-only)")
            if self._old_sys_path and self._old_os_env:
                res.append("")
                res.append("Custom environment stored:")
                res.append(f"    sys.path: {self._old_sys_path}")
                res.append(f"    os.environ: {self._old_os_env}")
        return "\n".join(res)

    @property
    def config_file(self):
        """The environment configuration file."""
        return self._config_file

    @property
    def config_dir(self):
        """The directory where the configuration files are stored."""
        return self._config_dir

    @property
    def data_dir(self):
        """The directory where the data files are stored."""
        return self._data_dir

    @property
    def initialised(self):
        return all(getattr(self, path, None) is not None and getattr(self, path, None).exists()
                   for path in self._paths.keys())

    @property
    def compiled(self):
        raise NotImplementedError("This property should be implemented in the subclass.")

    @property
    def ready(self):
        return self.initialised and self.compiled

    @property
    def temp_dir(self):
        if not self._temp_dir:
            self._temp_dir = tempfile.TemporaryDirectory()
        return FsPath(self._temp_dir.name)

    @temp_dir.setter
    def temp_dir(self, value):
        if value is None:
            if self._temp_dir:
                self._temp_dir.cleanup()
                self._temp_dir = None
            return
        if self._temp_dir:
            self._temp_dir.cleanup()
            self._temp_dir = None
        if isinstance(value, FsPath):
            value = value.resolve()
        elif isinstance(value, str):
            value = FsPath(value).resolve()
        else:
            raise TypeError("temp_dir must be a string or FsPath.")
        if not value.exists():
            raise FileNotFoundError(f"Provided temp_dir {value} does not exist!")
        self._temp_dir = tempfile.TemporaryDirectory(dir=value)

    @temp_dir.deleter
    def temp_dir(self):
        self.temp_dir = None

    def show(self):
        """Print the environment paths."""
        print(self)

    def save(self):
        data = {'paths': {}, 'read_only_paths': {}}
        for path in self._paths.keys():
            value = getattr(self, path, None)
            if value:
                value = FsPath(value).as_posix()
            data['paths'][path] = value
        for path in self._read_only_paths.keys():
            value = getattr(self, path, None)
            if value:
                value = FsPath(value).as_posix()
            data['read_only_paths'][path] = value
        # Check if anything changed
        with open(self._config_file, 'r') as fid:
            try:
                existing_data = json.load(fid)
            except json.JSONDecodeError:
                # File corrupted, need to save
                existing_data = {}
        if data == existing_data:
            return
        # If yes, save to file
        with open(self._config_file, 'w') as fid:
            json.dump(data, fid, indent=4)

    def load(self):
        if not self._config_file.exists():
            with open(self._config_file, 'w') as fid:
                json.dump({'paths': {}, 'read_only_paths': {}}, fid, indent=4)
        with open(self._config_file, 'r') as fid:
            try:
                data = json.load(fid)
            except json.JSONDecodeError:
                with open(self._config_file, 'w') as fid:
                    json.dump({'paths': {}, 'read_only_paths': {}}, fid, indent=4)
                return
        if 'paths' not in data or 'read_only_paths' not in data:
            return
        for key, value in data['paths'].items():
            setattr(self, key, FsPath(value) if value else None)
        for key, value in data['read_only_paths'].items():
            setattr(self, f'_{key}', FsPath(value) if value else None)

    def store_environment(self):
        self._old_sys_path = sys.path.copy()
        self._old_os_env = os.environ.copy()

    def restore_environment(self):
        if self._old_sys_path:
            sys.path = self._old_sys_path
            self._old_sys_path = None
        if self._old_os_env:
            os.environ = self._old_os_env
            self._old_os_env = None

    def brute_force_path(self, path):
        num_parents = 0
        if str(path) in self._paths:
            num_parents = self._paths[path]
            path = getattr(self, path)
        if str(path) in self._read_only_paths:
            num_parents = self._read_only_paths[path]
            path = getattr(self, path)
        if path is None:
            return
        path = FsPath(path).resolve()
        if num_parents > 0:
            path = path.parents[num_parents-1]
        if not path.exists():
            raise FileNotFoundError(f"Could not find path {path}!")
        try:
            cmd = run(['tree', path.as_posix()], stdout=PIPE, stderr=PIPE)
        except FileNotFoundError:
            # No tree executable. Return as no more can be done. TODO: alternatives for mac and windows?
            return
        if cmd.returncode != 0:
            stderr = cmd.stderr.decode('UTF-8').strip().split('\n')
            raise RuntimeError(f"Could not resolve {path} tree!\nError given is:\n{stderr}")

    def __getattr__(self, key):
        if key in self._paths | self._read_only_paths:
            value = getattr(self, f'_{key}', None)
            if value:
                return FsPath(value)
        else:
            raise AttributeError(f"{self.__class__.__name__} has no attribute '{key}'")

    def __setattr__(self, key, value):
        if key in self._paths.keys():
            if value:
                value = FsPath(value)
                if not self._in_constructor:
                    self.brute_force_path(value)
            old_value = getattr(self, f'_{key}', None)
            if value != old_value:
                super().__setattr__(f'_{key}', value)
                if not self._in_constructor:
                    self.save()
        elif key.startswith('_') and key[1:] in self._read_only_paths.keys():
            # Read-only attribute can only be set internally
            old_value = getattr(self, f'_{key}', None)
            if value != old_value:
                super().__setattr__(key, value)
                if not self._in_constructor:
                    self.save()
        elif key in self._read_only_paths.keys():
            raise AttributeError(f"Attribute '{key}' of {self.__class__.__name__} "
                               + f"is read-only!")
        else:
            super().__setattr__(key, value)

    def __delattr__(self, item):
        if item in self._paths.keys():
            self.__setattr__(self, f'_{item}', None)
            self.save()

    def assert_environment_ready(self):
        if not self.initialised:
            raise RuntimeError(f"{self.__class__.__name__} not initialised! "
                             + f"Please set all paths in the environment before "
                             + f"starting the engine.")
        if not self.compiled:
            raise RuntimeError(f"{self.__class__.__name__} not compiled! "
                             + f"Please compile before starting the engine.")

    def assert_installed(self, program, *, program_name=None, version_cmd="--version", verbose=False):
        if program_name is None:
            program_name = program
        if not hasattr(self, f'_{program}_installed'):
            try:
                cmd = run(["which", program], stdout=PIPE, stderr=PIPE)
            except Exception as err:
                new_err = RuntimeError(f"Could not run `which {program}`!")
                new_err.__cause__ = err
                setattr(self, f'_{program}_installed', [new_err, None])
            else:
                if cmd.returncode != 0:
                    stderr = cmd.stderr.decode('UTF-8').strip().split('\n')
                    err = RuntimeError(f"Could not run `which {program}` (output "
                                    f"error code {cmd.returncode})!\nError given is:\n"
                                    f"{stderr}")
                    setattr(self, f'_{program}_installed', [err, None])
                else:
                    file = cmd.stdout.decode('UTF-8').strip().split('\n')[0]
                    try:
                        cmd = run([program, version_cmd], stdout=PIPE, stderr=PIPE)
                    except Exception as err:
                        new_err = RuntimeError(f"Could not run `{program} {version_cmd}`!")
                        new_err.__cause__ = err
                        setattr(self, f'_{program}_installed', [file, new_err])
                    else:
                        if cmd.returncode != 0:
                            stderr = cmd.stderr.decode('UTF-8').strip().split('\n')
                            err = RuntimeError(f"Could not run `{program} {version_cmd}` "
                                            f"(output error code {cmd.returncode})!\nError "
                                            f"given is:\n{stderr}")
                            setattr(self, f'_{program}_installed', [file, err])
                        else:
                            version = cmd.stdout.decode('UTF-8').strip().split('\n')[0]
                            setattr(self, f'_{program}_installed', [file, version])
        file, version = getattr(self, f'_{program}_installed')
        if isinstance(file, Exception):
            # Which failed!
            raise file
        elif isinstance(version, Exception):
            # Version command failed!
            if verbose:
                print(f"Found {program_name} in {file}", flush=True)
            raise version
        else:
            if verbose:
                print(f"Found {program_name} ({version}) in {file}", flush=True)
            return file, version


    def assert_gcc_installed(self, minimum_version=9, verbose=False):
        gcc = os.environ['CC'] or 'gcc'
        _, version = self.assert_installed(gcc, program_name='CC', version_cmd='-dumpversion',
                                           verbose=verbose)
        if int(version.split('.')[0]) < minimum_version:
            self._gcc_installed = False
            raise RuntimeError(f"Need gcc {minimum_version} or higher, but found gcc {version}!")

    def assert_gxx_installed(self, minimum_version=9, verbose=False):
        gxx = os.environ['CXX'] or 'g++'
        _, version = self.assert_installed(gxx, program_name='CXX', version_cmd='-dumpversion',
                                           verbose=verbose)
        if int(version.split('.')[0]) < minimum_version:
            self._gxx_installed = False
            raise RuntimeError(f"Need gxx {minimum_version} or higher, but found gxx {version}!")

    def assert_gfortran_installed(self, minimum_version=9, verbose=False):
        gfortran = os.environ['FC'] or 'gfortran'
        _, version = self.assert_installed(gfortran, program_name='FC', version_cmd='-dumpversion',
                                           verbose=verbose)
        if int(version.split('.')[0]) < minimum_version:
            self._gfortran_installed = False
            raise RuntimeError(f"Need gfortran {minimum_version} or higher, but found gfortran {version}!")
