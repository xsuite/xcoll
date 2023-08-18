# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2023.                 #
# ######################################### #

import numpy as np
import subprocess
from pathlib import Path
import time

import xobjects as xo


default_fluka_path = Path('/', 'eos', 'project-f', 'flukafiles', 'fluka-coupling', 'fluka4-3.3.x86-Linux-gfor9',
                          'bin', 'rfluka').resolve()
default_flukaserver_path = Path('/', 'eos', 'project-f', 'flukafiles', 'fluka-coupling', 'fluka_coupling',
                                'fluka', 'flukaserver').resolve()

network_file = "network.nfo"
fluka_log  = "server_output.log"
server_log = "server_output.log"


class FlukaEngine(xo.HybridClass):
    _xofields = {
        'network_port':    xo.Int64,
        'n_alloc':         xo.Int64
    }

    # The engine is a singleton
    def __new__(cls, *args, **kwargs):
        if not hasattr(cls, 'instance'):
            cls.instance = super().__new__(cls)
            cls.instance._initialised = False
        return cls.instance

    def __del__(self, *args, **kwargs):
        self.stop_server()

    def __init__(self, fluka=None, flukaserver=None, verbose=True, testing=False, **kwargs):
        if(self._initialised):
            return
        self._initialised = True
        if '_xobject' not in kwargs:

            # Initialise python fields
            self._verbose = verbose
            self._cwd = None
            self._network_nfo = None
            self._log = None
            self._log_fid = None
            self._server_process = None
            self.server_pid = None
            self._gfortran_installed = False
            self._flukaio_connected = False
            self._warning_given = False
            kwargs.setdefault('network_port', 0)
            kwargs.setdefault('n_alloc', 500000)

            # Get paths to executables
            if fluka is None:
                if testing:
                    self._fluka = "NEED_TO_MAKE_MOCKUP"
                else:
                    self._fluka = default_fluka_path
            elif testing:
                raise ValueError("Cannot specify fluka executable when 'testing=True'!")
            else:
                self._fluka = Path(fluka).resolve()
            if not self._fluka.exists():
                raise ValueError(f"Could not find fluka executable {self._fluka}!")
            if flukaserver is None:
                if testing:
                    self._flukaserver = "NEED_TO_MAKE_MOCKUP"
                else:
                    self._flukaserver = default_flukaserver_path
            elif testing:
                raise ValueError("Cannot specify flukaserver executable when 'testing=True'!")
            else:
                self._flukaserver = Path(flukaserver).resolve()
            if not self._flukaserver.exists():
                raise ValueError(f"Could not find flukaserver executable {self._flukaserver}!")

            self.test_gfortran()
            try:
                from .pyflukaf import pyfluka_init
                pyfluka_init(n_alloc=kwargs['n_alloc'])
            except ImportError as error:
                self._warn_pyfluka(error)

        super().__init__(**kwargs)


    def _warn_pyfluka(self, error):
        if not self._warning_given:
            print("Warning: Failed to import pyfluka (did you compile?). " \
                + "FlukaCollimators will be installed but are not trackable.\n")
            print(error)
            self._warning_given = True


    @classmethod
    def test_gfortran(cls, *args, **kwargs):
        cls(*args, **kwargs)
        this = cls.instance
        try:
            cmd = subprocess.run(["gfortran", "-dumpversion"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        except FileNotFoundError:
            this._gfortran_installed = False
            raise RuntimeError("Could not find gfortran installation! Need gfortran 9 or higher for fluka.")
        if cmd.returncode == 0:
            version = int(cmd.stdout.decode('UTF-8').strip().split('\n')[0])
            if version < 9:
                this._gfortran_installed = False
                raise RuntimeError(f"Need gfortran 9 or higher for fluka, but found gfortran {version}!")
            this._gfortran_installed = True
            if this._verbose:
                cmd2 = subprocess.run(["which", "gfortran"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                if cmd2.returncode == 0:
                    file = cmd2.stdout.decode('UTF-8').strip().split('\n')[0]
                    print(f"Found gfortran {version} in {file}")
                else:
                    stderr = cmd2.stderr.decode('UTF-8').strip().split('\n')
                    raise RuntimeError(f"Could not run 'which gfortran'! Error given is:\n{stderr}")
        else:
            stderr = cmd.stderr.decode('UTF-8').strip().split('\n')
            this._gfortran_installed = False
            raise RuntimeError(f"Could not run gfortran! Verify its installation.")


    @classmethod
    def start_server(cls, input_file, *args, **kwargs):
        cls(*args, **kwargs)
        this = cls.instance
        if this.is_running():
            print("Server already running.")
            return
        if not this._gfortran_installed:
            return

        # Check files
        input_file = Path(input_file).resolve()
        if not input_file.exists():
            raise ValueError(f"Input file {input_file.as_posix()} not found!")
        this._cwd = input_file.parent
        insertion_file = this._cwd / "insertion.txt"
        if not insertion_file.exists():
            raise ValueError(f"Insertion file {insertion_file.as_posix()} not found!")

        # Start with cleaning files
        files_to_delete = [network_file, fluka_log, server_log]
        for f in files_to_delete:
            if f is not None:
                f = this._cwd / f
                if f.exists():
                    f.unlink()

        # Declare network
        this._network_nfo = this._cwd / network_file
        cmd = subprocess.run(["hostname"], cwd=this._cwd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        if cmd.returncode == 0:
            host = cmd.stdout.decode('UTF-8').strip().split('\n')[0]
        else:
            stderr = cmd.stderr.decode('UTF-8').strip().split('\n')
            raise RuntimeError(f"Could not declare hostname! Error given is:\n{stderr}")
        with this._network_nfo.open('w') as fid:
            fid.write(f"{host}\n")

        # Start server
        this._log = this._cwd / server_log
        this._log_fid = this._log.open('w')
        this._server_process = subprocess.Popen([this._fluka, input_file, '-e', this._flukaserver, '-M', "1"], 
                                                 cwd=this._cwd, stdout=this._log_fid, stderr=this._log_fid)
        this.server_pid = this._server_process.pid
        time.sleep(1)
        if not this.is_running():
            raise ValueError(f"Could not start fluka server! See logfile {this._log}.")
        i = 0
        while True:
            time.sleep(2)
            i += 1
            if i%30 == 0:
                print("Network port not yet established. Waiting 30s.")
            with this._network_nfo.open('r') as fid:
                lines = fid.readlines()
                if len(lines) > 1:
                    this.network_port = int(lines[1].strip())
                    break
        print(f"Started fluka server on network port {this.network_port}.")
        try:
            from .pyflukaf import pyfluka_connect
            pyfluka_connect()
            this._flukaio_connected = True
        except ImportError as error:
            this._warn_pyfluka(error)


    @classmethod
    def stop_server(cls, *args, **kwargs):
        cls(*args, **kwargs)
        this = cls.instance
        # Stop flukaio conneciton
        if this._flukaio_connected:
            try:
                from .pyflukaf import pyfluka_close
                pyfluka_close()
            except ImportError as error:
                this._warn_pyfluka(error)
            this._flukaio_connected = False
        # If the Popen process is still running, terminate it
        if this._server_process is not None:
            if this._server_process.poll() is None:
                this._server_process.terminate()
            this._server_process = None
        this.server_pid = None
        # Close the file pointer to the log
        if this._log_fid is not None:
            if not this._log_fid.closed:
                this._log_fid.close()
            this._log_fid = None
        # Delete network file
        if this._network_nfo is not None and this._network_nfo.exists():
            this._network_nfo.unlink()


    @classmethod
    def is_running(cls, *args, **kwargs):
        cls(*args, **kwargs)
        this = cls.instance
        # Is the Popen process still running?
        if this._server_process is None or this._server_process.poll() is not None:
            this.stop_server()
            return False
        # Get username (need a whoami for the next command)
        cmd = subprocess.run(["whoami"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        if cmd.returncode == 0:
            whoami = cmd.stdout.decode('UTF-8').strip().split('\n')[0]
        else:
            stderr = cmd.stderr.decode('UTF-8').strip().split('\n')
            raise RuntimeError(f"Could not find username! Error given is:\n{stderr}")
        # Get fluka processes for this user
        cmd = subprocess.run(["ps", "-u", whoami], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        if cmd.returncode == 0:
            processes = cmd.stdout.decode('UTF-8').strip().split('\n')
        else:
            stderr = cmd.stderr.decode('UTF-8').strip().split('\n')
            raise RuntimeError(f"Could not list running processes! Error given is:\n{stderr}")
        processes = [proc for proc in processes if 'rfluka' in proc]
        if len(processes) == 0:
            # Could not find a running rfluka
            this.stop_server()
            return False
        elif len(processes) == 1 and str(this.server_pid) in processes[0] and 'defunct' not in processes[0]:
            return True
        elif np.any([str(this.server_pid) in proc for proc in processes if 'defunct' not in proc]):
            print("Warning: Found other instances of rfluka besides the current one!")
            return True
        else:
            print("Warning: Found other instances of rfluka but not the current one!")
            this.stop_server()
            return False




