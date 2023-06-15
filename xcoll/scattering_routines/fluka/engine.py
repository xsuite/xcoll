# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2023.                 #
# ######################################### #

import numpy as np
import subprocess
from pathlib import Path
import time

import xobjects as xo


class FLUKAEngine(xo.HybridClass):
    _xofields = {
        'network_port':    xo.Int64
    }

    # The engine is a singleton
    def __new__(cls, *args, **kwargs):
        if not hasattr(cls, 'instance'):
            cls.instance = super().__new__(cls)
            cls.instance._initialised = False
        return cls.instance


    def __init__(self, fluka=None, flukaserver=None, **kwargs):
        if(self._initialised):
            return
        self._initialised = True
        if '_xobject' not in kwargs:
            if fluka is None:
                self._fluka = Path('/', 'eos', 'project-f', 'flukafiles', 'fluka-coupling',
                                        'fluka4-3.3.x86-Linux-gfor9', 'bin', 'rfluka')
            else:
                self._fluka = Path(fluka)
            if flukaserver is None:
                self._flukaserver = Path('/', 'eos', 'project-f', 'flukafiles', 'fluka-coupling',
                                              'fluka_coupling', 'fluka', 'flukaserver')
            else:
                self._flukaserver = Path(flukaserver)
            self._cwd = None
            self._network_nfo = None
            self._log = None
            self._log_fid = None
            self._server_process = None
            self.server_pid = None
            kwargs.setdefault('network_port', 0)
            # test gfortran
        super().__init__(**kwargs)

    def __del__(self, *args, **kwargs):
        self.stop_server()


    @classmethod
    def start_server(cls, input_file, path=None, fluka=None, flukaserver=None):
        cls(fluka, flukaserver)
        this = cls.instance
        if this.is_running():
            print("Server already running.")
            return

        # Declare network
        this._cwd = Path.cwd() if path is None else Path(path)
        this._network_nfo = this._cwd / "network.nfo"
        cmd = subprocess.run(["hostname"], cwd=this._cwd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        if cmd.returncode == 0:
            host = cmd.stdout.decode('UTF-8').strip().split('\n')[0]
        else:
            stderr = cmd.stderr.decode('UTF-8').strip().split('\n')
            raise ValueError(f"Could not declare hostname! Error given is:\n{stderr}")
        with this._network_nfo.open('w') as fid:
            fid.write(f"{host}\n")

        # Start server
        this._log = this._cwd / "server_output.log"
        this._log_fid = this._log.open('w')
        this._server_process = subprocess.Popen([this._fluka, input_file, '-e', this._flukaserver, '-M', "1"], 
                                                 cwd=this._cwd, stdout=this._log_fid, stderr=this._log_fid)
        this.server_pid = this._server_process.pid
        time.sleep(1)
        if not this.is_running():
            raise ValueError(f"Could not start FLUKA server! See logfile {this._log}.")
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
        print(f"Started FLUKA server on network port {this.network_port}.")


    @classmethod
    def stop_server(cls, fluka=None, flukaserver=None):
        cls(fluka, flukaserver)
        this = cls.instance
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


    @classmethod
    def is_running(cls, fluka=None, flukaserver=None):
        cls(fluka, flukaserver)
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
            raise ValueError(f"Could not find username! Error given is:\n{stderr}")
        # Get FLUKA processes for this user
        cmd = subprocess.run(["ps", "-u", whoami], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        if cmd.returncode == 0:
            processes = cmd.stdout.decode('UTF-8').strip().split('\n')
        else:
            stderr = cmd.stderr.decode('UTF-8').strip().split('\n')
            raise ValueError(f"Could not list running processes! Error given is:\n{stderr}")
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
            return False




