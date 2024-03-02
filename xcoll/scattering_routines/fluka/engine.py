# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2023.                 #
# ######################################### #

import numpy as np
import subprocess
from pathlib import Path
import time

import xobjects as xo
import xpart as xp
import xpart.pdg as pdg
import xtrack as xt

from .reference_masses import source, masses


default_fluka_path = Path('/', 'eos', 'project-f', 'flukafiles', 'fluka-coupling', 'fluka4-3.4',
                          'bin', 'rfluka').resolve()
default_flukaserver_path = Path('/', 'eos', 'project-f', 'flukafiles', 'fluka-coupling', 'fluka_coupling',
                                'fluka', 'flukaserver').resolve()
network_file = "network.nfo"
fluka_log  = "server_output.log"
server_log = "server_output.log"


class FlukaEngine(xo.HybridClass):
    _xofields = {
        'network_port':    xo.Int64,
        'n_alloc':         xo.Int64,
        'timeout_sec':     xo.Int32,
        'particle_ref':    xp.Particles,
        'max_particle_id': xo.Int64
    }

    _extra_c_sources = [
        source
    ]

    # The engine is a singleton
    def __new__(cls, *args, **kwargs):
        if not hasattr(cls, 'instance'):
            cls.instance = super().__new__(cls)
            cls.instance._initialised = False
        return cls.instance

    def __del__(self, *args, **kwargs):
        self.stop_server(warn=False)

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
            self._tracking_init = False
            self._insertion = {}
            self._collimators = {}
            kwargs.setdefault('network_port', 0)
            kwargs.setdefault('n_alloc', 5000)
            kwargs.setdefault('timeout_sec', 36000) # 10 hours
            kwargs.setdefault('particle_ref', xp.Particles())
            kwargs.setdefault('max_particle_id', 0)

            # Get paths to executables
            if fluka is None:
                if testing:
                    self._fluka = "NEED_TO_MAKE_MOCKUP"
                else:
                    self._fluka = default_fluka_path
            elif testing:
                raise ValueError("Cannot specify fluka executable when `testing=True`!")
            else:
                self._fluka = Path(fluka).resolve()
            if flukaserver is None:
                if testing:
                    self._flukaserver = "NEED_TO_MAKE_MOCKUP"
                else:
                    self._flukaserver = default_flukaserver_path
            elif testing:
                raise ValueError("Cannot specify flukaserver executable when `testing=True`!")
            else:
                self._flukaserver = Path(flukaserver).resolve()
        super().__init__(**kwargs)


    def _warn_pyfluka(self, error):
        if not self._warning_given:
            print("Warning: Failed to import pyfluka (did you compile?). " \
                + "FlukaCollimators will be installed but are not trackable.\n", flush=True)
            print(error, flush=True)
            self._warning_given = True


    @classmethod
    def test_gfortran(cls, *args, **kwargs):
        cls(*args, **kwargs)
        this = cls.instance
        if not this._gfortran_installed:
            try:
                cmd = subprocess.run(["gfortran", "-dumpversion"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            except FileNotFoundError:
                this._gfortran_installed = False
                raise RuntimeError("Could not find gfortran installation! Need gfortran 9 or higher for fluka.")
            if cmd.returncode == 0:
                version = cmd.stdout.decode('UTF-8').strip().split('\n')[0]
                if int(version.split('.')[0]) < 9:
                    this._gfortran_installed = False
                    raise RuntimeError(f"Need gfortran 9 or higher for fluka, but found gfortran {version}!")
                this._gfortran_installed = True
                if this._verbose:
                    cmd2 = subprocess.run(["which", "gfortran"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                    if cmd2.returncode == 0:
                        file = cmd2.stdout.decode('UTF-8').strip().split('\n')[0]
                        print(f"Found gfortran {version} in {file}", flush=True)
                    else:
                        stderr = cmd2.stderr.decode('UTF-8').strip().split('\n')
                        raise RuntimeError(f"Could not run `which gfortran`!\nError given is:\n{stderr}")
            else:
                stderr = cmd.stderr.decode('UTF-8').strip().split('\n')
                this._gfortran_installed = False
                raise RuntimeError(f"Could not run gfortran! Verify its installation.\nError given is:\n{stderr}")


    @classmethod
    def start_server(cls, input_file, fluka_ids=None, *args, **kwargs):
        cls(*args, **kwargs)
        this = cls.instance
        if this.is_running():
            print("Server already running.", flush=True)
            return
        if not this._fluka.exists():
            raise ValueError(f"Could not find fluka executable {this._fluka}!")
        if not this._flukaserver.exists():
            raise ValueError(f"Could not find flukaserver executable {this._flukaserver}!")
        this.test_gfortran()
        try:
            from .pyflukaf import pyfluka_init
            pyfluka_init(n_alloc=this.n_alloc)
        except ImportError as error:
            this._warn_pyfluka(error)

        # Check files
        input_file = Path(input_file).resolve()
        if not input_file.exists():
            raise ValueError(f"Input file {input_file.as_posix()} not found!")
        this._cwd = input_file.parent
        insertion_file = this._cwd / "insertion.txt"
        if not insertion_file.exists():
            raise ValueError(f"Insertion file {insertion_file.as_posix()} not found!")

        # Start with cleaning files
        files_to_delete = [this._cwd / f for f in [network_file, fluka_log, server_log]]
        files_to_delete  = list(this._cwd.glob(f'ran{input_file.stem}*'))
        files_to_delete += list(this._cwd.glob(f'{input_file.stem}_*'))
        files_to_delete += list(this._cwd.glob(f'{input_file.stem}*.err'))
        files_to_delete += list(this._cwd.glob(f'{input_file.stem}*.log'))
        files_to_delete += list(this._cwd.glob(f'{input_file.stem}*.out'))
        files_to_delete += [this._cwd / f'fort.{n}' for n in [208, 251]]
        files_to_delete += [this._cwd / 'fluka_isotope.log']
        for f in files_to_delete:
            if f is not None:
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
                print("Network port not yet established. Waiting 30s.", flush=True)
            with this._network_nfo.open('r') as fid:
                lines = fid.readlines()
                if len(lines) > 1:
                    this.network_port = int(lines[1].strip())
                    break
        print(f"Started fluka server on network port {this.network_port}. "
            + f"Connecting (timeout: {this.timeout_sec})...", flush=True, end='')
        try:
            from .pyflukaf import pyfluka_connect
            pyfluka_connect(this.timeout_sec)
            this._flukaio_connected = True
        except ImportError as error:
            this._warn_pyfluka(error)
        print(f"done.", flush=True)

        # Store collimator info    TODO: need to get rid of fort.3 dependence etc
        if fluka_ids is None:
            this._read_fort3(path=input_file.parent)
        else:
            this._collimators = {coll: {'fluka_id': i} for coll,i in fluka_ids.items()}
        this._read_insertion(path=input_file.parent)
        this._read_gaps(input_file)
        # TODO: check for consistency (no repeating ids etc)


    @classmethod
    def stop_server(cls, warn=True, *args, **kwargs):
        cls(*args, **kwargs)
        this = cls.instance
        # Stop flukaio connection
        if this._flukaio_connected:
            print(f"Closing fluka server connection...", flush=True, end='')
            try:
                from .pyflukaf import pyfluka_close
                pyfluka_close()
                this._flukaio_connected = False
            except ImportError as error:
                if warn:
                    this._warn_pyfluka(error)
            print(f"done.", flush=True)
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
            print("Warning: Found other instances of rfluka besides the current one!", flush=True)
            return True
        else:
            print("Warning: Found other instances of rfluka but not the current one!", flush=True)
            this.stop_server()
            return False


    @property
    def collimators(self):
        return self._collimators


    def _read_insertion(self, path):
        if self.collimators == {}:
            raise ValueError("No FLUKA collimator id's defined yet! Do this first.")
        with open(path / 'insertion.txt', 'r') as fid:
            for line in fid.readlines():
                fluka_id = int(line.split()[0])
                name = [name for name, data in self._collimators.items() if data['fluka_id']==fluka_id]
                if len(name)==0:
                    continue
                name = name[0]
                self._collimators[name]['inactive_length'] = float(line.split()[-1])
                self._collimators[name]['INROT'] = line.split()[1:-1]


    def _read_fort3(self, path):
        fort3 = path / 'fort.3'
        if not fort3.exists():
            raise ValueError("No fort.3 found, cannot deduce FLUKA collimator ids!")
        with fort3.open('r') as fid:
            lines = fid.readlines()
        idx_start = [i for i, line in enumerate(lines) if line.strip().startswith('FLUKA')][0] + 1
        idx_end = [i for i, line in enumerate(lines[idx_start:]) if line.strip().startswith('NEXT')][0] + idx_start
        self._collimators = {}
        for line in lines[idx_start:idx_end]:
            line = [l.strip() for l in line.split()]
            if line[0].startswith('/') or line[0].startswith('!'):
                continue
            self._collimators[line[0]] = {
                'fluka_id': int(line[2]),
                'length': float(line[3]),
                'exit_marker': line[1]
            }


    def _read_gaps(self, input_file):
        # TODO: very hacky
        with input_file.open('r') as fid:
            lines = fid.readlines()
        for coll in self._collimators:
            try:
                idx = [i for i, line in enumerate(lines) if coll.upper() in line][0]
                hgap = float([line for line in lines[idx:] if line.startswith('* hGap')][0].split()[3])*1.e-3
            except:
                print(f"Warning: {coll} not found. Aperture set to 0!")
                hgap = 0.0
            self._collimators[coll]['jaw'] = hgap


    def _compare_fluka_mass(self, part, keep_p0c_constant=True):
        if part._capacity > 1:
            raise ValueError("`particle_ref` has to be a single particle!")
        pdg_id = part.pdg_id[0]
        if pdg_id is None or pdg_id == 0:
            raise ValueError("Need to set pdg_id of reference particle!")
        mass = part.mass0
        if pdg_id in masses:
            mass_fluka = masses[pdg_id][-1]
            if abs(mass-mass_fluka) > 1.:
                old_energy0 = part.energy0[0]
                part.mass0  = mass_fluka
                if keep_p0c_constant:
                    part._update_refs(p0c=part.p0c[0])
                else:
                    part._update_refs(energy0=old_energy0)
                print(f"Warning: given mass of {mass}eV for "
                    + f"{pdg.get_name_from_pdg_id(pdg_id)} differs from FLUKA "
                    + f"mass of {mass_fluka}eV. Overwritten by the latter.")
        else:
            print(f"Warning: No FLUKA reference mass known for particle "
                + f"{pdg.get_name_from_pdg_id(pdg_id)}!\nIf the refence mass "
                + f"provided ({mass}eV) differs from the one used internally "
                + f"by FLUKA, differences in energy might be observed.\nOnce "
                + f"the FLUKA reference mass is known, contact the devs to "
                + f"input it in the code.")


    def set_particle_ref(self, particle_ref=None, line=None):
        if self.has_particle_ref or self.max_particle_id!=0:
            print("Reference particle already set!")
            return

        if particle_ref is None:
            if line is None or line.particle_ref is None:
                raise ValueError("Line has no reference particle! "
                               + "Please provide one with `particle_ref`.")
            particle_ref = line.particle_ref
        else:
            if particle_ref._capacity > 1:
                raise ValueError("`particle_ref` has to be a single particle!")
            if line is not None:
                if line.particle_ref is not None \
                and not xt.line._dicts_equal(line.particle_ref.to_dict(), particle_ref.to_dict()):
                    print("Warning: Found different reference particle in line! Overwritten.")
                line.particle_ref = particle_ref
        self._compare_fluka_mass(particle_ref)
        self.particle_ref = particle_ref


    @property
    def has_particle_ref(self):
        initial = xp.Particles().to_dict()
        current = self.particle_ref.to_dict()
        return not xt.line._dicts_equal(initial, current)


    def init_tracking(self, max_particle_id):
        if self.max_particle_id == 0:
            particle_ref = self.particle_ref
            _, A0, Z0, name = pdg.get_properties_from_pdg_id(particle_ref.pdg_id[0])
            # FLUKA expects units of MeV
            E0 = particle_ref.energy0[0] / 1.e6
            p0c = particle_ref.p0c[0] / 1.e6
            m0 = particle_ref.mass0 / 1.e6
            q0 = particle_ref.q0
            try:
                from .pyflukaf import pyfluka_init_max_uid, pyfluka_set_synch_part
            except ImportError as error:
                self._warn_pyfluka(error)
                return
            pyfluka_init_max_uid(max_particle_id)
            pyfluka_set_synch_part(E0, p0c, m0, A0, Z0, q0)
            print(f"Set max_particle_id to {max_particle_id}, "
                + f"and reference particle to {name} with mass {m0}MeV and energy {E0}MeV.")
            self.max_particle_id = max_particle_id


