# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

import numpy as np
from subprocess import run, PIPE, Popen
from time import sleep

import xobjects as xo
import xpart as xp
import xpart.pdg as pdg
import xtrack as xt

try:
    from xaux import ClassProperty, FsPath  # TODO: once xaux is in Xsuite keep only this
except ImportError:
    from ...xaux import ClassProperty, FsPath

from .reference_masses import source, masses
from .paths import flukafile_resolve
from .paths import fluka as default_fluka_path
from .paths import flukaserver as default_flukaserver_path
from ..engine import BaseEngine
from ...general import _pkg_root


network_file = "network.nfo"
fluka_log  = "fluka.log"
server_log = "rfluka.log"


class FlukaEngine(BaseEngine):
    _xofields = {**BaseEngine._xofields,
        '_network_port':    xo.Int64,
        '_timeout_sec':     xo.Int32,
        '_max_particle_id': xo.Int64
    }

    _int32 = True
    _uses_input_file = True
    _uses_run_folder = True

    _extra_c_sources = [source]

    def __init__(self, **kwargs):
        if not self._initialised and '_xobject' not in kwargs:
            # Set element classes dynamically
            from ...beam_elements import FlukaCollimator
            self.__class__._element_classes = (FlukaCollimator,)
            # Initialise fluka-only defaults
            self._network_nfo = None
            self._log = None
            self._log_fid = None
            self._server_process = None
            self.server_pid = None
            self._gfortran_installed = False
            self._flukaio_connected = False
            kwargs.setdefault('_network_port', 0)
            kwargs.setdefault('_timeout_sec', 36000) # 10 hours
            kwargs.setdefault('_max_particle_id', 0)
            self._fluka = None
            self._flukaserver = None
            kwargs.setdefault('fluka', default_fluka_path)
            kwargs.setdefault('flukaserver', default_flukaserver_path)
        super().__init__(**kwargs)

    def __del__(self, *args, **kwargs):
        self.stop(warn=False)
        try:
            super().__del__()
        except AttributeError:
            pass

    def _warn_pyfluka(self, error):
        if not self._warning_given:
            print("Warning: Failed to import pyfluka (did you compile?). " \
                + "FlukaCollimators will be installed but are not trackable.\n", flush=True)
            print(error, flush=True)
            self._warning_given = True


    # ======================
    # === New Properties ===
    # ======================

    @ClassProperty
    def network_port(cls):
        return cls.get_self()._network_port

    @ClassProperty
    def timeout_sec(cls):
        return cls.get_self()._timeout_sec

    @timeout_sec.setter
    def timeout_sec(cls, val):
        cls.get_self()._timeout_sec = val

    @ClassProperty
    def max_particle_id(cls):
        return cls.get_self()._max_particle_id

    @ClassProperty
    def fluka(cls):
        return cls.get_self()._fluka

    @fluka.setter
    def fluka(cls, val):
        self = cls.get_self()
        self._fluka = FsPath(val).expanduser().resolve()
        _fluka = flukafile_resolve(self._fluka)
        if _fluka is not None:
            self._fluka = _fluka
        else:
            raise ValueError(f"Could not find fluka executable {self._fluka}!")

    @ClassProperty
    def flukaserver(cls):
        return cls.get_self()._flukaserver

    @flukaserver.setter
    def flukaserver(cls, val):
        self = cls.get_self()
        self._flukaserver = FsPath(val).expanduser().resolve()
        _flukaserver = flukafile_resolve(self._flukaserver)
        if _flukaserver is not None:
            self._flukaserver = _flukaserver
        else:
            raise ValueError(f"Could not find fluka executable {self._flukaserver}!")


    # ======================
    # === Public Methods ===
    # ======================

    @classmethod
    def test_gfortran(cls, **kwargs):
        self = cls.get_self(**kwargs)
        if not self._gfortran_installed:
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
                if self.verbose:
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
    def generate_input_file(cls, *, line=None, elements=None, names=None, prototypes_file=None,
                            include_files=None, **kwargs):
        self = cls.get_self(**kwargs)
        if self._element_dict == {}:
            # This is for the case that the method is not called within start() but manually in advance
            self._get_elements(line=line, elements=elements, names=names)
            element_dict = self._element_dict
            self._element_dict = {}
        else:
            element_dict = self._element_dict
        # TODO: make prototypes file based on element_dict
        if prototypes_file is None:
            if self.verbose:
                print("Using default prototypes file.")
            prototypes_file = _pkg_root / 'scattering_routines' / 'fluka' / 'data' / 'prototypes.lbp'
        # TODO: make include files based on beam
        if include_files is None:
            if self.verbose:
                print("Using default include files.")
            include_files = [
                _pkg_root / 'scattering_routines' / 'fluka' / 'data' / 'include_settings_beam.inp',
                _pkg_root / 'scattering_routines' / 'fluka' / 'data' / 'include_settings_physics.inp',
                _pkg_root / 'scattering_routines' / 'fluka' / 'data' / 'include_custom_scoring.inp'
            ]

        from .fluka_input import create_fluka_input
        return create_fluka_input(element_dict=element_dict, prototypes_file=prototypes_file,
                                  include_files=include_files, verbose=self.verbose)


    @classmethod
    def start(cls, *, line=None, elements=None, names=None, cwd=None, seed=None, particle_ref=None,
               input_file=None, prototypes_file=None, include_files=None, touches=True, **kwargs):
        self = cls.get_self(**kwargs)
        cls.test_gfortran()
        if self.verbose:
            print("Starting FLUKA engine...", flush=True)
        super().start(line=line, elements=elements, names=names, cwd=cwd, seed=seed,
                      particle_ref=particle_ref, input_file=input_file,
                      prototypes_file=prototypes_file, include_files=include_files, **kwargs)

        from .fluka_input import verify_insertion_file
        verify_insertion_file(self.input_file[1], self._element_dict)

        self._create_touches(touches)

        try:
            from .pyflukaf import pyfluka_init
            debug_level = 2 if self.verbose else 0
            pyfluka_init(n_alloc=self._capacity, debug_level=debug_level)
        except ImportError as error:
            self._warn_pyfluka(error)
            self.stop()
            return

        self._declare_network()
        self._start_server()

        try:
            from .pyflukaf import pyfluka_connect
            pyfluka_connect(self.timeout_sec)
            self._flukaio_connected = True
        except ImportError as error:
            self._warn_pyfluka(error)
            self.stop()
            return

        if self.verbose:
            print(f"done.", flush=True)


    @classmethod
    def stop(cls, clean=False, **kwargs):
        self = cls.get_self(**kwargs)
        # Stop flukaio connection
        if self._flukaio_connected:
            if self.verbose:
                print(f"Closing fluka server connection... ", flush=True)
            try:
                from .pyflukaf import pyfluka_close
                pyfluka_close()
            except ImportError as error:
                self._warn_pyfluka(error)
            self._flukaio_connected = False
            if self.verbose:
                print(f"done.", flush=True)
        # If the Popen process is still running, terminate it
        if self._server_process is not None:
            while self._server_process.poll() is None:
                sleep(1)
            self._server_process.terminate()
            self._server_process = None
        self.server_pid = None
        # Close the file pointer to the log
        if self._log_fid is not None:
            if not self._log_fid.closed:
                self._log_fid.close()
            self._log_fid = None
        self._log = None
        # Delete network file
        if self._network_nfo is not None and self._network_nfo.exists():
            self._network_nfo.unlink()
            self._network_nfo = None
        self._network_port = 0
        self._max_particle_id = 0
        super().stop(clean=clean, **kwargs)
#        if this._old_cwd is not None:
#            os.chdir(this._old_cwd)
#            this._old_cwd = None
#        if clean:
#            this.clean_output_files(clean_all=True)
#        # Unassign the prototypes
#        for name, ee in this._collimator_dict.items():
#            ee.assembly.remove_element(name, force=False)
#        this._cwd = None
#        this._network_nfo = None
#        this._log = None
#        this._log_fid = None
#        this._server_process = None
#        this._input_file = None
#        this.server_pid = None
#        this._gfortran_installed = False
#        this._flukaio_connected = False
#        this._warning_given = False
#        this._tracking_initialised = False
#        this._insertion = {}
#        this._collimator_dict = {}
#        this.network_port = 0
#        this.max_particle_id =  0
#        this.seed = -1


    @classmethod
    def is_running(cls, **kwargs):
        base_is_running = super().is_running(**kwargs)
        if not base_is_running is None:
            return base_is_running
        self = cls.get_self(**kwargs)
        # Is the Popen process still running?
        if self._server_process is None or self._server_process.poll() is not None:
            self.stop()
            return False
        # Get username (need a whoami for the next command)
        cmd = run(["whoami"], stdout=PIPE, stderr=PIPE)
        if cmd.returncode == 0:
            whoami = cmd.stdout.decode('UTF-8').strip().split('\n')[0]
        else:
            stderr = cmd.stderr.decode('UTF-8').strip().split('\n')
            raise RuntimeError(f"Could not find username! Error given is:\n{stderr}")
        # Get fluka processes for this user
        cmd = run(["ps", "-u", whoami], stdout=PIPE, stderr=PIPE)
        if cmd.returncode == 0:
            processes = cmd.stdout.decode('UTF-8').strip().split('\n')
        else:
            stderr = cmd.stderr.decode('UTF-8').strip().split('\n')
            raise RuntimeError(f"Could not list running processes! Error given is:\n{stderr}")
        processes = [proc for proc in processes if 'rfluka' in proc]
        if len(processes) == 0:
            # Could not find a running rfluka
            self.stop()
            return False
        elif len(processes) == 1 and str(self.server_pid) in processes[0] and 'defunct' not in processes[0]:
            return True
        elif np.any([str(self.server_pid) in proc for proc in processes if 'defunct' not in proc]):
            print("Warning: Found other instances of rfluka besides the current one!", flush=True)
            return True
        else:
            print("Warning: Found other instances of rfluka but not the current one!", flush=True)
            self.stop()
            return False


    @classmethod
    def clean_output_files(cls, input_file=None, cwd=None, clean_all=False, **kwargs):
        self = cls.get_self(**kwargs)
        if input_file is None:
            input_file = self.input_file[0]
            if cwd is None and self._cwd is not None:
                cwd = self._cwd
        else:
            if hasattr(input_file, '__iter__'):
                if len(input_file) == 1:
                    input_file = input_file[0]
                else:
                    raise ValueError("Only one input file can be given!")
            input_file = FsPath(input_file)
            if cwd is None:
                cwd = input_file.parent
# TODO check this logic
        if isinstance(cwd, FsPath):
            files_to_delete = [cwd / f for f in [network_file, fluka_log, server_log, 'fluka_isotope.log',
                                                 'fort.208', 'fort.251']]
            if input_file is not None:
                files_to_delete += list(cwd.glob(f'ran{input_file.stem}*'))
                files_to_delete += list(cwd.glob(f'{input_file.stem}*'))
            if not clean_all:
                files_to_delete  = [file for file in files_to_delete if FsPath(file).resolve() != input_file.resolve()]
            if clean_all:
                files_to_delete += [cwd / f for f in ['insertion.txt', 'new_collgaps.dat', 'relcol.dat']]
            for f in files_to_delete:
                if f is not None and f.exists():
                    f.unlink()
            if clean_all and cwd != FsPath.cwd():
                cwd.rmdir()


    @classmethod
    def init_tracking(cls, max_particle_id, **kwargs):
        self = cls.get_self(**kwargs)
        if self._tracking_initialised == False:
            self.assert_particle_ref()
            _, A0, Z0, name = pdg.get_properties_from_pdg_id(self.particle_ref.pdg_id[0])
            # FLUKA expects units of MeV
            E0 = self.particle_ref.energy0[0] / 1.e6
            p0c = self.particle_ref.p0c[0] / 1.e6
            m0 = self.particle_ref.mass0 / 1.e6
            q0 = self.particle_ref.q0
            try:
                from .pyflukaf import pyfluka_init_max_uid, pyfluka_set_synch_part
            except ImportError as error:
                self._warn_pyfluka(error)
                return
            if self.verbose:
                print(f"Setting max_particle_id to {max_particle_id}, "
                    + f"and reference particle to {name} with mass {m0}MeV "
                    + f"and energy {E0}MeV.", flush=True)
            pyfluka_init_max_uid(max_particle_id)
            pyfluka_set_synch_part(E0, p0c, m0, A0, Z0, q0)
            self._max_particle_id = max_particle_id
            self._tracking_initialised = True


    # =======================
    # === Private Methods ===
    # =======================


    # Expand the Base method to include the FLUKA reference mass check
    def _use_particle_ref(self, particle_ref=None, keep_p0c_constant=True):
        super()._use_particle_ref(particle_ref=particle_ref)
        part = self.particle_ref
        mass = part.mass0
        pdg_id = part.pdg_id[0]
        if pdg_id in masses:
            mass_fluka = masses[pdg_id][-1]
            if abs(mass-mass_fluka) > 1.:    # The mass differs more than 1eV from the FLUKA reference mass
                old_energy0 = part.energy0[0]
                part.mass0  = mass_fluka
                if keep_p0c_constant:
                    part._update_refs(p0c=part.p0c[0])
                else:
                    part._update_refs(energy0=old_energy0)
                assert np.isclose(part.energy0[0]**2, part.p0c[0]**2 + part.mass0**2)
                assert np.isclose(part.mass0, mass_fluka)
                print(f"Warning: given mass of {mass}eV for "
                    + f"{pdg.get_name_from_pdg_id(pdg_id)} differs from FLUKA "
                    + f"mass of {mass_fluka}eV.\nReference particle mass is "
                    + f"overwritten by the latter.")
        else:
            print(f"Warning: No FLUKA reference mass known for particle "
                + f"{pdg.get_name_from_pdg_id(pdg_id)}!\nIf the reference mass "
                + f"provided ({mass}eV) differs from the one used internally "
                + f"by FLUKA, differences in energy might be observed.\nOnce "
                + f"the FLUKA reference mass is known, contact the devs to "
                + f"input it in the code.")


    def _match_input_file(self, *, line=None, elements=None, names=None):
        # Read the elements in the input file and compare to the elements in the engine,
        # overwriting parameters where necessary
        from .fluka_input import get_collimators_from_input_file
        input_dict = get_collimators_from_input_file(self.input_file[0])
        for name in input_dict:
            if name not in self._element_dict:
                raise ValueError(f"Element {name} in input file not found in engine!")
        for name, el in self._element_dict.items():
            if name not in input_dict:
                print(f"Warning: FlukaCollimator {name} not in FLUKA input file! "
                    + f"Maybe it was fully open. Deactivated")
                el.active = False
                continue
            from ...beam_elements import FlukaCollimator
            if not isinstance(el, FlukaCollimator):
                raise ValueError("Element {name} is not a FLUKA element!")
            el.fluka_id = input_dict[name]['fluka_id']
            el.length_front = (input_dict[name]['length'] - el.length)/2
            el.length_back = (input_dict[name]['length'] - el.length)/2
            jaw = input_dict[name]['jaw']
            if not hasattr(jaw, '__iter__'):
                jaw = [jaw, -jaw]
            if jaw[0] is None and jaw[1] is None:
                el.jaw = None
            else:
                if jaw[0] is None:
                    if el.side != 'right':
                        print(f"Warning: {name} is right-sided in the input file, but not "
                            + "in the line! Overwritten by the former.")
                        el.side = 'right'
                elif el.jaw_L is None or not np.isclose(el.jaw_L, jaw[0], atol=1e-9):
                    print(f"Warning: Jaw_L of {name} differs from input file ({el.jaw_L} "
                        + f"vs {jaw[0]})! Overwritten.")
                    el.jaw_L = jaw[0]
                if jaw[1] is None:
                    if el.side != 'left':
                        print(f"Warning: {name} is left-sided in the input file, but not "
                            + "in the line! Overwritten by the former.")
                        el.side = 'left'
                elif el.jaw_R is None or not np.isclose(el.jaw_R, jaw[1], atol=1e-9):
                    print(f"Warning: Jaw_R of {name} differs from input file ({el.jaw_R} "
                        + f"vs {jaw[1]})! Overwritten.")
                    el.jaw_R = jaw[1]
            self._collimator_dict[name] = el
# TODO: check if self._collimator_dict is correct
        # TODO: tilts!!


    def _create_touches(self, touches):
        # Create touches file (relcol.dat)
        # First line is the number of collimators, second line is the IDs (no newline at end)
        if touches is True:
            touches = FsPath('relcol.dat').resolve()
            with touches.open('w') as fid:
                fid.write(f'{len(self._element_dict.keys())}\n')
                for _, el in self._element_dict.items():
                    fid.write(f'{el.fluka_id} ')
            self._input_file.append(touches)
        # Check if touches is a list of collimator names
        elif touches is not None and isinstance(touches, list):
            relcol = FsPath('relcol.dat').resolve()
            with relcol.open('w') as fid:
                fid.write(f'{len(touches)}\n')
                for touch in touches:
                    if touch not in self._element_dict:
                        raise ValueError(f"Collimator {touch} not in collimator dict!")
                    else:
                        fid.write(f'{self._element_dict[touch].fluka_id} ')
            self._input_file.append(relcol)
        # Check if touches is not wrongly set
        elif touches is not None and not touches is False:
            raise NotImplementedError("Only True/False or a list of collimstors is allowed "
                                    + "for `touches` for now.")


    def _declare_network(self):
        self._network_nfo = self._cwd / network_file
        cmd = run(["hostname"], cwd=self._cwd, stdout=PIPE, stderr=PIPE)
        if cmd.returncode == 0:
            host = cmd.stdout.decode('UTF-8').strip().split('\n')[0]
        else:
            stderr = cmd.stderr.decode('UTF-8').strip().split('\n')
            raise RuntimeError(f"Could not declare hostname! Error given is:\n{stderr}")
        with self._network_nfo.open('w') as fid:
            fid.write(f"{host}\n")


    def _start_server(self):
        self._log = self._cwd / server_log
        self._log_fid = self._log.open('w')
        self._server_process = Popen([self.fluka.as_posix(), self.input_file[0], '-e', self.flukaserver.as_posix(), '-M', "1"], 
                                    cwd=self._cwd, stdout=self._log_fid, stderr=self._log_fid)
        self.server_pid = self._server_process.pid
        sleep(1)
        if not self.is_running():
            raise RuntimeError(f"Could not start fluka server! See logfile {self._log}.")
        i = 0
        while True:
            sleep(2)
            i += 2
            if i%30 == 0 and self.verbose:
                print("Network port not yet established (or missing executable permission on "
                    + "flukaserver). Waiting 30s.", flush=True)
            with self._network_nfo.open('r') as fid:
                lines = fid.readlines()
                if len(lines) > 1:
                    self._network_port = int(lines[1].strip())
                    break
                if self._server_process.poll() is not None:
                    raise RuntimeError("The flukaserver died. Please check the FLUKA output.")
        if self.verbose:
            print(f"Started fluka server on network port {self.network_port}. "
                + f"Connecting (timeout: {self.timeout_sec})... ", flush=True, end='')
