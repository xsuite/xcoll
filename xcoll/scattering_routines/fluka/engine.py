# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

import numpy as np
from numbers import Number
from subprocess import run, PIPE, Popen
from time import sleep
 
import xobjects as xo
import xtrack.particles.pdg as pdg

try:
    from xaux import FsPath  # TODO: once xaux is in Xsuite keep only this
except (ImportError, ModuleNotFoundError):
    from ...xaux import FsPath

from .reference_masses import source, fluka_masses
from .environment import format_fluka_float
from .prototype import FlukaPrototype, FlukaAssembly
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
    _num_input_files = 3
    _uses_run_folder = True

    _extra_c_sources = [source]

    def __init__(self, **kwargs):
        # Set element classes dynamically
        from ...beam_elements import FlukaCollimator, FlukaCrystal
        self.__class__._element_classes = (FlukaCollimator, FlukaCrystal)
        # Initialise fluka-only defaults
        self._log = None
        self._log_fid = None
        self.server_pid = None
        self._network_nfo = None
        self._server_process = None
        self._flukaio_connected = False
        kwargs.setdefault('_timeout_sec', 36000)
        kwargs.setdefault('_network_port', -1)
        kwargs.setdefault('_max_particle_id', -1)
        super().__init__(**kwargs)


    # ======================
    # === New Properties ===
    # ======================

    @property
    def network_port(self):
        return self._network_port

    @property
    def timeout_sec(self):
        return self._timeout_sec

    @timeout_sec.setter
    def timeout_sec(self, val):
        if val is None:
            val = 36000
        if not isinstance(Number) or val <= 60:
            raise ValueError("`capacity` has to be an integer and larger than 60!")
        self._timeout_sec = val

    @property
    def max_particle_id(self):
        return self._max_particle_id


    # ============================
    # === Overwrite Properties ===
    # ============================

    @property
    def seed(self):
        return BaseEngine.seed.fget(self)

    @seed.setter
    def seed(self, val):
        BaseEngine.seed.fset(self, val)
        if self._seed >= 9.e8:   # FLUKA only accepts 8 digits
            self._seed = int(self._seed / 4.77219)

    @seed.deleter
    def seed(self):
        return BaseEngine.seed.fdel(self)


    # ======================
    # === Public methods ===
    # ======================

    def init_tracking(self, max_particle_id):
        if self._tracking_initialised == False:
            self.assert_particle_ref()
            _, A0, Z0, name = pdg.get_properties_from_pdg_id(self.particle_ref.pdg_id[0], long_name=True)
            # FLUKA expects units of MeV
            E0 = self.particle_ref.energy0[0] / 1.e6
            p0c = self.particle_ref.p0c[0] / 1.e6
            m0 = self.particle_ref.mass0 / 1.e6
            q0 = abs(self.particle_ref.q0)
            try:
                from pyflukaf import pyfluka_init_max_uid, pyfluka_set_synch_part
            except (ModuleNotFoundError, ImportError) as error:
                self._warn(error)
            self._print(f"Setting max_particle_id to {max_particle_id}, "
                      + f"and reference particle to {name} with mass {m0} MeV "
                      + f"and momentum {p0c} MeV.")
            pyfluka_init_max_uid(max_particle_id)
            pyfluka_set_synch_part(E0, p0c, m0, A0, Z0, q0)
            self._max_particle_id = max_particle_id - 1 # FLUKA starts counting from 1, Python from 0
            self._tracking_initialised = True


    def view(self):
        if self.input_file is None:
            return
        xc.fluka.environment.run_flair(self.input_file[0])


    # =================================
    # === Base methods to overwrite ===
    # =================================

    def _set_engine_properties(self, **kwargs):
        kwargs = super()._set_engine_properties(**kwargs)
        self._set_property('timeout_sec', kwargs)
        return kwargs

    def _generate_input_file(self, *, prototypes_file=None, include_files=[], **kwargs):
        from .fluka_input import create_fluka_input
        self._deactivate_unused_assemblies()
        input_file, kwargs = create_fluka_input(element_dict=self._element_dict,
                                particle_ref=self.particle_ref, prototypes_file=prototypes_file,
                                include_files=include_files, verbose=self.verbose, **kwargs)
        self._set_seed_in_input_file(input_file)
        # The only thing left in kwargs are parameters to start the engine
        return input_file, kwargs


    def _pre_start(self, **kwargs):
        import xcoll as xc
        xc.fluka.environment.assert_gfortran_installed()
        xc.fluka.environment.set_fluka_environment()
        return kwargs


    def _start_engine(self, touches=True, fortran_debug_level=0, **kwargs):
        from .fluka_input import verify_insertion_file
        verify_insertion_file(self.input_file[1], self._element_dict)
        self._create_touches(touches)
        self._init_fortran(fortran_debug_level)
        self._declare_network()
        self._start_server()


    def _stop_engine(self, **kwargs):
        self._stop_fortran()
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
        self._network_port = -1
        self._max_particle_id = -1
        return kwargs


    def _is_running(self):
        # Is the Popen process still running?
        if self._server_process is None:
            return False
        elif self._server_process.poll() is not None:
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
            self._print("Warning: Found other instances of rfluka besides the current one!")
            return True
        else:
            self._print("Warning: Found other instances of rfluka but not the current one!")
            self.stop()
            return False


    def _match_input_file(self):
        # Read the elements in the input file and compare to the elements in the engine,
        # overwriting parameters where necessary
        from .fluka_input import get_collimators_from_input_file
        input_dict = get_collimators_from_input_file(self.input_file[0])
        for name in input_dict:
            if name not in self._element_dict:
                raise ValueError(f"Element {name} in input file not found in engine!")
        for name, ee in self._element_dict.items():
            if name not in input_dict:
                self._print(f"Warning: FlukaCollimator {name} not in FLUKA input file! "
                          + f"Maybe it was fully open. Deactivated")
                self._deactivate_element(ee)
                continue
            self._assert_element(ee)
            ee.fluka_id = input_dict[name]['fluka_id']
            ee.length_front = (input_dict[name]['length'] - ee.length)/2
            ee.length_back = (input_dict[name]['length'] - ee.length)/2
            jaw = input_dict[name]['jaw']
            if jaw is not None and not hasattr(jaw, '__iter__'):
                jaw = [jaw, -jaw]
            if jaw is None or (jaw[0] is None and jaw[1] is None):
                ee.jaw = None
            else:
                if ee.assembly.is_crystal:
                    if jaw[0] is None:
                        if ee.side != 'right':
                            self._print(f"Warning: {name} is right-sided in the input file, "
                                    + "but not in the line! Overwritten by the former.")
                            ee.side = 'right'
                    elif not np.isclose(ee.jaw, jaw[0], atol=1e-9):
                        self._print(f"Warning: Jaw of {name} differs from input file "
                                + f"({ee.jaw} vs {jaw[0]})! Overwritten.")
                        ee.jaw = jaw[0]
                    if jaw[1] is None:
                        if ee.side != 'left':
                            self._print(f"Warning: {name} is left-sided in the input file, "
                                    + "but not in the line! Overwritten by the former.")
                            ee.side = 'left'
                    elif not np.isclose(ee.jaw, jaw[1], atol=1e-9):
                        self._print(f"Warning: Jaw of {name} differs from input file "
                                + f"({ee.jaw} vs {jaw[1]})! Overwritten.")
                        ee.jaw = jaw[1]
                else:
                    if jaw[0] is None:
                        if ee.side != 'right':
                            self._print(f"Warning: {name} is right-sided in the input file, "
                                    + "but not in the line! Overwritten by the former.")
                            ee.side = 'right'
                    elif ee.jaw_L is None or not np.isclose(ee.jaw_L, jaw[0], atol=1e-9):
                        self._print(f"Warning: Jaw_L of {name} differs from input file "
                                + f"({ee.jaw_L} vs {jaw[0]})! Overwritten.")
                        ee.jaw_L = jaw[0]
                    if jaw[1] is None:
                        if ee.side != 'left':
                            self._print(f"Warning: {name} is left-sided in the input file, "
                                    + f"but not in the line! Overwritten by the former.")
                            ee.side = 'left'
                    elif ee.jaw_R is None or not np.isclose(ee.jaw_R, jaw[1], atol=1e-9):
                        self._print(f"Warning: Jaw_R of {name} differs from input file "
                                + f"({ee.jaw_R} vs {jaw[1]})! Overwritten.")
                        ee.jaw_R = jaw[1]
        # TODO: tilts!!


    def _get_input_files_to_clean(self, input_file, cwd, **kwargs):
        if cwd is not None:
            files_to_delete = ['prototypes.lbp', 'assignmat.inp', 'linebuilder.log',
                               'new_collgaps.dat']
            files_to_delete = [cwd / f for f in files_to_delete]
            files_to_delete += list(cwd.glob(f'include_*.inp'))
            if input_file is not None:
                if not hasattr(input_file, '__iter__') or isinstance(input_file, str):
                    input_file = [input_file]
                files_to_delete += list(cwd.glob(f'{input_file[0].stem}*'))
                files_to_delete += self._all_input_files(input_file)
            return files_to_delete
        else:
            return []

    def _get_output_files_to_clean(self, input_file, cwd, **kwargs):
        if cwd is not None:
            files_to_delete = [network_file, fluka_log, server_log,
                            'fluka_isotope.log', 'fort.208', 'fort.251']
            files_to_delete = [cwd / f for f in files_to_delete]
            if input_file is not None:
                if not hasattr(input_file, '__iter__') or isinstance(input_file, str):
                    input_file = [input_file]
                files_to_delete += list(cwd.glob(f'ran{input_file[0].stem}*'))
                files_to_delete += list(cwd.glob(f'{input_file[0].stem}*'))

            # Do not delete the extra files generated with the input file (they are deleted with clean_input_files)
            _input_files = self._clean_input_files(input_file, cwd, clean_all=False)
            files_to_delete  = [file for file in files_to_delete
                                if FsPath(file).resolve() not in _input_files]

            if input_file is not None:
                # Do not delete the input file itself
                files_to_delete  = [file for file in files_to_delete
                                    if FsPath(file).resolve() != input_file[0].resolve()]
            if input_file is not None:
                files_to_delete += self._all_input_files(input_file)
            return files_to_delete
        else:
            return []

    def _all_input_files(self, input_file):
        if not hasattr(input_file, '__iter__') or isinstance(input_file, str):
            input_file = [input_file]
        for ff in ['insertion.txt','relcol.dat']:
            if ff not in [fff.name for fff in input_file]:
                input_file.append(input_file[0].parent / ff)
        return input_file


    # Expand the Base method to include the FLUKA reference mass check
    def _use_particle_ref(self, particle_ref=None, keep_p0c_constant=True):
        super()._use_particle_ref(particle_ref=particle_ref)
        part = self.particle_ref
        mass = part.mass0
        pdg_id = part.pdg_id[0]
        if abs(pdg_id) in fluka_masses:
            mass_fluka = fluka_masses[abs(pdg_id)][-1]
            if abs(mass-mass_fluka) > 0.01:    # The mass differs more than 0.01eV from the FLUKA reference mass
                old_energy0 = part.energy0[0]
                part.mass0  = mass_fluka
                if keep_p0c_constant:
                    part._update_refs(p0c=part.p0c[0])
                else:
                    part._update_refs(energy0=old_energy0)
                assert np.isclose(part.energy0[0]**2, part.p0c[0]**2 + part.mass0**2)
                assert np.isclose(part.mass0, mass_fluka)
                self._print(f"Warning: given mass of {mass} eV for "
                          + f"{pdg.get_name_from_pdg_id(pdg_id)} differs from FLUKA "
                          + f"mass of {mass_fluka} eV.\nReference particle mass is "
                          + f"overwritten by the latter.")
        else:
            self._print(f"Warning: No FLUKA reference mass known for particle "
                      + f"{pdg.get_name_from_pdg_id(pdg_id)}!\nIf the reference mass "
                      + f"provided ({mass} eV) differs from the one used internally "
                      + f"by FLUKA, differences in energy might be observed.\nOnce "
                      + f"the FLUKA reference mass is known, contact the devs to "
                      + f"input it in the code.")


    # =======================
    # === Private Methods ===
    # =======================

    def _deactivate_unused_assemblies(self):
        # Deactivate all assemblies that are not used in the input file
        for prototype in [*FlukaPrototype._assigned_registry.values(),
                          *FlukaAssembly._assigned_registry.values()]:
            for ee in prototype.elements:
                if ee.name not in self._element_dict:
                    self._deactivate_element(ee)

    def _init_fortran(self, fortran_debug_level=0):
        try:
            from pyflukaf import pyfluka_init
            pyfluka_init(n_alloc=self._capacity, debug_level=fortran_debug_level)
        except (ModuleNotFoundError, ImportError) as error:
            self._warn(error)


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
        import xcoll as xc
        log = self._cwd / server_log
        self._log = log
        self._log_fid = self._log.open('w')
        self._server_process = Popen([xc.fluka.environment.fluka.as_posix(),
                                      self.input_file[0].as_posix(), '-e',
                                      xc.fluka.environment.flukaserver.as_posix(), '-M', "1"],
                                     cwd=self._cwd, stdout=self._log_fid, stderr=self._log_fid)
        self.server_pid = self._server_process.pid
        sleep(1)
        if not self.is_running():
            raise RuntimeError(f"Could not start fluka server! See logfile {log}.")
        i = 0
        while True:
            sleep(2)
            i += 2
            if i%30 == 0:
                self._print("Network port not yet established (or missing executable "
                          + "permission on flukaserver). Waiting 30s.")
            with self._network_nfo.open('r') as fid:
                lines = fid.readlines()
                if len(lines) > 1:
                    self._network_port = int(lines[1].strip())
                    break
                if self._server_process.poll() is not None:
                    self.stop()
                    raise RuntimeError("The flukaserver died. Please check the FLUKA output.")
        self._print(f"Started fluka server on network port {self.network_port}. "
                  + f"Connecting (timeout: {self.timeout_sec})...   ", end='')
        try:
            from pyflukaf import pyfluka_connect
            pyfluka_connect(self.timeout_sec)
            self._flukaio_connected = True
        except (ModuleNotFoundError, ImportError) as error:
            self._warn(error)
        self._print(f"Done.")


    def _stop_fortran(self):
        if self._flukaio_connected:
            self._flukaio_connected = False
            self._print(f"Closing fluka server connection...   ", end='')
            try:
                from pyflukaf import pyfluka_close
                pyfluka_close()
            except (ModuleNotFoundError, ImportError) as error:
                self._warn(error)
            self._print(f"Done.")


    def _set_seed_in_input_file(self, input_file):
        input_file = FsPath(input_file[0])
        with input_file.open('r') as fid:
            lines = fid.readlines()
        for i, line in enumerate(lines):
            if 'RANDOMIZ' in line:
                lines[i] = f"RANDOMIZe        1.0{format_fluka_float(self.seed)}\n"
                break
        with input_file.open('w') as fid:
            fid.writelines(lines)


    def _create_touches(self, touches=None):
        # Create touches file (relcol.dat)
        # First line is the number of collimators, second line is the IDs (no newline at end)
        if touches is True:
            touches = list(self._element_dict.keys())
        # Check if touches is a list of collimator names
        if touches is not None and hasattr(touches, '__iter__') \
        and not isinstance(touches, str):
            relcol = FsPath('relcol.dat').resolve()
            with relcol.open('w') as fid:
                fid.write(f'{len(touches)}\n')
                for touch in touches:
                    if touch not in self._element_dict:
                        raise ValueError(f"Collimator {touch} not in collimator dict, "
                                        + "but asked to write FLUKA touches!")
                    else:
                        fid.write(f'{self._element_dict[touch].fluka_id} ')
            self._input_file.append(relcol)
        # Check if touches is not wrongly set
        elif touches is not None and not touches is False:
            raise NotImplementedError("Only True/False or a list of collimator names "
                                    + "is allowed for `touches` for now.")
