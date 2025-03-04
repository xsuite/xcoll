# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

import os
import time
import numpy as np
from subprocess import run, PIPE, Popen
from pathlib import Path
from time import sleep
import shutil

import xobjects as xo
import xpart as xp
import xpart.pdg as pdg
import xtrack as xt

from .reference_masses import source, masses
from .paths import flukafile_resolve
from .paths import fluka as default_fluka_path
from .paths import flukaserver as default_flukaserver_path
from ...general import _pkg_root


network_file = "network.nfo"
fluka_log  = "fluka.log"
server_log = "rfluka.log"


class FlukaEngine(xo.HybridClass):
    _xofields = {
        'network_port':    xo.Int64,
        '_capacity':       xo.Int64,
        'timeout_sec':     xo.Int32,
        'particle_ref':    xp.Particles,
        'max_particle_id': xo.Int64,
        'seed':            xo.Int64
    }

    _extra_c_sources = [
        source
    ]

    # The engine is a singleton
    def __new__(cls, **kwargs):
        if not hasattr(cls, 'instance'):
            cls.instance = super().__new__(cls)
            cls.instance._initialised = False
        return cls.instance

    def __del__(self, *args, **kwargs):
        self.stop(warn=False)
        try:
            super().__del__()
        except AttributeError:
            pass

    def __init__(self, fluka=None, flukaserver=None, verbose=None, testing=False, **kwargs):
        # Apply init variables
        if fluka is None:
            if not hasattr(self, '_fluka'):
                self._fluka = "NEED_TO_MAKE_MOCKUP"  if testing else default_fluka_path  # TODO
        elif testing:
            raise ValueError("Cannot specify fluka executable when `testing=True`!")
        else:
            self._fluka = Path(fluka).resolve()
        if flukaserver is None:
            if not hasattr(self, '_flukaserver'):
                self._flukaserver = "NEED_TO_MAKE_MOCKUP" if testing else default_flukaserver_path  # TODO
        elif testing:
            raise ValueError("Cannot specify flukaserver executable when `testing=True`!")
        else:
            self._flukaserver = Path(flukaserver).resolve()
        if verbose is None:
            if not hasattr(self, '_verbose'):
                self._verbose = True
        else:
            self._verbose = verbose

        # Apply kwargs
        if(self._initialised):
            for kk, vv in kwargs.items():
                if not hasattr(self, kk):
                    raise ValueError(f"Invalid attribute {kk} for FlukaEngine!")
                setattr(self, kk, vv)
            return

        if '_xobject' not in kwargs:
            # Initialise defaults
            self._cwd = None
            self._old_cwd = None
            self._network_nfo = None
            self._log = None
            self._log_fid = None
            self._server_process = None
            self._input_file = None
            self.server_pid = None
            self._gfortran_installed = False
            self._flukaio_connected = False
            self._warning_given = False
            self._tracking_initialised = False
            self._insertion = {}
            self._collimator_dict = {}
            kwargs.setdefault('network_port', 0)
            kwargs.setdefault('_capacity', 5000)
            kwargs.setdefault('timeout_sec', 36000) # 10 hours
            kwargs.setdefault('particle_ref', xp.Particles())
            kwargs.setdefault('max_particle_id', 0)
            kwargs.setdefault('seed', -1)

        super().__init__(**kwargs)
        self._initialised = True


    def _warn_pyfluka(self, error):
        if not self._warning_given:
            print("Warning: Failed to import pyfluka (did you compile?). " \
                + "FlukaCollimators will be installed but are not trackable.\n", flush=True)
            print(error, flush=True)
            self._warning_given = True


    @classmethod
    def test_gfortran(cls, **kwargs):
        cls(**kwargs)
        this = cls.instance
        if not this._gfortran_installed:
            try:
                cmd = run(["gfortran", "-dumpversion"], stdout=PIPE, stderr=PIPE)
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
                    cmd2 = run(["which", "gfortran"], stdout=PIPE, stderr=PIPE)
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
    def start(cls, *, input_file=None, line=None, elements=None, names=None, cwd=None,
              prototypes_file=None, include_files=[], debug_level=0, touches=True,
              reference_particle=None, p0c=None, **kwargs):
        from .fluka_input import create_fluka_input, get_collimators_from_input_file, \
                                 verify_insertion_file, _beam_include_file
        cls(**kwargs)
        this = cls.instance
        if this.is_running():
            print("Server already running.", flush=True)
            return

        _fluka = flukafile_resolve(this._fluka)
        if _fluka is not None:
            this._fluka = _fluka
        else:
            raise ValueError(f"Could not find fluka executable {this._fluka}!")
        _flukaserver = flukafile_resolve(this._flukaserver)
        if _flukaserver is not None:
            this._flukaserver = _flukaserver
        else:
            raise ValueError(f"Could not find fluka executable {this._flukaserver}!")

        this.test_gfortran()
        this._starting_server = True

        # Check files
        if input_file is not None:
            input_file = Path(input_file).expanduser().resolve()
        if prototypes_file is not None:
            prototypes_file = Path(prototypes_file).expanduser().resolve()
        if include_files is not None:
            include_files = [Path(ff).expanduser().resolve() for ff in include_files]
        if cwd is not None:
            cwd = Path(cwd).expanduser().resolve()
        else:
            # TODO: use xaux.ranID here
            import base64
            ran = base64.urlsafe_b64encode(os.urandom(8)).decode('utf-8')
            ran_str = ''.join(c if c.isalnum() else 'X' for c in ran)
            cwd = Path.cwd() / f'fluka_run_{ran_str}'
        this._cwd = cwd
        cwd.mkdir(parents=True, exist_ok=True)
        this._old_cwd = Path.cwd()
        os.chdir(cwd)

        # Create input file
        if input_file is None:
            if not include_files:
                include_files.append(_beam_include_file(this.particle_ref.pdg_id[0]))
            print("Using default include files.")

            include_files.append(_pkg_root / 'scattering_routines' / 'fluka' / 'data' / 'include_settings_physics.inp')
            include_files.append(_pkg_root / 'scattering_routines' / 'fluka' / 'data' / 'include_custom_scoring.inp')

            input_file, elements, names = create_fluka_input(line=line, elements=elements, names=names,
                                            prototypes_file=prototypes_file,
                                            include_files=include_files)

        if not input_file.exists():
            raise ValueError(f"Input file {input_file.as_posix()} not found!")
        if input_file.parent != Path.cwd():
            shutil.copy(input_file, Path.cwd())
            input_file = Path.cwd() / input_file.name
        this._input_file = input_file

        # Check insertion file
        insertion_file = input_file.parent / "insertion.txt"
        if not insertion_file.exists():
            raise ValueError(f"Insertion file {insertion_file} not found!")
        if insertion_file.parent != Path.cwd():
            shutil.copy(insertion_file, Path.cwd())
            insertion_file = Path.cwd() / insertion_file.name

        # Match collimators
        collimator_dict = get_collimators_from_input_file(input_file)
        verify_insertion_file(insertion_file, collimator_dict)
        this._match_collimators_to_engine(collimator_dict, line=line, elements=elements, names=names)

        # Create touches file (relcol.dat)
        # # First line is the number of collimators, second line is the IDs (no newline at end)
        if touches is True:
            touches = Path('relcol.dat').resolve()
            with touches.open('w') as fid:
                fid.write(f'{len(this._collimator_dict.keys())}\n')
                for _, el in this._collimator_dict.items():
                    fid.write(f'{el.fluka_id} ')
        # Check if touches is a list of collimator names
        elif touches is not None and isinstance(touches, list):
            relcol = Path('relcol.dat').resolve()
            with relcol.open('w') as fid:
                fid.write(f'{len(touches)}\n')
                for touch in touches:
                    if touch not in this._collimator_dict.keys():
                        raise ValueError(f"Collimator {touch} not in collimator dict!")
                    else:
                        fid.write(f'{this._collimator_dict[touch].fluka_id} ')
        # Check if touches is not wrongly set
        elif touches is not None and not touches is False:
            raise NotImplementedError("Only True/False or a list of collimstors is allowed "
                                    + "for `touches` for now.")

        cls.clean_output_files()

        try:
            from .pyflukaf import pyfluka_init
            pyfluka_init(n_alloc=this._capacity, debug_level=debug_level)
        except ImportError as error:
            this._warn_pyfluka(error)
            this.stop(warn=False)
            return
        # Set reference particle
        if not this._has_particle_ref():
            if reference_particle is not None:
                if isinstance(reference_particle, xp.Particles):
                    if reference_particle.pdg_id == 0:
                        raise ValueError("The given `reference_particle` has no valid pdg_id!")
                    this.set_particle_ref(particle_ref=reference_particle)
                elif p0c is not None:
                    this.set_particle_ref(particle_ref=xp.Particles.reference_from_pdg_id(
                                          reference_particle, p0c=p0c))
                else:
                    raise ValueError("When providing `reference_particle`, it should be an "
                                     "xp.Particles object or a PDG ID. In the latter case, "
                                     "provide `p0c` as well.")
            elif line is not None:
                if line.particle_ref is None:
                    print("The given line has no reference particle. Don't forget to set it later.")
                else:
                    if line.particle_ref.pdg_id == 0:
                        print("The reference particle of the given line has no valid pdg_id, and can "
                            + "hence not be used.\nDon't forget to set the reference particle later.")
                    else:
                        this.set_particle_ref(line=line)
            else:
                print("No line given, so reference particle not yet set. Don't forget to set it later.")
        # Declare network
        this._network_nfo = this._cwd / network_file
        cmd = run(["hostname"], cwd=this._cwd, stdout=PIPE, stderr=PIPE)
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
        this._server_process = Popen([this._fluka, input_file, '-e', this._flukaserver, '-M', "1"], 
                                     cwd=this._cwd, stdout=this._log_fid, stderr=this._log_fid)
        this.server_pid = this._server_process.pid
        sleep(1)
        this._starting_server = False
        if not this.is_running():
            raise ValueError(f"Could not start fluka server! See logfile {this._log}.")
        i = 0
        while True:
            sleep(2)
            i += 2
            if i%30 == 0:
                print("Network port not yet established (or missing executable permission on "
                    + "flukaserver). Waiting 30s.", flush=True)
            with this._network_nfo.open('r') as fid:
                lines = fid.readlines()
                if len(lines) > 1:
                    this.network_port = int(lines[1].strip())
                    break
                if this._server_process.poll() is not None:
                    raise RuntimeError("The flukaserver died. Please check the FLUKA output.")
        print(f"Started fluka server on network port {this.network_port}. "
            + f"Connecting (timeout: {this.timeout_sec})... ", flush=True, end='')
        try:
            from .pyflukaf import pyfluka_connect
            pyfluka_connect(this.timeout_sec)
            this._flukaio_connected = True
        except ImportError as error:
            this._warn_pyfluka(error)
        print(f"done.", flush=True)


    @classmethod
    def stop(cls, warn=True, clean=False, **kwargs):
        cls(**kwargs)
        this = cls.instance
        # Stop flukaio connection
        if this._flukaio_connected:
            print(f"Closing fluka server connection... ", flush=True, end='')
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
            while this._server_process.poll() is None:
                sleep(1)
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
        if this._old_cwd is not None:
            os.chdir(this._old_cwd)
            this._old_cwd = None
        if clean:
            this.clean_output_files(clean_all=True)
        # Unassign the prototypes
        for name, ee in this._collimator_dict.items():
            ee.assembly.remove_element(name, force=False)
        this._cwd = None
        this._network_nfo = None
        this._log = None
        this._log_fid = None
        this._server_process = None
        this._input_file = None
        this.server_pid = None
        this._gfortran_installed = False
        this._flukaio_connected = False
        this._warning_given = False
        this._tracking_initialised = False
        this._insertion = {}
        this._collimator_dict = {}
        this.network_port = 0
        this.max_particle_id =  0
        this.seed = -1


    @classmethod
    def is_running(cls, **kwargs):
        cls(**kwargs)
        this = cls.instance
        if hasattr(this, '_starting_server') and this._starting_server:
            return False
        # Is the Popen process still running?
        if this._server_process is None or this._server_process.poll() is not None:
            this.stop()
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
            this.stop()
            return False
        elif len(processes) == 1 and str(this.server_pid) in processes[0] and 'defunct' not in processes[0]:
            return True
        elif np.any([str(this.server_pid) in proc for proc in processes if 'defunct' not in proc]):
            print("Warning: Found other instances of rfluka besides the current one!", flush=True)
            return True
        else:
            print("Warning: Found other instances of rfluka but not the current one!", flush=True)
            this.stop()
            return False


    @classmethod
    def clean_output_files(cls, input_file=None, cwd=None, clean_all=False, **kwargs):
        cls(**kwargs)
        this = cls.instance
        if input_file is None:
            input_file = this._input_file
            if cwd is None and this._cwd is not None:
                cwd = this._cwd
        else:
            input_file = Path(input_file)
            if cwd is None:
                cwd = input_file.parent
        if isinstance(cwd, Path):
            files_to_delete = [cwd / f for f in [network_file, fluka_log, server_log, 'fluka_isotope.log',
                                                 'fort.208', 'fort.251']]
            if input_file is not None:
                files_to_delete += list(cwd.glob(f'ran{input_file.stem}*'))
                files_to_delete += list(cwd.glob(f'{input_file.stem}*'))
                files_to_delete  = [file for file in files_to_delete
                                    if Path(file).resolve() != input_file.resolve()]
            if clean_all:
                files_to_delete += [cwd / f for f in ['insertion.txt', 'new_collgaps.dat',
                                                      'prototypes.lbp', 'relcol.dat']]
                files_to_delete += [this._input_file]
            for f in files_to_delete:
                if f is not None and f.exists():
                    f.unlink()
            if clean_all and cwd != Path.cwd():
                cwd.rmdir()


    def _match_collimators_to_engine(self, collimator_dict, *, line=None, elements=None, names=None):
        from ...beam_elements import FlukaCollimator
        if line is None:
            if elements is None or names is None:
                raise ValueError("Need to provide either `line` or `elements` and `names`.")
        else:
            elements, names = line.get_elements_of_type(FlukaCollimator)
        if not hasattr(elements, '__iter__') or isinstance(elements, str):
            elements = [elements]
        if not hasattr(names, '__iter__') or isinstance(names, str):
            names = [names]
        assert len(elements) == len(names)
        self._collimator_dict = {}
        for el, name in zip(elements, names):
            if name not in collimator_dict:
                if el.active:
                    print(f"Warning: FlukaCollimator {name} not in FLUKA input file! "
                        + f"Maybe it was fully open. Deactivated")
                    el.active = False
                continue
            el.fluka_id = collimator_dict[name]['fluka_id']
            el.length_front = (collimator_dict[name]['length'] - el.length)/2
            el.length_back = (collimator_dict[name]['length'] - el.length)/2
            jaw = collimator_dict[name]['jaw']
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
        # TODO: tilts!!


    @classmethod
    def init_tracking(cls, max_particle_id, **kwargs):
        cls(**kwargs)
        this = cls.instance
        if this._tracking_initialised == False:
            particle_ref = this.particle_ref
            _, A0, Z0, name = pdg.get_properties_from_pdg_id(particle_ref.pdg_id[0])
            # FLUKA expects units of MeV
            E0 = particle_ref.energy0[0] / 1.e6
            p0c = particle_ref.p0c[0] / 1.e6
            m0 = particle_ref.mass0 / 1.e6
            q0 = particle_ref.q0
            try:
                from .pyflukaf import pyfluka_init_max_uid, pyfluka_set_synch_part
            except ImportError as error:
                this._warn_pyfluka(error)
                return
            print(f"Setting max_particle_id to {max_particle_id}, "
                + f"and reference particle to {name} with mass {m0}MeV "
                + f"and energy {E0}MeV.", flush=True)
            pyfluka_init_max_uid(max_particle_id)
            pyfluka_set_synch_part(E0, p0c, m0, A0, Z0, q0)
            this.max_particle_id = max_particle_id
            this._tracking_initialised = True


    @classmethod
    def set_particle_ref(cls, particle_ref=None, line=None, **kwargs):
        cls(**kwargs)
        this = cls.instance
        # if this._has_particle_ref():
        #     print("Reference particle already set!")
        #     return
        overwrite_particle_ref_in_line = False
        if particle_ref is None:
            if line is None or line.particle_ref is None:
                raise ValueError("Line has no reference particle! "
                               + "Please provide one with `particle_ref`.")
            if line.particle_ref.pdg_id == 0:
                raise ValueError("`line.particle_ref` needs to have a valid pdg_id")
            particle_ref = line.particle_ref
        else:
            if particle_ref._capacity > 1:
                raise ValueError("`particle_ref` has to be a single particle!")
            if particle_ref.pdg_id == 0:
                raise ValueError("`particle_ref` needs to have a valid pdg_id")
            if line is not None:
                if line.particle_ref is not None \
                and not xt.line._dicts_equal(line.particle_ref.to_dict(), particle_ref.to_dict()):
                    overwrite_particle_ref_in_line = True
        this._compare_fluka_mass(particle_ref)
        this.particle_ref = particle_ref
        if overwrite_particle_ref_in_line:
            print("Warning: Found different reference particle in line! Overwritten.")
            line.particle_ref = particle_ref

    def _has_particle_ref(self):
        initial = xp.Particles().to_dict()
        current = self.particle_ref.to_dict()
        return not xt.line._dicts_equal(initial, current)

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

