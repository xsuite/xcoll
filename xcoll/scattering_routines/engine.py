# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

import numpy as np
from pathlib import Path
from numbers import Number

import xobjects as xo
import xtrack as xt
import xtrack.particles.pdg as pdg

from .geometry import XcollGeometry
from .physics_settings import PhysicsSettingsHelper
from ..materials import Material
from ..interaction_record import InteractionRecord
try:
    # TODO: once xaux is in Xsuite keep only this
    from xaux import FsPath, ranID
except (ImportError, ModuleNotFoundError):
    from ..xaux import FsPath, ranID


class BaseEngine(xo.HybridClass):
    _xofields = {
        '_particle_ref':      xt.Particles._XoStruct,
        '_seed':              xo.UInt64,
        '_capacity':          xo.Int64,
        '_relative_capacity': xo.Int64,
    }

    _int32 = False
    _only_protons = False
    _element_classes = None
    _uses_input_file = False
    _num_input_files = 1
    _uses_run_folder = False

    _depends_on = [Material, InteractionRecord, xt.RandomUniform,
                   xt.RandomExponential, xt.RandomNormal, xt.RandomRutherford,
                   xt.Drift, XcollGeometry]

    def __init__(self, **kwargs):
        if self._element_classes is None:
            raise NotImplementedError(f"{self.__class__.__name__} needs to define `_element_classes`!")
        # Initialise defaults
        self._cwd = None
        self._line = None
        self._masses = None
        self._verbose = False
        self._input_file = None
        self._environment = None
        self._element_dict = {}
        self._warning_given = False
        self._element_index = 0
        self._tracking_initialised = False
        self._deactivated_elements = {}
        kwargs.setdefault('_particle_ref', xt.Particles())
        kwargs.setdefault('_seed', 0)
        kwargs.setdefault('_capacity', 0)
        super().__init__(**kwargs)
        self._physics_settings = PhysicsSettingsHelper(self)

    def __del__(self, *args, **kwargs):
        self.stop()

    def _warn(self, error=None):
        if not self._warning_given:
            print(f"Warning: Failed to import {self.__class__.__name__} environment "
                + f"(did you compile?).\n{self.name.capitalize()} elements can be installed "
                + f"but are not trackable.", flush=True)
            self._warning_given = True
        self.stop()
        if error:
            raise error

    def _print(self, *args, **kwargs):
        if self.verbose:
            kwargs.setdefault('flush', True)
            print(*args, **kwargs)


    # ==================
    # === Properties ===
    # ==================

    @property
    def environment(self):
        return self._environment

    @property
    def name(self):
        return self.__class__.__name__.replace('Engine', '').lower()

    @property
    def cwd(self):
        return self._cwd

    @property
    def verbose(self):
        return self._verbose

    @verbose.setter
    def verbose(self, val):
        self._verbose = val

    @property
    def line(self):
        return self._line

    @line.setter
    def line(self, val):
        if not val is None and not isinstance(val, xt.Line):
            self.stop()
            raise ValueError("`line` has to be an xt.Line object!")
        self._line = val

    @line.deleter
    def line(self):
        self._line = None

    @property
    def particle_ref(self):
        initial = xt.Particles().to_dict()
        current = self._particle_ref.to_dict()
        if xt.line._dicts_equal(initial, current):
            return None
        else:
            return self._particle_ref

    @particle_ref.setter
    def particle_ref(self, val):
        if val is None:
            self._particle_ref = xt.Particles()
        else:
            if isinstance(val, xt.line.LineParticleRef):
                val = val._resolved
            if not isinstance(val, xt.Particles):
                self.stop()
                raise ValueError("`particle_ref` has to be an xt.Particles object!")
            if val._capacity > 1:
                self.stop()
                raise ValueError("`particle_ref` has to be a single particle!")
            pdg_id = val.pdg_id[0]
            if pdg_id == 0:
                if self._only_protons:
                    pdg_id = pdg.get_pdg_id_from_name('proton')
                else:
                    self.stop()
                    raise ValueError(f"{self.__class__.__name__} allows the use of particles "
                                   + f"different than protons. Hence, `particle_ref` "
                                   + f"needs to have a valid pdg_id.")
            elif self._only_protons and pdg_id != pdg.get_pdg_id_from_name('proton'):
                self.stop()
                raise ValueError("{self.__class__.__name__} only supports protons!")
            self._particle_ref = val
            self._particle_ref.pdg_id[0] = pdg_id
        self._physics_settings.update()

    @particle_ref.deleter
    def particle_ref(self):
        self.particle_ref = None

    @property
    def capacity(self):
        if self._capacity == 0:
            return None
        else:
            return int(self._capacity)

    @capacity.setter
    def capacity(self, val):
        if val is None:
            val = 0
        if not isinstance(val, Number) or val < 0:
            self.stop()
            raise ValueError("`capacity` has to be a positive integer!")
        self._capacity = int(val)

    @capacity.deleter
    def capacity(self):
        self.capacity = None

    @property
    def relative_capacity(self):
        if self._relative_capacity == 0:
            return None
        else:
            return int(self._relative_capacity)

    @relative_capacity.setter
    def relative_capacity(self, val):
        if val is None:
            val = 2
        if not isinstance(val, Number) or val < 0:
            self.stop()
            raise ValueError("`relative_capacity` has to be a positive integer!")
        if val <= 1:
            self.stop()
            raise ValueError("`relative_capacity` has to be larger than 1!")
        self._relative_capacity = int(val)

    @relative_capacity.deleter
    def relative_capacity(self):
        self.relative_capacity = None

    @property
    def seed(self):
        if self._seed == 0:
            return None
        else:
            return self._seed

    @seed.setter
    def seed(self, val):
        if val is None:
            val = 0
        if not isinstance(val, Number) or val < 0:
            self.stop()
            raise ValueError("`seed` has to be a positive integer!")
        val = int(val)
        if self._int32:
            new_val = np.uint32(val)
        else:
            new_val = np.uint64(val)
        if new_val != val:
            self._print(f"Warning: type change for seed {val}. Using {new_val}.")
        self._seed = new_val

    @seed.deleter
    def seed(self):
        self.seed = None

    @property
    def input_file(self):
        if self._uses_input_file:
            return self._input_file

    @property
    def element_dict(self):
        return self._element_dict

    def __getattr__(self, name):
        if name != '_physics_settings' and hasattr(self, '_physics_settings') and name in self._physics_settings.all_flags:
            return getattr(self._physics_settings, name)
        raise AttributeError(f"'{self.__class__.__name__}' object has no attribute '{name}'")

    def __setattr__(self, name, value):
        if hasattr(self, '_physics_settings') and name in self._physics_settings.all_flags:
            return setattr(self._physics_settings, name, value)
        return super().__setattr__(name, value)

    # ======================
    # === Public Methods ===
    # ======================

    def physics_settings(self):
        return self._physics_settings.show()

    def start(self, *, clean=True, input_file=None, **kwargs):
        if not self.environment:
            self.stop()
            raise RuntimeError(f"{self.name.capitalize()} environment not set up! "
                             + f"Do not manually create an instance of the engine.")
        self.environment.assert_environment_ready()
        if self.is_running():
            self._print("Engine already running.")
            return

        # Clean up any leftover failed runs
        self.stop(clean=clean)

        print(f"Starting {self.__class__.__name__}...   ", flush=True, end='')
        if self.verbose:
            print("", flush=True)

        kwargs = self._pre_start(**kwargs)

        # This needs to be set in the ChildEngine, either in _start_engine() or at the start of tracking
        self._tracking_initialised = False

        # Resolve input file path before defining cwd
        if input_file:
            if not isinstance(input_file, (str,Path)):
                self.stop()
                raise ValueError("`input_file` has to be a string or Path!")
            if self._num_input_files > 1:
                input_file = [FsPath(f).expanduser().resolve() for f in input_file]
                for f in input_file:
                    if not f.exists():
                        self.stop()
                        raise ValueError(f"Input file {f} does not exist!")
            else:
                input_file = FsPath(input_file).expanduser().resolve()
                if not input_file.exists():
                    self.stop()
                    raise ValueError(f"Input file {input_file} does not exist!")

        # Set all engine properties that have a setter (this will remove these properties from the kwargs)
        self._starting_or_stopping = True # We need this to allow changing the element settings which otherwise are locked
        kwargs = self._set_engine_properties(**kwargs)

        # Do some preparations before creating the input file if needed
        kwargs = self._pre_input(**kwargs)

        # Create input file if needed (this will remove the kwargs relevant to the input file and physics)
        kwargs = self._use_input_file(input_file, **kwargs)
        if clean:
            self.clean_input_files(clean_all=False)
        self._starting_or_stopping = False

        # Start the engine in the ChildEngine
        self._start_engine(**kwargs)

        # Done starting
        if self.verbose:
            print(f"{self.__class__.__name__} started.", flush=True)
        else:
            print(f"Done.", flush=True)

    def stop(self, clean=False, **kwargs):
        kwargs = self._stop_engine(**kwargs)
        if clean:
            self.clean(clean_all=True, **kwargs)
        self._starting_or_stopping = True # We need this to allow changing the element settings which otherwise are locked
        self._restore_engine_properties(clean=clean)
        self._starting_or_stopping = False
        self._warning_given = False
        self._tracking_initialised = False

    def is_running(self):
        if hasattr(self, '_starting_or_stopping') and self._starting_or_stopping:
            # We need this to allow changing the element settings which otherwise are locked
            return False
        # If we get here, we cannot say if the engine is running or not and we need an
        # implementation in the child class
        return self._is_running()


    def generate_input_file(self, *, clean=True, filename=None, **kwargs):
        '''This method manually generates an input file without starting the engine'''
        if not self._uses_input_file:
            self.stop()
            raise ValueError(f"{self.__class__.__name__} does not use input files!")
        if self._element_dict:
            self.stop()
            raise ValueError("Elements already assigned to engine (cannot regenerate input "
                           + "file after starting engine)!")

        # Set all engine properties that have a setter (this will remove these properties from the kwargs)
        cwd = kwargs.get('cwd')
        if cwd and FsPath(cwd).resolve() == FsPath.cwd().resolve():
            print("Warning: Cannot use current working directory as input folder. "
                + "Generating temporary folder instead.")
            kwargs.pop('cwd')
        self._starting_or_stopping = True # We need this to allow changing the element settings which otherwise are locked
        kwargs = self._set_engine_properties(**kwargs)

        # Do some preparations before creating the input file if needed
        kwargs = self._pre_input(**kwargs)

        # Create input file
        input_file, _ = self._generate_input_file(**kwargs)
        if not hasattr(input_file, '__iter__') or isinstance(input_file, str):
            # Some engines might create multiple input files (like Fluka)
            input_file = [input_file]

        # Move input file to desired location
        if filename is None:
            if input_file[0].parent != FsPath.cwd():
                new_input_file = [input_file[0].rename(Path.cwd() / input_file[0].name)]
            else:
                new_input_file = [input_file[0]]
        else:
            new_input_file = [input_file[0].rename(filename)]
        for file in input_file[1:]:
            if file.parent != new_input_file[0].parent:
                new_input_file.append(file.rename(new_input_file[0].parent / file.name))
            else:
                new_input_file.append(file)

        # Clean up
        if clean:
            if input_file[0].parent == new_input_file[0].parent:
                self.clean_input_files(clean_all=False, input_file=input_file)
            else:
                self.clean_input_files(clean_all=True, input_file=input_file)
        self._restore_engine_properties(clean=clean)
        self._starting_or_stopping = False

        return new_input_file[0] if self._num_input_files==1 else new_input_file


    def assert_particle_ref(self):
        if self.particle_ref is None:
            self.stop()
            raise ValueError(f"{self.__class__.__name__} reference particle not set!")

    def assert_ready_to_track_or_skip(self, coll, particles, _necessary_attributes=[], keep_p0c_constant=True):
        '''Asserts that the engine is ready to track the given element and particles.'''
        self._assert_element(coll)

        missing_attributes = False
        for attr in _necessary_attributes:
            if not hasattr(coll, attr) or not getattr(coll, attr):
                missing_attributes = True

        if not coll.active or not coll._tracking or not coll.jaw or missing_attributes:
            return False

        npart = particles._num_active_particles
        if npart == 0:
            return False
        if not isinstance(particles._buffer.context, xo.ContextCpu):
            self.stop()
            raise ValueError(f"{self.__class__.__name__} only supports CPU contexts!")

        assert self.environment.compiled
        if not self.is_running():
            self.stop()
            raise ValueError(f"{self.__class__.__name__} not yet running!\nPlease do this "
                           + f"first, by calling xcoll.{self.__class__.__name__}.start().")

        self.assert_particle_ref()
        if abs(particles.mass0 - self.particle_ref.mass0) > 1e-3:
            # This should only happen at the start of the tracking.
            pdg_id = self.particle_ref.pdg_id[0]
            if self._masses is not None and abs(pdg_id) in self._masses:
                # The reference particle in the engine was updated to match the scattering code.
                # Match the particles object to it.
                new_mass = self.particle_ref.mass0
                old_energy0 = particles.energy0[0]
                particles.mass0 = new_mass
                if keep_p0c_constant:
                    particles._update_refs(p0c=particles.p0c[0])
                else:
                    particles._update_refs(energy0=old_energy0)
                assert np.isclose(particles.energy0[0]**2, particles.p0c[0]**2 + particles.mass0**2)
                assert np.isclose(particles.mass0, new_mass)
                # No need to update chi, delta, etc as they depend on the RATIO m0/m and this is unchanged.
                self._print(f"Warning: reference mass in particles object differ"
                          + f"from reference mass in engine. Overwritten by the latter.")
            else:
                # The reference particle in the engine was changed unintentionally
                self.stop()
                raise ValueError(f"Error in reference mass of `particles`: not in sync with "
                            + f"{self.name} reference particle!\nRebuild the particles object "
                            + f"using the {self.__class__.__name__} reference particle.")
        if abs(particles.q0 - self.particle_ref.q0) > 1e-3:
            self.stop()
            raise ValueError(f"Error in reference charge of `particles`: not in sync with "
                           + f"{self.name} reference particle!\nRebuild the particles object "
                           + f"using the {self.__class__.__name__} reference particle.")
        if not self._only_protons:
            if np.any([pdg_id == 0 for pdg_id in particles.pdg_id]):
                self.stop()
                raise ValueError("Some particles are missing the pdg_id!")
            if particles._num_active_particles + particles._num_lost_particles == particles._capacity \
            and not np.any(particles.particle_id != particles.parent_particle_id):
                # Only raise this error at the start, e.g. when no secondaries are present yet.
                # It will get caught later during tracking, which will provide a more logical error.
                self.stop()
                raise ValueError("Particles capacity equal to size! Please provide extra capacity "
                               + "for secondaries.")
        return True


    def clean(self, **kwargs):
        '''Cleans all input and output files generated by the engine.
        Use clean_all (default: False) to also delete generated input files.'''
        self.clean_input_files(**kwargs)
        self.clean_output_files(**kwargs)

    def clean_input_files(self, clean_all=False, **kwargs):
        kwargs = self._get_input_cwd_for_cleaning(**kwargs)
        files_to_delete = self._get_input_files_to_clean(**kwargs)
        if not clean_all:
            files_to_delete = [f for f in files_to_delete
                    if f not in self._all_input_files(kwargs['input_file'])]
        for f in files_to_delete:
            if f is not None and f.exists():
                f.unlink()

    def clean_output_files(self, clean_all=False, **kwargs):
        kwargs = self._get_input_cwd_for_cleaning(**kwargs)
        files_to_delete = self._get_output_files_to_clean(**kwargs)
        if not clean_all:
            files_to_delete = [f for f in files_to_delete
                    if f not in self._all_input_files(kwargs['input_file'])]
        for f in files_to_delete:
            if f is not None and f.exists():
                f.unlink()


    # =======================
    # === Private Methods ===
    # =======================

    # For all the engine fields, they can either be set in advance on the engine,
    # or they can be set when the engine is started. In the latter case, the values
    # are temporary and the original will be restored when the engine is stopped.

    def _set_property(self, prop, kwargs):
        val = kwargs.pop(prop, None)
        if val is not None:
            # We only need to update the property when it is not None
            setattr(self, f'_old_{prop}', getattr(self, prop))
            setattr(self, prop, val)

    def _set_engine_properties(self, **kwargs):
        # We need to set the following properties first as they are needed by the others
        self._set_property('verbose', kwargs)
        self._set_property('line', kwargs)
        # The following properties have a specific logic
        self._use_seed(kwargs.pop('seed', None))
        self._use_particle_ref(kwargs.pop('particle_ref', None))
        self._sync_line_particle_ref()
        self._get_elements(kwargs.pop('elements', None), kwargs.pop('names', None))
        self._set_cwd(kwargs.pop('cwd', None))
        # Now we can set the rest of the properties
        self._set_property('capacity', kwargs)
        self._set_property('relative_capacity', kwargs)
        for ff in self._physics_settings.all_flags:
            self._set_property(ff, kwargs)
        return kwargs

    def _restore_engine_properties(self, clean=False):
        # Reset particle_ref in the line
        if hasattr(self, '_old_line_particle_ref'):
            self.line.particle_ref = self._old_line_particle_ref
            del self._old_line_particle_ref
        # The following properties have a specific logic
        self._reactivate_elements()
        self._reset_cwd(clean=clean)
        # Reset all other properties
        self_attributes = self.__dict__.copy()
        for kk, vv in self_attributes.items():
            if kk.startswith('_old_'):
                prop = kk[5:]
                setattr(self, prop, vv)
                delattr(self, kk)
        self._input_file = None
        self._element_dict = {}

    def _use_seed(self, seed=None):
        if seed is None:
            if self.seed is None:
                rng = np.random.default_rng()
                if self._int32:
                    self.seed = rng.integers(0, int(2**31)) # 32-bit signed
                else:
                    self.seed = rng.integers(0, int(2**63)) # 64-bit signed
        else:
            self._old_seed = self.seed
            self.seed = seed
        self._print(f"Using seed {self.seed}.")

    def _use_particle_ref(self, particle_ref=None, keep_p0c_constant=True):
        # Prefer: provided particle_ref > existing particle_ref > particle_ref from line
        if particle_ref is not None:
            self._old_particle_ref = self.particle_ref
            self.particle_ref = particle_ref
        elif self.particle_ref is None:
            if self.line is None or not hasattr(self.line, 'particle_ref') \
            or self.line.particle_ref is None:
                self.stop()
                raise ValueError("Need to provide either a line with a reference "
                               + "particle, or `particle_ref`.")
            self._old_particle_ref = self.particle_ref
            self.particle_ref = self.line.particle_ref
        self._print(f"Using {pdg.get_name_from_pdg_id(self.particle_ref.pdg_id[0])} "
                  + f"with momentum {self.particle_ref.p0c[0]/1.e9:.1f} GeV.")
        if self._masses is not None:
            mass = self.particle_ref.mass0
            pdg_id = self.particle_ref.pdg_id[0]
            if abs(pdg_id) in self._masses:
                new_mass = self._masses[abs(pdg_id)]
                if abs(mass-new_mass)/mass > 1.e-12:
                    old_energy0 = self.particle_ref.energy0[0]
                    self.particle_ref.mass0  = new_mass
                    if keep_p0c_constant:
                        self.particle_ref._update_refs(p0c=self.particle_ref.p0c[0])
                    else:
                        self.particle_ref._update_refs(energy0=old_energy0)
                    assert np.isclose(self.particle_ref.energy0[0]**2, self.particle_ref.p0c[0]**2 + self.particle_ref.mass0**2)
                    assert np.isclose(self.particle_ref.mass0, new_mass)
                    self._print(f"Warning: given mass of {mass} eV for "
                            + f"{pdg.get_name_from_pdg_id(pdg_id)} differs from {self.name} "
                            + f"mass of {new_mass} eV.\nReference particle mass is "
                            + f"overwritten by the latter.")
            else:
                self._print(f"Warning: No {self.name} reference mass known for particle "
                        + f"{pdg.get_name_from_pdg_id(pdg_id)}!\nIf the reference mass "
                        + f"provided ({mass} eV) differs from the one used internally "
                        + f"by {self.name}, differences in energy might be observed.\nOnce "
                        + f"the {self.name} reference mass is known, contact the devs to "
                        + f"input it in the code.")

    def _sync_line_particle_ref(self):
        if self.line is None:
            return
        if self.line.particle_ref is not None \
        and not xt.line._dicts_equal(self.line.particle_ref.to_dict(),
                                     self.particle_ref.to_dict()):
            self._print("Found different reference particle in line. Temporarily overwritten.")
            val = self.line.particle_ref
            if isinstance(val, xt.line.LineParticleRef):
                val = val._resolved
            self._old_line_particle_ref = val
            self.line.particle_ref = self.particle_ref

    def _get_new_element_name(self):
        name = f"{self.name}_el_{self._element_index}"
        self._element_index += 1
        return name

    def _assert_element(self, element):
        if not isinstance(element, self._element_classes):
            self.stop()
            raise ValueError(f"Element {element.name} is not a "
                            + ", or a ".join([c.__name__ for c in self._element_classes])
                            + ".")

    def _get_elements(self, elements=None, names=None):
        if elements is not None and (not hasattr(elements, '__iter__') or isinstance(elements, str)):
            elements = [elements]
        if names is not None and (not hasattr(names, '__iter__') or isinstance(names, str)):
            names = [names]
        if self.line is None:
            if elements is None:
                self.stop()
                raise ValueError("Need to provide either `line` or `elements`.")
            if names is None:
                names = []
                for ee in elements:
                    if hasattr(ee, 'name') and ee.name:
                        names.append(ee.name)
                    else:
                        name = self._get_new_element_name()
                        names.append(name)
                        ee.name = name
            elif len(names) == len(elements):
                for ee, name in zip(elements, names):
                    if hasattr(ee, 'name'):
                        if ee.name != name:
                            self._print(f"Warning: Element name {ee.name} changed to {name}.")
                            ee.name = name
            else:
                self.stop()
                raise ValueError("Length of `elements` and `names` doesn't match.")
        else:
            if elements is not None:
                self.stop()
                raise ValueError("Cannot provide both `line` and `elements`.")
            if names is None:
                elements, names = self.line.get_elements_of_type(self._element_classes)
                # This method can return duplicate elements if there are subclasses (like Geant4Collimator
                # and Geant4CollimatorTip) as elements of the child type will also be found by searching
                # for the parent type.
                d = {}
                for nn, ee in zip(names, elements):
                    d[nn] = ee # last one wins, not important here
                names, elements = list(d.keys()), list(d.values())
            else:
                elements = [self.line[name] for name in names]
        this_names = []
        this_elements = []
        for ee, name in zip(elements, names):
            self._assert_element(ee)
            if ee.jaw is None:
                self._print(f"Warning: Jaw not set for {name}. Ignoring.")
                self._deactivate_element(ee)
            elif not ee.active:
                self._print(f"Warning: Element {name} is not active. Ignoring.")
                self._deactivate_element(ee)
            else:
                this_names.append(name)
                this_elements.append(ee)
        if len(this_elements) == 0:
            self.stop()
            raise ValueError(f"No active {self.name} elements found!")
        if len(set(this_names)) != len(this_names):
            self.stop()
            raise ValueError(f"Duplicate names found in {self.name} elements: {this_names}. "
                           + f"Please provide unique names for each element.")
        self._element_dict = dict(zip(this_names, this_elements))

    def _deactivate_element(self, el):
        self._deactivated_elements[el.name] = [el, el.active or True]
        if hasattr(el, 'active'):
            el.active = False
        self._remove_element(el)

    def _reactivate_elements(self):
        names = list(self._deactivated_elements.keys())
        for name in names:
            ee, was_active = self._deactivated_elements.pop(name)
            self._restore_element(ee)
            if hasattr(ee, 'active'):
                ee.active = was_active

    def _set_cwd(self, cwd=None):
        if self._uses_run_folder:
            if cwd is not None:
                cwd = FsPath(cwd).expanduser().resolve()
                if cwd.exists() and cwd != FsPath.cwd():
                    # Check if the folder already exists
                    # If the specified folder is the current one, we do not need to rename
                    i = 0
                    while (cwd.parent / f'{cwd.name}_{i:0>4}').exists():
                        i += 1
                        if i > 9999:
                            self.stop()
                            raise ValueError(f"Too many folders with the same "
                                           + f"name {cwd}!")
                    cwd = cwd.parent / f'{cwd.name}_{i:0>4}'
            else:
                ran_str = ranID(only_alphanumeric=True)
                cwd = FsPath.cwd() / f'{self.name}_run_{ran_str}'
            self._cwd = cwd
            cwd.mkdir(parents=True, exist_ok=True)
            # self._old_cwd = FsPath.cwd()
            # os.chdir(cwd)

    def _reset_cwd(self, clean=False):
        if clean and self.cwd is not None and self.cwd != FsPath.cwd():
            try:
                # Only works if the folder is empty
                self.cwd.rmdir()
            except:
                self._print(f"Warning: Failed to remove temporary folder "
                            + f"{self.cwd}.")
        self._cwd = None

    def _use_input_file(self, input_file=None, **kwargs):
        if self._uses_input_file:
            if input_file is None:
                input_file, kwargs = self._generate_input_file(**kwargs)
            if not hasattr(input_file, '__iter__') or isinstance(input_file, (str,Path)):
                # Some engines might need multiple input files (like Fluka)
                input_file = [input_file]
            input_file = [FsPath(f).expanduser().resolve() for f in input_file]
            new_files = []
            for file in input_file:
                if not file.exists():
                    self.stop()
                    raise ValueError(f"Input file {file.as_posix()} not found!")
                if file.parent != self.cwd and self._uses_run_folder:
                    if self.cwd is None:
                        self.stop()
                        raise ValueError("Cannot copy input file to working directory: "
                                       + "working directory not set!")
                    file.copy_to(self.cwd, method='mount')
                    new_files.append(self.cwd / file.name)
                else:
                    new_files.append(file)
            self._input_file = new_files[0] if self._num_input_files==1 else new_files
            self._match_input_file()
        return kwargs

    def _get_input_cwd_for_cleaning(self, **kwargs):
        # Get the input file and the cwd
        input_file = kwargs.get('input_file', None)
        cwd = kwargs.get('cwd', None)
        if input_file is None:
            if self.input_file is not None:
                input_file = self.input_file
            if cwd is None and self.cwd is not None:
                cwd = self.cwd
        else:
            if not hasattr(input_file, '__iter__') or isinstance(input_file, (str,Path)):
                input_file = [input_file]
            input_file = [FsPath(ff) for ff in input_file]
            if cwd is None:
                cwd = input_file[0].parent
            if len(input_file) == 1:
                input_file = input_file[0]
        kwargs['input_file'] = input_file
        if self._uses_run_folder:
            kwargs['cwd'] = cwd
        else:
            kwargs['cwd'] = FsPath.cwd()
        return kwargs


    def _mask_particle_return_types(self, pdg_id, q_new):
        if self.return_all:
            # Allow everything and exclude
            mask_new = np.ones_like(pdg_id, dtype=bool)
        else:
            # Allow nothing and include
            mask_new = np.zeros_like(pdg_id, dtype=bool)

        # General categories
        mask_new[pdg_id > 1000000000] = self.return_ions
        # PDG ID of mesons: from .*0XX. where X != 0 and . is any digit
        mask_new[(pdg_id > 0) & (pdg_id // 10 % 10 != 0) & (pdg_id // 100 % 10 != 0)
                              & (pdg_id // 1000 % 10 == 0)] = self.return_other_mesons
        mask_new[(pdg_id < 0) & (-pdg_id // 10 % 10 != 0) & (-pdg_id // 100 % 10 != 0)
                              & (-pdg_id // 1000 % 10 == 0)] = self.return_other_mesons
        # PDG ID of baryons: from XXX. where X != 0 and . is any digit
        mask_new[(pdg_id > 1000) & (pdg_id < 9000) & (pdg_id // 10 % 10 != 0) & (pdg_id // 100 % 10 != 0)
                                 & (pdg_id // 1000 % 10 != 0)] = self.return_other_baryons    # PDG ID of from XX0X is a diquark
        mask_new[(pdg_id < -1000) & (pdg_id > -9000) & (-pdg_id // 10 % 10 != 0) & (-pdg_id // 100 % 10 != 0)
                                  & (-pdg_id // 1000 % 10 != 0)] = self.return_other_baryons

        if not self.return_neutral:
            # General modifier, has to be before more specific return types,
            # as other neutral particles might have been specifically activated.
            mask_new[np.abs(q_new) < 1.e-12] = False

        mask_new[pdg_id == 22] = self.return_photons
        mask_new[(pdg_id == 11) | (pdg_id == -11)] = self.return_electrons
        mask_new[(pdg_id == 12) | (pdg_id == -12)] = self.return_electrons and self.return_neutrinos
        mask_new[(pdg_id == 13) | (pdg_id == -13)] = self.return_muons
        mask_new[(pdg_id == 14) | (pdg_id == -14)] = self.return_muons and self.return_neutrinos
        mask_new[(pdg_id == 15) | (pdg_id == -15)] = self.return_tauons
        mask_new[(pdg_id == 16) | (pdg_id == -16)] = self.return_tauons and self.return_neutrinos
        mask_new[(pdg_id == 211) | (pdg_id == -211)] = self.return_pions
        mask_new[(pdg_id == 111)] = self.return_pions and self.return_neutral
        mask_new[(pdg_id == 321) | (pdg_id == -321)] = self.return_kaons
        mask_new[(pdg_id == 130) | (pdg_id == -130)] = self.return_kaons and self.return_neutral
        mask_new[(pdg_id == 310) | (pdg_id == -310)] = self.return_kaons and self.return_neutral
        mask_new[(pdg_id == 311) | (pdg_id == -311)] = self.return_kaons and self.return_neutral
        mask_new[(pdg_id == 2212) | (pdg_id == -2212)] = self.return_protons
        mask_new[(pdg_id == 2112) | (pdg_id == -2112)] = self.return_neutrons
        return mask_new


    # =================================================
    # === Methods to be overwritten by child engine ===
    # =================================================

    def _start_engine(self, **kwargs):
        raise NotImplementedError(f"Need to implement `_start_engine` for "
                                   + f"{self.__class__.__name__}!")

    def _stop_engine(self, **kwargs):
        raise NotImplementedError(f"Need to implement `_stop_engine` for "
                                   + f"{self.__class__.__name__}!")

    def _is_running(self, **kwargs):
        raise NotImplementedError(f"Need to implement `_is_running` for "
                                   + f"{self.__class__.__name__}!")

    def _match_input_file(self, **kwargs):
        if self._uses_input_file:
            raise NotImplementedError(f"Need to implement `_match_input_file` for "
                                   + f"{self.__class__.__name__}!")

    def _generate_input_file(self, **kwargs):
        if self._uses_input_file:
            raise NotImplementedError("Need to implement `_generate_input_file` for "
                                   + f"{self.__class__.__name__}!")


    # ============================================================
    # === Methods to be optionally overwritten by child engine ===
    # ============================================================

    def _get_input_files_to_clean(self, input_file=None, cwd=None, clean_all=False, **kwargs):
        return []

    def _get_output_files_to_clean(self, input_file=None, cwd=None, clean_all=False, **kwargs):
        return []

    def _all_input_files(self, input_file=None):
        if self._uses_input_file:
            if self._num_input_files == 1:
                return [self._input_file]
            else:
                raise NotImplementedError(f"Need to implement `_all_input_files` for "
                                    + f"{self.__class__.__name__}!")
        else:
            return []

    def _remove_element(self, el):
        pass

    def _restore_element(self, el):
        pass

    def _pre_start(self, **kwargs):
        return kwargs

    def _pre_input(self, **kwargs):
        # Do some preparations before creating the input file if needed
        return kwargs
