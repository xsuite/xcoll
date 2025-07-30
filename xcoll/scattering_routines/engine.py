# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

import os
import numpy as np
from numbers import Number
from functools import wraps

import xobjects as xo
import xtrack as xt
import xtrack.particles.pdg as pdg

try:
    # TODO: once xaux is in Xsuite keep only this
    from xaux import FsPath, ranID
except (ImportError, ModuleNotFoundError):
    from ..xaux import FsPath, ranID


class BaseEngine(xo.HybridClass):
    _xofields = {
        '_particle_ref': xt.Particles._XoStruct,
        '_seed':         xo.UInt64,
        '_capacity':     xo.Int64,
    }

    _int32 = False
    _only_protons = False
    _element_classes = None
    _uses_input_file = False
    _num_input_files = 1
    _uses_run_folder = False

    def __init__(self, **kwargs):
        if self._element_classes is None:
            raise NotImplementedError(f"{self.__class__.__name__} needs to define `_element_classes`!")
        # Initialise defaults
        self._cwd = None
        self._line = None
        self._verbose = False
        self._input_file = None
        self._element_dict = {}
        self._warning_given = False
        self._environment = None
        self._element_index = 0
        self._tracking_initialised = False
        self._deactivated_elements = {}
        kwargs.setdefault('_particle_ref', xt.Particles())
        kwargs.setdefault('_seed', 0)
        kwargs.setdefault('_capacity', 0)
        super().__init__(**kwargs)

    def __del__(self, *args, **kwargs):
        self.stop(warn=False)

    def _warn(self, error=None):
        if not self._warning_given:
            print(f"Warning: Failed to import {self.__class__.__name__} environment "
                + f"(did you compile?).\n{self.name.capitalize()} elements can be installed "
                + f"but are not trackable.", flush=True)
            if error:
                print(f"Error: {error}", flush=True)
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
            if not isinstance(val, xt.Particles):
                raise ValueError("`particle_ref` has to be an xt.Particles object!")
            if val._capacity > 1:
                raise ValueError("`particle_ref` has to be a single particle!")
            pdg_id = val.pdg_id[0]
            if pdg_id == 0:
                if self._only_protons:
                    pdg_id = pdg.get_pdg_id_from_name('proton')
                else:
                    raise ValueError(f"{self.__class__.__name__} allows the use of particles "
                                   + f"different than protons. Hence, `particle_ref` "
                                   + f"needs to have a valid pdg_id.")
            elif self._only_protons and pdg_id != pdg.get_pdg_id_from_name('proton'):
                raise ValueError("{self.__class__.__name__} only supports protons!")
            self._particle_ref = val
            self._particle_ref.pdg_id[0] = pdg_id

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
            raise ValueError("`capacity` has to be a positive integer!")
        self._capacity = int(val)

    @capacity.deleter
    def capacity(self):
        self.capacity = None

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


    # ======================
    # === Public Methods ===
    # ======================

    def start(self, *, clean=True, input_file=None, **kwargs):
        if not self.environment:
            raise RuntimeError(f"{self.name.capitalize()} environment not set up! "
                             + f"Do not manually create an instance of the engine.")
        if not self.environment.initialised:
            raise RuntimeError(f"{self.name.capitalize()} environment not initialised! "
                             + f"Please set all paths in the environment before "
                             + f"starting the engine.")
        if not self.environment.compiled:
            raise RuntimeError(f"{self.name.capitalize()} interface not compiled! "
                             + f"Please compile before starting the engine.")
        if self.is_running():
            self._print("Engine already running.")
            return

        # Clean up any leftover failed runs
        self.stop(clean=clean)

        print(f"Starting {self.__class__.__name__}...   ", flush=True, end='')
        kwargs = self._pre_start(**kwargs)

        # This needs to be set in the ChildEngine, either in _start_engine() or at the start of tracking
        self._tracking_initialised = False

        # Set all engine properties that have a setter (this will remove these properties from the kwargs)
        kwargs = self._set_engine_properties(**kwargs)

        # Create input file if needed (this will remove the kwargs relevant to the input file and physics)
        kwargs = self._use_input_file(input_file, **kwargs)
        if clean:
            self.clean_input_files(clean_all=False)
        self._preparing_input = False

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
        self._restore_engine_properties(clean=clean)
        self._warning_given = False
        self._tracking_initialised = False

    def is_running(self):
        if hasattr(self, '_preparing_input') and self._preparing_input:
            # We need this to allow changing the element settings which otherwise are locked
            return False
        # If we get here, we cannot say if the engine is running or not and we need an
        # implementation in the child class
        return self._is_running()


    def generate_input_file(self, *, clean=True, filename=None, **kwargs):
        # This method manually generates an input file without starting the engine
        if not self._uses_input_file:
            raise ValueError(f"{self.__class__.__name__} does not use input files!")
        if self._element_dict:
            raise ValueError("Elements already assigned to engine (cannot regenerate input "
                           + "file after starting engine)!")

        # Set all engine properties that have a setter (this will remove these properties from the kwargs)
        cwd = kwargs.get('cwd')
        if cwd and FsPath(cwd).resolve() == FsPath.cwd().resolve():
            print("Warning: Cannot use current working directory as input folder.")
            kwargs.pop('cwd')
        kwargs = self._set_engine_properties(**kwargs)

        # Create input file
        input_file, _ = self._generate_input_file(**kwargs)
        if not hasattr(input_file, '__iter__') or isinstance(input_file, str):
            # Some engines might create multiple input files (like Fluka)
            input_file = [input_file]
        if filename is None:
            if hasattr(self, '_old_cwd') and self._old_cwd is not None:
                new_input_file = [input_file[0].rename(self._old_cwd / input_file[0].name)]
            else:
                new_input_file = [input_file[0]]
        else:
            new_input_file = [input_file[0].rename(filename)]
        for i, file in enumerate(input_file[1:]):
            new_input_file.append(file.rename(new_input_file[0].parent / file.name))

        if clean:
            self.clean_input_files(clean_all=True, input_file=input_file)
        self._restore_engine_properties(clean=clean)

        return new_input_file


    def assert_particle_ref(self):
        if self.particle_ref is None:
            raise ValueError(f"{self.__class__.__name__} reference particle not set!")

    def assert_ready_to_track_or_skip(self, coll, particles, _necessary_attributes=[]):
        self._assert_element(coll)

        missing_attributes = False
        for attr in _necessary_attributes:
            if not hasattr(coll, attr) or not getattr(coll, attr):
                missing_attributes = True

        if not coll.active or not coll._tracking or not coll.jaw or missing_attributes:
            coll._equivalent_drift.track(particles)
            return True

        npart = particles._num_active_particles
        if npart == 0:
            return True
        if not isinstance(particles._buffer.context, xo.ContextCpu):
            raise ValueError(f"{self.__class__.__name__} only supports CPU contexts!")

        assert self.environment.compiled
        if not self.is_running():
            raise ValueError(f"{self.__class__.__name__} not yet running!\nPlease do this "
                           + f"first, by calling xcoll.{self.__class__.__name__}.start().")

        self.assert_particle_ref()
        if abs(particles.mass0 - self.particle_ref.mass0) > 1e-3:
            raise ValueError(f"Error in reference mass of `particles`: not in sync with "
                           + f"{self.name} reference particle!\nRebuild the particles object "
                           + f"using the {self.__class__.__name__} reference particle.")
        if abs(particles.q0 - self.particle_ref.q0) > 1e-3:
            raise ValueError(f"Error in reference charge of `particles`: not in sync with "
                           + f"{self.name} reference particle!\nRebuild the particles object "
                           + f"using the {self.__class__.__name__} reference particle.")
        if not self._only_protons and np.any([pdg_id == 0 for pdg_id in particles.pdg_id]):
            raise ValueError("Some particles are missing the pdg_id!")
        return False


    def clean(self, **kwargs):
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
        self._preparing_input = True  # We need this to allow changing the element settings which otherwise are locked
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
        return kwargs

    def _restore_engine_properties(self, clean=False):
        self._preparing_input = False
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
                if self._int32:
                    self.seed = np.random.randint(0, int(2**32))
                else:
                    self.seed = np.random.randint(0, int(2**64))
        else:
            self._old_seed = self.seed
            self.seed = seed
        self._print(f"Using seed {self.seed}.")

    def _use_particle_ref(self, particle_ref=None):
        # Prefer: provided particle_ref > existing particle_ref > particle_ref from line
        if particle_ref is not None:
            self._old_particle_ref = self.particle_ref
            self.particle_ref = particle_ref
        elif self.particle_ref is None:
            if self.line is None or not hasattr(self.line, 'particle_ref') \
            or self.line.particle_ref is None:
                raise ValueError("Need to provide either a line with a reference "
                               + "particle, or `particle_ref`.")
            self._old_particle_ref = self.particle_ref
            self.particle_ref = self.line.particle_ref
        self._print(f"Using {pdg.get_name_from_pdg_id(self.particle_ref.pdg_id[0])} "
                  + f"with momentum {self.particle_ref.p0c[0]/1.e9:.1f} GeV.")

    def _sync_line_particle_ref(self):
        if self.line is None:
            return
        if self.line.particle_ref is not None \
        and not xt.line._dicts_equal(self.line.particle_ref.to_dict(),
                                     self.particle_ref.to_dict()):
            self._print("Found different reference particle in line. Temporarily overwritten.")
            self._old_line_particle_ref = self.line.particle_ref
            self.line.particle_ref = self.particle_ref

    def _get_new_element_name(self):
        name = f"{self.name}_el_{self._element_index}"
        self._element_index += 1
        return name

    def _deactivate_element(self, el):
        self._deactivated_elements[el.name] = [el, el.active or True]
        if hasattr(el, 'active'):
            el.active = False
        self._remove_element(el)

    def _assert_element(self, element):
        if not isinstance(element, self._element_classes):
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
                raise ValueError("Need to provide either `line` or `elements`.")
            if names is None:
                names = []
                for ee in elements:
                    if hasattr(ee, 'name'):
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
                raise ValueError("Length of `elements` and `names` doesn't match.")
        else:
            if elements is not None:
                raise ValueError("Cannot provide both `line` and `elements`.")
            if names is None:
                elements, names = self.line.get_elements_of_type(self._element_classes)
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
            raise ValueError(f"No active {self.name} elements found!")
        self._element_dict = dict(zip(this_names, this_elements))

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
                            raise ValueError(f"Too many folders with the same "
                                           + f"name {cwd}!")
                    cwd = cwd.parent / f'{cwd.name}_{i:0>4}'
            else:
                ran_str = ranID(only_alphanumeric=True)
                cwd = FsPath.cwd() / f'{self.name}_run_{ran_str}'
            self._cwd = cwd
            cwd.mkdir(parents=True, exist_ok=True)
            self._old_cwd = FsPath.cwd()
            os.chdir(cwd)

    def _reset_cwd(self, clean=False):
        if hasattr(self, '_old_cwd'):
            if self._old_cwd is not None:
                os.chdir(self._old_cwd)
                if clean and self._cwd is not None and self._cwd != FsPath.cwd():
                    try:
                        self._cwd.rmdir()
                    except:
                        self._print(f"Warning: Failed to remove temporary folder "
                                  + f"{self._cwd}.")
            del self._old_cwd
        self._cwd = None


    def _use_input_file(self, input_file=None, **kwargs):
        if self._uses_input_file:
            if input_file is None:
                input_file, kwargs = self._generate_input_file(**kwargs)
            if not hasattr(input_file, '__iter__') or isinstance(input_file, str):
                # Some engines might need multiple input files (like Fluka)
                input_file = [input_file]
            input_file = [FsPath(f).expanduser().resolve() for f in input_file]
            new_files = []
            for file in input_file:
                if not file.exists():
                    raise ValueError(f"Input file {file.as_posix()} not found!")
                if file.parent != FsPath.cwd() and self._uses_run_folder:
                    file.copy_to(FsPath.cwd(), method='mount')
                    new_files.append(FsPath.cwd() / file.name)
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
            if cwd is None and self._cwd is not None:
                cwd = self._cwd
        else:
            if not hasattr(input_file, '__iter__') or isinstance(input_file, str):
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
