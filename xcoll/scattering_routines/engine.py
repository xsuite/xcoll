# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

import os
import numpy as np

import xobjects as xo
from xobjects.hybrid_class import MetaHybridClass
import xtrack as xt
import xtrack.particles.pdg as pdg

try:
    # TODO: once xaux is in Xsuite keep only this
    from xaux import ClassProperty, ClassPropertyMeta, FsPath, singleton, ranID
except (ImportError, ModuleNotFoundError):
    from ..xaux import ClassProperty, ClassPropertyMeta, FsPath, singleton, ranID


class BaseEngineMeta(MetaHybridClass, ClassPropertyMeta):
    def __new__(cls, name, bases, data):
        new_class = MetaHybridClass.__new__(cls, name, bases, data)
        return ClassPropertyMeta.__new__(cls, name, bases, data, new_class)


@singleton(allow_underscore_vars_in_init=False)
class BaseEngine(xo.HybridClass, metaclass=BaseEngineMeta):
    _xofields = {
        '_particle_ref': xt.Particles._XoStruct,
        '_seed':         xo.UInt64,
        '_capacity':     xo.Int64,
    }

    _int32 = False
    _only_protons = False
    _element_classes = None
    _uses_input_file = False
    _uses_run_folder = False

    def __init__(self, **kwargs):
        if self._element_classes is None:
            raise NotImplementedError(f"{self.__class__.__name__} needs to define `_element_classes`!")
        if np.any([key[0] != '_' for key in self._xofields.keys()]):
            raise ValueError(f"All fields in `{self.__class__.__name__}._xofields` have "
                            + f"to start with an underscore! This is to ensure to work "
                            + f"correctly with `ClassProperty`.")
        # Initialise defaults
        self._cwd = None
        self._line = None
        self._verbose = False
        self._input_file = None
        self._element_dict = {}
        self._warning_given = False
        self._tracking_initialised = False
        kwargs.setdefault('_particle_ref', xt.Particles())
        kwargs.setdefault('_seed', 0)
        kwargs.setdefault('_capacity', 0)
        super().__init__(**{key: value for key, value in kwargs.items()
                            if key in self._xofields.keys() or key == '_xobject'})

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

    @ClassProperty
    def name(cls):
        return cls.__name__.replace('Engine', '').lower()

    @ClassProperty
    def verbose(cls):
        return cls.get_self()._verbose

    @verbose.setter
    def verbose(cls, val):
        cls.get_self()._verbose = val

    @ClassProperty
    def line(cls):
        return cls.get_self()._line

    @line.setter
    def line(cls, val):
        if not val is None and not isinstance(val, xt.Line):
            raise ValueError("`line` has to be an xt.Line object!")
        cls.get_self()._line = val

    @line.deleter
    def line(cls):
        cls.get_self()._line = None

    @ClassProperty
    def particle_ref(cls):
        self = cls.get_self()
        initial = xt.Particles().to_dict()
        current = self._particle_ref.to_dict()
        if xt.line._dicts_equal(initial, current):
            return None
        else:
            return self._particle_ref

    @particle_ref.setter
    def particle_ref(cls, val):
        self = cls.get_self()
        if val is None:
            self._particle_ref = xt.Particles()
        else:
            if not isinstance(val, xt.Particles):
                raise ValueError("`particle_ref` has to be an xt.Particles object!")
            if val._capacity > 1:
                raise ValueError("`particle_ref` has to be a single particle!")
            pdg_id = val.pdg_id[0]
            if pdg_id == 0:
                if cls._only_protons:
                    pdg_id = pdg.get_pdg_id_from_name('proton')
                else:
                    raise ValueError(f"{cls.__name__} allows the use of particles "
                                   + f"different than protons. Hence, `particle_ref` "
                                   + f"needs to have a valid pdg_id.")
            elif cls._only_protons and pdg_id != pdg.get_pdg_id_from_name('proton'):
                raise ValueError("{cls.__name__} only supports protons!")
            self._particle_ref = val
            self._particle_ref.pdg_id[0] = pdg_id

    @particle_ref.deleter
    def particle_ref(cls):
        cls.particle_ref = None

    @ClassProperty
    def capacity(cls):
        self = cls.get_self()
        if self._capacity == 0:
            return None
        else:
            return int(self._capacity)

    @capacity.setter
    def capacity(cls, val):
        if val is None:
            val = 0
        cls.get_self()._capacity = int(val)

    @capacity.deleter
    def capacity(cls):
        cls.get_self().capacity = None

    @ClassProperty
    def seed(cls):
        self = cls.get_self()
        if self._seed == 0:
            return None
        else:
            return self._seed

    @seed.setter
    def seed(cls, val):
        self = cls.get_self()
        if val is None:
            val = 0
        val = int(val)
        if cls._int32:
            new_val = np.uint32(val)
        else:
            new_val = np.uint64(val)
        if new_val != val:
            self._print(f"Warning: type change for seed {val}. Using {new_val}.")
        self._seed = new_val

    @seed.deleter
    def seed(cls):
        cls.seed = None

    @ClassProperty
    def input_file(cls):
        if cls._uses_input_file:
            return cls.get_self()._input_file


    # ======================
    # === Public Methods ===
    # ======================

    @classmethod
    def start(cls, *, line=None, elements=None, names=None, cwd=None, seed=None,
              particle_ref=None, input_file=None, clean=True, **kwargs):
        self = cls.get_self(**kwargs)
        kwargs, _ = cls.filter_kwargs(**kwargs)
        for key in kwargs.keys():
            if key.startswith('_'):
                raise ValueError(f"Unknown keyword argument '{key}'!")

        if self.is_running():
            self._print("Engine already running.")
            return

        # Clean up any leftover failed runs
        cls.stop(clean=clean)

        if self.verbose:
            print(f"Starting {cls.__name__}...", flush=True)
        else:
            print(f"Starting {cls.__name__}...   ", flush=True, end='')
        self._pre_start(**kwargs)

        # This needs to be set in the ChildEngine, either in _start_engine() or at the start of tracking
        self._tracking_initialised = False

        self._prepare_input(line=line, elements=elements, names=names, cwd=cwd, seed=seed,
                            particle_ref=particle_ref, **kwargs)
        self._use_input_file(input_file, **kwargs)
        if clean:
            self.clean_input_files(clean_all=False)
        self._preparing_input = False
        self._start_engine(**kwargs)
        if self.verbose:
            print(f"{cls.__name__} started.", flush=True)
        else:
            print(f"Done.", flush=True)


    @classmethod
    def stop(cls, clean=False, **kwargs):
        self = cls.get_self(**kwargs)
        kwargs, _ = cls.filter_kwargs(**kwargs)
        self._stop_engine(**kwargs)
        if clean:
            self.clean(clean_all=True, **kwargs)
        self._restore_input(clean=clean)
        self._warning_given = False
        self._tracking_initialised = False


    @classmethod
    def assert_particle_ref(cls):
        if cls.get_self().particle_ref is None:
            raise ValueError(f"{cls.__name__} reference particle not set!")


    @classmethod
    def generate_input_file(cls, *, line=None, elements=None, names=None, cwd=None, seed=None,
                            particle_ref=None, filename=None, **kwargs):
        # This method manually generates an input file without starting the engine
        self = cls.get_self(**kwargs)
        kwargs, _ = cls.filter_kwargs(**kwargs)
        if not self._uses_input_file:
            raise ValueError(f"{cls.__name__} does not use input files!")
        if self._element_dict:
            raise ValueError("Elements already assigned to engine (cannot regenerate input "
                           + "file after starting engine)!")
        self._prepare_input(line=line, elements=elements, names=names, cwd=cwd, seed=seed,
                            particle_ref=particle_ref, **kwargs)
        input_file = self._generate_input_file(**kwargs)
        if not hasattr(input_file, '__iter__') or isinstance(input_file, str):
            # Some engines might need multiple input files (like Fluka)
            input_file = [input_file]
        if filename is None:
            if hasattr(self, '_old_cwd') and self._old_cwd is not None:
                input_file[0] = input_file[0].rename(self._old_cwd / input_file[0].name)
        else:
            input_file[0] = input_file[0].rename(filename)
        for i, file in enumerate(input_file[1:]):
            input_file[i+1] = file.rename(input_file[0].parent / file.name)
        if len(input_file) == 1:
            input_file = input_file[0]
        self.clean(clean_all=True)
        self._restore_input(clean=True)
        return input_file


    @classmethod
    def is_running(cls, **kwargs):
        self = cls.get_self() # Do not pass kwargs! We are dealing with the setters manually
        if hasattr(self, '_preparing_input') and self._preparing_input:
            # We need this to allow changing the element settings which otherwise are locked
            return False
        # If we get here, we cannot say if the engine is running or not and we need an
        # implementation in the child class
        return self._is_running(**kwargs)


    @classmethod
    def clean(cls, **kwargs):
        cls.clean_input_files(**kwargs)
        cls.clean_output_files(**kwargs)

    @classmethod
    def clean_input_files(cls, **kwargs):
        self = cls.get_self()
        kwargs = self._get_input_cwd_for_cleaning(**kwargs)
        self._clean_input_files(**kwargs)

    @classmethod
    def clean_output_files(cls, **kwargs):
        self = cls.get_self()
        kwargs = self._get_input_cwd_for_cleaning(**kwargs)
        self._clean_output_files(**kwargs)


    # =======================
    # === Private Methods ===
    # =======================

    def _prepare_input(self, *, line=None, elements=None, names=None, cwd=None, seed=None,
              particle_ref=None, **kwargs):
        self._preparing_input = True  # We need this to allow changing the element settings which otherwise are locked
        self._use_seed(seed)
        self._use_line(line)
        self._use_particle_ref(particle_ref)
        self._sync_line_particle_ref()
        self._get_elements(elements, names)
        self._set_cwd(cwd)

    def _restore_input(self, clean=False):
        self._preparing_input = False
        self._reset_seed()
        self._reset_particle_ref() # This has to be done before resetting the line
        self._reset_line()
        self._reset_cwd(clean=clean)
        self._input_file = None
        self._element_dict = {}

    # For all the following fields, they can either be set in advance on the engine,
    # or they can be set when the engine is started. In the latter case, the values
    # are temporary and the original will be restored when the engine is stopped.

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

    def _reset_seed(self):
        if hasattr(self, '_old_seed'):
            self.seed = self._old_seed
            del self._old_seed

    def _use_line(self, line=None):
        if line is not None:
            self._old_line = self.line
            self.line = line

    def _reset_line(self):
        if hasattr(self, '_old_line'):
            self.line = self._old_line
            del self._old_line

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

    def _reset_particle_ref(self):
        if hasattr(self, '_old_particle_ref'):
            self.particle_ref = self._old_particle_ref
            del self._old_particle_ref
        if hasattr(self, '_old_line_particle_ref'):
            self.line.particle_ref = self._old_line_particle_ref
            del self._old_line_particle_ref

    def _sync_line_particle_ref(self):
        if self.line is None:
            return
        if self.line.particle_ref is not None \
        and not xt.line._dicts_equal(self.line.particle_ref.to_dict(),
                                     self.particle_ref.to_dict()):
            self._print("Found different reference particle in line. Temporarily overwritten.")
            self._old_line_particle_ref = self.line.particle_ref
            self.line.particle_ref = self.particle_ref

    def _get_elements(self,elements=None, names=None):
        if elements is not None and (not hasattr(elements, '__iter__') or isinstance(elements, str)):
            elements = [elements]
        if names is not None and (not hasattr(names, '__iter__') or isinstance(names, str)):
            names = [names]
        if self.line is None:
            if elements is None:
                raise ValueError("Need to provide either `line` or `elements`.")
            if names is None:
                names = [f"{self.name}_el_{i}" for i, _ in enumerate(elements)]
            if len(names) != len(elements):
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
            if not isinstance(ee, self._element_classes):
                raise ValueError(f"Element {ee} is not a "
                                + ", or a ".join([c.__name__ for c in self._element_classes])
                                + ".")
            if ee.jaw is None:
                self._print(f"Warning: Jaw not set for {name}. Ignoring.")
                self._remove_element(name, ee)
            elif not ee.active:
                self._print(f"Warning: Element {name} is not active. Ignoring.")
                self._remove_element(name, ee)
            else:
                this_names.append(name)
                this_elements.append(ee)
        if len(this_elements) == 0:
            raise ValueError(f"No active {self.name} elements found!")
        self._element_dict = dict(zip(this_names, this_elements))

    def _set_cwd(self, cwd=None):
        if self._uses_run_folder:
            if cwd is not None:
                cwd = FsPath(cwd).expanduser().resolve()
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
                input_file = self._generate_input_file(**kwargs)
            if not hasattr(input_file, '__iter__') or isinstance(input_file, str):
                # Some engines might need multiple input files (like Fluka)
                input_file = [input_file]
            input_file = [FsPath(f).expanduser().resolve() for f in input_file]
            new_files = []
            for file in input_file:
                if not file.exists():
                    raise ValueError(f"Input file {file.as_posix()} not found!")
                if file.parent != FsPath.cwd() and self._uses_run_folder:
                    file.copy_to(FsPath.cwd())
                    new_files.append(FsPath.cwd() / file.name)
                else:
                    new_files.append(file)
            self._input_file = new_files[0] if len(new_files)==1 else new_files
            self._match_input_file()


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

    def _clean_input_files(self, input_file=None, cwd=None, clean_all=False, **kwargs):
        pass

    def _clean_output_files(self, input_file=None, cwd=None, clean_all=False, **kwargs):
        pass

    def _remove_element(self, name=None, el=None, **kwargs):
        pass

    def _pre_start(self, **kwargs):
        pass
