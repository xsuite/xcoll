# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

import os
import numpy as np
import shutil

import xobjects as xo
import xpart as xp
import xtrack as xt
import xtrack.particles.pdg as pdg

try:
    # TODO: once xaux is in Xsuite keep only this
    from xaux import ClassProperty, ClassPropertyMeta, FsPath, singleton
except ImportError:
    from ..xaux import ClassProperty, ClassPropertyMeta, FsPath, singleton


class BaseEngineMeta(xo.hybrid_class.MetaHybridClass, ClassPropertyMeta):
    pass

@singleton
class BaseEngine(xo.HybridClass, metaclass=BaseEngineMeta):
    _xofields = {
        '_particle_ref': xp.Particles,
        '_seed':         xo.UInt64,
        '_capacity':     xo.Int64,
    }

    _int32 = False
    _element_classes = None
    _only_protons = False
    _uses_input_file = False
    _uses_run_folder = False

    def __init__(self, **kwargs):
        if not self._initialised:
            if np.any([key[0] != '_' for key in self._xofields.keys()]):
                raise ValueError(f"All fields in `{self.__class__.__name__}._xofields` have "
                               + f"to start with an underscore! This is to ensure to work "
                               + f"correctly with `ClassProperty`.")
            if '_xobject' not in kwargs:
                # Initialise defaults
                self._cwd = None
                self._line = None
                self._verbose = False
                self._input_file = None
                self._element_dict = {}
                self._warning_given = False
                self._tracking_initialised = False
                kwargs.setdefault('_particle_ref', xp.Particles())
                kwargs.setdefault('_seed', 0)
                kwargs.setdefault('_capacity', 0)
            filtered_kwargs = {}
            remaining_kwargs = {}
            for key, value in kwargs.items():
                if key in self._xofields.keys() or key == '_xobject':
                    filtered_kwargs[key] = value
                else:
                    remaining_kwargs[key] = value
            super().__init__(**filtered_kwargs)
            kwargs = remaining_kwargs
            self._initialised = True
        # Apply kwargs
        for kk, vv in kwargs.items():
            if not hasattr(self.__class__, kk):
                raise ValueError(f"Invalid attribute {kk} for {self.__class__.__name__}!")
            setattr(self, kk, vv)

    def __del__(self, *args, **kwargs):
        self.stop(warn=False)


    def _warn(self, error=None):
        if not self._warning_given:
            print(f"Warning: Failed to import {self.__class__.__name__} environment "
                + f" (did you compile?).\n{self.__class__.__name__.replace('Engine', '')} "
                + f"elements will be installed but are not trackable.\n", flush=True)
            if error:
                print(error, flush=True)
            self._warning_given = True


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
        initial = xp.Particles().to_dict()
        current = self._particle_ref.to_dict()
        if xt.line._dicts_equal(initial, current):
            return None
        else:
            return self._particle_ref

    @particle_ref.setter
    def particle_ref(cls, val):
        self = cls.get_self()
        if val is None:
            self._particle_ref = xp.Particles()
        else:
            if not isinstance(val, xp.Particles):
                raise ValueError("`particle_ref` has to be an xp.Particles object!")
            if val._capacity > 1:
                raise ValueError("`particle_ref` has to be a single particle!")
            if val.pdg_id[0] == 0:
                if cls._only_protons:
                    val.pdg_id[0] = pdg.get_pdg_id_from_name('proton')
                else:
                    raise ValueError("`particle_ref` needs to have a valid pdg_id")
            elif cls._only_protons and val.pdg_id[0] != pdg.get_pdg_id_from_name('proton'):
                raise ValueError("{cls.__name__} only supports protons!")
            self._particle_ref = val

    @particle_ref.deleter
    def particle_ref(cls):
        cls.get_self()._particle_ref = xp.Particles()

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
        raise ValueError("Not allowed.")

    @ClassProperty
    def seed(cls):
        self = cls.get_self()
        if self._seed == 0:
            return None
        else:
            return self._seed

    @seed.setter
    def seed(cls, val):
        if val is None:
            val = 0
        val = int(val)
        if cls._int32:
            new_val = np.uint32(val)
        else:
            new_val = np.uint64(val)
        if new_val != val:
            print(f"Warning: type change for seed {val}. Using {new_val}.")
        cls.get_self()._seed = new_val

    @seed.deleter
    def seed(cls):
        cls.get_self()._seed = 0

    @ClassProperty
    def input_file(cls):
        return cls.get_self()._input_file


    # ======================
    # === Public Methods ===
    # ======================


    @classmethod
    def start(cls, *, line=None, elements=None, names=None, cwd=None, seed=None,
              particle_ref=None, input_file=None, **kwargs):
        self = cls.get_self(**kwargs)
        if self.is_running() is None:
            raise NotImplementedError(f"Need to implement `is_running` for {cls.__name__}!")
        elif self.is_running() is True:
            if self.verbose:
                print("Engine already running.", flush=True)
            return

        self._starting_engine = True  # We need this to allow changing the element settings which otherwise are locked
        self._use_seed(seed)
        self._use_line(line)
        self._use_particle_ref(particle_ref)
        self._sync_line_particle_ref()
        self._get_elements(line=line, elements=elements, names=names)
        self._set_cwd(cwd=cwd)
        self._use_input_file(input_file=input_file, **kwargs)
        self.clean_output_files()
        self._starting_engine = False


    @classmethod
    def stop(cls, clean=False, **kwargs):
        self = cls.get_self(**kwargs)
        if hasattr(self, '_old_seed'):
            self.seed = self._old_seed
            del self._old_seed
        if hasattr(self, '_old_line'):
            self.line = self._old_line
            del self._old_line
        if hasattr(self, '_old_particle_ref'):
            self.particle_ref = self._old_particle_ref
            del self._old_particle_ref
            self._sync_line_particle_ref()
        if hasattr(self, '_old_cwd') and self._old_cwd is not None:
            os.chdir(self._old_cwd)
            del self._old_cwd
        if clean:
            self.clean_output_files(clean_all=True)
        self._cwd = None
        self._input_file = None
        self._element_dict = {}
        self._warning_given = False
        self._tracking_initialised = False


    @classmethod
    def assert_particle_ref(cls, **kwargs):
        if cls.get_self(**kwargs).particle_ref is None:
            raise ValueError(f"{cls.__name__} reference particle not set!")


    # =======================
    # === Private Methods ===
    # =======================

    # For all the following fields, they can either be set in advance on the engine,
    # or they can be set when the engine is started. In the latter case, the values
    # are temporary and the original will be restored when the engine is stopped.

    def _use_line(self, line=None):
        self._old_line = self.line
        self.line = line

    def _use_seed(self, seed=None):
        self._old_seed = self.seed
        if seed is not None:
            self.seed = seed
        else:
            if self.seed is None:
                if self._int32:
                    self.seed = np.random.randint(0, int(2**32))
                else:
                    self.seed = np.random.randint(0, int(2**64))
        if self.verbose:
            print(f"Using seed {self.seed}.")

    def _use_particle_ref(self, particle_ref=None):
        # Prefer: provided particle_ref > existing particle_ref > particle_ref from line
        self._old_particle_ref = self.particle_ref
        if particle_ref is not None:
            self.particle_ref = particle_ref
        elif self.particle_ref is None:
            if self.line is None or not hasattr(self.line, 'particle_ref') \
            or self.line.particle_ref is None:
                raise ValueError("Need to provide either a line with a reference "
                               + "particle, or `particle_ref`.")
            self.particle_ref = self.line.particle_ref
        if self.verbose:
            print(f"Using {xp.get_name_from_pdg_id(self.particle_ref.pdg_id[0])}.")

    def _sync_line_particle_ref(self):
        if self.line is None:
            return
        if self.line.particle_ref is not None \
        and not xt.line._dicts_equal(self.line.particle_ref.to_dict(),
                                     self.particle_ref.to_dict()):
            overwrite_particle_ref_in_line = True
        if overwrite_particle_ref_in_line:
            print("Warning: Found different reference particle in line! Temporarily overwritten.")
            self.line.particle_ref = self.particle_ref

    def _get_elements(self, line=None, elements=None, names=None):
        if self._element_classes is None:
            raise NotImplementedError(f"{self.__class__.__name__} needs to define `_element_classes`!")
        if line is None:
            if elements is None:
                raise ValueError("Need to provide either `line` or `elements`.")
            if not hasattr(elements, '__iter__') or isinstance(elements, str):
                elements = [elements]
            if names is None:
                names = [f"{self.__class__.name}_el_{i}" for i, _ in enumerate(elements)]
            else:
                if not hasattr(names, '__iter__') or isinstance(names, str):
                    names = [names]
                if len(names) != len(elements):
                    raise ValueError("Length of `elements` and `names` doesn't match.")
        else:
            if elements is not None:
                raise ValueError("Cannot provide both `line` and `elements`.")
            if names is None:
                elements, names = line.get_elements_of_type(self._element_classes)
            else:
                if not hasattr(names, '__iter__') or isinstance(names, str):
                    names = [names]
                elements = [line[name] for name in names]
        elements = [el for el in elements if el.active and el.jaw is not None]
        names    = [name for el, name in zip(elements,names) if el.active and el.jaw is not None]
        for el in elements:
            if not isinstance(el, self._element_classes):
                raise ValueError(f"Element {el} is not a "
                                + ", or a ".join([c.__name__ for c in self._element_classes])
                                + ".")
        if len(elements) == 0:
            raise ValueError(f"No active {self.name} elements found!")
        self._element_dict = dict(zip(names, elements))

    def _set_cwd(self, cwd):
        if self._uses_run_folder:
            if cwd is not None:
                cwd = FsPath(cwd).expanduser().resolve()
            else:
                # TODO: use xaux.ranID here
                import base64
                ran = base64.urlsafe_b64encode(os.urandom(8)).decode('utf-8')
                ran_str = ''.join(c if c.isalnum() else 'X' for c in ran)
                cwd = FsPath.cwd() / f'{self.name}_run_{ran_str}'
            self._cwd = cwd
            cwd.mkdir(parents=True, exist_ok=True)
            self._old_cwd = FsPath.cwd()
            os.chdir(cwd)

    def _use_input_file(self, input_file=None, **kwargs):
        if self._uses_input_file:
            if input_file is None:
                if not hasattr(self, 'generate_input_file'):
                    raise NotImplementedError("Need to implement `generate_input_file`"
                                              "for {cls.__name__}!")
                input_file = self.generate_input_file(**kwargs)
            if not hasattr(input_file, '__iter__') or isinstance(input_file, str):
                # Some engines might need multiple input files (like Fluka)
                input_file = [input_file]
            input_file = [FsPath(f).expanduser().resolve() for f in input_file]
            new_files = []
            for file in input_file:
                if not file.exists():
                    raise ValueError(f"Input file {file.as_posix()} not found!")
                if file.parent != FsPath.cwd():
                    shutil.copy(file, FsPath.cwd())
                    new_files.append(FsPath.cwd() / file.name)
                else:
                    new_files.append(file)
            self._input_file = new_files[0] if len(new_files)==1 else new_files
            self._match_input_file()


    # ===============================
    # === Methods to be inherited ===
    # ===============================

    @classmethod
    def is_running(cls, **kwargs):
        self = cls.get_self(**kwargs)
        if hasattr(self, '_starting_engine') and self._starting_engine:
            # We need this to allow changing the element settings which otherwise are locked
            return False
        # If we get here, we cannot say if the engine is running or not and we need an
        # implementation in the child class
        return None

    @classmethod
    def clean_output_files(cls, **kwargs):
        pass

    def _match_input_file(self, **kwargs):
        raise NotImplementedError("Need to implement `_match_input_file` for {cls.__name__}!")

    @classmethod
    def generate_input_file(cls, **kwargs):
        raise NotImplementedError("Need to implement `generate_input_file` for {cls.__name__}!")
