# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

import numpy as np
import os
from pathlib import Path
import shutil

import xobjects as xo
import xpart as xp
import xtrack as xt

from create_dat import create_dat_file
from ...beam_elements.k2 import _K2Collimator

class K2Engine(xo.HybridClass):
    _xofields = {
        '_capacity':              xo.Int64,
        'random_generator_seed':  xo.Int64
    }

    # The engine is a singleton
    def __new__(cls, **kwargs):
        if not hasattr(cls, 'instance'):
            cls.instance = super().__new__(cls)
            cls.instance._initialised = False
        return cls.instance


    def __init__(self, **kwargs):
        if(self._initialised):
            return
        self._initialised = True
        self._warning_given = False
        kwargs.setdefault('_capacity', 50000)
        kwargs.setdefault('random_generator_seed', None)  # Allow seed to be set to None to get default:
        if kwargs['random_generator_seed'] is None:
            kwargs['random_generator_seed'] = np.random.randint(1, 10000000)
        if '_xobject' not in kwargs:
            # Initialise defaults  
            self._cwd = None               # TODO: this is needed right
            self._input_file = None
            self._warning_given = False
            self._collimator_dict = {}
        super().__init__(**kwargs)
        K2Engine.reset()

    def _warn(self, init=False):
        if init == True and not self.instance._warning_given:
            print("Warning: Failed to run pyk2_init.\n")
            self._warning_given = True
        if not self.instance._warning_given:
            print("Warning: Failed to import pyK2 (did you compile?).\n"
                + "_K2collimators will be installed but are not trackable.")
            self._warning_given = True


    @classmethod
    def reset(cls, seed=None):
        if seed is not None:
            cls.instance.random_generator_seed = seed
        try:
            from .pyk2f import pyk2_init
        except ImportError:
            cls.instance._warn()
        else:
            pyk2_init(n_alloc=cls.instance._capacity, random_generator_seed=cls.instance.random_generator_seed)

    @classmethod
    def start(cls, *, input_file=None, line=None, elements=None, seed=None, names=None, cwd=None,
              prototypes_file=None, include_files=None, debug_level=0, **kwargs):
        
        cls(**kwargs)
        this = cls.instance
        # should we check if not this._k2.exists()? 
        if seed is not None:
            cls.instance.random_generator_seed = seed
        try:
            from .pyk2f import pyk2_init
        except ImportError:
            cls.instance._warn()
        else:
            pyk2_init(n_alloc=cls.instance._capacity, random_generator_seed=cls.instance.random_generator_seed)
        
        if cwd is not None:
            cwd = Path(cwd).resolve()
            cwd.mkdir(parents=True, exist_ok=True)
            this._old_cwd = Path.cwd()
            os.chdir(cwd)
        else:
            cwd = Path.cwd()

        this._cwd = cwd

        if input_file is not None:
            input_file = Path(input_file).resolve()
        else:
            if line is None:
                raise ValueError("Either 'input_file' or 'line' must be given.") # must we not always need line
            # if prototypes_file is None:
            #     print("Using default prototypes file.")
            #     prototypes_file = _pkg_root / 'scattering_routines' / 'k2' / 'prototypes.dat'
            input_file = create_dat_file(line, cwd=None, filename=None) # TODO: make
        if not input_file.exists():
            raise ValueError(f"Input file {input_file.as_posix()} not found!")
        if input_file.parent != Path.cwd():
            shutil.copy(input_file, Path.cwd())
            input_file = Path.cwd() / input_file.name
        this._input_file = input_file

        try:
            # TODO: need to find a way to get optics without line
            # num_coll from line? input file?
            pyk2_init(n_alloc=this.n_alloc, colldb_input_fname, random_generator_seed=seed, num_coll, \
                      betax, betay, alphax, alphay, orbx, orby, orbxp, orbyp, gamma, emit)
        except ImportError as error:
            this._warn(init=True)
            return

        # Match collimators 
        if line is None:
            if elements is None and names is None:
                raise ValueError("Either 'line' or 'elements' and 'names' must be given.")
        else:
            elements, names = line.get_elements_of_type(_K2Collimator)
        this._match_collimators_to_engine(elements, names)

        # reference particle
        initial = xp.Particles().to_dict()
        current = line.particle_ref.to_dict()
        if xt.line._dicts_equal(initial, current):
            if line is not None:
                if line.particle_ref is None:
                    print("Warning: This line does not have any reference particle. Set it later.")
                else:
                    this.set_particle_ref(line=line)
            else:
                print("No line given. Set reference particle later.")
    
    @classmethod
    def set_particle_ref(cls, particle_ref=None,line=None,**kwargs): # TODO: is this reallllly needed
        cls(**kwargs)
        this = cls.instance
        overwrite_particle_ref = False

        if particle_ref is None:
            if line is None or line.particle_ref is None:
                raise ValueError("Line has no reference particle. Please provide one with '.particle_ref'.")
            if line.particle_ref.pdg_id == 0:
                raise ValueError("'line.particle_ref' needs a valid pdg_id.")
            particle_ref = line.particle_ref
        else:  
            if particle_ref._capacity > 1: 
                raise ValueError("`particle_ref` has to be a single particle!")
            if particle_ref.pdg_id == 0:
                raise ValueError("'particle_ref' needs a valid pdg_id.")
            if line is not None and not \
                xt.line._dicts_equal(line.particle_ref.to_dict(), particle_ref.to_dict()):
                overwrite_particle_ref = True
        this._particle_ref = particle_ref
        if overwrite_particle_ref:
            print("Warning: Different reference particles found in line! Overwritten.")
            line.particle_ref = particle_ref

    def _match_collimators_to_engine(self, elements, names):
        if not hasattr(elements, '__iter__') or isinstance(elements, str):
            elements = [elements]
        if not hasattr(names, '__iter__') or isinstance(names, str):
            names = [names]
        assert len(elements) == len(names), "Number of elements and names must match!"

        for name in names:
            if name not in self._collimator_dict:
                raise ValueError(f"K2 Collimator '{name}' not found in input file!")
        for el, name in zip(elements, names):
            offset = self._collimator_dict[name]['offset']
            # TODO check if right sided is onesided etc etc match between line and input file; how :)) 
            
        

    