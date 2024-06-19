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
    def start(cls, *, line=None, seed=None, cwd=None):
        
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

        if line is None:
            raise ValueError("'line' must be given.")
        file = create_dat_file(line, cwd=None, filename=None) # TODO: make
        if not file.exists():
            raise ValueError(f"File {file.as_posix()} not found!")
        if file.parent != Path.cwd():
            shutil.copy(file, Path.cwd())
            file = Path.cwd() / file.name
        try:
            # TODO orbit ?? + emittance is from colldb or per collimator not from line
            tw = line.twiss()
            elements, names = line.get_elements_of_type(_K2Collimator)
            if not hasattr(elements, '__iter__') or isinstance(elements, str):
                elements = [elements]

            pyk2_init(n_alloc=this.n_alloc, file=file, random_generator_seed=seed, num_coll=len(elements), \
                      betax=tw.betx, betay=tw.bety, alphax=tw.alpx, alphay=tw.alpy, orbx=orbx, \
                      orby, orbxp, orbyp, gamma=tw.gamma, emit)

        except ImportError as error:
            this._warn(init=True)
            return      
  
    # TODO: Would be nice; but works different from FLUKA 
    def is_running(self):
        pass
    def __setattribute__(self, name, value):
        # if name in ['gap', 'gap_L', 'gap_R', 'jaw', 'jaw_L', 'jaw_R', 'jaw_LU', 'jaw_LD', 'jaw_RU', 
        # 'jaw_RD', 'tilt', 'tilt_L', 'tilt_R']:
        if K2Engine.is_running():
            raise ValueError('Engine is running; K2Collimator is frozen.')
        super().__setattribute__(name, value)
    