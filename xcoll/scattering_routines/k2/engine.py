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

from ...sixtrack_input import create_dat_file


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
            print("Warning: Failed to run pyk2_init (did you compile?) \n")
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
        from ...beam_elements.k2 import _K2Collimator # should this be here..?
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
        
        # prepare for init
        tw = line.twiss()
        elements, names = line.get_elements_of_type(_K2Collimator)
        if not hasattr(elements, '__iter__') or isinstance(elements, str):
            elements = [elements]
        if not hasattr(names, '__iter__') or isinstance(names,str):
            names = [names]
        nemitt_x = None
        nemitt_y = None
        for el in elements:
            if hasattr(el, 'optics') and el.optics is not None:
                if nemitt_x is None:
                    nemitt_x = el.nemitt_x
                if nemitt_y is None:
                    nemitt_y = el.nemitt_y
                if not np.isclose(el.nemitt_x, nemitt_x) \
                or not np.isclose(el.nemitt_x, nemitt_x):
                    raise ValueError("Not all collimators have the same "
                                   + "emittance. This is not supported.")
                else: 
                    emittance = np.sqrt(nemitt_x**2 + nemitt_y**2)

        # TODO: make
        file = create_dat_file(line=line, cwd=cwd, elements=elements, names=names)
        if not file.exists():
            raise ValueError(f"File {file.as_posix()} not found!")
        if file.parent != Path.cwd():
            shutil.copy(file, Path.cwd())
            file = Path.cwd() / file.name

        try:
            pyk2_init(n_alloc=this.n_alloc, file=file, random_generator_seed=seed, \
                     num_coll=len(elements), betax=tw.betx, betay=tw.bety, alphax=tw.alpx, \
                     alphay=tw.alpy, orbx=tw.x, orby=tw.y, orbxp=tw.xp, orbyp=tw.yp, \
                     gamma=(np.sqrt(tw.gamx**2 + tw.gamy**2)), emit=emittance)
        except ValueError:
           this._warn(init=True)
           return      
  
    # TODO: Would be nice; but works different from FLUKA 
    def is_running(self):
        if self._initialised:
            return True
        else:
            return False

    def __setattribute__(self, name, value):
        # if name in ['gap', 'gap_L', 'gap_R', 'jaw', 'jaw_L', 'jaw_R', 'jaw_LU', 'jaw_LD', 'jaw_RU', 
        # 'jaw_RD', 'tilt', 'tilt_L', 'tilt_R']:
        if K2Engine.is_running():
            raise ValueError('Engine is running; K2Collimator is frozen.')
        super().__setattribute__(name, value)
    