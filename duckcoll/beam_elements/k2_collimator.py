import numpy as np

from xcoll import BaseCollimator, Material, CrystalMaterial
from ..scattering_routines.k2 import track

import xobjects as xo


# TODO: remove dx, dy, offset, tilt, as this should only be in colldb (and here only the jaw positions)
class K2Collimator(BaseCollimator):
    _xofields = { **BaseCollimator._xofields,
        'onesided':   xo.Int8,
        'material':   Material,
        '_tracking':  xo.Int8
    }

    isthick = True
    iscollective = True
    behaves_like_drift = True
    skip_in_loss_location_refinement = True

    _skip_in_to_dict       = BaseCollimator._skip_in_to_dict
    _store_in_to_dict      = BaseCollimator._store_in_to_dict
    _internal_record_class = BaseCollimator._internal_record_class


    def __init__(self, **kwargs):
        if '_xobject' not in kwargs:
            kwargs.setdefault('onesided', False)
            kwargs.setdefault('material', None)
            if kwargs['material'] is None:
                raise ValueError("Need to provide a material to the collimator!")
            kwargs.setdefault('_tracking', True)
        super().__init__(**kwargs)
        if '_xobject' not in kwargs:
            from ..scattering_routines.k2.engine import K2Engine
            K2Engine()  # initialise the engine if it does not exist yet


    def track(self, particles):  # TODO: write impacts
        track(self, particles)
        return



class K2Crystal(BaseCollimator):
    _xofields = { **BaseCollimator._xofields,
        'align_angle': xo.Float64,  #  = - sqrt(eps/beta)*alpha*nsigma
        'bend':        xo.Float64,
        'xdim':        xo.Float64,
        'ydim':        xo.Float64,
        'thick':       xo.Float64,
        'crytilt':     xo.Float64,
        'miscut':      xo.Float64,
        'orient':      xo.Float64,
        'onesided':    xo.Int8,
        'material':    CrystalMaterial,
        '_tracking':  xo.Int8
    }

    isthick = True
    iscollective = True
    behaves_like_drift = True
    skip_in_loss_location_refinement = True

    _skip_in_to_dict       = BaseCollimator._skip_in_to_dict
    _store_in_to_dict      = BaseCollimator._store_in_to_dict
    _internal_record_class = BaseCollimator._internal_record_class


    def __init__(self, **kwargs):
        if '_xobject' not in kwargs:
            kwargs.setdefault('onesided', False)
            if kwargs['material'] is None:
                raise ValueError("Need to provide a material to the collimator!")
            kwargs.setdefault('bend', 0)
            kwargs.setdefault('xdim', 0)
            kwargs.setdefault('ydim', 0)
            kwargs.setdefault('thick', 0)
            kwargs.setdefault('crytilt', 0)
            kwargs.setdefault('miscut', 0)
            kwargs.setdefault('orient', 0)
            kwargs.setdefault('_tracking', True)
        super().__init__(**kwargs)
        if '_xobject' not in kwargs:
            from ..scattering_routines.k2.engine import K2Engine
            K2Engine()  # initialise the engine if it does not exist yet


    def track(self, particles):  # TODO: write impacts
        track(self, particles)
        return

