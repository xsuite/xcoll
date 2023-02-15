# copyright ############################### #
# This file is part of the Xcoll Package.  #
# Copyright (c) CERN, 2023.                 #
# ######################################### #

import numpy as np

import xobjects as xo
import xtrack as xt

from ..tables import CollimatorImpacts
from ..general import _pkg_root

# class MetaCollimator(xt.base_element.MetaBeamElement, ABCMeta):
#     pass



class InvalidCollimator(xt.BeamElement):
    _xofields = {
        'length': xo.Float64
    }

    isthick = True
    behaves_like_drift = True
    skip_in_loss_location_refinement = True

    # InvalidCollimator catches unallowed cases, like backtracking through a collimator
    _extra_c_sources = [
        _pkg_root.joinpath('headers','particle_states.h'),
        _pkg_root.joinpath('beam_elements','collimators_src','invalid_collimator.h')
    ]

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def get_backtrack_element(self, _context=None, _buffer=None, _offset=None):
        return self.__class__(length=-self.length,
                              _context=_context, _buffer=_buffer, _offset=_offset)



class BaseCollimator(xt.BeamElement):#, metaclass=MetaCollimator):
    _xofields = {
        'inactive_front': xo.Float64,
        'active_length': xo.Float64,
        'inactive_back': xo.Float64,
        'jaw_F_L': xo.Float64,
        'jaw_F_R': xo.Float64,
        'jaw_B_L': xo.Float64,
        'jaw_B_R': xo.Float64,
        'dx': xo.Float64,
        'dy': xo.Float64,
        'cos_z': xo.Float64,
        'sin_z': xo.Float64,
        '_active': xo.Int8
    }

    isthick = True
    behaves_like_drift = True
    skip_in_loss_location_refinement = True
    
    _skip_in_to_dict  = ['_active', 'cos_z', 'sin_z']
    _store_in_to_dict = ['is_active', 'angle']
    _internal_record_class = CollimatorImpacts

    _extra_c_sources = [
        _pkg_root.joinpath('beam_elements','collimators_src','base_collimator.h')
    ]

    _depends_on = [InvalidCollimator, xt.Drift, xt.XYShift, xt.SRotation]

    def __init__(self, **kwargs):
        # TODO: quick hack to avoid instantiation; did not manage to get it to work correclty with ABC
        if self.__class__.__name__ == 'BaseCollimator':
            raise Exception("Abstract class 'BaseCollimator' cannot be instantiated!")
        if '_xobject' not in kwargs:
            kwargs.setdefault('jaw_F_L', 1)
            kwargs.setdefault('jaw_F_R', -1)
            kwargs.setdefault('jaw_B_L', 1)
            kwargs.setdefault('jaw_B_R', -1)
            kwargs.setdefault('inactive_front', 0)
            kwargs.setdefault('inactive_back', 0)
            kwargs.setdefault('dx', 0)
            kwargs.setdefault('dy', 0)
            angle = kwargs.pop('angle', 0)
            anglerad = angle / 180. * np.pi
            kwargs['cos_z'] = np.cos(anglerad)
            kwargs['sin_z'] = np.sin(anglerad)
            is_active = kwargs.pop('is_active', True)
            is_active = 1 if is_active == True else is_active
            is_active = 0 if is_active == False else is_active
            kwargs['_active'] = is_active
        super().__init__(**kwargs)


    @property
    def angle(self):
        return np.arctan2(self.sin_z, self.cos_z) * (180.0 / np.pi)

    @angle.setter
    def angle(self, angle):
        anglerad = angle / 180. * np.pi
        self.cos_z = np.cos(anglerad)
        self.sin_z = np.sin(anglerad)

    @property
    def is_active(self):
        return True if self._active == 1 else False

    @is_active.setter
    def is_active(self, is_active):
        is_active = 1 if is_active == True else is_active
        is_active = 0 if is_active == False else is_active
        self._active = is_active

    @property
    def length(self):
        return (self.inactive_front + self.active_length + self.inactive_back)
