# copyright ############################### #
# This file is part of the Xcoll Package.  #
# Copyright (c) CERN, 2023.                 #
# ######################################### #

import numpy as np

import xobjects as xo
import xtrack as xt

from ..colldb import _get_LR, _set_LR, _get_LRUD, _set_LRUD
from ..tables import CollimatorImpacts
from ..general import _pkg_root


class InvalidCollimator(xt.BeamElement):
    _xofields = {
        'length': xo.Float64
    }

    isthick = True
    behaves_like_drift = True
    skip_in_loss_location_refinement = True
    allow_backtrack = True

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



class BaseCollimator(xt.BeamElement):
    _xofields = {
        'inactive_front': xo.Float64,  # Drift before jaws
        'active_length':  xo.Float64,  # Length of jaws
        'inactive_back':  xo.Float64,  # Drift after jaws
        'jaw_LU':         xo.Float64,  # jaw left upstream
        'jaw_RU':         xo.Float64,  # jaw right upstream
        'jaw_LD':         xo.Float64,  # jaw left downstream
        'jaw_RD':         xo.Float64,  # jaw right downstream
        'ref_xU':         xo.Float64,  # upstream center of collimator reference frame
        'ref_yU':         xo.Float64,
        'ref_xD':         xo.Float64,  # downstream center of collimator reference frame
        'ref_yD':         xo.Float64,
        'sin_zL':         xo.Float64,  # angle of left jaw
        'cos_zL':         xo.Float64,
        'sin_zR':         xo.Float64,  # angle of right jaw
        'cos_zR':         xo.Float64,
        '_active':        xo.Int8
    }

    isthick = True
    behaves_like_drift = True
    skip_in_loss_location_refinement = True
    
    _skip_in_to_dict  = ['_active', 'jaw_LU', 'jaw_RU', 'jaw_LD', 'jaw_RD',
                         'ref_xU', 'ref_yU', 'ref_xD', 'ref_yD',
                         'sin_zL', 'cos_zL', 'sin_zR', 'cos_zR']
    _store_in_to_dict = ['is_active', 'angle', 'jaw', 'reference_center']
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
            _set_LRUD(kwargs, 'jaw', kwargs.pop('jaw', 1),
                      neg=True, name=self.__class__.__name__)
            _set_LRUD(kwargs, 'reference_center',
                      kwargs.pop('reference_center', 0), name=self.__class__.__name__)
            kwargs.setdefault('inactive_front', 0)
            kwargs.setdefault('inactive_back', 0)
            kwargs['sin_zL'], kwargs['cos_zL'], kwargs['sin_zR'], kwargs['cos_zR'] = _angle_setter(angle)
            is_active = kwargs.pop('is_active', True)
            is_active = 1 if is_active == True else is_active
            is_active = 0 if is_active == False else is_active
            kwargs['_active'] = is_active
        super().__init__(**kwargs)


    @property
    def jaw(self):
        return _get_LRUD(self, 'jaw', neg=True)

    @jaw.setter
    def jaw(self, val):
        _set_LRUD(self, 'jaw', val, neg=True)

    @property
    def angle(self):
        angle_L = round(np.arctan2(self.sin_zL, self.cos_zL) * (180.0 / np.pi), 9)
        angle_R = round(np.arctan2(self.sin_zR, self.cos_zR) * (180.0 / np.pi), 9)
        return angle_L if angle_L==angle_R else [angle_L, angle_R]

    def _angle_setter(self, val):
        if not hasattr(val, '__iter__'):
            val = [val]
        if isinstance(val, str):
            raise ValueError(f"Error in setting angle: not a number!")
        elif len(val) == 2:
            anglerad_L = val[0] / 180. * np.pi
            anglerad_R = val[1] / 180. * np.pi
        elif len(val) == 1:
            anglerad_L = val[0] / 180. * np.pi
            anglerad_R = val[0] / 180. * np.pi
        else:
            raise ValueError(f"Error in setting angle: must have one or two (L, R) values!")
        return np.sin(anglerad_L), np.cos(anglerad_L), np.sin(anglerad_R), np.cos(anglerad_R)

    @angle.setter
    def angle(self, angle):
        self.sin_zL, self.cos_zL, self.sin_zR, self.cos_zR = _angle_setter(angle)

    @property
    def reference_center(self):
        return _get_LRUD(self, 'ref', name_LU='_xU', name_RU='_yU', name_LD='_xD', name_RD='_yD')

    @reference_center.setter
    def reference_center(self, ref):
        _set_LRUD(self, 'ref', ref, name_LU='_xU', name_RU='_yU', name_LD='_xD', name_RD='_yD')

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
