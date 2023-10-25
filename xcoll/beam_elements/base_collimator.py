# copyright ############################### #
# This file is part of the Xcoll Package.  #
# Copyright (c) CERN, 2023.                 #
# ######################################### #

import numpy as np

import xobjects as xo
import xtrack as xt

from ..collimator_settings import _get_LR, _set_LR, _get_LRUD, _set_LRUD
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
        'jaw_L':          xo.Float64,  # left jaw (distance to ref)
        'jaw_R':          xo.Float64,  # right jaw
        'ref_x':          xo.Float64,  # center of collimator reference frame
        'ref_y':          xo.Float64,
        'sin_zL':         xo.Float64,  # angle of left jaw
        'cos_zL':         xo.Float64,
        'sin_zR':         xo.Float64,  # angle of right jaw
        'cos_zR':         xo.Float64,
        'sin_yL':         xo.Float64,  # tilt of left jaw (around jaw midpoint)
        'cos_yL':         xo.Float64,
        'tan_yL':         xo.Float64,
        'sin_yR':         xo.Float64,  # tilt of right jaw (around jaw midpoint)
        'cos_yR':         xo.Float64,
        'tan_yR':         xo.Float64,
        '_side':          xo.Int8,     # is it a onesided collimator?
        'active':         xo.Int8
    }

    isthick = True
    behaves_like_drift = True
    skip_in_loss_location_refinement = True
    
    _skip_in_to_dict  = ['jaw_L', 'jaw_R', 'ref_x', 'ref_y',
                         'sin_yL', 'cos_yL', 'tan_yL', 'sin_yR', 'cos_yR', 'tan_yR',
                         'sin_zL', 'cos_zL', 'sin_zR', 'cos_zR', '_side']
    _store_in_to_dict = ['angle', 'tilt', 'jaw', 'reference_center', 'side']
    # Extra fields (only in Python): angle_L, angle_R, tilt_L, tilt_R, jaw_LU, jaw_LD, jaw_RU, jaw_RD
    _internal_record_class = CollimatorImpacts

    _extra_c_sources = [
        _pkg_root.joinpath('beam_elements','collimators_src','base_collimator.h')
    ]

    _depends_on = [InvalidCollimator, xt.Drift, xt.XYShift, xt.SRotation, xt.YRotation]


    def __init__(self, **kwargs):
        # TODO: quick hack to avoid instantiation; did not manage to get it to work correclty with ABC
        if self.__class__.__name__ == 'BaseCollimator':
            raise Exception("Abstract class 'BaseCollimator' cannot be instantiated!")

        if '_xobject' not in kwargs:
            # Set jaw
            if 'jaw' in kwargs:
                if 'jaw_L' in kwargs or 'jaw_R' in kwargs:
                    raise ValuError("Cannot use both 'jaw' and 'jaw_L/R'!")
                _set_LR(kwargs, 'jaw', kwargs.pop('jaw'), neg=True)
            else:
                kwargs.setdefault('jaw_L', 1)
                kwargs.setdefault('jaw_R', -1)

            # Set reference_center
            if 'reference_center' in kwargs:
                if 'ref_x' in kwargs or 'ref_y' in kwargs:
                    raise ValuError("Cannot use both 'reference_center' and 'ref_x/y'!")
                _set_LR(kwargs, 'ref', kwargs.pop('reference_center'), name_L='_x', name_R='_y')
            else:
                kwargs.setdefault('ref_x', 0)
                kwargs.setdefault('ref_y', 0)

            # Set side
            kwargs['_side'] = _side_setter(kwargs.pop('side', 'both'))

            # Set angle
            if 'angle' in kwargs:
                if 'angle_L' in kwargs or 'angle_R' in kwargs:
                    raise ValuError("Cannot use both 'angle' and 'angle_L/R'!")
                kwargs['sin_zL'], kwargs['cos_zL'], _, kwargs['sin_zR'], kwargs['cos_zR'], _ \
                            = _angle_setter(kwargs.pop('angle', 0))
            else:
                anglerad_L = kwargs.pop('angle_L', 0) / 180. * np.pi
                kwargs['sin_zL'] = np.sin(anglerad_L)
                kwargs['cos_zL'] = np.cos(anglerad_L)
                anglerad_R = kwargs.pop('angle_R', 0) / 180. * np.pi
                kwargs['sin_zR'] = np.sin(anglerad_R)
                kwargs['cos_zR'] = np.cos(anglerad_R)

            # Set tilt
            if 'tilt' in kwargs:
                if 'tilt_L' in kwargs or 'tilt_R' in kwargs:
                    raise ValuError("Cannot use both 'tilt' and 'tilt_L/R'!")
                kwargs['sin_yL'], kwargs['cos_yL'], kwargs['tan_yL'], \
                kwargs['sin_yR'], kwargs['cos_yR'], kwargs['tan_yR'] \
                                = _angle_setter(kwargs.pop('tilt', 0), rad=True)
            else:
#                 anglerad_L = kwargs.pop('tilt_L', 0) / 180. * np.pi
                anglerad_L = kwargs.pop('tilt_L', 0)
                kwargs['sin_yL'] = np.sin(anglerad_L)
                kwargs['cos_yL'] = np.cos(anglerad_L)
                kwargs['tan_yL'] = np.cos(anglerad_L)
#                 anglerad_R = kwargs.pop('tilt_R', 0) / 180. * np.pi
                anglerad_R = kwargs.pop('tilt_R', 0)
                kwargs['sin_yR'] = np.sin(anglerad_R)
                kwargs['cos_yR'] = np.cos(anglerad_R)
                kwargs['tan_yR'] = np.cos(anglerad_R)

            # Set lengths
            if 'length' in kwargs.keys():
                if 'active_length' in kwargs.keys():
                    if 'inactive_front' in kwargs.keys():
                        if 'inactive_back' in kwargs.keys():
                            raise ValueError("Too many length variables used at initialisation!")
                        kwargs['inactive_back'] = kwargs.pop('length') - kwargs['inactive_front'] - kwargs['active_length']
                    elif 'inactive_back' in kwargs.keys():
                        if 'inactive_front' in kwargs.keys():
                            raise ValueError("Too many length variables used at initialisation!")
                        kwargs['inactive_front'] = kwargs.pop('length') - kwargs['inactive_back'] - kwargs['active_length']
                    else:
                        kwargs['inactive_front'] = (kwargs['length']- kwargs['active_length'])/2
                        kwargs['inactive_back'] = (kwargs.pop('length') - kwargs['active_length'])/2
                elif 'inactive_front' in kwargs.keys():
                    if 'inactive_back' in kwargs.keys():
                        kwargs['active_length'] = kwargs.pop('length') - kwargs['inactive_back'] - kwargs['inactive_front']
                    else:
                        kwargs['active_length'] = kwargs.pop('length') - kwargs['inactive_front']
                elif 'inactive_back' in kwargs.keys():
                    kwargs['active_length'] = kwargs.pop('length') - kwargs['inactive_back']
                else:
                    kwargs['active_length'] = kwargs.pop('length')
            kwargs.setdefault('inactive_front', 0)
            kwargs.setdefault('inactive_back', 0)
            kwargs.setdefault('active_length', 0)
            kwargs.setdefault('active', True)

        super().__init__(**kwargs)


    @property
    def jaw(self):
        return _get_LR(self, 'jaw', neg=True)

    @jaw.setter
    def jaw(self, val):
        _set_LR(self, 'jaw', val, neg=True)

    @property
    def jaw_LU(self):
        return self.jaw_L - self.sin_yL*self.active_length/2

    @jaw_LU.setter   # This assumes you keep jaw_LD fixed, hence both jaw_L and the tilt change
    def jaw_LU(self, jaw_LU):
        jaw_LD = self.jaw_LD
        self.jaw_L  = (jaw_LU+jaw_LD)/2
        self.sin_yL = (jaw_LD-jaw_LU)/self.active_length
        self.cos_yL = np.sqrt(1-self.sin_yL**2)
        self.tan_yL = self.sin_yL / self.cos_yL

    @property
    def jaw_LD(self):
        return self.jaw_L + self.sin_yL*self.active_length/2

    @jaw_LD.setter   # This assumes you keep jaw_LU fixed, hence both jaw_L and the tilt change
    def jaw_LD(self, jaw_LD):
        jaw_LU = self.jaw_LU
        self.jaw_L  = (jaw_LU+jaw_LD)/2
        self.sin_yL = (jaw_LD-jaw_LU)/self.active_length
        self.cos_yL = np.sqrt(1-self.sin_yL**2)
        self.tan_yL = self.sin_yL / self.cos_yL

    @property
    def jaw_RU(self):
        return self.jaw_R - self.sin_yR*self.active_length/2

    @jaw_RU.setter   # This assumes you keep jaw_RD fixed, hence both jaw_R and the tilt change
    def jaw_RU(self, jaw_RU):
        jaw_RD = self.jaw_RD
        self.jaw_R  = (jaw_RU+jaw_RD)/2
        self.sin_yR = (jaw_RD-jaw_RU)/self.active_length
        self.cos_yR = np.sqrt(1-self.sin_yR**2)
        self.tan_yR = self.sin_yR / self.cos_yR

    @property
    def jaw_RD(self):
        return self.jaw_R + self.sin_yR*self.active_length/2

    @jaw_RD.setter   # This assumes you keep jaw_RU fixed, hence both jaw_R and the tilt change
    def jaw_RD(self, jaw_RD):
        jaw_RU = self.jaw_RU
        self.jaw_R  = (jaw_RU+jaw_RD)/2
        self.sin_yR = (jaw_RD-jaw_RU)/self.active_length
        self.cos_yR = np.sqrt(1-self.sin_yR**2)
        self.tan_yR = self.sin_yR / self.cos_yR

    @property
    def angle_L(self):
        return round(np.arctan2(self.sin_zL, self.cos_zL) * (180.0 / np.pi), 10)

    @angle_L.setter
    def angle_L(self, angle_L):
        anglerad_L = angle_L / 180. * np.pi
        self.sin_zL = np.sin(anglerad_L)
        self.cos_zL = np.cos(anglerad_L)

    @property
    def angle_R(self):
        return round(np.arctan2(self.sin_zR, self.cos_zR) * (180.0 / np.pi), 10)

    @angle_R.setter
    def angle_R(self, angle_R):
        anglerad_R = angle_R / 180. * np.pi
        self.sin_zR = np.sin(anglerad_R)
        self.cos_zR = np.cos(anglerad_R)

    @property
    def angle(self):
        return self.angle_L if self.angle_L==self.angle_R else [self.angle_L, self.angle_R]

    @angle.setter
    def angle(self, angle):
        self.sin_zL, self.cos_zL, _, self.sin_zR, self.cos_zR, _ = _angle_setter(angle)

    @property
    def reference_center(self):
        return [self.ref_x, self.ref_y]

    @reference_center.setter
    def reference_center(self, ref):
        _set_LR(self, 'ref', ref, name_L='_x', name_R='_y')

    @property
    def length(self):
        return (self.inactive_front + self.active_length + self.inactive_back)

    @length.setter
    def length(self, val):
        raise ValueError("The parameter 'length' can only be set at initialisation. " 
                       + "Use 'active_length', 'inactive_front', and/or 'inactive_back'.")

    @property
    def side(self):
        if self._side == 0:
            return 'both'
        elif self._side == 1:
            return 'left'
        elif self._side == 2:
            return 'right'

    @side.setter
    def side(self, side):
        self._side = _side_setter(side)

    # TODO: tilts are in rad! Do we want that? It's a bit inconsistent with angle which is in deg...
    # ==============================================================================================
    @property
    def tilt_L(self):
#         return round(np.arctan2(self.sin_yL, self.cos_yL) * (180.0 / np.pi), 10)
        return round(np.arctan2(self.sin_yL, self.cos_yL), 10)

    @tilt_L.setter
    def tilt_L(self, tilt_L):
#         anglerad_L = tilt_L / 180. * np.pi
        anglerad_L = tilt_L
        self.sin_yL = np.sin(anglerad_L)
        self.cos_yL = np.cos(anglerad_L)
        self.tan_yL = np.tan(anglerad_L)

    @property
    def tilt_R(self):
#         return round(np.arctan2(self.sin_yR, self.cos_yR) * (180.0 / np.pi), 10)
        return round(np.arctan2(self.sin_yR, self.cos_yR), 10)

    @tilt_R.setter
    def tilt_R(self, tilt_R):
#         anglerad_R = tilt_R / 180. * np.pi
        anglerad_R = tilt_R
        self.sin_yR = np.sin(anglerad_R)
        self.cos_yR = np.cos(anglerad_R)
        self.tan_yR = np.tan(anglerad_R)

    @property
    def tilt(self):
        return self.tilt_L if self.tilt_L==self.tilt_R else [self.tilt_L, self.tilt_R]

    @tilt.setter
    def tilt(self, tilt):
        self.sin_yL, self.cos_yL, self.tan_yL, self.sin_yR, self.cos_yR, self.tan_yR = _angle_setter(tilt, rad=True)

    def jaw_func(self, pos):
        positions = ['LU', 'RU', 'LD', 'RD']
        if not pos in positions:
            raise ValueError(f"Parameter {pos} needs to be one of {positions}!")
        point_x = getattr(self, 'ref_x')
        point_y = getattr(self, 'ref_y')
        sinz    = getattr(self, 'sin_z' + pos[0])
        cosz    = getattr(self, 'cos_z' + pos[0])
        # Shift to the jaw, whose location is given as the shortest distance:
        point_x += getattr(self, 'jaw_' + pos) * cosz
        point_y += getattr(self, 'jaw_' + pos) * sinz
        return lambda t: (point_x - t*sinz, point_y + t*cosz)


def _side_setter(val):
    if isinstance(val, str):
        if val.lower() == 'both' or val == '+-' or val == '-+':
            return 0
        elif val.lower() == 'left' or val.lower() == 'l' or val == '+':
            return 1
        elif val.lower() == 'right' or val.lower() == 'r' or val == '-':
            return 2
    raise ValueError(f"Unkown setting {val} for 'side'! Choose from "
                   + f"('left', 'L', '+'), ('right', 'R', '-'), or ('both', '+-').")

def _angle_setter(val, rad=False):
    if not hasattr(val, '__iter__'):
        val = [val]
    conversion = 1 if rad else np.pi / 180.
    if isinstance(val, str):
        raise ValueError(f"Error in setting angle: not a number!")
    elif len(val) == 2:
        anglerad_L = val[0] * conversion
        anglerad_R = val[1] * conversion
    elif len(val) == 1:
        anglerad_L = val[0] * conversion
        anglerad_R = val[0] * conversion
    else:
        raise ValueError(f"Error in setting angle: must have one or two (L, R) values!")
    return np.sin(anglerad_L), np.cos(anglerad_L), np.tan(anglerad_L), \
           np.sin(anglerad_R), np.cos(anglerad_R), np.tan(anglerad_R)


