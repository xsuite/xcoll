# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2023.                 #
# ######################################### #

import numpy as np

import xobjects as xo
import xtrack as xt

from ..collimator_settings import _get_LR, _set_LR, _get_LRUD, _set_LRUD
from ..impacts import CollimatorImpacts
from ..general import _pkg_root


class InvalidXcoll(xt.BeamElement):
    _xofields = {
        'length': xo.Float64
    }

    isthick = True
    behaves_like_drift = True
#     allow_track = False   # Need to wait for xtrack release to implement
    skip_in_loss_location_refinement = True
    allow_backtrack = True

    # InvalidXcoll catches unallowed cases, like backtracking through a collimator
    _extra_c_sources = [
        _pkg_root.joinpath('headers','particle_states.h'),
        _pkg_root.joinpath('headers','checks.h'),
        _pkg_root.joinpath('beam_elements','collimators_src','invalid_collimator.h')
    ]

    _depends_on = [xt.RandomRutherford]  # Needed for checks

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def get_backtrack_element(self, _context=None, _buffer=None, _offset=None):
        return self.__class__(length=-self.length,
                              _context=_context, _buffer=_buffer, _offset=_offset)


class BaseBlock(xt.BeamElement):
    _xofields = {
        'length': xo.Float64
    }

    isthick = True
    behaves_like_drift = True
#     allow_track = False   # Need to wait for xtrack release to implement
    skip_in_loss_location_refinement = True

    _extra_c_sources = [
        _pkg_root.joinpath('beam_elements','collimators_src','base_block.h')
    ]
    _depends_on = [InvalidXcoll]
    _internal_record_class = CollimatorImpacts

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def get_backtrack_element(self, _context=None, _buffer=None, _offset=None):
        return InvalidXcoll(length=-self.length,
                              _context=_context, _buffer=_buffer, _offset=_offset)


class BaseCollimator(xt.BeamElement):
    _xofields = {
        'length': xo.Float64,  # Length of jaws
        'jaw_L':  xo.Float64,  # left jaw (distance to ref)
        'jaw_R':  xo.Float64,  # right jaw
        'sin_zL': xo.Float64,  # angle of left jaw
        'cos_zL': xo.Float64,
        'sin_zR': xo.Float64,  # angle of right jaw
        'cos_zR': xo.Float64,
        'sin_yL': xo.Float64,  # tilt of left jaw (around jaw midpoint)
        'cos_yL': xo.Float64,
        'tan_yL': xo.Float64,
        'sin_yR': xo.Float64,  # tilt of right jaw (around jaw midpoint)
        'cos_yR': xo.Float64,
        'tan_yR': xo.Float64,
        '_side':  xo.Int8,     # is it a onesided collimator?
        'active': xo.Int8
    }

    isthick = True
    behaves_like_drift = True
#     allow_track = False   # Need to wait for xtrack release to implement
    skip_in_loss_location_refinement = True

    _skip_in_to_dict  = ['jaw_L', 'jaw_R', 'sin_yL', 'cos_yL', 'tan_yL', 'sin_yR',
                         'cos_yR', 'tan_yR', 'sin_zL', 'cos_zL', 'sin_zR', 'cos_zR', '_side']
    _store_in_to_dict = ['angle', 'tilt', 'jaw', 'side']
    # Extra fields (only in Python): angle_L, angle_R, tilt_L, tilt_R, jaw_LU, jaw_LD, jaw_RU, jaw_RD

    _extra_c_sources = [
        _pkg_root.joinpath('beam_elements','collimators_src','base_collimator.h')
    ]
    _depends_on = [InvalidXcoll, xt.Drift, xt.XYShift, xt.SRotation, xt.YRotation]
    _internal_record_class = CollimatorImpacts


    def __init__(self, **kwargs):
        # TODO: quick hack to avoid instantiation; did not manage to get it to work correclty with ABC
        if self.__class__.__name__ == 'BaseCollimator':
            raise Exception("Abstract class 'BaseCollimator' cannot be instantiated!")

        if '_xobject' not in kwargs:
            # Set jaw
            if 'jaw' in kwargs:
                if 'jaw_L' in kwargs or 'jaw_R' in kwargs:
                    raise ValueError("Cannot use both 'jaw' and 'jaw_L/R'!")
                _set_LR(kwargs, 'jaw', kwargs.pop('jaw'), neg=True)
            else:
                kwargs.setdefault('jaw_L', 1)
                kwargs.setdefault('jaw_R', -1)

            # Set side
            kwargs['_side'] = _side_setter(kwargs.pop('side', 'both'))

            # Set angle
            if 'angle' in kwargs:
                if 'angle_L' in kwargs or 'angle_R' in kwargs:
                    raise ValueError("Cannot use both 'angle' and 'angle_L/R'!")
                kwargs['sin_zL'], kwargs['cos_zL'], _, kwargs['sin_zR'], kwargs['cos_zR'], _ \
                            = _angle_setter(kwargs.pop('angle', 0))
            else:
                anglerad_L = kwargs.pop('angle_L', 0) / 180. * np.pi
                kwargs['sin_zL'] = np.sin(anglerad_L)
                kwargs['cos_zL'] = np.cos(anglerad_L)
                anglerad_R = kwargs.pop('angle_R', 0) / 180. * np.pi
                kwargs['sin_zR'] = np.sin(anglerad_R)
                kwargs['cos_zR'] = np.cos(anglerad_R)
                if abs(anglerad_L - anglerad_R) >= np.pi/2.:
                    raise ValueError(f"Angles of both jaws differ more than 90 degrees.")

            # Set tilt
            if 'tilt' in kwargs:
                if 'tilt_L' in kwargs or 'tilt_R' in kwargs:
                    raise ValueError("Cannot use both 'tilt' and 'tilt_L/R'!")
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
                if abs(anglerad_L - anglerad_R) >= np.pi/2.:
                    raise ValueError(f"Tilts of both jaws differ more than 90 degrees.")

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
        return self.jaw_L - self.sin_yL*self.length/2

    @jaw_LU.setter   # This assumes you keep jaw_LD fixed, hence both jaw_L and the tilt change
    def jaw_LU(self, jaw_LU):
        jaw_LD = self.jaw_LD
        self.jaw_L  = (jaw_LU+jaw_LD)/2
        self.sin_yL = (jaw_LD-jaw_LU)/self.length
        self.cos_yL = np.sqrt(1-self.sin_yL**2)
        self.tan_yL = self.sin_yL / self.cos_yL

    @property
    def jaw_LD(self):
        return self.jaw_L + self.sin_yL*self.length/2

    @jaw_LD.setter   # This assumes you keep jaw_LU fixed, hence both jaw_L and the tilt change
    def jaw_LD(self, jaw_LD):
        jaw_LU = self.jaw_LU
        self.jaw_L  = (jaw_LU+jaw_LD)/2
        self.sin_yL = (jaw_LD-jaw_LU)/self.length
        self.cos_yL = np.sqrt(1-self.sin_yL**2)
        self.tan_yL = self.sin_yL / self.cos_yL

    @property
    def jaw_RU(self):
        return self.jaw_R - self.sin_yR*self.length/2

    @jaw_RU.setter   # This assumes you keep jaw_RD fixed, hence both jaw_R and the tilt change
    def jaw_RU(self, jaw_RU):
        jaw_RD = self.jaw_RD
        self.jaw_R  = (jaw_RU+jaw_RD)/2
        self.sin_yR = (jaw_RD-jaw_RU)/self.length
        self.cos_yR = np.sqrt(1-self.sin_yR**2)
        self.tan_yR = self.sin_yR / self.cos_yR

    @property
    def jaw_RD(self):
        return self.jaw_R + self.sin_yR*self.length/2

    @jaw_RD.setter   # This assumes you keep jaw_RU fixed, hence both jaw_R and the tilt change
    def jaw_RD(self, jaw_RD):
        jaw_RU = self.jaw_RU
        self.jaw_R  = (jaw_RU+jaw_RD)/2
        self.sin_yR = (jaw_RD-jaw_RU)/self.length
        self.cos_yR = np.sqrt(1-self.sin_yR**2)
        self.tan_yR = self.sin_yR / self.cos_yR

    @property
    def angle_L(self):
        return round(np.arctan2(self.sin_zL, self.cos_zL) * (180.0 / np.pi), 10)

    @angle_L.setter
    def angle_L(self, angle_L):
        if abs(angle_L - self.angle_R) >= 90.:
            raise ValueError(f"Angles of both jaws differ more than 90 degrees.")
        anglerad_L = angle_L / 180. * np.pi
        self.sin_zL = np.sin(anglerad_L)
        self.cos_zL = np.cos(anglerad_L)

    @property
    def angle_R(self):
        return round(np.arctan2(self.sin_zR, self.cos_zR) * (180.0 / np.pi), 10)

    @angle_R.setter
    def angle_R(self, angle_R):
        if abs(self.angle_L - angle_R) >= 90.:
            raise ValueError(f"Angles of both jaws differ more than 90 degrees.")
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
        if abs(tilt_L - self.tilt_R) >= np.pi/2.:
            raise ValueError(f"Tilts of both jaws differ more than 90 degrees.")
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
        if abs(self.tilt_L - tilt_R) >= np.pi/2.:
            raise ValueError(f"Tilts of both jaws differ more than 90 degrees.")
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
        if pos[0] == 'L': 
            other_pos = 'R' 
        else: 
            other_pos = 'L'
        point_x = ((getattr(self, 'jaw_' + pos[0]) * getattr(self, 'cos_z' + pos[0]) 
                    + getattr(self, 'jaw_' + other_pos) * getattr(self, 'cos_z' + other_pos))/2)
        point_y = ((getattr(self, 'jaw_' + pos[0]) * getattr(self, 'sin_z' + pos[0]) 
                    + getattr(self, 'jaw_' + other_pos) * getattr(self, 'sin_z' + other_pos))/2)
        if not pos in positions:
            raise ValueError(f"Parameter {pos} needs to be one of {positions}!")
        sinz    = getattr(self, 'sin_z' + pos[0])
        cosz    = getattr(self, 'cos_z' + pos[0])
        # Shift to the jaw, whose location is given as the shortest distance:
        point_x += getattr(self, 'jaw_' + pos) * cosz
        point_y += getattr(self, 'jaw_' + pos) * sinz
        return lambda t: (point_x - t*sinz, point_y + t*cosz)

    @property
    def active_length(self):
        raise ValueError("`active_length`is deprecated. Please use `length`.")

    @property
    def inactive_front(self):
        raise ValueError("`inactive_front`is deprecated. Collimators now only "
                       + "contain their active length (implemented as `length`).")

    @property
    def inactive_back(self):
        raise ValueError("`inactive_back`is deprecated. Collimators now only "
                       + "contain their active length (implemented as `length`).")


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
    if abs(anglerad_L - anglerad_R) >= np.pi/2.:
        raise ValueError(f"Angles of both jaws differ more than 90 degrees.")
    return np.sin(anglerad_L), np.cos(anglerad_L), np.tan(anglerad_L), \
           np.sin(anglerad_R), np.cos(anglerad_R), np.tan(anglerad_R)


