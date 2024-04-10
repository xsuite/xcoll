# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2023.                 #
# ######################################### #

import numpy as np

import xobjects as xo
import xtrack as xt

from ..impacts import CollimatorImpacts
from ..general import _pkg_root


class InvalidXcoll(xt.BeamElement):
    _xofields = {
        'length': xo.Float64
    }

    isthick = True
    behaves_like_drift = True
    allow_track = False
    skip_in_loss_location_refinement = True
    allow_backtrack = True

    # InvalidXcoll catches unallowed cases, like backtracking through a collimator
    _extra_c_sources = [
        _pkg_root.joinpath('headers','particle_states.h'),
        _pkg_root.joinpath('headers','checks.h')
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
    allow_track = False
    behaves_like_drift = True
    skip_in_loss_location_refinement = True

    _depends_on = [InvalidXcoll]

    _internal_record_class = CollimatorImpacts

    # This is an abstract class and cannot be instantiated
    def __new__(cls, *args, **kwargs):
        if cls == BaseBlock:
            raise Exception("Abstract class `BaseBlock` cannot be instantiated!")
        instance = super().__new__(cls)
        return instance

    def enable_scattering(self):
        if hasattr(self, '_tracking'):
            self._tracking = True

    def disable_scattering(self):
        if hasattr(self, '_tracking'):
            self._tracking = False

    def get_backtrack_element(self, _context=None, _buffer=None, _offset=None):
        return InvalidXcoll(length=-self.length,
                              _context=_context, _buffer=_buffer, _offset=_offset)


class BaseCollimator(xt.BeamElement):
    _xofields = {
        'length':  xo.Float64,
        # Collimator angle
        '_sin_zL':        xo.Float64,
        '_cos_zL':        xo.Float64,
        '_sin_zR':        xo.Float64,
        '_cos_zR':        xo.Float64,
        '_sin_zDiff':     xo.Float64, # Angle of right jaw: difference with respect to angle of left jaw
        '_cos_zDiff':     xo.Float64,
        '_jaws_parallel': xo.Int8,
        # Jaw corners (this is the x-coordinate in the rotated frame)
        '_jaw_LU': xo.Float64,  # left upstream
        '_jaw_RU': xo.Float64,  # right upstream
        '_jaw_LD': xo.Float64,  # left downstream
        '_jaw_RD': xo.Float64,  # right downstream
        # Tilts (superfluous but added to speed up calculations)
        '_sin_yL': xo.Float64,
        '_cos_yL': xo.Float64,
        '_tan_yL': xo.Float64,
        '_sin_yR': xo.Float64,
        '_cos_yR': xo.Float64,
        '_tan_yR': xo.Float64,
        # Other
        '_side':   xo.Int8,
        'active':  xo.Int8,
        'record_touches':      xo.Int8,
        'record_interactions': xo.Int8
    }

    isthick = True
    allow_track = False
    behaves_like_drift = True
    skip_in_loss_location_refinement = True

    # _skip_in_to_dict  = ['_jaw_LU', '_jaw_RU', '_jaw_LD', '_jaw_RD', '_sin_yL', '_cos_yL', '_tan_yL',
    #                      '_sin_yR', '_cos_yR', '_tan_yR', '_sin_zL', '_cos_zL', '_sin_zR', '_cos_zR', '_side']
    _skip_in_to_dict  = [f for f in _xofields if f.startswith('_')]
    _store_in_to_dict = ['angle', 'jaw', 'tilt', 'side']
    # Extra fields (only in Python): angle_L, angle_R, jaw_L, jaw_R, tilt_L, tilt_R

    _depends_on = [InvalidXcoll, xt.Drift, xt.XYShift, xt.SRotation, xt.YRotation]

    _extra_c_sources = [
        _pkg_root.joinpath('beam_elements','collimators_src','collimator_geometry.h')
    ]

    _internal_record_class = CollimatorImpacts

    # This is an abstract class and cannot be instantiated
    def __new__(cls, *args, **kwargs):
        if cls == BaseCollimator:
            raise Exception("Abstract class `BaseCollimator` cannot be instantiated!")
        instance = super().__new__(cls)
        return instance

    def __init__(self, **kwargs):
        to_assign = {}
        if '_xobject' not in kwargs:
            # Set angle
            if 'angle' in kwargs:
                for key in ['angle_L', 'angle_R']:
                    if key in kwargs:
                        raise ValueError(f"Cannot use both `angle` and `{key}`!")
                to_assign['angle'] = kwargs.pop('angle')
            else:
                to_assign['angle_L'] = kwargs.pop('angle_L', 0)
                to_assign['angle_R'] = kwargs.pop('angle_R', 0)

            # Set jaw
            if 'jaw' in kwargs:
                for key in ['jaw_L', 'jaw_R', 'jaw_LU', 'jaw_LD', 'jaw_RU', 'jaw_RD']:
                    if key in kwargs:
                        raise ValueError(f"Cannot use both `jaw and `{key}`!")
                to_assign['jaw'] = kwargs.pop('jaw')
            elif 'jaw_L' in kwargs or 'jaw_R' in kwargs:
                for key in ['jaw_LU', 'jaw_LD', 'jaw_RU', 'jaw_RD']:
                    if key in kwargs:
                        raise ValueError(f"Cannot use `jaw_L` or `jaw_R` together with `{key}`!")
                to_assign['jaw_L'] = kwargs.pop('jaw_L', 1)
                to_assign['jaw_R'] = kwargs.pop('jaw_R', -1)
            elif 'jaw_LU' in kwargs or 'jaw_LD' in kwargs or 'jaw_RU' in kwargs or 'jaw_RD' in kwargs:
                for key in ['tilt', 'tilt_L', 'tilt_R']:
                    if key in kwargs:
                        raise ValueError(f"Cannot use both `jaw_LU` etc with `{key}`!")
                to_assign['jaw_LU'] = kwargs.pop('jaw_LU', 1)
                to_assign['jaw_RU'] = kwargs.pop('jaw_RU', -1)
                to_assign['jaw_LD'] = kwargs.pop('jaw_LD', 1)
                to_assign['jaw_RD'] = kwargs.pop('jaw_RD', -1)
            else:
                to_assign['jaw'] = 1

            # Set tilt
            if 'tilt' in kwargs:
                for key in ['tilt_L', 'tilt_R']:
                    if key in kwargs:
                        raise ValueError(f"Cannot use both `tilt` and `{key}`!")
                to_assign['tilt'] = kwargs.pop('tilt')
            else:
                to_assign['tilt_L'] = kwargs.pop('tilt_L', 0)
                to_assign['tilt_R'] = kwargs.pop('tilt_R', 0)

            # Set others
            to_assign['side'] = kwargs.pop('side', 'both')
            kwargs.setdefault('active', True)
            kwargs.setdefault('record_touches', False)
            kwargs.setdefault('record_interactions', False)
            kwargs.setdefault('_jaw_LU', 1)
            kwargs.setdefault('_jaw_RU', -1)
            kwargs.setdefault('_jaw_LD', 1)
            kwargs.setdefault('_jaw_RD', -1)

        super().__init__(**kwargs)
        for key, val in to_assign.items():
            setattr(self, key, val)
        self._verify_consistency()


    # Main collimator angle
    # =====================

    @property
    def angle(self):
        return self.angle_L if self.angle_L==self.angle_R else [self.angle_L, self.angle_R]

    @angle.setter
    def angle(self, val):
        if not hasattr(val, '__iter__'):
            self.angle_L = val
            self.angle_R = val
        elif len(val) == 1:
            self.angle_L = val[0]
            self.angle_R = val[0]
        elif len(val) == 2:
            self.angle_L = val[0]
            self.angle_R = val[1]
        else:
            raise ValueError(f"The attribute `angle` should be of the form LR or [L, R] "
                           + f"but got {val}.")

    @property
    def angle_L(self):
        return round(np.rad2deg(np.arctan2(self._sin_zL, self._cos_zL)), 10)

    @angle_L.setter
    def angle_L(self, angle_L):
        self._sin_zL = np.sin(np.deg2rad(angle_L))
        self._cos_zL = np.cos(np.deg2rad(angle_L))
        if np.isclose(self.angle_R, angle_L):
            self._jaws_parallel = True
            self._sin_zDiff = 0.
            self._cos_zDiff = 1.
        else:
            self._jaws_parallel = False
            self._sin_zDiff = np.sin(np.deg2rad(self.angle_R - angle_L))
            self._cos_zDiff = np.cos(np.deg2rad(self.angle_R - angle_L))

    @property
    def angle_R(self):
        return round(np.rad2deg(np.arctan2(self._sin_zL, self._cos_zL)), 10)

    @angle_R.setter
    def angle_R(self, angle_R):
        self._sin_zR = np.sin(np.deg2rad(angle_R))
        self._cos_zR = np.cos(np.deg2rad(angle_R))
        if np.isclose(self.angle_L, angle_R):
            self._jaws_parallel = True
            self._sin_zDiff = 0.
            self._cos_zDiff = 1.
        else:
            self._jaws_parallel = False
            self._sin_zDiff = np.sin(np.deg2rad(angle_R - self.angle_L))
            self._cos_zDiff = np.cos(np.deg2rad(angle_R - self.angle_L))


    # Jaw attributes
    # ==============

    @property
    def jaw(self):
        if self.tilt_L == 0 and self.tilt_R == 0:
            return [self._jaw_LU, self._jaw_RU]
        else:
            return [[self._jaw_LU, self._jaw_RU], [self._jaw_LD, self._jaw_RD]]

    @jaw.setter   # Keeps the tilts unless all 4 corners are specified
    def jaw(self, val):
        if not hasattr(val, '__iter__'):
            self.jaw_L = val
            self.jaw_R = -val
            return
        elif len(val) == 1:
            self.jaw_L = val[0]
            self.jaw_R = -val[0]
            return
        elif len(val) == 2:
            if hasattr(val[0], '__iter__'):
                if hasattr(val[1], '__iter__') and len(val[0]) == 2 and len(val[1]) == 2:
                    self.jaw_LU = val[0][0]
                    self.jaw_RU = val[0][1]
                    self.jaw_LD = val[1][0]
                    self.jaw_RD = val[1][1]
                    return
            else:
                self.jaw_L = val[0]
                self.jaw_R = val[1]
                return
        # If we got here, val is incompatible
        raise ValueError(f"The attribute `jaw` should be of the form [L, R] or "
                       + f"[[LU, RU], [LD, RD], but got {val}.")

    @property
    def jaw_L(self):
        return (self._jaw_LU + self._jaw_LD) / 2

    @jaw_L.setter   # This moves both jaw_LU and jaw_LD in parallel
    def jaw_L(self, val):
        diff = val - self.jaw_L
        self._jaw_LU += diff
        self._jaw_LD += diff

    @property
    def jaw_R(self):
        return (self._jaw_RU + self._jaw_RD) / 2

    @jaw_R.setter   # This moves both jaw_RU and jaw_RD in parallel
    def jaw_R(self, val):
        diff = val - self.jaw_R
        self._jaw_RU += diff
        self._jaw_RD += diff

    @property
    def jaw_LU(self):
        return self._jaw_LU

    @jaw_LU.setter   # This assumes jaw_LD remains fixed, hence both jaw_L and the tilt change
    def jaw_LU(self, jaw_LU):
        self._jaw_LU = jaw_LU
        self._update_tilt()

    @property
    def jaw_LD(self):
        return self._jaw_LD

    @jaw_LD.setter   # This assumes jaw_LU remains fixed, hence both jaw_L and the tilt change
    def jaw_LD(self, jaw_LD):
        self._jaw_LD = jaw_LD
        self._update_tilt()

    @property
    def jaw_RU(self):
        return self._jaw_RU

    @jaw_RU.setter   # This assumes jaw_RD remains fixed, hence both jaw_R and the tilt change
    def jaw_RU(self, jaw_RU):
        self._jaw_RU = jaw_RU
        self._update_tilt()

    @property
    def jaw_RD(self):
        return self._jaw_RD

    @jaw_RD.setter   # This assumes jaw_RU remains fixed, hence both jaw_R and the tilt change
    def jaw_RD(self, jaw_RD):
        self._jaw_RD = jaw_RD
        self._update_tilt()


    # Tilt attributes
    # ===============

    # TODO: tilts are in rad! Do we want that? It's a bit inconsistent with angle which is in deg...

    @property
    def tilt(self):
        return self.tilt_L if self.tilt_L==self.tilt_R else [self.tilt_L, self.tilt_R]

    @tilt.setter
    def tilt(self, val):
        if not hasattr(val, '__iter__'):
            self.tilt_L = val
            self.tilt_R = val
        elif len(val) == 1:
            self.tilt_L = val[0]
            self.tilt_R = val[0]
        elif len(val) == 2:
            self.tilt_L = val[0]
            self.tilt_R = val[1]
        else:
            raise ValueError

    @property
    def tilt_L(self):
        return round(np.arctan2(self._sin_yL, self._cos_yL), 10)

    @tilt_L.setter   # This assumes jaw_L remains fixed (hence jaw_LU and jaw_LD change)
    def tilt_L(self, tilt_L):
        self._sin_yL = np.sin(tilt_L)
        self._cos_yL = np.cos(tilt_L)
        self._tan_yL = np.tan(tilt_L)
        jaw_L = self.jaw_L
        self._jaw_LD = jaw_L + self._sin_yL * self.length / 2.
        self._jaw_LU = jaw_L - self._sin_yL * self.length / 2.

    @property
    def tilt_R(self):
        return round(np.arctan2(self._sin_yR, self._cos_yR), 10)

    @tilt_R.setter   # This assumes jaw_R remains fixed (hence jaw_RU and jaw_RD change)
    def tilt_R(self, tilt_R):
        self._sin_yR = np.sin(tilt_R)
        self._cos_yR = np.cos(tilt_R)
        self._tan_yR = np.tan(tilt_R)
        jaw_R = self.jaw_R
        self._jaw_RD = jaw_R + self._sin_yR * self.length / 2.
        self._jaw_RU = jaw_R - self._sin_yR * self.length / 2.

    def _update_tilt(self):
        self._sin_yL = (self.jaw_LD - self.jaw_LU) / self.length
        self._cos_yL = np.sqrt(1 - self._sin_yL**2)
        self._tan_yL = self._sin_yL / self._cos_yL
        self._sin_yR = (self.jaw_RD - self.jaw_RU) / self.length
        self._cos_yR = np.sqrt(1 - self._sin_yR**2)
        self._tan_yR = self._sin_yR / self._cos_yR


    # Other attributes
    # ================

    @property
    def side(self):
        if self._side == 0:
            return 'both'
        elif self._side == 1:
            return 'left'
        elif self._side == -1:
            return 'right'

    @side.setter
    def side(self, val):
        if isinstance(val, str):
            if val.lower() == 'both' or val == '+-' or val == '-+':
                self._side = 0
                return
            elif val.lower() == 'left' or val.lower() == 'l' or val == '+':
                self._side = 1
                return
            elif val.lower() == 'right' or val.lower() == 'r' or val == '-':
                self._side = -1
                return
        raise ValueError(f"Unkown setting {val} for 'side'! Choose from "
                       + f"('left', 'L', '+'), ('right', 'R', '-'), or ('both', '+-').")

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


    # Methods
    # =======

    def enable_scattering(self):
        if hasattr(self, '_tracking'):
            self._tracking = True

    def disable_scattering(self):
        if hasattr(self, '_tracking'):
            self._tracking = False

    def _verify_consistency(self):
        # Verify angles
        if abs(self.angle_L - self.angle_R) >= 90.:
            raise ValueError("Angles of both jaws differ more than 90 degrees!")
        ang = abs(np.arccos(self._cos_zL))
        ang = np.pi - ang if ang > np.pi/2 else ang
        assert np.isclose(ang, abs(np.arcsin(self._sin_zL)))
        ang = abs(np.arccos(self._cos_zR))
        ang = np.pi - ang if ang > np.pi/2 else ang
        assert np.isclose(ang, abs(np.arcsin(self._sin_zR)))
        if np.isclose(self.angle_L, self.angle_R):
            assert self._jaws_parallel == True
            assert np.isclose(self._sin_zL, self._sin_zR)
            assert np.isclose(self._cos_zL, self._cos_zR)
            assert np.isclose(self._sin_zDiff, 0.)
            assert np.isclose(self._cos_zDiff, 1.)
        else:
            assert self._jaws_parallel == False
            assert np.isclose(self._sin_zDiff, self._cos_zL*self._sin_zR - self._sin_zL*self._cos_zR)
            assert np.isclose(self._cos_zDiff, self._cos_zL*self._cos_zR + self._sin_zL*self._sin_zR)
        if abs(self.tilt_L - self.tilt_R) >= 90.:
            raise ValueError("Tilts of both jaws differ more than 90 degrees!")
        ang = abs(np.arccos(self._cos_yL))
        ang = np.pi - ang if ang > np.pi/2 else ang
        assert np.isclose(ang, abs(np.arcsin(self._sin_yL)))
        assert np.isclose(self._sin_yL/self._cos_yL, self._tan_yL)
        ang = abs(np.arccos(self._cos_yR))
        ang = np.pi - ang if ang > np.pi/2 else ang
        assert np.isclose(ang, abs(np.arcsin(self._sin_yR)))
        assert np.isclose(self._sin_yR/self._cos_yR, self._tan_yR)

        # Verify bools
        assert self._side in [-1, 1, 0]
        assert isinstance(self._jaws_parallel, bool) or self._jaws_parallel in [0, 1]
        assert isinstance(self.active, bool) or self.active in [0, 1]
        assert isinstance(self.record_touches, bool) or self.record_touches in [0, 1]
        assert isinstance(self.record_interactions, bool) or self.record_interactions in [0, 1]

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

