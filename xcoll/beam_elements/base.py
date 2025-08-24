# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

import numpy as np

import xobjects as xo
import xtrack as xt

from ..interaction_record import InteractionRecord
from ..general import _pkg_root
from ..headers.particle_states import particle_states_src


OPEN_JAW = 3
OPEN_GAP = 999


class BaseXcoll(xt.BeamElement):
    _xofields = {
        'length': xo.Float64
    }

    isthick = True
    needs_rng = False
    allow_track = False
    behaves_like_drift = True
    allow_rot_and_shift = False
    allow_loss_refinement = True
    skip_in_loss_location_refinement = True

    _extra_c_sources = [
        particle_states_src,
        _pkg_root.joinpath('headers','checks.h')
    ]

    _depends_on = [xt.RandomRutherford]  # Needed for checks

    def __init__(self, **kwargs):
        super().__init__(**kwargs)


class BaseBlock(xt.BeamElement):
    _xofields = {
        'length':                xo.Float64,
        'active':                xo.Int8,
        '_record_interactions':  xo.Int8
    }

    isthick = True
    needs_rng = False
    allow_track = False
    allow_double_sided = True
    behaves_like_drift = True
    allow_rot_and_shift = False
    allow_loss_refinement = True
    skip_in_loss_location_refinement = True

    _noexpr_fields = {'name'}
    _skip_in_to_dict  = ['_record_interactions']
    _store_in_to_dict = ['name', 'record_impacts', 'record_exits', 'record_scatterings']

    _depends_on = [BaseXcoll]

    _internal_record_class = InteractionRecord

    # This is an abstract class and cannot be instantiated
    def __new__(cls, *args, **kwargs):
        if cls is BaseBlock:
            raise Exception("Abstract class `BaseBlock` cannot be instantiated!")
        instance = super().__new__(cls)
        return instance

    def __init__(self, **kwargs):
        to_assign = {}
        if '_xobject' not in kwargs:
            # Set name (useful for bookkeeping like in FLUKA)
            to_assign['name'] = kwargs.pop('name', None)
            # Set active
            kwargs.setdefault('active', True)
            to_assign['record_impacts'] = kwargs.pop('record_impacts', False)
            to_assign['record_exits'] = kwargs.pop('record_exits', False)
            to_assign['record_scatterings'] = kwargs.pop('record_scatterings', False)
        super().__init__(**kwargs)
        # Careful: non-xofields are not passed correctly between copy's / to_dict. This messes with flags etc..
        # We also have to manually initialise them for xobject generation
        for key, val in to_assign.items():
            setattr(self, key, val)
        BaseBlock._verify_consistency(self)

    def copy(self, **kwargs):
        obj = super().copy(**kwargs)
        obj.name = self.name
        return obj

    @property
    def name(self):
        if not hasattr(self, '_name'):
            self._name = None
        return self._name

    @name.setter
    def name(self, val):
        self._name = val

    def enable_scattering(self):
        if hasattr(self, '_tracking'):
            if hasattr(self, 'optics') and self.optics is None and \
            (hasattr(self, '_gap_L_set_manually') and self._gap_L_set_manually() \
            or hasattr(self, '_gap_R_set_manually') and self._gap_R_set_manually()):
                raise ValueError("Gap set but optics not yet assigned! "
                               + "Cannot enable scattering.")
            self._tracking = True

    def disable_scattering(self):
        if hasattr(self, '_tracking'):
            self._tracking = False

    @property
    def record_impacts(self):
        return bool(self._record_interactions % 2)

    @record_impacts.setter
    def record_impacts(self, val):
        if not isinstance(val, bool):
            raise ValueError("`record_impacts` must be a boolean value.")
        if val and not self.record_impacts:
            self._record_interactions += 1
        elif not val and self.record_impacts:
            self._record_interactions -= 1

    @property
    def record_exits(self):
        return bool((self._record_interactions >> 1) % 2)

    @record_exits.setter
    def record_exits(self, val):
        if not isinstance(val, bool):
            raise ValueError("`record_exits` must be a boolean value.")
        if val and not self.record_exits:
            self._record_interactions += 2
        elif not val and self.record_exits:
            self._record_interactions -= 2

    @property
    def record_scatterings(self):
        return bool((self._record_interactions >> 2) % 2)

    @record_scatterings.setter
    def record_scatterings(self, val):
        if not isinstance(val, bool):
            raise ValueError("`record_scatterings` must be a boolean value.")
        if val and not self.record_scatterings:
            self._record_interactions += 4
        elif not val and self.record_scatterings:
            self._record_interactions -= 4

    def _verify_consistency(self):
        assert isinstance(self.active, bool) or self.active in [0, 1]
        assert self._record_interactions in list(range(8))


class BaseCollimator(BaseBlock):
    _xofields = BaseBlock._xofields | {
        # Collimator angle
        '_sin_zL':        xo.Float64,
        '_cos_zL':        xo.Float64,
        '_sin_zR':        xo.Float64,
        '_cos_zR':        xo.Float64,
        '_sin_zDiff':     xo.Float64, # Angle of right jaw: difference with respect to angle of left jaw
        '_cos_zDiff':     xo.Float64,
        '_jaws_parallel': xo.Int8,
        # Jaw corners (this is the x-coordinate in the rotated frame)
        '_jaw_LU':        xo.Float64,  # left upstream
        '_jaw_RU':        xo.Float64,  # right upstream
        '_jaw_LD':        xo.Float64,  # left downstream
        '_jaw_RD':        xo.Float64,  # right downstream
        # Tilts (superfluous but added to speed up calculations)
        '_sin_yL':        xo.Float64,
        '_cos_yL':        xo.Float64,
        '_tan_yL':        xo.Float64,
        '_sin_yR':        xo.Float64,
        '_cos_yR':        xo.Float64,
        '_tan_yR':        xo.Float64,
        # Other
        '_side':          xo.Int8,
        # These are not used in C, but need to be an xofield to get them in the to_dict:
        '_align':         xo.Int8,
        '_gap_L':         xo.Float64,
        '_gap_R':         xo.Float64,
        '_nemitt_x':      xo.Float64,
        '_nemitt_y':      xo.Float64
    }

    isthick = True
    needs_rng = False
    allow_track = False
    allow_double_sided = True
    behaves_like_drift = True
    allow_rot_and_shift = False
    allow_loss_refinement = True
    skip_in_loss_location_refinement = True

    _noexpr_fields = {'align', 'side', 'name'}
    _skip_in_to_dict  = [*BaseBlock._skip_in_to_dict,
                         *[f for f in _xofields if f.startswith('_')]]
    _store_in_to_dict = [*BaseBlock._store_in_to_dict, 'angle', 'jaw', 'tilt', 'gap',
                         'side', 'align', 'emittance']

    _depends_on = [BaseBlock]

    _internal_record_class = BaseBlock._internal_record_class


    # This is an abstract class and cannot be instantiated
    def __new__(cls, *args, **kwargs):
        if cls is BaseCollimator:
            raise Exception("Abstract class `BaseCollimator` cannot be instantiated!")
        instance = super().__new__(cls)
        return instance

    def __init__(self, **kwargs):
        to_assign = {}
        if '_xobject' not in kwargs:
            # Set side
            to_assign['side'] = kwargs.pop('side', 'both')

            # Set angle
            if 'angle' in kwargs:
                for key in ['angle_L', 'angle_R']:
                    if key in kwargs:
                        raise ValueError(f"Cannot use both `angle` and `{key}`!")
                to_assign['angle'] = kwargs.pop('angle')
            else:
                to_assign['angle_L'] = kwargs.pop('angle_L', 0)
                to_assign['angle_R'] = kwargs.pop('angle_R', 0)

            # We do not allow any combination of jaw_ and gap_ attributes
            # (except when jaw=..., gap=None or jaw=None, gap=... is used, as this is how the colldb installs it)
            kwargs = {kk: vv for kk, vv in kwargs.items() if not vv is None}

            # Set jaw
            if 'jaw' in kwargs:
                for key in ['jaw_L', 'jaw_R', 'jaw_LU', 'jaw_LD', 'jaw_RU', 'jaw_RD', 'gap', 'gap_L', 'gap_R']:
                    if key in kwargs:
                        raise ValueError(f"Cannot use both `jaw` and `{key}`!")
                if hasattr(kwargs['jaw'], '__iter__') and hasattr(kwargs['jaw'][0], '__iter__'):
                    for key in ['tilt', 'tilt_L', 'tilt_R']:
                        if key in kwargs:
                            raise ValueError(f"Cannot specify jaw corners and `{key}`!")
                to_assign['jaw'] = kwargs.pop('jaw')
            elif 'jaw_L' in kwargs or 'jaw_R' in kwargs:
                for key in ['jaw_LU', 'jaw_LD', 'jaw_RU', 'jaw_RD', 'gap', 'gap_L', 'gap_R']:
                    if key in kwargs:
                        raise ValueError(f"Cannot use `jaw_L` or `jaw_R` together with `{key}`!")
                to_assign['jaw_L'] = kwargs.pop('jaw_L', None)
                to_assign['jaw_R'] = kwargs.pop('jaw_R', None)
            elif 'jaw_LU' in kwargs or 'jaw_LD' in kwargs or 'jaw_RU' in kwargs or 'jaw_RD' in kwargs:
                for key in ['tilt', 'tilt_L', 'tilt_R', 'gap', 'gap_L', 'gap_R']:
                    if key in kwargs:
                        raise ValueError(f"Cannot use both `jaw_LU` etc with `{key}`!")
                to_assign['jaw_LU'] = kwargs.pop('jaw_LU', None)
                to_assign['jaw_RU'] = kwargs.pop('jaw_RU', None)
                to_assign['jaw_LD'] = kwargs.pop('jaw_LD', None)
                to_assign['jaw_RD'] = kwargs.pop('jaw_RD', None)
            kwargs.setdefault('_jaw_LU', OPEN_JAW)   # Important that these are initialised, in 
            kwargs.setdefault('_jaw_RU', -OPEN_JAW)  # order to keep a tilt (given together with
            kwargs.setdefault('_jaw_LD', OPEN_JAW)   # a gap) when optics are not known yet.
            kwargs.setdefault('_jaw_RD', -OPEN_JAW)

            # Set tilt
            if 'tilt' in kwargs:
                for key in ['tilt_L', 'tilt_R']:
                    if key in kwargs:
                        raise ValueError(f"Cannot use both `tilt` and `{key}`!")
                to_assign['tilt'] = kwargs.pop('tilt')
            elif 'tilt_L' in kwargs or 'tilt_R' in kwargs:
                to_assign['tilt_L'] = kwargs.pop('tilt_L', None)
                to_assign['tilt_R'] = kwargs.pop('tilt_R', None)
            kwargs.setdefault('_sin_yL', 0)
            kwargs.setdefault('_cos_yL', 1)
            kwargs.setdefault('_sin_yR', 0)
            kwargs.setdefault('_cos_yR', 1)

            # Set gap
            if 'gap' in kwargs:
                for key in ['gap_L', 'gap_R']:
                    if key in kwargs:
                        raise ValueError(f"Cannot use both `gap` and `{key}`!")
                to_assign['gap'] = kwargs.pop('gap')
            elif 'gap_L' in kwargs or 'gap_R' in kwargs:
                to_assign['gap_L'] = kwargs.pop('gap_L', None)
                to_assign['gap_R'] = kwargs.pop('gap_R', None)
            kwargs.setdefault('_gap_L', OPEN_GAP)
            kwargs.setdefault('_gap_R', -OPEN_GAP)

            # Set others
            to_assign['align'] = kwargs.pop('align', 'upstream')
            to_assign['emittance'] = kwargs.pop('emittance', None)

        super().__init__(**kwargs)
        # Careful: non-xofields are not passed correctly between copy's / to_dict. This messes with flags etc..
        # We also have to manually initialise them for xobject generation
        if not hasattr(self, '_optics'):
            self._optics = None
        for key, val in to_assign.items():
            setattr(self, key, val)
        self._update_tilts()
        BaseCollimator._verify_consistency(self)


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
        self._apply_optics()

    @property
    def angle_R(self):
        return round(np.rad2deg(np.arctan2(self._sin_zR, self._cos_zR)), 10)

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
        self._apply_optics()


    # Jaw attributes
    # ==============

    @property
    def jaw(self):
        if self.jaw_L is None and self.jaw_R is None:
            return None
        elif (self.side == 'left'  and self.tilt_L == 0) \
        or   (self.side == 'right' and self.tilt_R == 0) \
        or   (self.tilt_L == 0 and self.tilt_R == 0):
            return [self.jaw_L, self.jaw_R]
        else:
            return [[self.jaw_LU, self.jaw_LD], [self.jaw_RU, self.jaw_RD]]

    @jaw.setter   # Keeps the tilts unless all 4 corners are specified
    def jaw(self, val):
        if not hasattr(val, '__iter__') or len(val) == 1:
            val = val[0] if hasattr(val, '__iter__') else val
            if self.side == 'left':
                self.jaw_L = val
            elif self.side == 'right':
                # self.jaw_R = -val if val is not None else None
                self.jaw_R = val
            else:
                self.jaw_L = val
                self.jaw_R = -val if val is not None else None
            return
        elif len(val) == 2:
            if hasattr(val[0], '__iter__'):
                if hasattr(val[1], '__iter__') and len(val[0]) == 2 and len(val[1]) == 2:
                    self.jaw_LU = val[0][0]
                    self.jaw_LD = val[0][1]
                    self.jaw_RU = val[1][0]
                    self.jaw_RD = val[1][1]
                    return
            else:
                self.jaw_L = val[0]
                self.jaw_R = val[1]
                return
        # If we got here, val is incompatible
        raise ValueError(f"The attribute `jaw` should be of the form [L, R] or "
                       + f"[[LU, LD], [RU, RD], but got {val}.")

    @property
    def jaw_L(self):
        jaw_L = (self._jaw_LU + self._jaw_LD) / 2
        if not np.isclose(jaw_L, OPEN_JAW, atol=1.e-10):  # open position
            return jaw_L

    @jaw_L.setter   # This moves both jaw_LU and jaw_LD in parallel
    def jaw_L(self, val):
        if self.side == 'right' and val is not None:
            val = None
            print("Warning: Ignored value for jaw_L (right-sided collimator).")
        if val is None:
            val = OPEN_JAW
            self._gap_L = OPEN_GAP
        diff = val - (self._jaw_LU + self._jaw_LD) / 2
        self._jaw_LU += diff
        self._jaw_LD += diff
        self._update_gaps(only_L=True)

    @property
    def jaw_R(self):
        jaw_R = (self._jaw_RU + self._jaw_RD) / 2
        if not np.isclose(jaw_R, -OPEN_JAW, atol=1.e-10):  # open position
            return jaw_R

    @jaw_R.setter   # This moves both jaw_RU and jaw_RD in parallel
    def jaw_R(self, val):
        if self.side == 'left' and val is not None:
            val = None
            print("Warning: Ignored value for jaw_R (left-sided collimator).")
        if val is None:
            val = -OPEN_JAW
            self._gap_R = -OPEN_GAP
        diff = val - (self._jaw_RU + self._jaw_RD) / 2
        self._jaw_RU += diff
        self._jaw_RD += diff
        self._update_gaps(only_R=True)

    @property
    def jaw_LU(self):
        if not np.isclose((self._jaw_LU + self._jaw_LD) / 2,
                          OPEN_JAW, atol=1.e-10):  # open position
            return self._jaw_LU

    @jaw_LU.setter   # This assumes jaw_LD remains fixed, hence both jaw_L and the tilt change
    def jaw_LU(self, val):
        if self.side == 'right':
            if val is not None:
                print("Warning: Ignored value for jaw_LU (right-sided collimator).")
            return
        if val is None:
            raise ValueError("Cannot set corner to None! Use open_jaws() or set jaw_L to None.")
        self._jaw_LU = val
        self._update_tilts()   # Extra, to update tilts which are also in C for efficiency
        self._update_gaps(only_L=True)

    @property
    def jaw_LD(self):
        if not np.isclose((self._jaw_LU + self._jaw_LD) / 2,
                          OPEN_JAW, atol=1.e-10):  # open position
            return self._jaw_LD

    @jaw_LD.setter   # This assumes jaw_LU remains fixed, hence both jaw_L and the tilt change
    def jaw_LD(self, val):
        if self.side == 'right':
            if val is not None:
                print("Warning: Ignored value for jaw_LD (right-sided collimator).")
            return
        if val is None:
            raise ValueError("Cannot set corner to None! Use open_jaws() or set jaw_L to None.")
        self._jaw_LD = val
        self._update_tilts()   # Extra, to update tilts which are also in C for efficiency
        self._update_gaps(only_L=True)

    @property
    def jaw_RU(self):
        if not np.isclose((self._jaw_RU + self._jaw_RD) / 2,
                          -OPEN_JAW, atol=1.e-10):  # open position
            return self._jaw_RU

    @jaw_RU.setter   # This assumes jaw_RD remains fixed, hence both jaw_R and the tilt change
    def jaw_RU(self, val):
        if self.side == 'left':
            if val is not None:
                print("Warning: Ignored value for jaw_RU (left-sided collimator).")
            return
        if val is None:
            raise ValueError("Cannot set corner to None! Use open_jaws() or set jaw_R to None.")
        self._jaw_RU = val
        self._update_tilts()   # Extra, to update tilts which are also in C for efficiency
        self._update_gaps(only_R=True)

    @property
    def jaw_RD(self):
        if not np.isclose((self._jaw_RU + self._jaw_RD) / 2,
                          -OPEN_JAW, atol=1.e-10):  # open position
            return self._jaw_RD

    @jaw_RD.setter   # This assumes jaw_RU remains fixed, hence both jaw_R and the tilt change
    def jaw_RD(self, val):
        if self.side == 'left':
            if val is not None:
                print("Warning: Ignored value for jaw_RD (left-sided collimator).")
            return
        if val is None:
            raise ValueError("Cannot set corner to None! Use open_jaws() or set jaw_R to None.")
        self._jaw_RD = val
        self._update_tilts()   # Extra, to update tilts which are also in C for efficiency
        self._update_gaps(only_R=True)

    @property
    def jaw_s_LU(self):
        return self.length/2 * (1 - self._cos_yL)

    @property
    def jaw_s_LD(self):
        return self.length/2 * (1 + self._cos_yL)

    @property
    def jaw_s_RU(self):
        return self.length/2 * (1 - self._cos_yR)

    @property
    def jaw_s_RD(self):
        return self.length/2 * (1 + self._cos_yR)

    def open_jaws(self, keep_tilts=False):
        self.jaw_L = None
        self.jaw_R = None
        if not keep_tilts:
            self.tilt = 0

    def _update_tilts(self):
        if self.side != 'right':
            if self.length > 0:
                self._sin_yL = (self._jaw_LD - self._jaw_LU) / self.length
                self._cos_yL = np.sqrt(1 - self._sin_yL**2)
                self._tan_yL = self._sin_yL / self._cos_yL
            else:
                self._sin_yL = 0.
                self._cos_yL = 1.
                self._tan_yL = 0.
        if self.side != 'left':
            if self.length > 0:
                self._sin_yR = (self._jaw_RD - self._jaw_RU) / self.length
                self._cos_yR = np.sqrt(1 - self._sin_yR**2)
                self._tan_yR = self._sin_yR / self._cos_yR
            else:
                self._sin_yR = 0.
                self._cos_yR = 1.
                self._tan_yR = 0.

    def _update_gaps(self, only_L=False, only_R=False):
        # If we had set a value for the gap manually, this needs to be updated
        # as well after setting the jaw
        if self._gap_L_set_manually() and not only_R:
            self._gap_L = self.gap_L
        if self._gap_R_set_manually() and not only_L:
            self._gap_R = self.gap_R


    # Tilt attributes
    # ===============

    # TODO: tilts are in rad! Do we want that? It's a bit inconsistent with angle which is in deg...

    @property
    def tilt(self):
        if self.side == 'left':
            return self.tilt_L
        elif self.side == 'right':
            return self.tilt_R
        else:
            return self.tilt_L if self.tilt_L==self.tilt_R else [self.tilt_L, self.tilt_R]

    @tilt.setter
    def tilt(self, val):
        if not hasattr(val, '__iter__') or len(val) == 1:
            val = val[0] if hasattr(val, '__iter__') else val
            if self.side == 'left':
                self.tilt_L = val
            elif self.side == 'right':
                self.tilt_R = val
            else:
                self.tilt_L = val
                self.tilt_R = val
        elif len(val) == 2:
            self.tilt_L = val[0]
            self.tilt_R = val[1]
        else:
            raise ValueError(f"The attribute `tilt` should be of the form LR or [L, R] ")

    @property
    def tilt_L(self):
        if self.side != 'right':
            return round(np.arctan2(self._sin_yL, self._cos_yL), 10)

    @tilt_L.setter   # This assumes jaw_L remains fixed (hence jaw_LU and jaw_LD change)
    def tilt_L(self, val):
        if self.side == 'right' and val != 0:
            val = 0
            print("Warning: Ignored value for tilt_L (right-sided collimator).")
        if val != 0:
            print("Warning: Setting a tilt does not preserve the hierarchy, as there "
                + "will always be one corner that tightens (the tilt is applied at "
                + "the centre of the jaw).")
            if val > np.pi/2 or val < -np.pi/2:
                raise ValueError("Tilts larger than 90 degrees are not supported.")
        self._sin_yL = np.sin(val)
        self._cos_yL = np.cos(val)
        self._tan_yL = np.tan(val)
        jaw_L = (self._jaw_LU + self._jaw_LD) / 2
        self._jaw_LD = jaw_L + self._sin_yL * self.length / 2.
        self._jaw_LU = jaw_L - self._sin_yL * self.length / 2.

    @property
    def tilt_R(self):
        if self.side != 'left':
            return round(np.arctan2(self._sin_yR, self._cos_yR), 10)

    @tilt_R.setter   # This assumes jaw_R remains fixed (hence jaw_RU and jaw_RD change)
    def tilt_R(self, val):
        if self.side == 'left' and val != 0:
            val = 0
            print("Warning: Ignored value for tilt_R (left-sided collimator).")
        if val != 0:
            print("Warning: Setting a tilt does not preserve the hierarchy, as there "
                + "will always be one corner that tightens (the tilt is applied at "
                + "the centre of the jaw).")
            if val > np.pi/2 or val < -np.pi/2:
                raise ValueError("Tilts larger than 90 degrees are not supported.")
        self._sin_yR = np.sin(val)
        self._cos_yR = np.cos(val)
        self._tan_yR = np.tan(val)
        jaw_R = (self._jaw_RU + self._jaw_RD) / 2
        self._jaw_RD = jaw_R + self._sin_yR * self.length / 2.
        self._jaw_RU = jaw_R - self._sin_yR * self.length / 2.


    # Optics
    # ======

    @property
    def optics(self):
        return self._optics

    def optics_ready(self):
        return self.emittance is not None and self.optics is not None

    def assign_optics(self, *, nemitt_x=None, nemitt_y=None, beta_gamma_rel=None, name=None, twiss=None,
                      twiss_upstream=None, twiss_downstream=None):
        if nemitt_x is None:
            if self.nemitt_x is None:
                raise ValueError("Need to provide `nemitt_x`.")
        else:
            self.nemitt_x = nemitt_x
        if nemitt_y is None:
            if self.nemitt_y is None:
                raise ValueError("Need to provide `nemitt_y`.")
        else:
            self.nemitt_y = nemitt_y
        if beta_gamma_rel is None:
            raise ValueError("Need to provide `beta_gamma_rel`.")
        if twiss is None:
            if twiss_upstream is None or twiss_downstream is None:
                raise ValueError("Use either `twiss` or `twiss_upstream` and `twiss_downstream`.")
            if name is None:
                if len(twiss_upstream.name) > 1 or len(twiss_downstream.name) > 1:
                    raise ValueError("Need to provide `name` or twisses that are a single row each.")
                tw_up   = twiss_upstream
                tw_down = twiss_downstream
            else:
                tw_up   = twiss_upstream.rows[name]
                tw_down = twiss_downstream.rows[name]
        elif twiss_downstream is not None or twiss_downstream is not None:
            raise ValueError("Use either `twiss` or `twiss_upstream` and `twiss_downstream`.")
        elif name is None:
            raise ValueError("When using `twiss`, need to provide the name as well.")
        else:
            tw_up   = twiss.rows[name]
            tw_down = twiss.rows[twiss.rows.indices[[name]]+1]
        if not np.isclose(tw_up.s[0] + self.length, tw_down.s[0]):
            raise ValueError(f"Downstream twiss not compatible with length {self.length}m.")
        self._optics = {
            'upstream': tw_up,
            'downstream': tw_down,
            'beta_gamma_rel': beta_gamma_rel
        }
        self._apply_optics()

    @property
    def nemitt_x(self):
        if self._nemitt_x == 0:
            return None
        return self._nemitt_x

    @nemitt_x.setter
    def nemitt_x(self, val):
        if val is None:
            val = 0
        elif val <= 0:
            raise ValueError(f"The field `nemitt_x` should be positive, but got {val}.")
        self._nemitt_x = val
        self._apply_optics()

    @property
    def gemitt_x(self):
        if self.nemitt_x is not None and self.optics_ready():
            return self.nemitt_x / self.optics['beta_gamma_rel']

    @property
    def nemitt_y(self):
        if self._nemitt_y == 0:
            return None
        return self._nemitt_y

    @nemitt_y.setter
    def nemitt_y(self, val):
        if val is None:
            val = 0
        elif val <= 0:
            raise ValueError(f"The field `nemitt_y` should be positive, but got {val}.")
        self._nemitt_y = val
        self._apply_optics()

    @property
    def gemitt_y(self):
        if self.nemitt_y is not None and self.optics_ready():
            return self.nemitt_y / self.optics['beta_gamma_rel']

    @property
    def emittance(self):
        if self.nemitt_x is not None and self.nemitt_y is not None:
            if np.isclose(self.nemitt_x, self.nemitt_y):
                return self.nemitt_x
            else:
                return [self.nemitt_x, self.nemitt_y]

    @emittance.setter
    def emittance(self, val):
        if val is None:
            self._nemitt_x = 0
            self._nemitt_y = 0
        else:
            if not hasattr(val, '__iter__'):
                val = [val]
            if len(val) == 1:
                val = [val[0], val[0]]
            assert len(val) == 2
            if val[0] <= 0 or val[1] <= 0:
                raise ValueError(f"The field `emittance` should be positive, but got {val}.")
            self._nemitt_x = val[0]
            self._nemitt_y = val[1]
            self._apply_optics()

    @property
    def sigma(self):
        if self.optics_ready():
            betx = self.optics[self.align]['betx'][0]
            bety = self.optics[self.align]['bety'][0]
            sigma_x = np.sqrt(betx*self.nemitt_x/self.optics['beta_gamma_rel'])
            sigma_y = np.sqrt(bety*self.nemitt_y/self.optics['beta_gamma_rel'])
            if hasattr(self, '_cos_zL'):
                sigma_L = np.sqrt((sigma_x*self._cos_zL)**2 + (sigma_y*self._sin_zL)**2)
                sigma_R = np.sqrt((sigma_x*self._cos_zR)**2 + (sigma_y*self._sin_zR)**2)
                return [sigma_L, sigma_R], [sigma_x, sigma_y]
            else:  # crystal
                sigma = np.sqrt((sigma_x*self._cos_z)**2 + (sigma_y*self._sin_z)**2)
                return sigma, [sigma_x, sigma_y]

    @property
    def co(self):
        if self.optics_ready():
            x = self.optics[self.align]['x'][0]
            y = self.optics[self.align]['y'][0]
            if hasattr(self, '_cos_zL'):
                co_L = x*self._cos_zL + y*self._sin_zL
                co_R = x*self._cos_zR + y*self._sin_zR
                return [co_L, co_R], [x, y]
            else:  # crystal
                co = x*self._cos_z + y*self._sin_z
                return co, [x, y]

    @property
    def co_perp(self):
        if self.optics_ready():
            x = self.optics[self.align]['x'][0]
            y = self.optics[self.align]['y'][0]
            if hasattr(self, '_cos_zL'):
                co_L = -x*self._sin_zL + y*self._cos_zL
                co_R = -x*self._sin_zR + y*self._cos_zR
                return [co_L, co_R], [x, y]
            else:  # crystal
                co = -x*self._sin_z + y*self._cos_z
                return co, [x, y]

    @property
    def divergence(self):
        if self.optics_ready():
            alfx = self.optics[self.align]['alfx'][0]
            alfy = self.optics[self.align]['alfy'][0]
            betx = self.optics[self.align]['betx'][0]
            bety = self.optics[self.align]['bety'][0]
            divx = -np.sqrt(self.gemitt_x/betx)*alfx
            divy = -np.sqrt(self.gemitt_y/bety)*alfy
            if hasattr(self, '_cos_zL'):
                if self.side != 'right':
                    return divx if abs(self.angle_L) < 1e-6 else divy
                else:
                    return divx if abs(self.angle_R) < 1e-6 else divy
            else:
                return divx if abs(self.angle) < 1e-6 else divy

    @property
    def align(self):
        if self._align == 0:
            return 'upstream'
        elif self._align == 1:
            return 'downstream'
        else:
            raise ValueError(f"The attribute `align` can only be 'upstream' or "
                            +f"'downstream', but stored as {self._align}.")

    @align.setter
    def align(self, val):
        if val == 'upstream':
            self._align = 0
        elif val == 'downstream':
            self._align = 1
        else:
            raise ValueError(f"The attribute `align` can only be 'upstream' or "
                            +f"'downstream', but got {val}.")
        self._apply_optics()


    # Gap attributes
    # ==============

    @property
    def gap(self):
        if self.gap_L is None and self.gap_R is None:
            return None
        elif self.gap_R is not None and self.gap_L == -self.gap_R:
            return self.gap_L
        else:
            return [self.gap_L, self.gap_R]

    @gap.setter
    def gap(self, val):
        if not hasattr(val, '__iter__') or len(val) == 1:
            val = val[0] if hasattr(val, '__iter__') else val
            if self.side == 'left':
                self.gap_L = val
            elif self.side == 'right':
                # self.gap_R = -val if val is not None else None
                self.gap_R = val
            else:
                self.gap_L = val
                self.gap_R = -val if val is not None else None
            return
        elif len(val) == 2:
            if not hasattr(val[0], '__iter__') \
            and not hasattr(val[1], '__iter__'):
                if val[0] is not None and val[1] is not None:
                    if val[0] <= val[1]:
                        raise ValueError(f"The attribute `gap_L` should be larger "
                                       + f"than `gap_R` but got {val}.")
                self.gap_L = val[0]
                self.gap_R = val[1]
                return
        # If we got here, val is incompatible
        raise ValueError(f"The attribute `gap` should be of the form `gap` or "
                       + f"`[gap_L, gap_R]`, but got {val}.")

    @property
    def gap_L(self):
        if self.side != 'right':
            if self.optics_ready() and self.jaw_L is not None:
                return round((self.jaw_L - self.co[0][0])/self.sigma[0][0], 6)
            elif self._gap_L_set_manually():
                return self._gap_L

    @gap_L.setter
    def gap_L(self, val):
        if val is None:
            val = OPEN_GAP
            self.jaw_L = None
        elif val <= 0:
            raise ValueError(f"The field `gap_L` should be positive, but got {val}.")
        self._gap_L = val
        self._apply_optics(only_L=True)

    @property
    def gap_R(self):
        if self.side != 'left':
            if self.optics_ready() and self.jaw_R is not None:
                return round((self.jaw_R - self.co[0][1])/self.sigma[0][1], 6)
            elif self._gap_R_set_manually():
                return self._gap_R

    @gap_R.setter
    def gap_R(self, val):
        if val is None:
            val = -OPEN_GAP
            self.jaw_R = None
        elif val >= 0:
            raise ValueError(f"The field `gap_R` should be negative, but got {val}.")
        self._gap_R = val
        self._apply_optics(only_R=True)

    @property
    def gap_LU(self):
        if self.gap_L is not None and self.optics_ready():
            return round(self._gap_L - self._sin_yL * self.length / 2. / self.sigma[0][0], 6)

    @property
    def gap_LD(self):
        if self.gap_L is not None and self.optics_ready():
            return round(self._gap_L + self._sin_yL * self.length / 2. / self.sigma[0][0], 6)

    @property
    def gap_RU(self):
        if self.gap_R is not None and self.optics_ready():
            return round(self._gap_R - self._sin_yR * self.length / 2. / self.sigma[0][1], 6)

    @property
    def gap_RD(self):
        if self.gap_R is not None and self.optics_ready():
            return round(self._gap_R + self._sin_yR * self.length / 2. / self.sigma[0][1], 6)

    def _gap_L_set_manually(self):
        return not np.isclose(self._gap_L, OPEN_GAP)

    def _gap_R_set_manually(self):
        return not np.isclose(self._gap_R, -OPEN_GAP)

    def _apply_optics(self, only_L=False, only_R=False):
        if self.optics_ready():
            # Only if we have set a value for the gap manually, this needs to be updated
            if self._gap_L_set_manually() and not only_R:
                self.jaw_L = self._gap_L * self.sigma[0][0] + self.co[0][0]
            if self._gap_R_set_manually() and not only_L:
                self.jaw_R = self._gap_R * self.sigma[0][1] + self.co[0][1]


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
                self.gap_R = None
                return
            elif val.lower() == 'right' or val.lower() == 'r' or val == '-':
                self._side = -1
                self.gap_L = None
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

    def generate_pencil(self, num_particles, *, side='+-', pencil_spread=1e-6,
                        impact_parameter=0, sigma_z=7.61e-2, twiss=None, longitudinal=None,
                        longitudinal_betatron_cut=None, tw=None, **kwargs):
        if not hasattr(self, '_line') or not hasattr(self, '_name'):
            raise ValueError("Collimator is missing a pointer to the line. Install collimators "
                           + "with `line.collimators.install()` (or use "
                           + "`xcoll.initial_distribution.generate_pencil_on_collimator()`).")
        from xcoll.initial_distribution import generate_pencil_on_collimator
        return generate_pencil_on_collimator(line=self._line, name=self._name, side=side,
                        num_particles=num_particles, pencil_spread=pencil_spread, tw=tw,
                        impact_parameter=impact_parameter, sigma_z=sigma_z, twiss=twiss,
                        longitudinal=longitudinal, longitudinal_betatron_cut=longitudinal_betatron_cut,
                        **kwargs)

    def generate_delta(self, *, plane, position_mm, nemitt_x, nemitt_y, betatron_cut=0,
                       match_at_front=True, twiss=None):
        if not hasattr(self, '_line') or not hasattr(self, '_name'):
            raise ValueError("Collimator is missing a pointer to the line. Install collimators "
                           + "with `line.collimators.install()` (or use "
                           + "`xcoll.initial_distribution.generate_delta_from_dispersion()`).")
        from xcoll.initial_distribution import generate_delta_from_dispersion
        return generate_delta_from_dispersion(line=self._line, at_element=self._name, plane=plane,
                        position_mm=position_mm, nemitt_x=nemitt_x, nemitt_y=nemitt_y, twiss=twiss,
                        betatron_cut=betatron_cut, match_at_front=match_at_front)

    def _verify_consistency(self):
        BaseBlock._verify_consistency(self)
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
        if self.side != 'right':
            ang = abs(np.arccos(self._cos_yL))
            ang = np.pi - ang if ang > np.pi/2 else ang
            assert np.isclose(ang, abs(np.arcsin(self._sin_yL)))
            assert np.isclose(self._sin_yL/self._cos_yL, self._tan_yL)
        if self.side != 'left':
            ang = abs(np.arccos(self._cos_yR))
            ang = np.pi - ang if ang > np.pi/2 else ang
            assert np.isclose(ang, abs(np.arcsin(self._sin_yR)))
            assert np.isclose(self._sin_yR/self._cos_yR, self._tan_yR)

        # Verify bools
        if '_side' in self._xofields:  # Not the case e.g. for FlukaCollimator
            assert self._side in [-1, 1, 0]
        assert isinstance(self._jaws_parallel, bool) or self._jaws_parallel in [0, 1]

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


class BaseCrystal(BaseBlock):
    _xofields = BaseBlock._xofields | {
        # Collimator angle
        '_sin_z':             xo.Float64,
        '_cos_z':             xo.Float64,
        # Jaw corners (this is the x-coordinate in the rotated frame)
        '_jaw_U':             xo.Float64,
        # Tilts (not superfluous)
        '_sin_y':             xo.Float64,
        '_cos_y':             xo.Float64,
        '_tan_y':             xo.Float64,
        # Other
        '_side':              xo.Int8,
        # These are not used in C, but need to be an xofield to get them in the to_dict:
        '_align':             xo.Int8,
        '_gap':               xo.Float64,
        '_nemitt_x':          xo.Float64,
        '_nemitt_y':          xo.Float64,
        # Crystal specific
        '_bending_radius':    xo.Float64,
        '_bending_angle':     xo.Float64,
        '_width':              xo.Float64,
        '_height':             xo.Float64
        # 'thick':              xo.Float64
    }


    isthick = True
    needs_rng = False
    allow_track = False
    allow_double_sided = False
    behaves_like_drift = True
    allow_rot_and_shift = False
    allow_loss_refinement = True
    skip_in_loss_location_refinement = True

    _noexpr_fields    = {'align', 'side', 'name'}
    _skip_in_to_dict  = [*BaseBlock._skip_in_to_dict, *[f for f in _xofields if f.startswith('_')]]
    _store_in_to_dict = [*BaseBlock._store_in_to_dict, 'angle', 'jaw', 'tilt', 'gap', 'side', 'align',
                         'emittance', 'width', 'height', 'bending_radius', 'bending_angle']

    _depends_on = [BaseCollimator]

    _internal_record_class = BaseBlock._internal_record_class

    # This is an abstract class and cannot be instantiated
    def __new__(cls, *args, **kwargs):
        if cls is BaseCrystal:
            raise Exception("Abstract class `BaseCrystal` cannot be instantiated!")
        instance = super().__new__(cls)
        return instance

    def __init__(self, **kwargs):
        to_assign = {}
        if '_xobject' not in kwargs:
            # Set side
            to_assign['side'] = kwargs.pop('side', 'left')

            # Set angle
            to_assign['angle'] = kwargs.pop('angle', 0)

            # We do not allow any combination of jaw_ and gap_ attributes
            # (except when jaw=..., gap=None or jaw=None, gap=... is used, as this is how the colldb installs it)
            kwargs = {kk: vv for kk, vv in kwargs.items() if not vv is None}


            # Set jaw
            if 'jaw' in kwargs:
                for key in ['jaw_U', 'jaw_D', 'gap']:
                    if key in kwargs:
                        raise ValueError(f"Cannot use both `jaw` and `{key}`!")
                to_assign['jaw'] = kwargs.pop('jaw')
            elif 'jaw_D' in kwargs:
                for key in ['tilt', 'gap']:
                    if key in kwargs:
                        raise ValueError(f"Cannot use both `jaw_D` with `{key}`!")
                if not 'jaw_U' in kwargs:
                    raise ValueError("Need to provide `jaw_U` when setting `jaw_D`!")
                to_assign['jaw_U'] = kwargs.pop('jaw_U')
                to_assign['jaw_D'] = kwargs.pop('jaw_D')
            elif 'jaw_U' in kwargs:
                if 'gap' in kwargs:
                    raise ValueError(f"Cannot use both `jaw_U` and `gap`!")
                to_assign['jaw_U'] = kwargs.pop('jaw_U')
            # TODO: correct sign if right-sided
            kwargs.setdefault('_jaw_U', OPEN_JAW)

            # Set gap
            if 'gap' in kwargs:
                for key in ['jaw', 'jaw_U', 'jaw_D']:
                    if key in kwargs:
                        raise ValueError(f"Cannot use both `gap` and `{key}`!")
                to_assign['gap'] = kwargs.pop('gap')
            # TODO: correct sign if right-sided
            kwargs.setdefault('_gap', OPEN_GAP)

            # Set tilt
            if 'tilt' in kwargs:
                to_assign['tilt'] = kwargs.pop('tilt')

            # Set others
            to_assign['align'] = kwargs.pop('align', 'upstream')
            to_assign['emittance'] = kwargs.pop('emittance', None)
            kwargs.setdefault('active', True)
            kwargs.setdefault('_sin_y', 0)
            kwargs.setdefault('_cos_y', 1)

            # Set crystal specific
            if 'bending_angle' in kwargs:
                if 'bending_radius' in kwargs:
                    raise ValueError("Need to choose between 'bending_radius' and 'bending_angle'!")
                to_assign['bending_angle'] = kwargs.pop('bending_angle')
            else:
                to_assign['bending_radius'] = kwargs.pop('bending_radius', 1)
            to_assign['width'] = kwargs.pop('width', 1)
            to_assign['height'] = kwargs.pop('height', 1)

        super().__init__(**kwargs)
        # Careful: non-xofields are not passed correctly between copy's / to_dict. This messes with flags etc..
        # We also have to manually initialise them for xobject generation
        if not hasattr(self, '_optics'):
            self._optics = None
        for key, val in to_assign.items():
            setattr(self, key, val)
        if self.side == 'right':
            if np.isclose(self._jaw_U, OPEN_JAW):
                self._jaw_U *= -1
            if np.isclose(self._gap, OPEN_GAP):
                self._gap *= -1
        BaseCrystal._verify_consistency(self)


    # Main crystal angle
    # ==================

    @property
    def angle(self):
        return round(np.rad2deg(np.arctan2(self._sin_z, self._cos_z)), 10)

    @angle.setter
    def angle(self, val):
        self._sin_z = np.sin(np.deg2rad(val))
        self._cos_z = np.cos(np.deg2rad(val))
        self._apply_optics()


    # Jaw attributes
    # ==============

    @property
    def jaw(self):
        return self.jaw_U

    @jaw.setter
    def jaw(self, val):
        if val is None:
            if self.side == 'left':
                val = OPEN_JAW
            elif self.side == 'right':
                val = -OPEN_JAW
            else:
                raise ValueError("Cannot determine side. Something is wrong with the collimator!")
        self.jaw_U = val

    @property
    def jaw_U(self):
        if not np.isclose(abs(self._jaw_U), OPEN_JAW, atol=1.e-10):  # open position
            return self._jaw_U

    @jaw_U.setter
    def jaw_U(self, val):
        if val is None:
            raise ValueError("Cannot set corner to None! Use open_jaws() or set jaw to None.")
        self._jaw_U = val
        self._update_gaps()

    @property
    def jaw_D(self):
        if not np.isclose(abs(self._jaw_U), OPEN_JAW, atol=1.e-10):  # open position
            length = self.length
            if (self.side == 'left' and self.bending_radius < 0) \
            or (self.side == 'right' and self.bending_radius > 0):
                # Correction for inner corner point
                length -= self.width*np.sin(abs(self._bending_angle))
            shift = np.tan(self._bending_angle/2)*self._cos_y + self._sin_y
            return self._jaw_U + length*shift

    @jaw_D.setter
    def jaw_D(self, val):
        if val is None:
            self.tilt = 0
        else:
            shift = (val - self._jaw_U )/self.length * np.cos(self._bending_angle/2)
            self._sin_y  = shift*np.cos(self._bending_angle/2)
            self._sin_y -= np.sin(self._bending_angle/2)*np.sqrt(1 - shift**2)
            self._cos_y  = np.sqrt(1 - self._sin_y**2)
            self._tan_y = self._sin_y / self._cos_y
        self._update_gaps()

    def open_jaws(self, keep_tilts=False):
        self.jaw = None
        if not keep_tilts:
            self.tilt = 0

    def _update_gaps(self):
        # If we had set a value for the gap manually, this needs to be updated
        # as well after setting the jaw
        if self._gap_set_manually():
            self._gap = self.gap


    # Tilt attributes
    # ===============

    # TODO: tilts are in rad! Do we want that? It's a bit inconsistent with angle which is in deg...

    @property
    def tilt(self):
        return round(np.arctan2(self._sin_y, self._cos_y), 10)

    @tilt.setter   # This assumes jaw_U remains fixed (hence jaw_D changes)
    def tilt(self, val):
        if self.side == 'left':
            if val < min(0, self.bending_angle/2):
                print("Warning: Setting a negative tilt does not preserve the hierarchy, as the "
                    + "crystal tightens towards the beam.")
        elif self.side == 'right':
            if val > min(0, -self.bending_angle/2):
                print("Warning: Setting a positive tilt does not preserve the hierarchy, as the "
                    + "crystal tightens towards the beam.")
        if val > np.pi/2 or val < -np.pi/2:
            raise ValueError("Tilts larger than 90 degrees are not supported.")
        self._sin_y = np.sin(val)
        self._cos_y = np.cos(val)
        self._tan_y = np.tan(val)


    # Optics
    # ======

    @property
    def optics(self):
        return self._optics

    def optics_ready(self):
        return BaseCollimator.optics_ready(self)

    def assign_optics(self, *, nemitt_x=None, nemitt_y=None, beta_gamma_rel=None, name=None, twiss=None,
                      twiss_upstream=None, twiss_downstream=None):
        return BaseCollimator.assign_optics(self, nemitt_x=nemitt_x, nemitt_y=nemitt_y,
                                            beta_gamma_rel=beta_gamma_rel, name=name, twiss=twiss,
                                            twiss_upstream=twiss_upstream, twiss_downstream=twiss_downstream)

    @property
    def nemitt_x(self):
        return BaseCollimator.nemitt_x.fget(self)

    @nemitt_x.setter
    def nemitt_x(self, val):
        BaseCollimator.nemitt_x.fset(self, val)

    @property
    def gemitt_x(self):
        return BaseCollimator.gemitt_x.fget(self)

    @property
    def nemitt_y(self):
        return BaseCollimator.nemitt_y.fget(self)

    @nemitt_y.setter
    def nemitt_y(self, val):
        BaseCollimator.nemitt_y.fset(self, val)

    @property
    def gemitt_y(self):
        return BaseCollimator.gemitt_y.fget(self)

    @property
    def emittance(self):
        return BaseCollimator.emittance.fget(self)

    @emittance.setter
    def emittance(self, val):
        BaseCollimator.emittance.fset(self, val)

    @property
    def sigma(self):
        return BaseCollimator.sigma.fget(self)

    @property
    def co(self):
        return BaseCollimator.co.fget(self)

    @property
    def divergence(self):
        return BaseCollimator.divergence.fget(self)

    @property
    def align(self):
        return BaseCollimator.align.fget(self)

    @align.setter
    def align(self, val):
        if val != 'upstream':
            raise NotImplementedError("Crystals cannot be aligned to the downstream optics!")
        BaseCollimator.align.fset(self, val)

    def align_to_beam_divergence(self):
        if not self.optics_ready():
            raise ValueError("Optics not assigned! Cannot align to beam divergence.")
        if self.gap is None:
            raise ValueError("Need to set `gap` to align to beam divergence.")
        self.tilt = self.divergence * self.gap


    # Gap attributes
    # ==============

    @property
    def gap(self):
        if self.optics_ready() and self.jaw_U is not None:
            return round((self.jaw_U - self.co[0])/self.sigma[0], 6)
        elif not self._gap_set_manually():
            return None
        else:
            return self._gap

    # Gap is always positive, irrespective of the side
    @gap.setter
    def gap(self, val):
        if val is None:
            val = OPEN_GAP
            self.jaw = None
        if hasattr(val, '__iter__'):
            raise ValueError("The attribute `gap` should be a single value, not a list.")
        if val <= 0:
            raise ValueError(f"The field `gap` should be positive, but got {val}.")
        self._gap = val
        self._apply_optics()

    def _gap_set_manually(self):
        return not np.isclose(self._gap, OPEN_GAP)

    def _apply_optics(self):
        if self.optics_ready():
            # Only if we have set a value for the gap manually, this needs to be updated
            if self._gap_set_manually():
                self.jaw_U = self._gap * self.sigma[0] + self.co[0]


    # Other attributes
    # ================

    @property
    def bending_radius(self):
        return self._bending_radius

    @bending_radius.setter
    def bending_radius(self, bending_radius):
        bending_angle = np.arcsin(self.length/bending_radius)
        if abs(bending_angle) > np.pi/2:
            raise ValueError("Bending angle cannot be larger than 90 degrees!")
        self._bending_radius = bending_radius
        self._bending_angle = bending_angle

    @property
    def bending_angle(self):
        return self._bending_angle

    @bending_angle.setter
    def bending_angle(self, bending_angle):
        if abs(bending_angle) > np.pi/2:
            raise ValueError("Bending angle cannot be larger than 90 degrees!")
        self._bending_angle = bending_angle
        self._bending_radius = self.length / np.sin(bending_angle)

    @property
    def width(self):
        return self._width

    @width.setter
    def width(self, val):
        if val <= 0:
            raise ValueError(f"The field `width` should be positive, but got {val}.")
        self._width = val

    @property
    def height(self):
        return self._height

    @height.setter
    def height(self, val):
        if val <= 0:
            raise ValueError(f"The field `height` should be positive, but got {val}.")
        self._height = val

    @property
    def side(self):
        return BaseCollimator.side.fget(self)

    @side.setter
    def side(self, val):
        temp = self._side
        BaseCollimator.side.fset(self, val)
        if self._side == 0:
            self._side = temp
            raise ValueError("Crystal cannot be two-sided! Please set `side` "
                           + "to 'left' or 'right'.")


    # Methods
    # =======

    def generate_pencil(self, **kwargs):
        return BaseCollimator.generate_pencil(self, **kwargs)

    def _verify_consistency(self):
        BaseBlock._verify_consistency(self)
        # Verify angles
        ang = abs(np.arccos(self._cos_z))
        ang = np.pi - ang if ang > np.pi/2 else ang
        assert np.isclose(ang, abs(np.arcsin(self._sin_z)))
        ang = abs(np.arccos(self._cos_y))
        ang = np.pi - ang if ang > np.pi/2 else ang
        assert np.isclose(ang, abs(np.arcsin(self._sin_y)))
        assert np.isclose(self._sin_y/self._cos_y, self._tan_y)
        # Verify bools
        if '_side' in self._xofields:
            assert self._side in [-1, 0, 1]
        # Crystal specific
        if '_bending_radius' in self._xofields and '_bending_angle' in self._xofields:
            assert isinstance(self._bending_radius, float) and not np.isclose(self._bending_radius, 0)
            assert isinstance(self._bending_angle, float) and abs(self._bending_angle) <= np.pi/2
            assert np.isclose(self._bending_angle, np.arcsin(self.length/self._bending_radius))

