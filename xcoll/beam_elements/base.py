# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

import numpy as np

import xobjects as xo
import xtrack as xt

from ..interaction_record import InteractionRecord
from ..general import _pkg_root


OPEN_JAW = 3.
OPEN_GAP = 999.


class InvalidXcoll(xt.BeamElement):
    _xofields = {
        'length': xo.Float64
    }

    isthick = True
    behaves_like_drift = True
    allow_track = False
    skip_in_loss_location_refinement = True
    allow_loss_refinement = True

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
        'length':             xo.Float64,
        'record_touches':     xo.Int8,
        'record_scatterings': xo.Int8
    }

    isthick = True
    allow_track = False
    behaves_like_drift = True
    skip_in_loss_location_refinement = True

    _depends_on = [InvalidXcoll]

    _internal_record_class = InteractionRecord

    # This is an abstract class and cannot be instantiated
    def __new__(cls, *args, **kwargs):
        if cls == BaseBlock:
            raise Exception("Abstract class `BaseBlock` cannot be instantiated!")
        instance = super().__new__(cls)
        return instance

    def __init__(self, **kwargs):
        if '_xobject' not in kwargs:
            kwargs.pop('use_prebuilt_kernels', None)
        super().__init__(**kwargs)

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
        'active':         xo.Int8,
        'record_touches':     xo.Int8,
        'record_scatterings': xo.Int8,
        # These are not used in C, but need to be an xofield to get them in the to_dict:
        '_align':         xo.Int8,
        '_gap_L':         xo.Float64,
        '_gap_R':         xo.Float64,
        '_nemitt_x':      xo.Float64,
        '_nemitt_y':      xo.Float64
    }

    isthick = True
    allow_track = False
    behaves_like_drift = True
    skip_in_loss_location_refinement = True

    _skip_in_to_dict  = [f for f in _xofields if f.startswith('_')]
    _store_in_to_dict = ['angle', 'jaw', 'tilt', 'gap', 'side', 'align', 'emittance']

    _depends_on = [InvalidXcoll, xt.Drift, xt.XYShift, xt.SRotation, xt.YRotation]

    _extra_c_sources = [
        _pkg_root.joinpath('beam_elements','collimators_src','collimator_geometry.h')
    ]

    _internal_record_class = InteractionRecord


    # This is an abstract class and cannot be instantiated
    def __new__(cls, *args, **kwargs):
        if cls == BaseCollimator:
            raise Exception("Abstract class `BaseCollimator` cannot be instantiated!")
        instance = super().__new__(cls)
        return instance

    def __init__(self, **kwargs):
        to_assign = {}
        if '_xobject' not in kwargs:
            kwargs.pop('use_prebuilt_kernels', None)
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

            # Set jaw
            if 'jaw' in kwargs:
                for key in ['jaw_L', 'jaw_R', 'jaw_LU', 'jaw_LD', 'jaw_RU', 'jaw_RD', 'gap', 'gap_L', 'gap_R']:
                    if key in kwargs:
                        raise ValueError(f"Cannot use both `jaw` and `{key}`!")
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
            else:
                to_assign['tilt_L'] = kwargs.pop('tilt_L', 0)
                to_assign['tilt_R'] = kwargs.pop('tilt_R', 0)

            # Set gap
            if 'gap' in kwargs:
                for key in ['jaw', 'jaw_L', 'jaw_R', 'jaw_LU', 'jaw_LD', 'jaw_RU', 'jaw_RD', 'gap_L', 'gap_R']:
                    if key in kwargs:
                        raise ValueError(f"Cannot use both `gap` and `{key}`!")
                to_assign['gap'] = kwargs.pop('gap')
            elif 'gap_L' in kwargs or 'gap_R' in kwargs:
                for key in ['jaw', 'jaw_L', 'jaw_R', 'jaw_LU', 'jaw_LD', 'jaw_RU', 'jaw_RD', 'gap']:
                    if key in kwargs:
                        raise ValueError(f"Cannot use both `gap` and `{key}`!")
                to_assign['gap_L'] = kwargs.pop('gap_L', None)
                to_assign['gap_R'] = kwargs.pop('gap_R', None)
            kwargs.setdefault('_gap_L', OPEN_GAP)
            kwargs.setdefault('_gap_R', -OPEN_GAP)

            # Set others
            to_assign['align'] = kwargs.pop('align', 'upstream')
            to_assign['emittance'] = kwargs.pop('emittance', 0)
            kwargs.setdefault('active', True)
            kwargs.setdefault('record_touches', False)
            kwargs.setdefault('record_scatterings', False)

        super().__init__(**kwargs)
        # Careful: non-xo fields are not passed correctly between copy's / to_dict. This messes with flags etc..
        # We also have to manually initialise them for xobject generation
        if not hasattr(self, '_optics'):
            self._optics = None
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
        self._apply_optics()

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
            return [[self.jaw_LU, self.jaw_RU], [self.jaw_LD, self.jaw_RD]]

    @jaw.setter   # Keeps the tilts unless all 4 corners are specified
    def jaw(self, val):
        if not hasattr(val, '__iter__') or len(val) == 1:
            val = val[0] if hasattr(val, '__iter__') else val
            if self.side == 'left':
                self.jaw_L = val
            elif self.side == 'right':
                self.jaw_R = -val if val is not None else None
            else:
                self.jaw_L = val
                self.jaw_R = -val if val is not None else None
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
        self._update_gaps()

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
            self._gap_R = OPEN_GAP
        diff = val - (self._jaw_RU + self._jaw_RD) / 2
        self._jaw_RU += diff
        self._jaw_RD += diff
        self._update_gaps()

    @property
    def jaw_LU(self):
        if not self.jaw_L is None:
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
        self._update_gaps()

    @property
    def jaw_LD(self):
        if not self.jaw_L is None:
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
        self._update_gaps()

    @property
    def jaw_RU(self):
        if not self.jaw_R is None:
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
        self._update_gaps()

    @property
    def jaw_RD(self):
        if not self.jaw_R is None:
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
        self._update_gaps()

    def open_jaws(self, keep_tilts=False):
        self.jaw_L = None
        self.jaw_R = None
        if not keep_tilts:
            self.tilt = 0

    def _update_tilts(self):
        if self.side != 'right':
            self._sin_yL = (self.jaw_LD - self.jaw_LU) / self.length
            self._cos_yL = np.sqrt(1 - self._sin_yL**2)
            self._tan_yL = self._sin_yL / self._cos_yL
        if self.side != 'left':
            self._sin_yR = (self.jaw_RD - self.jaw_RU) / self.length
            self._cos_yR = np.sqrt(1 - self._sin_yR**2)
            self._tan_yR = self._sin_yR / self._cos_yR

    def _update_gaps(self):
        # If we had set a value for the gap manually, this needs to be updated
        # as well after setting the jaw
        if self._gap_L_set_manually():
            self._gap_L = self.gap_L
        if self._gap_R_set_manually():
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
            raise ValueError

    @property
    def tilt_L(self):
        if self.side != 'right':
            return round(np.arctan2(self._sin_yL, self._cos_yL), 10)

    @tilt_L.setter   # This assumes jaw_L remains fixed (hence jaw_LU and jaw_LD change)
    def tilt_L(self, val):
        if self.side == 'right' and val != 0:
            val = 0
            print("Warning: Ignored value for tilt_L (right-sided collimator).")
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
        from xcoll import element_classes
        if not isinstance(self, element_classes):
            raise ValueError("Please install collimator before assigning optics.")
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
            tw_down = twiss.rows[twiss.mask[[name]]+1]
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
        self._nemitt_x = val
        self._apply_optics()

    @property
    def nemitt_y(self):
        if self._nemitt_y == 0:
            return None
        return self._nemitt_y

    @nemitt_y.setter
    def nemitt_y(self, val):
        self._nemitt_y = val
        self._apply_optics()

    @property
    def emittance(self):
        if self.nemitt_x is not None and self.nemitt_y is not None:
            if np.isclose(self.nemitt_x, self.nemitt_y):
                return self.nemitt_x
            else:
                return [self.nemitt_x, self.nemitt_y]

    @emittance.setter
    def emittance(self, val):
        if not hasattr(val, '__iter__'):
            val = [val]
        if len(val) == 1:
            val = [val[0], val[0]]
        assert len(val) == 2
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
            sigma_L = np.sqrt((sigma_x*self._cos_zL)**2 + (sigma_y*self._sin_zL)**2)
            sigma_R = np.sqrt((sigma_x*self._cos_zR)**2 + (sigma_y*self._sin_zR)**2)
            return [sigma_L, sigma_R], [sigma_x, sigma_y]

    @property
    def co(self):
        if self.optics_ready():
            x = self.optics[self.align]['x'][0]
            y = self.optics[self.align]['y'][0]
            co_L = x*self._cos_zL + y*self._sin_zL
            co_R = x*self._cos_zR + y*self._sin_zR
            return [co_L, co_R], [x, y]
    
    @property
    def divergence(self):
        if self.optics_ready():
            alfx = self.optics[self.align]['alfx'][0]
            alfy = self.optics[self.align]['alfy'][0]
            betx = self.optics[self.align]['betx'][0]
            bety = self.optics[self.align]['bety'][0]
            divx = -np.sqrt(self.nemitt_x/self.optics['beta_gamma_rel']/betx)*alfx
            divy = -np.sqrt(self.nemitt_y/self.optics['beta_gamma_rel']/bety)*alfy
            if self.side != 'right':
                return divx if abs(self.angle_L) < 1e-6 else divy
            else:
                return divx if abs(self.angle_R) < 1e-6 else divy

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
                self.gap_R = -val if val is not None else None
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
            elif np.isclose(self._gap_L, OPEN_GAP):
                return None
            else:
                return self._gap_L

    @gap_L.setter
    def gap_L(self, val):
        if val is None:
            val = OPEN_GAP
            self.jaw_L = None
        self._gap_L = val
        self._apply_optics()

    @property
    def gap_R(self):
        if self.side != 'left':
            if self.optics_ready() and self.jaw_R is not None:
                return round((self.jaw_R - self.co[0][1])/self.sigma[0][1], 6)
            elif np.isclose(self._gap_R, -OPEN_GAP):
                return None
            else:
                return self._gap_R

    @gap_R.setter
    def gap_R(self, val):
        if val is None:
            val = -OPEN_GAP
            self.jaw_R = None
        self._gap_R = val
        self._apply_optics()

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

    def _apply_optics(self):
        if self.optics_ready():
            # Only if we have set a value for the gap manually, this needs to be updated
            if self._gap_L_set_manually():
                self.jaw_L = self._gap_L * self.sigma[0][0] + self.co[0][0]
            if self._gap_R_set_manually():
                self.jaw_R = self._gap_R * self.sigma[0][1] + self.co[0][1]
            if hasattr(self, 'align_angle'):
                if self.gap_L is not None: # Hack to do in this order because crystal is onesided; but side might not be assigned yet
                    self.align_angle = self.divergence * self.gap_L
                elif self.gap_R is not None:
                    self.align_angle = self.divergence * self.gap_R


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

    def enable_scattering(self):
        if hasattr(self, '_tracking'):
            if self.optics is None:
                raise ValueError("Optics not assigned! Cannot enable scattering.")
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
        if self.side == 'both' and abs(self.tilt_L - self.tilt_R) >= 90.:
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
        assert isinstance(self.record_scatterings, bool) or self.record_scatterings in [0, 1]

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

