# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

import numpy as np
import xobjects as xo
import xtrack as xt

from .base import BaseCollimator, BaseCrystal
from ..general import _pkg_root
from ..xaux import track_construction, super_if_being_constructed
from ..scattering_routines.geant4 import Geant4Engine, track_pre, track_core, track_post
from ..materials import _DEFAULT_MATERIAL, _resolve_material


@track_construction
class Geant4Collimator(BaseCollimator):
    _xofields = BaseCollimator._xofields | {
        'geant4_id': xo.String,
        'length_front': xo.Float64,   # Hard-coded to correct 250nm margin in BDSIM
        'length_back':  xo.Float64
    }

    isthick = True
    needs_cpu = True    # TODO: adapt track.py to work with particles on GPU context
    allow_track = True
    iscollective = True
    behaves_like_drift = True
    allow_rot_and_shift = False
    allow_loss_refinement = True
    skip_in_loss_location_refinement = True

    _depends_on = [BaseCollimator, Geant4Engine]

    _extra_c_sources = [
        _pkg_root.joinpath('beam_elements','elements_src','geant4_collimator.h')
    ]

    _noexpr_fields         = {*BaseCollimator._noexpr_fields, 'material'}
    _skip_in_to_dict       = BaseCollimator._skip_in_to_dict
    _store_in_to_dict      = [*BaseCollimator._store_in_to_dict, 'material']
    _internal_record_class = BaseCollimator._internal_record_class
    _allowed_fields_when_frozen = BaseCollimator._allowed_fields_when_frozen

    def __init__(self, **kwargs):
        import xcoll as xc
        if xc.geant4.engine.is_running():
            raise ValueError('Cannot create Geant4Collimator while engine is running.')
        to_assign = {}
        if '_xobject' not in kwargs:
            kwargs.setdefault('geant4_id', ''.ljust(16))
            to_assign['name'] = xc.geant4.engine._get_new_element_name()
            to_assign['material'] = kwargs.pop('material', None)
            kwargs['_material'] = _DEFAULT_MATERIAL
        super().__init__(**kwargs)
        for key, val in to_assign.items():
            setattr(self, key, val)
        if not hasattr(self, '_equivalent_drift'):
            self._equivalent_drift = xt.Drift(length=self.length)
            self._equivalent_drift.model = 'exact'
        self.length_front = 250e-9
        self.length_back = -250e-9

    @property
    def angle(self):
        return BaseCollimator.angle.fget(self)

    @angle.setter
    def angle(self, val):
        if hasattr(val, '__iter__') and len(val) == 2 and val[0] != val[1]:
            raise ValueError('The Geant4 scattering engine does not '
                           + 'support unequal jaw rotation angles')
        BaseCollimator.angle.fset(self, val)

    @property
    def material(self):
        if self._material != _DEFAULT_MATERIAL:
            return self._material

    @material.setter
    def material(self, material):
        material = _resolve_material(material, ref='geant4')
        if self.material != material:
            self._material = material

    def enable_scattering(self):
        import xcoll as xc
        xc.geant4.interface.assert_environment_ready()
        if not xc.geant4.engine.is_running():
            raise RuntimeError("Geant4 engine is not running.")
        super().enable_scattering()

    def track(self, part):
        if track_pre(self, part):
            super().track(part)
            track_core(self, part)
            track_post(self, part)
        else:
            self._drift(part)

    def _drift(self, particles, length=None):
        if length is None:
            length = self.length
        if length != self.length:
            old_length = self._equivalent_drift.length
            self._equivalent_drift.length = length
        self._equivalent_drift.track(particles)
        if length != self.length:
            self._equivalent_drift.length = old_length

    @super_if_being_constructed
    def __setattr__(self, name, value):
        import xcoll as xc
        if name not in self._allowed_fields_when_frozen \
        and xc.geant4.engine.is_running():
            raise ValueError('Engine is running; Geant4Collimator is frozen.')
        super().__setattr__(name, value)


@track_construction
class Geant4CollimatorTip(Geant4Collimator):
    _xofields = Geant4Collimator._xofields | {
        'tip_thickness': xo.Float64,
        '_tip_length': xo.Float64,  # m, length of the jaw face at the beam for a tapered jaw; 0 (default)
                                    # means a flat (untapered) jaw, as before. length is always the full
                                    # length of the jaw (see tip_length and taper_angle)
        'jaw_depth': xo.Float64,    # m, transverse depth of each jaw from its edge; 0 (default)
                                    # means a 2 m wide Geant4 box, as before
        'taper_depth': xo.Float64,  # m, depth from the jaw edge where the taper ends, the jaw keeping
                                    # a constant length beyond it; 0 (default) tapers the whole jaw
        'jaw_width': xo.Float64     # m, full width of each jaw across its plane of motion; 0 (default)
                                    # means BDSIM's default of 12 mm
    }

    isthick = True
    needs_cpu = True
    allow_track = True
    iscollective = True
    behaves_like_drift = True
    skip_in_loss_location_refinement = True

    _depends_on = [*Geant4Collimator._depends_on]

    _noexpr_fields         = {*Geant4Collimator._noexpr_fields, 'tip_material'}
    _skip_in_to_dict       = Geant4Collimator._skip_in_to_dict
    _store_in_to_dict      = [*Geant4Collimator._store_in_to_dict, 'tip_material', 'tip_length']
    _internal_record_class = Geant4Collimator._internal_record_class
    _allowed_fields_when_frozen = Geant4Collimator._allowed_fields_when_frozen

    _extra_c_sources = [
        _pkg_root.joinpath('beam_elements', 'elements_src', 'geant4_collimator_tip.h')
    ]

    def __init__(self, **kwargs):
        import xcoll as xc
        if xc.geant4.engine.is_running():
            raise ValueError('Cannot create Geant4CollimatorTip while engine is running.')
        to_assign = {}
        if '_xobject' not in kwargs:
            to_assign['tip_material'] = kwargs.pop('tip_material', None)
            kwargs['_tip_material'] = _DEFAULT_MATERIAL
            kwargs['_tip_length'] = _resolve_taper(kwargs)
        super().__init__(**kwargs)
        for key, val in to_assign.items():
            setattr(self, key, val)

    @property
    def tip_material(self):
        if self._tip_material != _DEFAULT_MATERIAL:
            return self._tip_material

    @tip_material.setter
    def tip_material(self, tip_material):
        tip_material = _resolve_material(tip_material, ref='geant4')
        if self.tip_material != tip_material:
            self._tip_material = tip_material

    @property
    def tip_length(self):
        """Length of the jaw face at the beam [m]; the element length if not tapered."""
        return self._tip_length if self._tip_length > 0 else self.length

    @tip_length.setter
    def tip_length(self, val):
        # A tip as long as the element is a flat (untapered) jaw
        self._tip_length = 0 if np.isclose(val, self.length, atol=1e-12, rtol=0) else val

    @property
    def taper_angle(self):
        """Taper angle of the jaw front and back faces [rad]; 0 if not tapered.

        The element length is the full length of the jaw. With a taper, each jaw is
        tip_length long at the beam and grows with depth, up to length where the taper
        ends (taper_depth, or jaw_depth if not set), keeping that length beyond. Hence
        length = tip_length + 2*depth/tan(taper_angle).
        """
        if self._tip_length <= 0:
            return 0.
        depth = _taper_depth(self.jaw_depth, self.taper_depth)
        if depth <= 0 or self.length <= self._tip_length:
            return np.nan  # Invalid geometry, rejected when starting the engine
        return np.arctan(2*depth/(self.length - self._tip_length))

    @taper_angle.setter
    def taper_angle(self, val):
        # The element length (the footprint of the collimator in the line) is kept,
        # and tip_length is adapted
        if val == 0:
            self._tip_length = 0
            return
        self._tip_length = _tip_length_from_angle(self.length, val,
                                    _taper_depth(self.jaw_depth, self.taper_depth))


def _taper_depth(jaw_depth, taper_depth):
    # Depth from the jaw edge over which a jaw is tapered
    return taper_depth if taper_depth > 0 else jaw_depth


def _taper_length(taper_angle, depth):
    # How much longer than its face at the beam a tapered jaw is
    if not 0 < taper_angle < np.pi/2:
        raise ValueError("The taper angle must be in (0, 90) degrees!")
    if depth <= 0:
        raise ValueError("A taper needs jaw_depth or taper_depth!")
    return 2*depth/np.tan(taper_angle)


def _tip_length_from_angle(length, taper_angle, depth):
    tip_length = length - _taper_length(taper_angle, depth)
    if tip_length <= 0:
        raise ValueError(f"The length ({length} m) is too short for the taper angle and "
                       + f"depth ({depth} m)!")
    return tip_length


def _resolve_taper(kwargs):
    # Work out the length of the jaw face at the beam (0 if not tapered) from two of
    # length, tip_length and the taper angle, as BDSIM does, with
    # length = tip_length + 2*depth/tan(taper_angle). If length is not given, it is
    # calculated from the other two and set in kwargs.
    angles = {kk: kwargs.pop(kk) for kk in ['taper_angle', 'taper_angle_deg', 'taper_angle_rad']
              if kk in kwargs}
    if len(angles) > 1:
        raise ValueError("Use only one of `taper_angle`, `taper_angle_deg` or `taper_angle_rad`.")
    taper_angle = np.deg2rad(angles['taper_angle_deg']) if 'taper_angle_deg' in angles \
                  else next(iter(angles.values()), 0)
    tip_length = kwargs.pop('tip_length', 0)
    length = kwargs.get('length', 0)
    if taper_angle == 0:
        if tip_length == 0 or np.isclose(tip_length, length, atol=1e-12, rtol=0):
            return 0  # Flat (untapered) jaw
        if tip_length < 0 or tip_length > length:
            raise ValueError("tip_length must be positive and shorter than length!")
        return tip_length
    depth = _taper_depth(kwargs.get('jaw_depth', 0), kwargs.get('taper_depth', 0))
    if length > 0 and tip_length > 0:
        expected = tip_length + _taper_length(taper_angle, depth)
        if not np.isclose(expected, length, atol=1e-6, rtol=0):
            raise ValueError(f"length ({length} m) does not match tip_length, the taper angle "
                           + f"and depth, which give {expected} m: specify only two of "
                           + "length, tip_length and the taper angle.")
        return tip_length
    if length > 0:
        return _tip_length_from_angle(length, taper_angle, depth)
    if tip_length <= 0:
        raise ValueError("A taper angle needs length or tip_length!")
    kwargs['length'] = tip_length + _taper_length(taper_angle, depth)
    return tip_length


@track_construction
class Geant4Crystal(BaseCrystal):

    needs_cpu = True

    def __init__(self, **kwargs):
        import xcoll as xc
        if xc.geant4.engine.is_running():
            raise ValueError('Cannot create Geant4Crystal while engine is running.')
        raise NotImplementedError("Geant4Crystal not yet implemented.")
