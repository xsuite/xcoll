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
    allow_track = True
    iscollective = True
    behaves_like_drift = True
    allow_rot_and_shift = False
    allow_loss_refinement = True
    skip_in_loss_location_refinement = True
    allow_no_prebuilt_kernel = True

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
        'taper_angle': xo.Float64,  # rad; 0 (default) means a flat (untapered) jaw, as before. With a
                                    # taper, length is that of the flat jaw face at the beam, and the
                                    # jaws are longer than the element (see taper_extension)
        'jaw_depth': xo.Float64,    # m, transverse depth of each jaw from its edge; 0 (default)
                                    # means a 2 m wide Geant4 box, as before
        'taper_depth': xo.Float64,  # m, depth from the jaw edge where the taper ends, the jaw keeping
                                    # a constant length beyond it; 0 (default) tapers the whole jaw
        'jaw_width': xo.Float64     # m, full width of each jaw across its plane of motion; 0 (default)
                                    # means BDSIM's default of 12 mm
    }

    isthick = True
    allow_track = True
    iscollective = True
    behaves_like_drift = True
    skip_in_loss_location_refinement = True
    allow_no_prebuilt_kernel = True

    _depends_on = [*Geant4Collimator._depends_on]

    _noexpr_fields         = {*Geant4Collimator._noexpr_fields, 'tip_material'}
    _skip_in_to_dict       = Geant4Collimator._skip_in_to_dict
    _store_in_to_dict      = [*Geant4Collimator._store_in_to_dict, 'tip_material']
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
            taper_angle_deg = kwargs.pop('taper_angle_deg', None)
            taper_angle_rad = kwargs.pop('taper_angle_rad', None)
            if taper_angle_deg is not None and taper_angle_rad is not None:
                raise ValueError("Use only one of `taper_angle_deg` or "
                                + "`taper_angle_rad`, not both.")
            elif taper_angle_deg is not None:
                kwargs['taper_angle'] = np.deg2rad(taper_angle_deg)
            elif taper_angle_rad is not None:
                kwargs['taper_angle'] = taper_angle_rad
            kwargs.setdefault('taper_angle', 0)
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
    def taper_extension(self):
        """How far a tapered jaw extends beyond each end of the element [m].

        The element length is that of the flat jaw face at the beam; from there
        each jaw grows with depth, up to length + 2*depth/tan(taper_angle) where
        the taper ends (taper_depth, or jaw_depth if not set or smaller), and
        keeps that length beyond. Geant4 tracks through that full length,
        centred on the element, so the lattice should leave this much room on
        either side.
        """
        if self.taper_angle <= 0:
            return 0.
        depth = self.jaw_depth
        if self.taper_depth > 0 and (depth <= 0 or self.taper_depth < depth):
            depth = self.taper_depth
        return depth / np.tan(self.taper_angle)


@track_construction
class Geant4Crystal(BaseCrystal):

    allow_no_prebuilt_kernel = True

    def __init__(self, **kwargs):
        import xcoll as xc
        if xc.geant4.engine.is_running():
            raise ValueError('Cannot create Geant4Crystal while engine is running.')
        raise NotImplementedError("Geant4Crystal not yet implemented.")
