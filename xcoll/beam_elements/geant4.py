# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

from contextlib import contextmanager

import xobjects as xo
import xtrack as xt

from .base import BaseCollimator, BaseCrystal
from ..general import _pkg_root
from ..scattering_routines.geant4 import Geant4Engine, track_pre, track_core, track_post
from ..materials import _DEFAULT_MATERIAL, _resolve_material


class Geant4Collimator(BaseCollimator):
    _xofields = BaseCollimator._xofields | {
        'geant4_id': xo.String,
        '_tracking': xo.Int8,
        '_acc_ionisation_loss':  xo.Float64,  # TODO: this is not very robust, for when a track is done with new particles etc
    }

    isthick = True
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
    _skip_in_to_dict       = [*BaseCollimator._skip_in_to_dict, '_tracking', '_acc_ionisation_loss']
    _store_in_to_dict      = [*BaseCollimator._store_in_to_dict, 'material']
    _internal_record_class = BaseCollimator._internal_record_class

    _allowed_fields_when_frozen = ['_tracking', '_acc_ionisation_loss']

    def __new__(cls, *args, **kwargs):
        with cls._in_constructor():
            self = super().__new__(cls, *args, **kwargs)
        return self

    def __init__(self, **kwargs):
        import xcoll as xc
        with self.__class__._in_constructor(self):
            to_assign = {}
            if '_xobject' not in kwargs:
                kwargs.setdefault('_tracking', True)
                kwargs.setdefault('_acc_ionisation_loss', -1.)
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

    @property
    def angle(self):
        return BaseCollimator.angle.fget(self)

    @angle.setter
    def angle(self, val):
        if hasattr(val, '__iter__'):
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
        xc.geant4.environment.assert_environment_ready()
        if not xc.geant4.engine.is_running():
            raise RuntimeError("Geant4 engine is not running.")
        super().enable_scattering()

    def track(self, part):
        if track_pre(self, part):
            # super().track(part)
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

    def __setattr__(self, name, value):
        import xcoll as xc
        if name not in self._allowed_fields_when_frozen \
        and xc.geant4.engine.is_running():
            raise ValueError('Engine is running; Geant4Collimator is frozen.')
        super().__setattr__(name, value)

    @classmethod
    @contextmanager
    def _in_constructor(cls, self=None):
        original_setattr = cls.__setattr__
        if self is not None:
            self._being_constructed_ = True
        def new_setattr(self, *args, **kwargs):
            return super().__setattr__( *args, **kwargs)
        cls.__setattr__ = new_setattr
        try:
            yield
        finally:
            cls.__setattr__ = original_setattr
            if self is not None:
                self._being_constructed_ = False

    def _being_constructed(self):
        if hasattr(self, '_being_constructed_'):
            return self._being_constructed_
        else:
            return False


class Geant4CollimatorTip(Geant4Collimator):

    _xofields = Geant4Collimator._xofields | {
        'tip_thickness': xo.Float64
    }

    isthick = True
    allow_track = True
    iscollective = True
    behaves_like_drift = True
    skip_in_loss_location_refinement = True

    _extra_c_sources = [
        _pkg_root.joinpath('beam_elements', 'elements_src', 'geant4_collimator_tip.h')
    ]

    _depends_on = [*Geant4Collimator._depends_on]

    _noexpr_fields         = {*Geant4Collimator._noexpr_fields, 'tip_material'}
    _skip_in_to_dict       = [*Geant4Collimator._skip_in_to_dict, '_tracking', '_acc_ionisation_loss']
    _store_in_to_dict      = [*Geant4Collimator._store_in_to_dict, 'tip_material']
    _internal_record_class = Geant4Collimator._internal_record_class

    def __new__(cls, *args, **kwargs):
        with cls._in_constructor():
            self = super().__new__(cls, *args, **kwargs)
        return self

    def __init__(self, **kwargs):
        with self.__class__._in_constructor():
            to_assign = {}
            if '_xobject' not in kwargs:
                to_assign['tip_material'] = kwargs.pop('tip_material', None)
                kwargs['_tip_material'] = _DEFAULT_MATERIAL
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


class Geant4Crystal(BaseCrystal):
    def __init__(self, **kwargs):
        raise NotImplementedError("Geant4Crystal not yet implemented.")
