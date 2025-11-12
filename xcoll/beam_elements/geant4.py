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
from ..materials import(Material, CrystalMaterial, RefMaterial, db as material_db,
                        _DEFAULT_MATERIAL, _DEFAULT_CRYSTALMATERIAL)


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
    _skip_in_to_dict       = [*BaseCollimator._skip_in_to_dict, '_material']
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
        if material is None:
            material = _DEFAULT_MATERIAL
        elif isinstance(material, dict):
            material = Material.from_dict(material)
        elif isinstance(material, str):
            material = material_db[material]
        elif isinstance(material, RefMaterial):
            if material.geant4_name is None:
                raise ValueError(f"RefMaterial {material} does not have a Geant4 name!")
        if not isinstance(material, Material) \
        or isinstance(material, CrystalMaterial):
            raise ValueError(f"Invalid material of type {type(material)}!")
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


class Geant4Crystal(BaseCrystal):
    def __init__(self, **kwargs):
        raise NotImplementedError("Geant4Crystal not yet implemented.")
