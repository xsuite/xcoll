# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

from contextlib import contextmanager

import xobjects as xo
import xtrack as xt

from .base import BaseCollimator, BaseCrystal
from ..general import _pkg_root
from ..scattering_routines.geometry import XcollGeometry
from ..scattering_routines.geant4 import Geant4Engine
from ..scattering_routines.geant4 import track as track_in_python
from ..materials import Material, CrystalMaterial, RefMaterial, db as material_db


class Geant4Collimator(BaseCollimator):
    _xofields = BaseCollimator._xofields | {
        'geant4_id': xo.Int16,
        '_material': xo.String,
        '_tracking': xo.Int8
    }

    isthick = True
    allow_track = True
    iscollective = True
    behaves_like_drift = True
    skip_in_loss_location_refinement = True

    _depends_on = [BaseCollimator, XcollGeometry, Geant4Engine]

    _extra_c_sources = [
        _pkg_root.joinpath('beam_elements','elements_src','geant4_collimator.h')
    ]

    _skip_in_to_dict       = [*BaseCollimator._skip_in_to_dict, '_material']
    _store_in_to_dict      = [*BaseCollimator._store_in_to_dict, 'material']
    _internal_record_class = BaseCollimator._internal_record_class

    _allowed_fields_when_frozen = ['_tracking']

    def __new__(cls, *args, **kwargs):
        with cls._in_constructor():
            self = super().__new__(cls, *args, **kwargs)
        return self

    def __init__(self, **kwargs):
        with self.__class__._in_constructor():
            to_assign = {}
            if '_xobject' not in kwargs:
                kwargs.setdefault('_tracking', True)
                to_assign['material'] = kwargs.pop('material', None)
                kwargs['_material'] = ''.ljust(16)
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
        return self._material

    @material.setter
    def material(self, material):
        if material is None:
            material = Material()
        elif isinstance(material, dict):
            material = Material.from_dict(material)
        elif isinstance(material, str):
            material = material_db[material]
        if not isinstance(material, Material):
            raise ValueError("Invalid material!")
        if self.material != material:
            self._material = material


    def track(self, part):
        if self._record_interactions > 0:
            part2 = part.copy()
            super().track(part2)
        track_in_python(self, part)


    def __setattr__(self, name, value):
        import xcoll as xc
        if name not in self._allowed_fields_when_frozen \
        and xc.geant4.engine.is_running():
            raise ValueError('Engine is running; Geant4Collimator is frozen.')
        super().__setattr__(name, value)

    @classmethod
    @contextmanager
    def _in_constructor(cls):
        original_setattr = cls.__setattr__
        def new_setattr(self, *args, **kwargs):
            return super().__setattr__( *args, **kwargs)
        cls.__setattr__ = new_setattr
        try:
            yield
        finally:
            cls.__setattr__ = original_setattr


class Geant4Crystal(BaseCrystal):
    def __init__(self, **kwargs):
        raise NotImplementedError("Geant4Crystal not yet implemented.")
