# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2023.                 #
# ######################################### #

import random
import string
from contextlib import contextmanager

import xobjects as xo
import xtrack as xt

from .base import BaseCollimator
from ..scattering_routines.geant4 import Geant4Engine, track
from ..scattering_routines.everest.materials import _SixTrack_to_xcoll, SixTrack_from_xcoll, \
                                            SixTrack_from_xcoll_crystal, Material, CrystalMaterial


def _new_id64(len=16):
    chars = string.ascii_letters + string.digits + '+/'
    return ''.join(random.choice(chars) for i in range(len))


class Geant4Collimator(BaseCollimator):
    _xofields = BaseCollimator._xofields | {
        'geant4_id': xo.String,
        '_material':  xo.String,
        '_tracking': xo.Int8
    }

    isthick = True
    allow_track = True
    iscollective = True
    behaves_like_drift = True
    skip_in_loss_location_refinement = True

    _skip_in_to_dict       = [*BaseCollimator._skip_in_to_dict, '_material']
    _store_in_to_dict      = [*BaseCollimator._store_in_to_dict, 'material']
    _internal_record_class = BaseCollimator._internal_record_class

    _depends_on = [BaseCollimator, Geant4Engine]

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
                kwargs['geant4_id'] = kwargs.get('geant4_id', _new_id64()).ljust(16)
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
    def material(self, val):
        # TODO: better material handling
        if isinstance(val, Material):
            self._material = SixTrack_from_xcoll(val)
        elif not isinstance(val, str):
            raise ValueError("Material should be an Everest `Material` or a string.")
        else:
            val = val.strip()
            if val.lower() == "va":
                raise ValueError("SixTrack material 'VA' not supported. Use a drift.")
            elif val.lower() == "bl":
                raise ValueError("SixTrack material 'BL' not supported. Use a BlackAbsorber.")
        geant4_materials = {
            'c':    'AC150GPH',
            'cu':   'Cu',
            'mogr': 'MG6403Fc',
            'mo': 'Mo',
            'cucd': 'CUDIAM75',
            'iner': 'INERM180'
        }
        if val.lower() not in geant4_materials:
            # TODO: need to check with BDSIM configuration file etc
            raise ValueError(f"Material {val} not yet supported.")
        self._material = geant4_materials[val.lower()]


    def track(self, part):
        track(self, part)


    def __setattr__(self, name, value):
        if name not in self._allowed_fields_when_frozen and Geant4Engine.is_running():
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
