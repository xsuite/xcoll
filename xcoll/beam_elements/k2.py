# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

import numpy as np

import xobjects as xo
import xpart as xp
import xtrack as xt

from .base import BaseCollimator, BaseCrystal, InvalidXcoll
from ..scattering_routines.k2 import K2Engine, track
from ..scattering_routines.everest.materials import _SixTrack_to_xcoll, SixTrack_from_xcoll, \
                                            SixTrack_from_xcoll_crystal, Material, CrystalMaterial
from ..general import _pkg_root


class _K2Collimator(BaseCollimator):
    _xofields = {**BaseCollimator._xofields,
        '_k2_id':    xo.Int32,
        '_material': xo.String,
        '_tracking': xo.Int8
    }

    isthick = True
    needs_rng = True
    allow_track = True
    iscollective = True
    behaves_like_drift = True
    skip_in_loss_location_refinement = True

    _skip_in_to_dict       = [*BaseCollimator._skip_in_to_dict, '_material']
    _store_in_to_dict      = [*BaseCollimator._store_in_to_dict, 'material']
    _internal_record_class = BaseCollimator._internal_record_class

    _depends_on = [BaseCollimator, K2Engine]


    def __init__(self, **kwargs):
        to_assign = {}
        if '_xobject' not in kwargs:
            to_assign['material'] = kwargs.pop('material', None)
            kwargs['_material'] = ''.ljust(4)
            kwargs.setdefault('_tracking', True)
            kwargs['_k2_id'] = -1
        super().__init__(**kwargs)
        for key, val in to_assign.items():
            setattr(self, key, val)
        if not hasattr(self, '_equivalent_drift'):
            self._equivalent_drift = xt.Drift(length=self.length)

    @property
    def material(self):
        return self._material

    @material.setter
    def material(self, val):
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
            elif val.lower() not in [mat.lower() for mat in _SixTrack_to_xcoll.keys()]:
                raise ValueError(f"Unknown material {val} (should be a SixTrack material code)!")
            correct_material_capitalisation = {mat.lower(): mat for mat in _SixTrack_to_xcoll.keys()}
            self._material = correct_material_capitalisation[val.lower()]


    def track(self, part):
        if self._tracking:
            track(self, part)


    def __setattribute__(self, name, value):
        if K2Engine.is_running():
            raise ValueError('K2Engine is running; _K2Collimator is frozen.')
        super().__setattribute__(name, value)


    def get_backtrack_element(self, _context=None, _buffer=None, _offset=None):
        return InvalidXcoll(length=-self.length, _context=_context,
                            _buffer=_buffer, _offset=_offset)


class _K2Crystal(BaseCrystal):
    _xofields = {**BaseCrystal._xofields,
        '_k2_id':    xo.Int8,
        'miscut':    xo.Float64,
        '_orient':   xo.Int8,
        '_material': xo.String,
        '_tracking': xo.Int8
    }

    isthick = True
    needs_rng = True
    allow_track = True
    iscollective = True
    behaves_like_drift = True
    skip_in_loss_location_refinement = True

    _skip_in_to_dict       = [*BaseCrystal._skip_in_to_dict, '_orient', '_material']
    _store_in_to_dict      = [*BaseCrystal._store_in_to_dict, 'lattice', 'material']
    _internal_record_class = BaseCrystal._internal_record_class

    _depends_on = [BaseCrystal, K2Engine]


    def __init__(self, **kwargs):
        to_assign = {}
        if '_xobject' not in kwargs:
            to_assign['material'] = kwargs.pop('material', None)
            kwargs['_material'] = ''.ljust(4)
            to_assign['lattice'] = kwargs.pop('lattice', 'strip')
            kwargs.setdefault('miscut', 0)
            kwargs.setdefault('_tracking', True)
            kwargs['_k2_id'] = -1
        super().__init__(**kwargs)
        for key, val in to_assign.items():
            setattr(self, key, val)
        if not hasattr(self, '_equivalent_drift'):
            self._equivalent_drift = xt.Drift(length=self.length)

    @property
    def lattice(self):
        if self._orient == 1:
            return 'strip'
        elif self._orient == 2:
            return 'quasi-mosaic'
        else:
            raise ValueError(f"Illegal value {self._orient} for '_orient'!")

    @lattice.setter
    def lattice(self, lattice):
        if lattice == 'strip' or lattice == '110' or lattice == 110:
            self._orient = 1
        elif lattice == 'quasi-mosaic' or lattice == '111' or lattice == 111:
            self._orient = 2
        else:
            raise ValueError(f"Illegal value {lattice} for 'lattice'! "
                            + "Only use 'strip' (110) or 'quasi-mosaic' (111).")

    @property
    def material(self):
        return self._material

    @material.setter
    def material(self, val):
        if isinstance(val, CrystalMaterial):
            self._material = SixTrack_from_xcoll_crystal(val)
        elif not isinstance(val, str):
            raise ValueError("Material should be an Everest `CrystalMaterial` or a string.")
        else:
            val = val.strip()
            if val.lower() == "va":
                raise ValueError("SixTrack material 'VA' not supported. Use a drift.")
            elif val.lower() == "bl":
                raise ValueError("SixTrack material 'BL' not supported. Use a BlackAbsorber.")
            elif val.lower() not in [mat.lower() for mat, vv in _SixTrack_to_xcoll.items() if len(vv) > 1]:
                raise ValueError(f"Unknown crystal material {val} (should be a SixTrack material code)!")
            correct_material_capitalisation = {mat.lower(): mat for mat in _SixTrack_to_xcoll.keys()}
            self._material = correct_material_capitalisation[val.lower()]


    def track(self, part):
        if self._tracking:
            track(self, part)


    def __setattribute__(self, name, value):
        if K2Engine.is_running():
            raise ValueError('K2Engine is running; _K2Crystal is frozen.')
        super().__setattribute__(name, value)


    def get_backtrack_element(self, _context=None, _buffer=None, _offset=None):
        return InvalidXcoll(length=-self.length, _context=_context,
                            _buffer=_buffer, _offset=_offset)

