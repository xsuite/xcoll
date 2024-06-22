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
from ..scattering_routines.everest.materials import SixTrack_to_xcoll
from ..general import _pkg_root


class _K2Collimator(BaseCollimator):
    _xofields = {**BaseCollimator._xofields,
        '_material':        xo.String,
        '_tracking':        xo.Int8
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
        super().__init__(**kwargs)
        for key, val in to_assign.items():
            setattr(self, key, val)


    @property
    def material(self):
        return self._material

    @material.setter
    def material(self, val):
        if not isinstance(val, str):
            raise ValueError("Material should be a string.")
        if val not in SixTrack_to_xcoll:
            raise ValueError(f"Unknown material {val} (should be a SixTrack material code)!")
        self._material = val.strip()


    @property
    def track(self, part):
        if self._tracking:
            track(self, part)


    def __setattribute__(self, name, value):
        if K2Engine._file_generated:
            # TODO: test if file exists or somethin
            raise ValueError('Input file generated; K2Collimator is frozen.')
        super().__setattribute__(name, value)


    def get_backtrack_element(self, _context=None, _buffer=None, _offset=None):
        return InvalidXcoll(length=-self.length, _context=_context,
                            _buffer=_buffer, _offset=_offset)


class _K2Crystal(BaseCrystal):
    _xofields = {**BaseCrystal._xofields,
        'miscut':             xo.Float64,
        '_orient':            xo.Int8,
        '_material':        xo.String,
        '_tracking':        xo.Int8
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
        super().__init__(**kwargs)
        for key, val in to_assign.items():
            setattr(self, key, val)

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


    def get_backtrack_element(self, _context=None, _buffer=None, _offset=None):
        return InvalidXcoll(length=-self.length, _context=_context,
                            _buffer=_buffer, _offset=_offset)

