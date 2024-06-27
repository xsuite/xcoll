# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

import xobjects as xo
import xpart as xp
import xtrack as xt

from .base import BaseCollimator
from ..scattering_routines.fluka import track, FlukaEngine
from ..scattering_routines.everest.materials import SixTrack_to_xcoll


class FlukaCollimator(BaseCollimator):
    _xofields = { **BaseCollimator._xofields,
        'fluka_id':            xo.Int16,    # Do not change! Should be 16 bit because of FlukaIO type
        'accumulated_energy':  xo.Float64,
        'length_front':        xo.Float64,
        'length_back':         xo.Float64,
        '_frozen':             xo.Int8,
        '_tracking':           xo.Int8,
        '_material':           xo.String
    }

    isthick = True
    allow_track = True
    iscollective = True
    behaves_like_drift = True
    skip_in_loss_location_refinement = True

    _skip_in_to_dict       = [*BaseCollimator._skip_in_to_dict, '_material']
    _store_in_to_dict      = [*BaseCollimator._store_in_to_dict, 'material']
    _internal_record_class = BaseCollimator._internal_record_class

    _depends_on = [BaseCollimator, FlukaEngine]

    def __init__(self, **kwargs):
        to_assign = {}
        if '_xobject' not in kwargs:
            to_assign['material'] = kwargs.pop('material', None)
            kwargs['_material'] = 'NO NAME'.ljust(55)  # Pre-allocate 64 byte using whitespace
            kwargs.setdefault('_frozen', False)
            kwargs.setdefault('_tracking', True)
        super().__init__(**kwargs)
        for key, val in to_assign.items():
            setattr(self, key, val)
        if not hasattr(self, '_equivalent_drift'):
            self._equivalent_drift = xt.Drift(length=self.length)

    def track(self, part):
        track(self, part)


    @property
    def material(self):
        return self._material.strip()

    @material.setter
    def material(self, val):
        if FlukaEngine.is_running():
            raise ValueError('Engine is running; FlukaCollimator is frozen.')
        if val is None:
            self._material = 'NO NAME'.ljust(55)
            return
        if not val in SixTrack_to_xcoll:
            raise ValueError(f'Unknown material: {val}')
        self._material = val


    def __setattribute__(self, name, value):
        if FlukaEngine.is_running():
            raise ValueError('Engine is running; FlukaCollimator is frozen.')
        super().__setattribute__(name, value)
