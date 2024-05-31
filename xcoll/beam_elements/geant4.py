# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2023.                 #
# ######################################### #

import numpy as np
import random
import string

import xobjects as xo
from xtrack import Drift

from .base import BaseCollimator
from ..scattering_routines.geant4 import Geant4Engine, track


def _new_id64(len=16):
    chars = string.ascii_letters + string.digits + '+/'
    return ''.join(random.choice(chars) for i in range(len))


class Geant4Collimator(BaseCollimator):
    _xofields = BaseCollimator._xofields | {
        '_id':       xo.String,
        'material':  xo.String,
        '_tracking': xo.Int8
    }

    isthick = True
    iscollective = True
    allow_track = True
    behaves_like_drift = True
    skip_in_loss_location_refinement = True

    _skip_in_to_dict       = [*BaseCollimator._skip_in_to_dict, '_id']
    _store_in_to_dict      = [*BaseCollimator._store_in_to_dict, 'collimator_id']
    _internal_record_class = BaseCollimator._internal_record_class

    _depends_on = [Drift, BaseCollimator, Geant4Engine]


    def __init__(self, **kwargs):
        if '_xobject' not in kwargs:
            kwargs.setdefault('collimator_id', _new_id64())
            kwargs['_id'] = kwargs.pop('collimator_id')
            kwargs.setdefault('_tracking', True)
            if kwargs.get('material') is None:
                raise ValueError("Need to provide a material to the collimator!")
        super().__init__(**kwargs)
        if '_xobject' not in kwargs:
            Geant4Engine().register_collimator(self)
            self.equivalent_drift = Drift(length=self.length)
        # TODO: is there a smarter way to hangle the drift tracking when inactive
        # TODO: should the inactive length before/after be handled here?

    # def __del__(self):
    #     # TODO: with a custom destructor, should all base class destructors be called explicitly?
    #     Geant4Engine().deregister_collimator(self)


    @property
    def collimator_id(self):
        return self._id


    def track(self, part):
        if self.active and self._tracking and Geant4Engine().connected:
            track(self.collimator_id, part)
        else:
            self.equivalent_drift.track(part)
