# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

from contextlib import contextmanager

import xobjects as xo
import xtrack as xt

from .base import BaseCollimator
from ..scattering_routines.fluka import track, FlukaEngine


class FlukaCollimator(BaseCollimator):
    _xofields = { **BaseCollimator._xofields,
        'fluka_id':            xo.Int16,    # Do not change! Should be 16 bit because of FlukaIO type
        'accumulated_energy':  xo.Float64,
        'length_front':        xo.Float64,
        'length_back':         xo.Float64,
        '_tracking':           xo.Int8
    }

    isthick = True
    allow_track = True
    iscollective = True
    behaves_like_drift = True
    skip_in_loss_location_refinement = True

    _skip_in_to_dict       = BaseCollimator._skip_in_to_dict
    _store_in_to_dict      = BaseCollimator._store_in_to_dict
    _internal_record_class = BaseCollimator._internal_record_class

    _depends_on = [BaseCollimator, FlukaEngine]

    _allowed_fields_when_frozen = ['_tracking']

    def __new__(cls, *args, **kwargs):
        with cls._in_constructor():
            self = super().__new__(cls, *args, **kwargs)
        return self

    def __init__(self, **kwargs):
        with self.__class__._in_constructor():
            if '_xobject' not in kwargs:
                kwargs.setdefault('_tracking', True)
            super().__init__(**kwargs)
            if not hasattr(self, '_equivalent_drift'):
                self._equivalent_drift = xt.Drift(length=self.length)

    def track(self, part):
        track(self, part)

    def __setattr__(self, name, value):
        if name not in self._allowed_fields_when_frozen and FlukaEngine.is_running():
            raise ValueError('Engine is running; FlukaCollimator is frozen.')
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
