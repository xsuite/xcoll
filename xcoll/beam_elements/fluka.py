# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

import xobjects as xo
import xpart as xp
import xtrack as xt

from .base import BaseCollimator
from ..scattering_routines.fluka import track, FlukaEngine


class FlukaCollimator(BaseCollimator):
    _xofields = { **BaseCollimator._xofields,
        'fluka_id':            xo.Int16,    # Do not change! Should be 16 bit because of FlukaIO type
        'accumulated_energy':  xo.Float64,
        'length_front':        xo.Float64,
        'length_back':         xo.Float64,
        '_frozen':              xo.Int8,
        '_tracking':           xo.Int8
    }

    isthick = True
    allow_track = True
    iscollective = True
    behaves_like_drift = True
    skip_in_loss_location_refinement = True

    _skip_in_to_dict       = [*BaseCollimator._skip_in_to_dict]
    _store_in_to_dict      = [*BaseCollimator._store_in_to_dict]
    _internal_record_class = BaseCollimator._internal_record_class

    _depends_on = [BaseCollimator, FlukaEngine]

    def __init__(self, **kwargs):
        if '_xobject' not in kwargs:
            kwargs.setdefault('_frozen', False)
            kwargs.setdefault('_tracking', True)
        super().__init__(**kwargs)

    def track(self, part):
        track(self, part)

    def __setattribute__(self, name, value):
        # if name in ['gap', 'gap_L', 'gap_R', 'jaw', 'jaw_L', 'jaw_R', 'jaw_LU', 'jaw_LD', 'jaw_RU', 'jaw_RD', 'tilt', 'tilt_L', 'tilt_R']:
        if FlukaEngine.is_running():
            raise ValueError('Engine is running; FlukaCollimator is frozen.')
        super().__setattribute__(name, value)

    # @BaseCollimator.gap.setter
    # def gap(self, value):
    #     if FlukaEngine.is_running():
    #         raise ValueError('Fluka collimators do not support gap changes')

    # @BaseCollimator.gap_L.setter
    # def gap_L(self, value):
    #     if FlukaEngine.is_running():
    #         raise ValueError('Fluka collimators do not support gap changes')

    # @BaseCollimator.gap_R.setter
    # def gap_R(self, value):
    #     if FlukaEngine.is_running():
    #         raise ValueError('Fluka collimators do not support gap changes')

    # @BaseCollimator.jaw.setter
    # def jaw(self, value):
    #     if FlukaEngine.is_running():
    #         raise ValueError('Fluka collimators do not support jaw changes')

    # @BaseCollimator.jaw_L.setter
    # def jaw_L(self, value):
    #     if FlukaEngine.is_running():
    #         raise ValueError('Fluka collimators do not support jaw changes')

    # @BaseCollimator.jaw_R.setter
    # def jaw_R(self, value):
    #     if FlukaEngine.is_running():
    #         raise ValueError('Fluka collimators do not support jaw changes')

    # @BaseCollimator.jaw_LU.setter
    # def jaw_LU(self, value):
    #     if FlukaEngine.is_running():
    #         raise ValueError('Fluka collimators do not support jaw changes')

    # @BaseCollimator.jaw_LD.setter
    # def jaw_LD(self, value):
    #     if FlukaEngine.is_running():
    #         raise ValueError('Fluka collimators do not support jaw changes')

    # @BaseCollimator.jaw_RU.setter
    # def jaw_RU(self, value):
    #     if FlukaEngine.is_running():
    #         raise ValueError('Fluka collimators do not support jaw changes')

    # @BaseCollimator.jaw_RD.setter
    # def jaw_RD(self, value):
    #     if FlukaEngine.is_running():
    #         raise ValueError('Fluka collimators do not support jaw changes')

    # @BaseCollimator.tilt.setter
    # def tilt(self, value):
    #     if FlukaEngine.is_running():
    #         raise ValueError('Fluka collimators do not support tilt changes')

    # @BaseCollimator.tilt_L.setter
    # def tilt_L(self, value):
    #     if FlukaEngine.is_running():
    #         raise ValueError('Fluka collimators do not support tilt changes')

    # @BaseCollimator.tilt_R.setter
    # def tilt_R(self, value):
    #     if FlukaEngine.is_running():
    #         raise ValueError('Fluka collimators do not support tilt changes')
