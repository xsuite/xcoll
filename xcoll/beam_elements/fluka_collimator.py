# copyright ############################### #
# This file is part of the Xcoll Package.  #
# Copyright (c) CERN, 2023.                 #
# ######################################### #

import xobjects as xo
import xpart as xp
import xtrack as xt

from .base_collimator import BaseCollimator, InvalidCollimator
from ..scattering_routines.fluka import track, FlukaEngine
from ..general import _pkg_root


class FlukaCollimator(BaseCollimator):
    _xofields = { **BaseCollimator._xofields,
        'fluka_id':  xo.Int64,
        '_tracking': xo.Int8
    }

    isthick = True
    iscollective = True
    behaves_like_drift = True
    skip_in_loss_location_refinement = True

    _skip_in_to_dict       = [*BaseCollimator._skip_in_to_dict]
    _store_in_to_dict      = [*BaseCollimator._store_in_to_dict]
    _internal_record_class = BaseCollimator._internal_record_class

    _depends_on = [BaseCollimator, FlukaEngine]

    _extra_c_sources = [
        _pkg_root.joinpath('beam_elements','collimators_src','everest_collimator.h')
    ]


    def __init__(self, **kwargs):
        if '_xobject' not in kwargs:
            kwargs.setdefault('_tracking', True)
        super().__init__(**kwargs)

    def track(self, part):
        track(self, part)
