# copyright ############################### #
# This file is part of the Xcoll Package.  #
# Copyright (c) CERN, 2023.                 #
# ######################################### #

import xtrack as xt
from .base_collimator import BaseCollimator
from ..general import _pkg_root


class BlackAbsorber(BaseCollimator):
    _extra_c_sources = [
        _pkg_root.joinpath('beam_elements/collimators_src/absorber.h')
    ]

    isthick = True
    behaves_like_drift = True
    skip_in_loss_location_refinement = True

    _depends_on = [xt.Drift, xt.SRotation, xt.XYShift]

    _skip_in_to_dict       = BaseCollimator._skip_in_to_dict
    _store_in_to_dict      = BaseCollimator._store_in_to_dict
    _internal_record_class = BaseCollimator._internal_record_class

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

