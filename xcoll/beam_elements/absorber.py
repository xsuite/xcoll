# copyright ############################### #
# This file is part of the Xcoll Package.  #
# Copyright (c) CERN, 2023.                 #
# ######################################### #

import xtrack as xt
import xobjects as xo
from .base_collimator import BaseCollimator
from ..general import _pkg_root


class BlackAbsorber(BaseCollimator):
    _xofields = { **BaseCollimator._xofields,
        '_tracking':        xo.Int8
    }

    _extra_c_sources = [
        _pkg_root.joinpath('beam_elements','collimators_src','absorber.h')
    ]

    isthick = True
    behaves_like_drift = True
    skip_in_loss_location_refinement = True

    _depends_on = [BaseCollimator, xt.Drift, xt.SRotation, xt.XYShift]

    _skip_in_to_dict       = BaseCollimator._skip_in_to_dict
    _store_in_to_dict      = BaseCollimator._store_in_to_dict
    _internal_record_class = BaseCollimator._internal_record_class

    def __init__(self, **kwargs):
        if '_xobject' not in kwargs:
            kwargs.setdefault('_tracking', True)
        super().__init__(**kwargs)

    def get_backtrack_element(self, _context=None, _buffer=None, _offset=None):
        return InvalidCollimator(length=-self.length, _context=_context, _buffer=_buffer, _offset=_offset)

