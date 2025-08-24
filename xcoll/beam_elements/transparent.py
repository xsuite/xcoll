# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

import xtrack as xt
import xobjects as xo
from .base import BaseCollimator, BaseCrystal
from ..scattering_routines.geometry import XcollGeometry
from ..general import _pkg_root


class TransparentCollimator(BaseCollimator):
    _xofields = BaseCollimator._xofields | {
        '_tracking':        xo.Int8
    }

    isthick = True
    needs_rng = False
    allow_track = True
    allow_double_sided = True
    behaves_like_drift = True
    allow_rot_and_shift = False
    allow_loss_refinement = True
    skip_in_loss_location_refinement = True

    _noexpr_fields         = BaseCollimator._noexpr_fields
    _skip_in_to_dict       = BaseCollimator._skip_in_to_dict
    _store_in_to_dict      = BaseCollimator._store_in_to_dict
    _internal_record_class = BaseCollimator._internal_record_class

    _depends_on = [BaseCollimator, XcollGeometry]

    _extra_c_sources = [
        _pkg_root.joinpath('beam_elements','elements_src','transparent_collimator.h')
    ]

    def __init__(self, **kwargs):
        if '_xobject' not in kwargs:
            kwargs.setdefault('_tracking', True)
        super().__init__(**kwargs)
        if not isinstance(self._context, xo.ContextCpu):
            raise ValueError('TransparentCollimator is currently not supported on GPU.')


class TransparentCrystal(BaseCrystal):
    _xofields = BaseCrystal._xofields | {
        '_tracking':        xo.Int8
    }

    isthick = True
    needs_rng = False
    allow_track = True
    allow_double_sided = False
    behaves_like_drift = True
    allow_rot_and_shift = False
    allow_loss_refinement = True
    skip_in_loss_location_refinement = True

    _noexpr_fields         = BaseCrystal._noexpr_fields
    _skip_in_to_dict       = BaseCrystal._skip_in_to_dict
    _store_in_to_dict      = BaseCrystal._store_in_to_dict
    _internal_record_class = BaseCrystal._internal_record_class

    _depends_on = [BaseCrystal, XcollGeometry]

    _extra_c_sources = [
        _pkg_root.joinpath('beam_elements','elements_src','transparent_crystal.h')
    ]

    def __init__(self, **kwargs):
        if '_xobject' not in kwargs:
            kwargs.setdefault('_tracking', True)
        super().__init__(**kwargs)
        if not isinstance(self._context, xo.ContextCpu):
            raise ValueError('TransparentCrystal is currently not supported on GPU.')
