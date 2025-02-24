# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

import xtrack as xt
import xobjects as xo
from .base import BaseCollimator, BaseCrystal, InvalidXcoll
from ..scattering_routines.geometry import XcollGeometry
from ..general import _pkg_root


class BlackAbsorber(BaseCollimator):
    _xofields = { **BaseCollimator._xofields,
        '_tracking':        xo.Int8
    }

    isthick = True
    allow_track = True
    behaves_like_drift = True
    allow_rot_and_shift = False
    skip_in_loss_location_refinement = True

    _depends_on = [BaseCollimator, XcollGeometry]

    _extra_c_sources = [
        _pkg_root.joinpath('beam_elements','elements_src','black_absorber.h')
    ]

    _skip_in_to_dict       = BaseCollimator._skip_in_to_dict
    _store_in_to_dict      = BaseCollimator._store_in_to_dict
    _internal_record_class = BaseCollimator._internal_record_class

    def __init__(self, **kwargs):
        if '_xobject' not in kwargs:
            kwargs.setdefault('_tracking', True)
        super().__init__(**kwargs)
        if not isinstance(self._context, xo.ContextCpu):
            raise ValueError('BlackAbsorber is currently not supported on GPU.')

    def get_backtrack_element(self, _context=None, _buffer=None, _offset=None):
        return InvalidXcoll(length=-self.length, _context=_context, _buffer=_buffer, _offset=_offset)


class BlackCrystal(BaseCrystal):
    _xofields = { **BaseCrystal._xofields,
        '_tracking':        xo.Int8
    }

    isthick = True
    allow_track = True
    behaves_like_drift = True
    skip_in_loss_location_refinement = True

    _depends_on = [BaseCrystal, XcollGeometry]

    _extra_c_sources = [
        _pkg_root.joinpath('beam_elements','elements_src','black_crystal.h')
    ]

    _skip_in_to_dict       = BaseCrystal._skip_in_to_dict
    _store_in_to_dict      = BaseCrystal._store_in_to_dict
    _internal_record_class = BaseCrystal._internal_record_class

    def __init__(self, **kwargs):
        if '_xobject' not in kwargs:
            kwargs.setdefault('_tracking', True)
        super().__init__(**kwargs)
        if not isinstance(self._context, xo.ContextCpu):
            raise ValueError('BlackCrystal is currently not supported on GPU.')

    def get_backtrack_element(self, _context=None, _buffer=None, _offset=None):
        return InvalidXcoll(length=-self.length, _context=_context, _buffer=_buffer, _offset=_offset)

