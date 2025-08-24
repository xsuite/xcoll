# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

import xtrack as xt
import xobjects as xo
from .base import BaseCollimator, BaseCrystal
from ..scattering_routines.geometry import XcollGeometry
from ..general import _pkg_root


class BlackAbsorber(BaseCollimator):
    _xofields = BaseCollimator._xofields | {
        '_tracking':        xo.Int8
    }

    allow_track      = True
    _depends_on      = [XcollGeometry]
    _extra_c_sources = ['#include <xcoll/beam_elements/elements_src/black_absorber.h>']

    def __init__(self, **kwargs):
        if '_xobject' not in kwargs:
            kwargs.setdefault('_tracking', True)
        super().__init__(**kwargs)
        if not isinstance(self._context, xo.ContextCpu):
            raise ValueError('BlackAbsorber is currently not supported on GPU.')


class BlackCrystal(BaseCrystal):
    _xofields = BaseCrystal._xofields | {
        '_tracking':        xo.Int8
    }

    allow_track      = True
    _depends_on      = [XcollGeometry]
    _extra_c_sources = ['#include <xcoll/beam_elements/elements_src/black_crystal.h>']

    def __init__(self, **kwargs):
        if '_xobject' not in kwargs:
            kwargs.setdefault('_tracking', True)
        super().__init__(**kwargs)
        if not isinstance(self._context, xo.ContextCpu):
            raise ValueError('BlackCrystal is currently not supported on GPU.')
