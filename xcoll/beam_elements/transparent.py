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

    allow_track      = True
    _depends_on      = [XcollGeometry, xt.RandomRutherford]
    _extra_c_sources = [
        xt._pkg_root / 'headers/checks.h',
        xt._pkg_root / 'headers/particle_states.h',
        xt._pkg_root / 'random/random_src/rutherford.h',
        _pkg_root / 'headers/particle_states.h',
        _pkg_root / 'headers/checks.h',
        _pkg_root / 'beam_elements/elements_src/transparent_collimator.h']

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

    allow_track      = True
    _depends_on      = [XcollGeometry, xt.RandomRutherford]
    _extra_c_sources = [
        xt._pkg_root / 'headers/checks.h',
        xt._pkg_root / 'headers/particle_states.h',
        xt._pkg_root / 'random/random_src/rutherford.h',
        _pkg_root / 'headers/checks.h',
        _pkg_root / 'headers/particle_states.h',
        _pkg_root / 'beam_elements/elements_src/transparent_crystal.h']

    def __init__(self, **kwargs):
        if '_xobject' not in kwargs:
            kwargs.setdefault('_tracking', True)
        super().__init__(**kwargs)
        if not isinstance(self._context, xo.ContextCpu):
            raise ValueError('TransparentCrystal is currently not supported on GPU.')
