# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

import xobjects as xo
from ...general import _pkg_root


class XcollGeometry(xo.HybridClass):
    _xofields = {}

    _depends_on = [xo.Float64, xo.Int64] # Hack: need something to depend on, otherwise the class is added twice in the cdefs during compilation

    _extra_c_sources = [
        _pkg_root.joinpath('scattering_routines','geometry','sort.h'),
        _pkg_root.joinpath('scattering_routines','geometry','segments.h'),
        _pkg_root.joinpath('scattering_routines','geometry','objects.h'),
        _pkg_root.joinpath('scattering_routines','geometry','methods.h'),
        _pkg_root.joinpath('scattering_routines','geometry','get_s.h')
    ]
