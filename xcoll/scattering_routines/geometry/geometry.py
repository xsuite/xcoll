# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

import xtrack as xt
from ...interaction_record import InteractionRecord
from ...general import _pkg_root

class XcollGeometry(xt.BeamElement):
    _xofields = {}

    allow_track = False

    _depends_on = [xt.Drift, xt.XYShift, xt.SRotation, xt.YRotation, InteractionRecord]

    _extra_c_sources = [
        _pkg_root.joinpath('scattering_routines','geometry','sort.h'),
        _pkg_root.joinpath('scattering_routines','geometry','segments.h'),
        _pkg_root.joinpath('scattering_routines','geometry','objects.h'),
        _pkg_root.joinpath('scattering_routines','geometry','methods.h'),
        _pkg_root.joinpath('scattering_routines','geometry','get_s.h'),
        _pkg_root.joinpath('scattering_routines','geometry','rotation.h'),
        _pkg_root.joinpath('scattering_routines','geometry','collimator_geometry.h'),
        _pkg_root.joinpath('scattering_routines','geometry','crystal_geometry.h')
    ]
