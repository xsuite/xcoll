# copyright ############################### #
# This file is part of the Xcoll package.   #
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
        _pkg_root.joinpath('geometry','old_geometry','sort.h'),
        _pkg_root.joinpath('geometry','old_geometry','segments.h'),
        _pkg_root.joinpath('geometry','old_geometry','objects.h'),
        _pkg_root.joinpath('geometry','old_geometry','methods.h'),
        _pkg_root.joinpath('geometry','old_geometry','get_s.h'),
        _pkg_root.joinpath('geometry','old_geometry','rotation.h'),
        _pkg_root.joinpath('geometry','old_geometry','collimator_geometry.h'),
        _pkg_root.joinpath('geometry','old_geometry','crystal_geometry.h')
    ]
