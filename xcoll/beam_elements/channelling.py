# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

import xobjects as xo

from ..general import _pkg_root
from .base import InvalidXcoll


class ChannellingDev(InvalidXcoll):
    _xofields = {
        'length': xo.Float64,
    }

    isthick = True
    behaves_like_drift = False
    allow_track = True
    needs_rng = False
    skip_in_loss_location_refinement = True
    allow_loss_refinement = False

    _depends_on = [InvalidXcoll]

    _extra_c_sources = [
        _pkg_root.joinpath('beam_elements','elements_src','elliptic_functions.h'),
        _pkg_root.joinpath('beam_elements','elements_src','channelling.h')
    ]

