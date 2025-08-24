# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

import xtrack as xt

from .materials import Material, CrystalMaterial
from ..geometry import XcollGeometry
from ...interaction_record import InteractionRecord
from ...general import _pkg_root


class EverestEngine(xt.BeamElement):
    _xofields = {}

    allow_track = False

    _depends_on = [Material, CrystalMaterial, InteractionRecord, xt.RandomUniform, xt.RandomExponential,
                   xt.RandomNormal, xt.RandomRutherford, xt.Drift, XcollGeometry]

    _extra_c_sources = [
        '#include <xcoll/scattering_routines/everest/amorphous.h>',
        '#include <xcoll/scattering_routines/everest/jaw.h>',
        '#include <xcoll/scattering_routines/everest/channelling.h>'
    ]
