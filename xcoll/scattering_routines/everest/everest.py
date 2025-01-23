# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

import xtrack as xt

from .materials import Material, CrystalMaterial
from ...geometry.old_geometry import XcollGeometry
from ...interaction_record import InteractionRecord
from ...general import _pkg_root


class EverestEngine(xt.BeamElement):
    _xofields = {}

    allow_track = False

    _depends_on = [Material, CrystalMaterial, InteractionRecord, xt.RandomUniform, xt.RandomExponential,
                   xt.RandomNormal, xt.RandomRutherford, xt.Drift, XcollGeometry]

    _extra_c_sources = [
        _pkg_root.joinpath('scattering_routines','geometry','rotation.h'),
        _pkg_root.joinpath('scattering_routines','everest','constants.h'),
        _pkg_root.joinpath('scattering_routines','everest','everest.h'),
        _pkg_root.joinpath('scattering_routines','everest','properties.h'),
        _pkg_root.joinpath('scattering_routines','everest','multiple_coulomb_scattering.h'),
        _pkg_root.joinpath('scattering_routines','everest','nuclear_interaction.h'),
        _pkg_root.joinpath('scattering_routines','everest','crystal_parameters.h'),
        _pkg_root.joinpath('scattering_routines','everest','amorphous.h'),
        _pkg_root.joinpath('scattering_routines','everest','jaw.h'),
        _pkg_root.joinpath('scattering_routines','everest','channeling.h')
    ]
