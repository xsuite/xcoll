# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2023.                 #
# ######################################### #

import xobjects as xo
import xtrack as xt

from .materials import Material, CrystalMaterial
from ...interaction_record import InteractionRecord
from ...general import _pkg_root


class EverestEngine(xo.HybridClass):
    _xofields = {}

    _depends_on = [Material, CrystalMaterial, InteractionRecord, xt.RandomUniform, xt.RandomExponential,
                   xt.RandomNormal, xt.RandomRutherford, xt.Drift]

    _extra_c_sources = [
        _pkg_root.joinpath('scattering_routines','everest','constants.h'),
        _pkg_root.joinpath('scattering_routines','everest','everest.h'),
        _pkg_root.joinpath('scattering_routines','everest','properties.h'),
#         _pkg_root.joinpath('scattering_routines','everest','scatter_init.h'),
        _pkg_root.joinpath('scattering_routines','everest','multiple_coulomb_scattering.h'),
        _pkg_root.joinpath('scattering_routines','everest','crystal_parameters.h'),
        _pkg_root.joinpath('scattering_routines','everest','amorphous.h'),
        _pkg_root.joinpath('scattering_routines','everest','jaw.h'),
        # _pkg_root.joinpath('scattering_routines','everest','scatter.h'),
        _pkg_root.joinpath('scattering_routines','everest','channeling.h'),
        _pkg_root.joinpath('scattering_routines','everest','crystal.h'),
        _pkg_root.joinpath('scattering_routines','everest','scatter_crystal.h')
    ]
