# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

import xcoll as xc


DEFAULT_XCOLL_ELEMENTS = [
    xc.BlackAbsorber,
    xc.BlackCrystal,
    xc.TransparentCollimator,
    xc.TransparentCrystal,
    xc.EverestBlock,
    xc.EverestCollimator,
    xc.EverestCrystal,
    xc.BlowUp,
    # xc.ParticleStatsMonitor,
    xc.EmittanceMonitor,
]

XCOLL_NON_TRACKING_ELEMENTS = []

if xc.fluka.interface.ready:
    XCOLL_NON_TRACKING_ELEMENTS += [xc.FlukaCollimator, xc.FlukaCrystal]
if xc.geant4.interface.ready:
    XCOLL_NON_TRACKING_ELEMENTS += [xc.Geant4Collimator, xc.Geant4CollimatorTip]

EXTRA_XCOLL_ELEMENTS = []
