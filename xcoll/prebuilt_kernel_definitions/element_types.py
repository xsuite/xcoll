# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

from ..beam_elements import *


DEFAULT_XCOLL_ELEMENTS = [
    BlackAbsorber,
    BlackCrystal,
    TransparentCollimator,
    TransparentCrystal,
    EverestBlock,
    EverestCollimator,
    EverestCrystal,
    FlukaCollimator,
    # FlukaCrystal,
    Geant4Collimator,
    # Geant4Crystal,
    BlowUp,
    # ParticleStatsMonitor,
    EmittanceMonitor
]
