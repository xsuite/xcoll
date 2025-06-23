# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

from .general import _pkg_root, __version__, citation

from .beam_elements import BlackAbsorber, BlackCrystal, TransparentCollimator, TransparentCrystal, \
                           EverestBlock, EverestCollimator, EverestCrystal, Geant4Collimator, BlowUp, \
                           EmittanceMonitor, collimator_classes, crystal_classes, element_classes
from .scattering_routines.everest import materials, Material, CrystalMaterial
from .scattering_routines.geant4 import Geant4Engine
from .colldb import CollimatorDatabase
from .interaction_record import InteractionRecord
from .rf_sweep import RFSweep
from .lossmap import LossMap, MultiLossMap


# print("If you use Xcoll in your simulations, please cite us :-)")
# print(citation)

