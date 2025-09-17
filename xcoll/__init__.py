# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

from .general import _pkg_root, __version__, citation

from .beam_elements import BlackAbsorber, BlackCrystal, TransparentCollimator, TransparentCrystal, \
                           EverestBlock, EverestCollimator, EverestCrystal, Geant4Collimator, Geant4Crystal, \
                           BlowUp, EmittanceMonitor, collimator_classes, crystal_classes, element_classes
from .scattering_routines.everest import materials, Material, CrystalMaterial
from .scattering_routines.geant4 import Geant4Engine
from .colldb import CollimatorDatabase
from .interaction_record import InteractionRecord
from .rf_sweep import RFSweep, prepare_rf_sweep
from .lossmap import LossMap, MultiLossMap
from .headers import particle_states

# Initialise Geant4 environment
from .scattering_routines.geant4.wrapper import Geant4Wrapper as _Geant4Wrapper
geant4 = _Geant4Wrapper()


# print("If you use Xcoll in your simulations, please cite us :-)")
# print(citation)

