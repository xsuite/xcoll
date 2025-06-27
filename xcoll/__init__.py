# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

from .general import _pkg_root, __version__, citation

from .beam_elements import BlackAbsorber, BlackCrystal, TransparentCollimator, TransparentCrystal, \
                           EverestBlock, EverestCollimator, EverestCrystal, BlowUp, EmittanceMonitor, \
                           FlukaCollimator, FlukaCrystal, collimator_classes, crystal_classes, element_classes
from .scattering_routines.everest import materials, Material, CrystalMaterial
from .scattering_routines.fluka import FlukaPrototype, FlukaAssembly, create_generic_assembly
from .colldb import CollimatorDatabase
from .interaction_record import InteractionRecord
from .rf_sweep import RFSweep
from .lossmap import LossMap, MultiLossMap
from .headers import particle_states

# Initialise FLUKA environment
from .scattering_routines.fluka.wrapper import FlukaWrapper as _FlukaWrapper
fluka = _FlukaWrapper()


# print("If you use Xcoll in your simulations, please cite us :-)")
# print(citation)
