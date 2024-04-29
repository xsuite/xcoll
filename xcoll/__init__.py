# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

from .general import _pkg_root, __version__, citation

from .beam_elements import BlackAbsorber, EverestBlock, EverestCollimator, EverestCrystal, element_classes
from .install import install_elements
from .line_tools import assign_optics_to_collimators, open_collimators, send_to_parking, enable_scattering, disable_scattering
from .scattering_routines.everest import materials, Material, CrystalMaterial
from .manager import CollimatorManager
from .colldb import CollimatorDatabase, load_SixTrack_colldb
from .rf_sweep import RFSweep
from .initial_distribution import *
from .lossmap import LossMap

# print("If you use Xcoll in your simulations, please cite us :-)")
# print(citation)

