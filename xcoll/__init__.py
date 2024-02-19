# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2023.                 #
# ######################################### #

from .general import _pkg_root, __version__, citation

from .beam_elements import BlackAbsorber, EverestBlock, EverestCollimator, EverestCrystal
from .scattering_routines.everest import materials, Material, CrystalMaterial
from .manager import CollimatorManager
from .colldb import CollimatorDatabase, load_SixTrack_colldb
from .rf_sweep import RFSweep
from .lossmap import LossMap

# print("If you use Xcoll in your simulations, please cite us :-)")
# print(citation)
