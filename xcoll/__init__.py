# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

from .general import _pkg_root, __version__, citation

from .beam_elements import BlackAbsorber, BlackCrystal, TransparentCollimator, TransparentCrystal, \
                           EverestBlock, EverestCollimator, EverestCrystal, BlowUp, EmittanceMonitor, \
                           ChannellingDev, collimator_classes, crystal_classes, element_classes, \
                           BentChannellingDev, \
                           BentChannellingDevM2V1o02, BentChannellingDevM2V1o04, BentChannellingDevM2V1o06, \
                           BentChannellingDevM2V1o08, BentChannellingDevM2V1o10, BentChannellingDevM2V1o12, \
                           BentChannellingDevM2V2o02, BentChannellingDevM2V2o04, BentChannellingDevM2V2o06, \
                           BentChannellingDevM2V2o08, BentChannellingDevM2V2o10, BentChannellingDevM2V2o12, \
                           BentChannellingDevM3V1o02, BentChannellingDevM3V1o04, BentChannellingDevM3V1o06, \
                           BentChannellingDevM3V1o08, BentChannellingDevM3V1o10, BentChannellingDevM3V1o12, \
                           BentChannellingDevM3V2o02, BentChannellingDevM3V2o04, BentChannellingDevM3V2o06, \
                           BentChannellingDevM3V2o08, BentChannellingDevM3V2o10, BentChannellingDevM3V2o12, \
                           BentChannellingDevM4V1o02, BentChannellingDevM4V1o04, BentChannellingDevM4V1o06, \
                           BentChannellingDevM4V1o08, BentChannellingDevM4V1o10, BentChannellingDevM4V1o12, \
                           BentChannellingDevM4V2o02, BentChannellingDevM4V2o04, BentChannellingDevM4V2o06, \
                           BentChannellingDevM4V2o08, BentChannellingDevM4V2o10, BentChannellingDevM4V2o12

 
from .scattering_routines.everest import materials, Material, CrystalMaterial
from .colldb import CollimatorDatabase
from .interaction_record import InteractionRecord
from .rf_sweep import RFSweep
from .lossmap import LossMap, MultiLossMap


# print("If you use Xcoll in your simulations, please cite us :-)")
# print(citation)

