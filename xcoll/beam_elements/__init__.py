# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

from .base import BaseBlock, BaseCollimator, BaseCrystal
from .absorber import BlackAbsorber, BlackCrystal
from .transparent import TransparentCollimator, TransparentCrystal
from .everest import EverestBlock, EverestCollimator, EverestCrystal
from .blowup import BlowUp
from .monitor import EmittanceMonitor
from .channelling import ChannellingDev
from .channelling import   BentChannellingDev, \
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

block_classes = tuple(v for v in globals().values()
                      if isinstance(v, type) and issubclass(v, BaseBlock) and v != BaseBlock
                      and v != BaseCollimator and v != BaseCrystal)
# Includes crystals
collimator_classes = tuple(v for v in globals().values()
                           if isinstance(v, type) and (issubclass(v, BaseCollimator) or issubclass(v, BaseCrystal))
                           and v != BaseCollimator and v != BaseCrystal)
crystal_classes = tuple(v for v in globals().values()
                        if isinstance(v, type) and issubclass(v, BaseCrystal) and v != BaseCrystal)

element_classes = block_classes + (BlowUp, EmittanceMonitor, ChannellingDev, BentChannellingDev,  BentChannellingDevM2V1o02, BentChannellingDevM2V1o04, BentChannellingDevM2V1o06,  BentChannellingDevM2V1o08, BentChannellingDevM2V1o10, BentChannellingDevM2V1o12, BentChannellingDevM2V2o02, BentChannellingDevM2V2o04, BentChannellingDevM2V2o06, BentChannellingDevM2V2o08, BentChannellingDevM2V2o10, BentChannellingDevM2V2o12, BentChannellingDevM3V1o02, BentChannellingDevM3V1o04, BentChannellingDevM3V1o06, BentChannellingDevM3V1o08, BentChannellingDevM3V1o10, BentChannellingDevM3V1o12,BentChannellingDevM3V2o02, BentChannellingDevM3V2o04, BentChannellingDevM3V2o06, BentChannellingDevM3V2o08, BentChannellingDevM3V2o10, BentChannellingDevM3V2o12, BentChannellingDevM4V1o02, BentChannellingDevM4V1o04, BentChannellingDevM4V1o06, BentChannellingDevM4V1o08, BentChannellingDevM4V1o10, BentChannellingDevM4V1o12, BentChannellingDevM4V2o02, BentChannellingDevM4V2o04, BentChannellingDevM4V2o06, BentChannellingDevM4V2o08, BentChannellingDevM4V2o10, BentChannellingDevM4V2o12)
