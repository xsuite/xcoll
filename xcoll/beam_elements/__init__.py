# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

from .base import BaseBlock, BaseCollimator, BaseCrystal
from .absorber import BlackAbsorber, BlackCrystal
from .everest import EverestBlock, EverestCollimator, EverestCrystal
from .blowup import BlowUp
from .monitor import EmittanceMonitor

block_classes = tuple(v for v in globals().values()
                      if isinstance(v, type) and issubclass(v, BaseBlock) and v != BaseBlock
                      and v != BaseCollimator and v != BaseCrystal)
# Includes crystals
collimator_classes = tuple(v for v in globals().values()
                           if isinstance(v, type) and (issubclass(v, BaseCollimator) or issubclass(v, BaseCrystal))
                           and v != BaseCollimator and v != BaseCrystal)
crystal_classes = tuple(v for v in globals().values()
                        if isinstance(v, type) and issubclass(v, BaseCrystal) and v != BaseCrystal)

element_classes = block_classes + (BlowUp, EmittanceMonitor)

# For compatibility with _K2Collimators
_all_block_classes      = block_classes
_all_collimator_classes = collimator_classes
_all_crystal_classes    = crystal_classes
