# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

from .base import BaseBlock, BaseCollimator, BaseCrystal
from .absorber import BlackAbsorber, BlackCrystal
from .everest import EverestBlock, EverestCollimator, EverestCrystal
from .geant4 import Geant4Collimator
from .blowup import BlowUp

block_classes = tuple(v for v in globals().values()
                      if isinstance(v, type) and issubclass(v, BaseBlock) and v != BaseBlock)
# Includes crystals
collimator_classes = tuple(v for v in globals().values()
                           if isinstance(v, type) and (issubclass(v, BaseCollimator) or issubclass(v, BaseCrystal))
                           and v != BaseCollimator and v != BaseCrystal)
crystal_classes = tuple(v for v in globals().values()
                        if isinstance(v, type) and issubclass(v, BaseCrystal) and v != BaseCrystal)
element_classes = block_classes + collimator_classes
