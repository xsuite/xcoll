# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

from .base import BaseBlock, BaseCollimator, BaseCrystal
from .absorber import BlackAbsorber, BlackCrystal
from .transparent import TransparentCollimator, TransparentCrystal
from .everest import EverestBlock, EverestCollimator, EverestCrystal
from .fluka import FlukaCollimator, FlukaCrystal
from .geant4 import Geant4Collimator, Geant4CollimatorTip, Geant4Crystal
from .blowup import BlowUp
from .monitor import ParticleStatsMonitor, EmittanceMonitor

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
