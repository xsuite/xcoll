# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

from .base import BaseBlock, BaseCollimator
from .absorber import BlackAbsorber
from .everest import EverestBlock, EverestCollimator, EverestCrystal

block_classes = tuple(v for v in globals().values()
                      if isinstance(v, type) and issubclass(v, BaseBlock) and v != BaseBlock)
collimator_classes = tuple(v for v in globals().values()
                           if isinstance(v, type) and issubclass(v, BaseCollimator) and v != BaseCollimator)
element_classes = block_classes + collimator_classes
