from .base import BaseBlock, BaseCollimator
from .absorber import BlackAbsorber
from .everest import EverestBlock, EverestCollimator, EverestCrystal

_all_collimator_types = {BlackAbsorber, EverestCollimator, EverestCrystal}
