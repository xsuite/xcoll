from .base_collimator import BaseCollimator
from .absorber import BlackAbsorber
from .everest_collimator import EverestCollimator, EverestCrystal
from .fluka_collimator import FlukaCollimator

_all_collimator_types = {BlackAbsorber, EverestCollimator, EverestCrystal, FlukaCollimator}
