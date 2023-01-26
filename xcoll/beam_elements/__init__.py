from .base_collimator import BaseCollimator
from .absorber import BlackAbsorber
from .everest_collimator import EverestCollimator, EverestCrystal
from .pyeverest_collimator import PyEverestCollimator, PyEverestCrystal
from .k2_collimator import K2Collimator, K2Crystal

_all_collimator_types = { BlackAbsorber, EverestCollimator, EverestCrystal, PyEverestCollimator, PyEverestCrystal, K2Collimator, K2Crystal }
