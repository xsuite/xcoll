from .base_collimator import BaseCollimator
from .absorber import BlackAbsorber
from .everest_collimator import EverestCollimator, EverestCrystal
from .pyeverest_collimator import PyEverestCollimator, PyEverestCrystal

_all_collimator_types = { BlackAbsorber, EverestCollimator, EverestCrystal, PyEverestCollimator, PyEverestCrystal }
