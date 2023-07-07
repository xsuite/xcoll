from .base_collimator import BaseCollimator
from .absorber import BlackAbsorber
from .everest_collimator import EverestCollimator, EverestCrystal
from .geant4_collimator import Geant4Collimator

_all_collimator_types = { BlackAbsorber, EverestCollimator, EverestCrystal, Geant4Collimator}
