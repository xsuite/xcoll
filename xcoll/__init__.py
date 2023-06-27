from .general import _pkg_root

from .beam_elements import BaseCollimator, BlackAbsorber, EverestCollimator, EverestCrystal
from .scattering_routines.everest import materials, Material, CrystalMaterial
from .manager import CollimatorManager
from .colldb import CollimatorDatabase, load_SixTrack_colldb

__version__ = '0.2.2'
