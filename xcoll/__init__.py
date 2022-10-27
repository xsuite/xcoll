from .general import _pkg_root

from .beam_elements import BlackAbsorber, K2Collimator, K2Crystal, K2Engine
from .manager import CollimatorManager
from .colldb import CollDB, load_SixTrack_colldb
from .beam_elements.k2 import materials
from .beam_elements.k2.materials import Material, CrystalMaterial

__version__ = '0.1.1'
