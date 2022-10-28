from .general import _pkg_root

from .beam_elements.collimators import BlackAbsorber
from .beam_elements.k2collimator import K2Collimator, K2Crystal
from .manager import CollimatorManager
from .colldb import CollDB, load_SixTrack_colldb
from .scattering_routines.k2 import materials, K2Engine
from .scattering_routines.k2.materials import Material, CrystalMaterial

__version__ = '0.1.1'
