from .general import _pkg_root

from .beam_elements import BlackAbsorber, K2Collimator, K2Crystal
from .scattering_routines.k2.engine import K2Engine
from .scattering_routines.k2 import materials
from .scattering_routines.k2.materials import Material, CrystalMaterial
from .manager import CollimatorManager
from .colldb import CollDB, load_SixTrack_colldb

__version__ = '0.1.2'
