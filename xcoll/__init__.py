from .general import _pkg_root

from .beam_elements.absorber import BlackAbsorber
from .beam_elements.everest_collimator import EverestCollimator, EverestCrystal
from .manager import CollimatorManager
from .colldb import CollDB, load_SixTrack_colldb
from .scattering_routines.everest import materials
from .scattering_routines.everest.materials import Material, CrystalMaterial

__version__ = '0.1.1'
