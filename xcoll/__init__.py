from .general import _pkg_root

from .beam_elements import BlackAbsorber, K2Collimator, K2Crystal, K2Engine
from .manager import CollimatorManager
from .colldb import CollDB, load_SixTrack_colldb

__version__ = '0.1.1'
