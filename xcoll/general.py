# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

from .xaux import FsPath as _FsPath

_pkg_root = _FsPath(__file__).parent.absolute()

citation = "F.F. Van der Veken, et al.: Recent Developments with the New Tools for Collimation Simulations in Xsuite, Proceedings of HB2023, Geneva, Switzerland."

# ======================
# Do not change
# ======================
__version__ = '0.12.5'
# ======================
