# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

try:
    from xaux import FsPath  # TODO: once xaux is in Xsuite keep only this
except (ImportError, ModuleNotFoundError):
    from .xaux import FsPath

_pkg_root = FsPath(__file__).parent.absolute()

citation = "F.F. Van der Veken, et al.: Recent Developments with the New Tools for Collimation Simulations in Xsuite, Proceedings of HB2023, Geneva, Switzerland."

# ======================
# Do not change
# ======================
__version__ = '0.11.6'
# ======================
