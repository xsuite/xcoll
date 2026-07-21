# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2026.                 #
# ######################################### #

from .general import _pkg_root


# Declare paths and libraries to xobjects
def get_build_info():
    return {
        "include_dirs": [_pkg_root.parent],
        "libraries": ["FlukaIO"],
        "library_dirs": [_pkg_root / "lib"],
    }
