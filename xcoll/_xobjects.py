# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2026.                 #
# ######################################### #

from ..general import _pkg_root


# Declare paths and libraries to xobjecst for prebuilding kernels
def get_build_info():
    root = _pkg_root.parent
    return {
        "include": [root.as_posix()],
        "lib": [(root / "xcoll" / "lib").as_posix()],
        "libraries": ["FlukaIO"],
    }
