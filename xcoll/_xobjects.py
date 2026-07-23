# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2026.                 #
# ######################################### #

import xcoll as xc


# Declare paths and libraries to xobjects
def get_build_info():
    libraries = set()
    library_dirs = set()
    if xc.fluka.interface.ready:
        libraries.add("FlukaIO")
        library_dirs.add(xc.fluka.interface.lib_dir)
    return {
        "include_dirs": [xc._pkg_root.parent],
        "libraries": list(libraries),
        "library_dirs": list(library_dirs),
    }
