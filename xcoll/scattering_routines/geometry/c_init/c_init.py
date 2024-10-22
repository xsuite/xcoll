# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

import xobjects as xo
from ....general import _pkg_root

XC_EPSILON = 1e-12

class GeomCInit(xo.Struct):
    _extra_c_sources = [
        "#ifndef XC_EPSILON\n#define XC_EPSILON " + f"{XC_EPSILON}" + "\n#endif",
        _pkg_root / 'scattering_routines' / 'geometry' / 'c_init' / 'defines.h',
        _pkg_root / 'scattering_routines' / 'geometry' / 'c_init' / 'sort.h',
        _pkg_root / 'scattering_routines' / 'geometry' / 'c_init' / 'methods.h',
        _pkg_root / 'scattering_routines' / 'geometry' / 'c_init' / 'find_root.h',
    ]

    # A Struct needs something to depend on, otherwise the class is added twice in the cdefs during compilation
    _depends_on = [xo.Float64]

