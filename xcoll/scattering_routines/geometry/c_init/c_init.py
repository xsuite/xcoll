# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

import xobjects as xo
from ....general import _pkg_root


XC_EPSILON = 1.e-15
XC_S_MAX = 1.e21


def xo_to_ctypes(args):
    if not hasattr(args, '__iter__') or isinstance(args, str):
        args = [args]
    return ", ".join([f"{arg.get_c_type()} {arg.name}" for arg in args])

def xo_to_cnames(args):
    if not hasattr(args, '__iter__') or isinstance(args, str):
        args = [args]
    return ", ".join([f"{arg.name}" for arg in args])


define_src = f"""
#ifndef XCOLL_GEOM_DEFINES_H
#define XCOLL_GEOM_DEFINES_H
#include <math.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>

#ifndef XC_EPSILON
#define XC_EPSILON {XC_EPSILON}
#endif

#ifndef XC_S_MAX
#define XC_S_MAX {XC_S_MAX}
#endif

#endif /* XCOLL_GEOM_DEFINES_H */
"""

class GeomCInit(xo.Struct):
    _extra_c_sources = [
        define_src,
        _pkg_root / 'scattering_routines' / 'geometry' / 'c_init' / 'sort.h',
        _pkg_root / 'scattering_routines' / 'geometry' / 'c_init' / 'methods.h',
        _pkg_root / 'scattering_routines' / 'geometry' / 'c_init' / 'find_root.h',
    ]

    # A Struct needs something to depend on, otherwise the class is added twice in the cdefs during compilation
    _depends_on = [xo.Float64]
