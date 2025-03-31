# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

import xobjects as xo


XC_GEOM_EPSILON = 1.e-15
XC_GEOM_S_MAX = 1.e21
XC_GEOM_ROOT_NEWTON_EPSILON = 1.e-10
XC_GEOM_ROOT_NEWTON_MAX_ITER = 100         # Maximum number of iterations in Newton's method
XC_GEOM_ROOT_NEWTON_DERIVATIVE_TOL = 1e-10 # Threshold for small derivative
XC_GEOM_ROOT_GRID_MAX_INTER = 10           # Maximum number of intervals for grid search
XC_GEOM_ROOT_GRID_POINTS = 1000            # Number of points to search in grid

define_src = f"""
#ifndef XCOLL_GEOM_DEFINES_H
#define XCOLL_GEOM_DEFINES_H
#include <math.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>

#ifndef XC_GEOM_EPSILON
#define XC_GEOM_EPSILON {XC_GEOM_EPSILON}
#endif

#ifndef XC_GEOM_S_MAX
#define XC_GEOM_S_MAX {XC_GEOM_S_MAX}
#endif

#ifndef XC_GEOM_ROOT_NEWTON_EPSILON
#define XC_GEOM_ROOT_NEWTON_EPSILON {XC_GEOM_ROOT_NEWTON_EPSILON}
#endif

#ifndef XC_GEOM_ROOT_NEWTON_MAX_ITER
#define XC_GEOM_ROOT_NEWTON_MAX_ITER {XC_GEOM_ROOT_NEWTON_MAX_ITER}
#endif

#ifndef XC_GEOM_ROOT_NEWTON_DERIVATIVE_TOL
#define XC_GEOM_ROOT_NEWTON_DERIVATIVE_TOL {XC_GEOM_ROOT_NEWTON_DERIVATIVE_TOL}
#endif

#ifndef XC_GEOM_ROOT_GRID_MAX_INTER
#define XC_GEOM_ROOT_GRID_MAX_INTER {XC_GEOM_ROOT_GRID_MAX_INTER}
#endif

#ifndef XC_GEOM_ROOT_GRID_POINTS
#define XC_GEOM_ROOT_GRID_POINTS {XC_GEOM_ROOT_GRID_POINTS}
#endif


#endif /* XCOLL_GEOM_DEFINES_H */
"""


def xo_to_ctypes(args):
    if not hasattr(args, '__iter__') or isinstance(args, str):
        args = [args]
    return ", ".join([f"{arg.get_c_type()} {arg.name}" for arg in args])

def xo_to_cnames(args):
    if not hasattr(args, '__iter__') or isinstance(args, str):
        args = [args]
    return ", ".join([f"{arg.name}" for arg in args])


class PyMethod:
    # Similar class as for the xt.BeamElement, but without the Metaclass magic
    # (and hence no need for PyMethodDescriptor)
    def __init__(self, kernel_name, element, element_name=None):
        self.kernel_name = kernel_name
        self.element = element
        self.element_name = element_name

    def __call__(self, *args, **kwargs):
        if len(args) > 0:
            raise ValueError("Kernel calling should be done with keyword arguments only!")
        instance = self.element
        context = instance._context
        # import pdb; pdb.set_trace()

        if instance.__class__._needs_compilation:
            # We don't have the HybridClass metaclass magic, so we need to manually replace ThisClass
            for ker in instance.__class__._kernels.values():
                for arg in ker.args:
                    if arg.atype == xo.ThisClass:
                        arg.atype = instance.__class__
            instance.__class__.compile_kernels(instance)#, save_source_as="temp.c")
            instance.__class__._needs_compilation = False
        kernel = context.kernels[self.kernel_name]
        if self.element_name:
            kwargs[self.element_name] = instance
        return kernel(**kwargs)
