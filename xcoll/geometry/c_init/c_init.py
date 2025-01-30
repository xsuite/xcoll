# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

import xobjects as xo
from ...general import _pkg_root


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
class McsLineParams(xo.Struct):
    Xo = xo.Float64
    Ax = xo.Float64
    s1 = xo.Float64
    x1 = xo.Float64
    s2 = xo.Float64
    x2 = xo.Float64

class PyMethod:
    # Similar class as for the xt.BeamElement, but without the Metaclass magic
    # (and hence no need for PyMethodDescriptor)
    def __init__(self, kernel_name, element, element_name=None):
        self.kernel_name = kernel_name
        self.element = element
        self.element_name = element_name

    def __call__(self, **kwargs):
        instance = self.element
        context = instance._context
        # import pdb; pdb.set_trace()

        if instance.__class__._needs_compilation:
            # We don't have the HybridClass metaclass magic, so we need to manually replace ThisClass
            for ker in instance.__class__._kernels.values():
                for arg in ker.args:
                    if arg.atype == xo.ThisClass:
                        arg.atype = instance.__class__
            instance.__class__.compile_kernels(instance, save_source_as="temp.c")
            instance.__class__._needs_compilation = False
        kernel = context.kernels[self.kernel_name]
        if self.element_name:
            kwargs[element_name] = instance

        return kernel(**kwargs)

class GeomCInit(xo.Struct):
    _extra_c_sources = [
        define_src,
        _pkg_root / 'geometry' / 'c_init' / 'sort.h',
        _pkg_root / 'geometry' / 'c_init' / 'methods.h',
        _pkg_root / 'geometry' / 'segments' / 'line.h',     # remove after testing is done, add the others as well
        _pkg_root / 'geometry' / 'segments' / 'halfopen_line.h',
        _pkg_root / 'geometry' / 'segments' / 'circular.h',
        _pkg_root / 'geometry' / 'c_init' / 'find_root.h',
    ]

    # A Struct needs something to depend on, otherwise the class is added twice in the cdefs during compilation
    _depends_on = [McsLineParams]

    _needs_compilation = True

    _kernels = {"grid_search_and_newton":
                xo.Kernel(
                    c_name="grid_search_and_newton",
                    args=[
                        xo.Arg(xo.Int8, "traj_id"),
                        xo.Arg(xo.Int8, "seg_id"),
                        xo.Arg(xo.Float64, "s_min"),
                        xo.Arg(xo.Float64, "s_max"),
                        xo.Arg(xo.Float64, name="roots", pointer=True),
                        xo.Arg(xo.Float64, "max_crossings"),
                        xo.Arg(McsLineParams, name="params"),
                        xo.Arg(xo.Int8, name="number_of_roots", pointer=True)
                    ],
                    ret=None)
                }
    def __getattr__(self, attr):
        if attr in self._kernels:
            return PyMethod(kernel_name=attr, element=self, element_name='shape')
        raise ValueError(f"Attribute {attr} not found in {self.__class__.__name__}")

# class McsHalfOpenLineParams(xo.Struct):
#     Xo = xo.Float64
#     Ax = xo.Float64
#     s1 = xo.Float64
#     x1 = xo.Float64
#     s2 = xo.Float64
#     x2 = xo.Float64

# class McsCircleParams(xo.Struct):
#     Xo = xo.Float64
#     Ax = xo.Float64
#     R  = xo.Float64
#     sC = xo.Float64
#     xC = xo.Float64
#     x  = xo.Float64    

# params = McsLineParams(Xo=0, Ax=0, s1=0, x1=0, s2=0, x2=0)
# max_crossings = 8
# roots = np.zeros(max_crossings, dtype=np.float64)
# number_of_roots = np.zeros(1, dtype=np.int8)
# GeomCInit.grid_search_and_newton(traj_id=0, seg_id=0, s_min=0, s_max=1, roots=roots, params=params,
#                                 max_crossings=max_crossings, number_of_roots=number_of_roots)



