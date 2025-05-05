# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

import numpy as np

import xobjects as xo
from ...general import _pkg_root

from .c_init import define_src, PyMethod


class BoundingBox(xo.Struct):
    rC = xo.Float64        # length of position vector to first vertex
    sin_tC = xo.Float64    # angle of position vector to first vertex, [radians]
    cos_tC = xo.Float64
    proj_l = xo.Float64    # projection of position vector on length: rC * (cos_t*cos_tC + sin_t*sin_tC)
    proj_w = xo.Float64    # projection of position vector on width:  rC * (cos_t*sin_tC - sin_t*cos_tC)
    l = xo.Float64         # length of the box
    w = xo.Float64         # width of the box
    sin_tb = xo.Float64    # orientation of the box (angle of length wrt horizontal)
    cos_tb = xo.Float64

    _kernels = {'overlaps': xo.Kernel(
                                c_name='BoundingBox_overlaps',
                                args=[xo.Arg(xo.ThisClass, name="b1"),
                                      xo.Arg(xo.ThisClass, name="b2")],
                                ret=xo.Arg(xo.Int8, name="overlaps")),
                'set_params': xo.Kernel(
                                c_name='BoundingBox_set_params',
                                args=[xo.Arg(xo.ThisClass, name="box"),
                                      xo.Arg(xo.Float64, name="rC"),
                                      xo.Arg(xo.Float64, name="sin_tC"),
                                      xo.Arg(xo.Float64, name="cos_tC"),
                                      xo.Arg(xo.Float64, name="l"),
                                      xo.Arg(xo.Float64, name="w"),
                                      xo.Arg(xo.Float64, name="sin_tb"),
                                      xo.Arg(xo.Float64, name="cos_tb")],
                                ret=None)}
    _needs_compilation = True
    _extra_c_sources = [
        define_src,
        _pkg_root / 'geometry' / 'c_init' / 'sort.h',
        _pkg_root / 'geometry' / 'c_init' / 'methods.h',
        _pkg_root / 'geometry' / 'c_init' / 'find_root.h',
        _pkg_root / 'geometry' / 'c_init' / 'bounding_box.h',
    ]

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        #self._set_params(**kwargs)

    def __str__(self):
        """Return str(self)."""
        return f"BoundingBox([{self.l:.3} x {self.w:.3}], {self.theta:.3} rad, ({self.s1:.3}, {self.x1:.3})"

    def __repr__(self):
        """Return repr(self)."""
        return f"<{str(self)} at {hex(id(self))}>"

    # Add kernel
    def __getattr__(self, attr):
        kernel_name = attr
        if kernel_name in self._kernels:
            return PyMethod(kernel_name=kernel_name, element=self, element_name='b1')
        raise ValueError(f"Attribute {attr} not found in {self.__class__.__name__}")

    @property
    def s1(self):
        return self.rC * self.cos_tC

    @property
    def x1(self):
        return self.rC * self.sin_tC

    @property
    def s2(self):
        return self.rC * self.cos_tC + self.l * self.cos_tb

    @property
    def x2(self):
        return self.rC * self.sin_tC + self.l * self.sin_tb

    @property
    def s3(self):
        return self.rC * self.cos_tC + self.l * self.cos_tb - self.w * self.sin_tb

    @property
    def x3(self):
        return self.rC * self.sin_tC + self.l * self.sin_tb + self.w * self.cos_tb

    @property
    def s4(self):
        return self.rC * self.cos_tC - self.w * self.sin_tb

    @property
    def x4(self):
        return self.rC * self.sin_tC + self.w * self.cos_tb

    @property
    def theta(self):
        return np.arctan2(self.sin_tb, self.cos_tb)

    @property
    def vertices(self):
        return (self.s1, self.x1), (self.s2, self.x2), (self.s3, self.x3), (self.s4, self.x4)
