# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

import numpy as np

import xobjects as xo

from ...general import _pkg_root
from ..c_init import BoundingBox
from ..c_init.c_init import define_src

class LineSegment(xo.Struct):
    """Line segment between two points (s1, x1) -- (s2, x2)"""
    s1 = xo.Float64
    x1 = xo.Float64
    s2 = xo.Float64
    x2 = xo.Float64
    box = BoundingBox

    _extra_c_sources = [define_src,
                        _pkg_root / 'geometry' / 'segments' / 'line.h']
    _kernels = {'update_box': xo.Kernel(
                                        c_name='LineSegment_update_box',
                                        args=[xo.Arg(xo.ThisClass, name="seg"),
                                              #xo.Arg(BoundingBox, name="box"),
                                              xo.Arg(xo.Float64, name="t1"),
                                              xo.Arg(xo.Float64, name="t2")],
                                        ret=None)}

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.box = BoundingBox()
        self.init_box(t1=0., t2=1.)

    def __str__(self):
        return f"LineSegment(({self.s1:.3}, {self.x1:.3}) -- ({self.s2:.3}, {self.x2:.3}))"

    def get_vertices(self):
        return (self.s1, self.x1), (self.s2, self.x2)

    def _translate_inplace(self, ds, dx):
        self.s1 += ds
        self.x1 += dx
        self.s2 += ds
        self.x2 += dx

    def _rotate_inplace(self, ps, px, angle):
        c = np.cos(angle)
        s = np.sin(angle)
        self._translate_inplace(-ps, -px)
        new_s1 = self.s1 * c - self.x1 * s
        new_x1 = self.s1 * s + self.x1 * c
        self.s1 = new_s1
        self.x1 = new_x1
        new_s2 = self.s2 * c - self.x2 * s
        new_x2 = self.s2 * s + self.x2 * c
        self.s2 = new_s2
        self.x2 = new_x2
        self._translate_inplace(ps, px)

    def init_box(self, t1, t2):
        if t1 >= t2:
            raise ValueError("t1 must be smaller than t2!")
        if t1 < 0 or t1 > 1:
            raise ValueError("t1 must be in [0, 1]!!")
        if t2 < 0 or t2 > 1:
            raise ValueError("t2 must be in [0, 1]!!")
        self.update_box(seg=self, t1=t1, t2=t2)
