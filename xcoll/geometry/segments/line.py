# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

import numpy as np

import xobjects as xo

from ...general import _pkg_root
from ..c_init import BoundingBox


class LineSegment(xo.Struct):
    """Line segment between two points (s1, x1) -- (s2, x2)"""
    s1 = xo.Float64
    x1 = xo.Float64
    s2 = xo.Float64
    x2 = xo.Float64
    _t1 = xo.Float64  # parameter along line for first point (default 0)
    _t2 = xo.Float64  # parameter along line for second point (default 1
    box = BoundingBox

    _extra_c_sources = [_pkg_root / 'geometry' / 'segments' / 'line.h']
    _kernels = {'init_bounding_box': xo.Kernel(
                                        c_name='LineSegment_init_bounding_box',
                                        args=[xo.Arg(xo.ThisClass, name="seg"),
                                                xo.Arg(BoundingBox, name="box"),
                                                xo.Arg(xo.Float64, name="t1"),
                                                xo.Arg(xo.Float64, name="t2")], # this is not parameters of mcs??
                                        ret=None)}

    def __init__(self, *args, **kwargs):
        t1 = kwargs.pop('t1', 0.)
        t2 = kwargs.pop('t2', 1.)
        super().__init__(*args, **kwargs)
        self._t1 = t1
        self._t2 = t2
        self.box = BoundingBox()
        self.init_bounding_box(box=self.box, t1=t1, t2=t2)

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

    @property
    def t1(self):
        return self._t1

    @t1.setter
    def t1(self, val):
        if val >= self._t2:
            raise ValueError("t1 must be smaller than t2!")
        if val < 0 or val > 1:
            raise ValueError("t1 must be in [0, 1]!")
        self._t1 = val
        self.init_bounding_box(box=self.box, t1=self._t1, t2=self._t2)

    @property
    def t2(self):
        return self._t2

    @t2.setter
    def t2(self, val):
        if val <= self._t1:
            raise ValueError("t2 must be larger than t1!")
        if val < 0 or val > 1:
            raise ValueError("t2 must be in [0, 1]!")
        self._t2 = val
        self.init_bounding_box(box=self.box, t1=self._t1, t2=self._t2)