# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

import numpy as np

import xobjects as xo

from ...general import _pkg_root
from ..c_init import BoundingBox
from ..trajectories import DriftTrajectory, CircularTrajectory, MultipleCoulombTrajectory

# Quartic vs Cubic Bezier

class BezierSegment(xo.Struct):
    """Bézier segment, defined by a start and end point P1 and P2, and two control points that define the curve"""
    _s1 = xo.Float64
    _x1 = xo.Float64
    _s2 = xo.Float64
    _x2 = xo.Float64
    _cs1 = xo.Float64
    _cx1 = xo.Float64
    _cs2 = xo.Float64
    _cx2 = xo.Float64
    _ts1 = xo.Float64 # First extremum in s
    _ts2 = xo.Float64 # Second extremum in s
    _tx1 = xo.Float64 # First extremum in x
    _tx2 = xo.Float64 # Second extremum in x
    _es1 = xo.Float64  # Value of first extremum in s
    _es2 = xo.Float64  # Value of second extremum in s
    _ex1 = xo.Float64  # Value of first extremum in x
    _ex2 = xo.Float64  # Value of second extremum in x
    box = BoundingBox

    _extra_c_sources = [_pkg_root / 'geometry' / 'segments' / 'bezier.h']

    _kernels = {'calculate_extrema': xo.Kernel(
                                c_name='BezierSegment_calculate_extrema',
                                args=[xo.Arg(xo.ThisClass, name="seg")],
                                ret=None),
                'init_bounding_box': xo.Kernel(
                        c_name='BezierSegment_init_bounding_box',
                        args=[xo.Arg(xo.ThisClass, name="seg"),
                              xo.Arg(xo.ThisClass,  name="box"),
                              xo.Arg(xo.Float64, name="t1"),
                              xo.Arg(xo.Float64, name="t2")], # this is not parameters of mcs??
                        ret=None)
                }   

    _max_crossings = {DriftTrajectory: 3, CircularTrajectory: 6, MultipleCoulombTrajectory: 6}

    def __init__(self, *, s1, x1, s2, x2, cs1, cx1, cs2, cx2, **kwargs):
        kwargs['_s1'] = s1
        kwargs['_x1'] = x1
        kwargs['_s2'] = s2
        kwargs['_x2'] = x2
        kwargs['_cs1'] = cs1
        kwargs['_cx1'] = cx1
        kwargs['_cs2'] = cs2
        kwargs['_cx2'] = cx2
        super().__init__(**kwargs)
        self.box = BoundingBox()
        self.init_bounding_box(box=self.box, t1=0., t2=1.)
        self.calculate_extrema()

    def __str__(self):
        return f"BezierSegment(({self.s1:.3}, {self.x1:.3})-c-({self.cs1:.3}, {self.cx1:.3}) -- " \
             + f"({self.cs2:.3}, {self.cx2:.3})-c-({self.s2:.3}, {self.x2:.3}))"

    def get_vertices(self):
        return (self.s1, self.x1), (self.s2, self.x2)

    def get_control_points(self):
        return (self.cs1, self.cx1), (self.cs2, self.cx2)

    def _translate_inplace(self, ds, dx):
        self.s1 += ds
        self.x1 += dx
        self.s2 += ds
        self.x2 += dx
        self.cs1 += ds
        self.cx1 += dx
        self.cs2 += ds
        self.cx2 += dx

    def _rotate_inplace(self, ps, px, angle):
        c = np.cos(angle)
        s = np.sin(angle)
        self._translate_inplace(-ps, -px)
        new_s1 = self.s1 * c - self.x1 * s
        new_x1 = self.s1 * s + self.x1 * c
        new_s2 = self.s2 * c - self.x2 * s
        new_x2 = self.s2 * s + self.x2 * c
        self.s1 = new_s1
        self.x1 = new_x1
        self.s2 = new_s2
        self.x2 = new_x2
        new_cs1 = self.cs1 * c - self.cx1 * s
        new_cx1 = self.cs1 * s + self.cx1 * c
        new_cs2 = self.cs2 * c - self.cx2 * s
        new_cx2 = self.cs2 * s + self.cx2 * c
        self.cs1 = new_cs1
        self.cx1 = new_cx1
        self.cs2 = new_cs2
        self.cx2 = new_cx2
        self._translate_inplace(ps, px)

    @property
    def s1(self):
        return self._s1

    @s1.setter
    def s1(self, value):
        self._s1 = value
        self.calculate_extrema()

    @property
    def x1(self):
        return self._x1

    @x1.setter
    def x1(self, value):
        self._x1 = value
        self.calculate_extrema()

    @property
    def s2(self):
        return self._s2

    @s2.setter
    def s2(self, value):
        self._s2 = value
        self.calculate_extrema()

    @property
    def x2(self):
        return self._x2

    @x2.setter
    def x2(self, value):
        self._x2 = value
        self.calculate_extrema()

    @property
    def cs1(self):
        return self._cs1

    @cs1.setter
    def cs1(self, value):
        self._cs1 = value
        self.calculate_extrema()

    @property
    def cx1(self):
        return self._cx1

    @cx1.setter
    def cx1(self, value):
        self._cx1 = value
        self.calculate_extrema()

    @property
    def cs2(self):
        return self._cs2

    @cs2.setter
    def cs2(self, value):
        self._cs2 = value
        self.calculate_extrema()

    @property
    def cx2(self):
        return self._cx2

    @cx2.setter
    def cx2(self, value):
        self._cx2 = value
        self.calculate_extrema()
