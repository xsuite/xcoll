# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

import numpy as np

import xobjects as xo

from ...general import _pkg_root
from ..c_init import BoundingBox


class HalfOpenLineSegment(xo.Struct):
    """Half-open line segment from a point (s, x) to +/-infinity along an angle t"""
    s1 = xo.Float64
    x1 = xo.Float64
    sin_t1 = xo.Float64 # angle (wrt s-axis) towards inf
    cos_t1 = xo.Float64
    box = BoundingBox

    _extra_c_sources = [_pkg_root / 'geometry' / 'segments' / 'halfopen_line.h']
    _kernels = {'init_bounding_box': xo.Kernel(
                                        c_name='HalfOpenLineSegment_init_bounding_box',
                                        args=[xo.Arg(xo.ThisClass, name="seg"),
                                              xo.Arg(BoundingBox, name="box"),
                                              xo.Arg(xo.Float64, name="t1"),
                                              xo.Arg(xo.Float64, name="t2")], # this is not parameters of mcs??
                                        ret=None)}
    def __init__(self, *args, **kwargs):
        if not ('theta1' in kwargs or 'sin_t1' in kwargs or 'cos_t1' in kwargs):
            raise ValueError("At least one of 'theta1', 'sin_t1', or 'cos_t1' must be provided!")
        if 'theta1' in kwargs:
            theta1 = kwargs.pop('theta1')
            kwargs['sin_t1'] = np.sin(theta1)
            kwargs['cos_t1'] = np.cos(theta1)
        elif 'sin_t1' in kwargs and 'cos_t1' not in kwargs:
            kwargs['cos_t1'] = np.sqrt(1 - kwargs['sin_t1']**2)
        elif 'cos_t1' in kwargs and 'sin_t1' not in kwargs:
            kwargs['sin_t1'] = np.sqrt(1 - kwargs['cos_t1']**2)
        t1     = kwargs.pop('t1', 0.)
        t2     = kwargs.pop('t2', 10.)
        super().__init__(*args, **kwargs)
        self.box = BoundingBox()
        self.init_bounding_box(box=self.box, t1=t1, t2=t2)

    def __str__(self):
        return f"HalfOpenLineSegment(({self.s1:.3}, {self.x1:.3}) -- " \
            + f"{np.rad2deg(self.theta1):.0f}" + u'\xb0' + " * inf)"

    def get_vertices(self):
        return ((self.s1, self.x1),)

    def _translate_inplace(self, ds, dx):
        self.s1 += ds
        self.x1 += dx

    def _rotate_inplace(self, ps, px, angle):
        c = np.cos(angle)
        s = np.sin(angle)
        self._translate_inplace(-ps, -px)
        new_s = self.s1 * c - self.x1 * s
        new_x = self.s1 * s + self.x1 * c
        self.s1 = new_s
        self.x1 = new_x
        self.t += angle
        while self.t < -np.pi:
            self.t += 2*np.pi
        while self.t < -np.pi:
            self.t += 2*np.pi
        self._translate_inplace(ps, px)

    @property
    def theta1(self):
        return self.round(np.arctan2(self.sin_t1, self.cos_t1))

    @theta1.setter
    def theta1(self, value):
        while value < -np.pi:
            value += 2*np.pi
        while value > np.pi:
            value -= 2*np.pi
        self.sin_t1 = np.sin(value)
        self.cos_t1 = np.cos(value)
