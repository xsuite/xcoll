# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

import numpy as np

import xobjects as xo

from ...general import _pkg_root
from ..c_init import GeomCInit


class HalfOpenLineSegment(xo.Struct):
    """Half-open line segment from a point (s, x) to +/-infinity along an angle t"""
    s1 = xo.Float64
    x1 = xo.Float64
    sin_t1 = xo.Float64 # angle (wrt s-axis) towards inf
    cos_t1 = xo.Float64

    _depends_on = [GeomCInit]
    _extra_c_sources = [_pkg_root / 'geometry' / 'segments' / 'halfopen_line.h']

    def __init__(self, *args, **kwargs):
        theta1 = kwargs.pop('theta1', None)
        super().__init__(*args, **kwargs)
        if theta1 is not None:
            self.theta1 = theta1

    def __str__(self):
        return f"HalfOpenLineSegment(({self.s:.3}, {self.x:.3}) -- " \
            + f"{np.rad2deg(self.t):.0f}" + u'\xb0' + " * inf)"

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

    def get_vertices(self):
        return ((self.s, self.x),)

    def _translate_inplace(self, ds, dx):
        self.s += ds
        self.x += dx

    def _rotate_inplace(self, ps, px, angle):
        c = np.cos(angle)
        s = np.sin(angle)
        self._translate_inplace(-ps, -px)
        new_s = self.s * c - self.x * s
        new_x = self.s * s + self.x * c
        self.s = new_s
        self.x = new_x
        self.t += angle
        while self.t < -np.pi:
            self.t += 2*np.pi
        while self.t < -np.pi:
            self.t += 2*np.pi
        self._translate_inplace(ps, px)
