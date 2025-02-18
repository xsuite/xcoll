# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

import numpy as np

import xobjects as xo

from ...general import _pkg_root
from ..c_init import GeomCInit


class CircularSegment(xo.Struct):
    """Circular arc segment, defined by a centre and radius, and the starting/end angles defined anti-clockwise"""
    R  = xo.Float64
    sR = xo.Float64  # s-coordinate of the centre
    xR = xo.Float64  # x-coordinate of the centre
    theta1 = xo.Float64  # Starting angle
    theta2 = xo.Float64  # Ending angle

    _depends_on = [GeomCInit]
    _extra_c_sources = [_pkg_root / 'geometry' / 'segments' / 'circular.h']

    def __init__(self, *args, **kwargs):
        if 'theta1' in kwargs:
            while kwargs['theta1'] < -np.pi:
                kwargs['theta1'] += 2*np.pi
            while kwargs['theta1'] > np.pi:
                kwargs['theta1'] -= 2*np.pi
        if 'theta2' in kwargs:
            while kwargs['theta2'] < -np.pi:
                kwargs['theta2'] += 2*np.pi
            while kwargs['theta2'] > np.pi:
                kwargs['theta2'] -= 2*np.pi
        super().__init__(*args, **kwargs)

    def __str__(self):
        p1, p2 = self.get_vertices()
        return f"CircularSegment(({p1[0]:.3}, {p1[1]:.3})-b-({np.rad2deg(self.theta1):.0f}" + u'\xb0' \
             + f":{np.rad2deg(self.theta2):.0f}" + u'\xb0' + f":{self.R:.3})-b-({p2[0]:.3}, {p2[1]:.3}))"

    def get_vertices(self):
        s1 = self.round(self.s + self.R*np.cos(self.theta1))
        x1 = self.round(self.x + self.R*np.sin(self.theta1))
        s2 = self.round(self.s + self.R*np.cos(self.theta2))
        x2 = self.round(self.x + self.R*np.sin(self.theta2))
        return (s1, x1), (s2, x2)

    def _translate_inplace(self, ds, dx):
        self.sC += ds
        self.xC += dx

    def _rotate_inplace(self, ps, px, angle):
        c = np.cos(angle)
        s = np.sin(angle)
        self._translate_inplace(-ps, -px)
        new_sC = self.s * c - self.x * s
        new_xC = self.s * s + self.x * c
        self.s = new_sC
        self.x = new_xC
        self.theta1 += angle
        self.theta2 += angle
        while self.theta1 < -np.pi:
            self.theta1 += 2*np.pi
        while self.theta1 < -np.pi:
            self.theta1 += 2*np.pi
        while self.theta2 < -np.pi:
            self.theta2 += 2*np.pi
        while self.theta2 < -np.pi:
            self.theta2 += 2*np.pi
        self._translate_inplace(ps, px)

