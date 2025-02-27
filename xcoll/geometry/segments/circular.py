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
    _theta1 = xo.Float64  # Starting angle
    _theta2 = xo.Float64  # Ending angle

    _depends_on = [GeomCInit]
    _extra_c_sources = [_pkg_root / 'geometry' / 'segments' / 'circular.h']

    def __init__(self, *args, **kwargs):
        if 'R' not in kwargs:
            raise ValueError("Radius must be provided")
        if kwargs['R'] < 0:
            raise ValueError("Radius must be positive")
        theta1 = kwargs.pop('theta1', -np.pi)
        theta2 = kwargs.pop('theta2', np.pi)
        super().__init__(*args, **kwargs)
        self.set_angles(theta1, theta2)

    def __str__(self):
        p1, p2 = self.get_vertices()
        return f"CircularSegment(({p1[0]:.3}, {p1[1]:.3})-b-({np.rad2deg(self.theta1):.0f}" + u'\xb0' \
             + f":{np.rad2deg(self.theta2):.0f}" + u'\xb0' + f":{self.R:.3})-b-({p2[0]:.3}, {p2[1]:.3}))"

    def get_vertices(self):
        return (self.s1, self.x1), (self.s2, self.x2)

    def get_control_points(self):
        return (self.sR, self.xR),

    def _translate_inplace(self, ds, dx):
        self.sR += ds
        self.xR += dx

    def _rotate_inplace(self, ps, px, angle):
        c = np.cos(angle)
        s = np.sin(angle)
        self._translate_inplace(-ps, -px)
        new_sR = self.sR * c - self.xR * s
        new_xR = self.sR * s + self.xR * c
        self.sR = new_sR
        self.xR = new_xR
        self.set_angles(self.theta1 + angle, self.theta2 + angle)
        self._translate_inplace(ps, px)

    @property
    def s1(self):
        return self.round(self.sR + self.R*np.cos(self.theta1))

    @property
    def x1(self):
        return self.round(self.xR + self.R*np.sin(self.theta1))

    @property
    def s2(self):
        return self.round(self.sR + self.R*np.cos(self.theta2))

    @property
    def x2(self):
        return self.round(self.xR + self.R*np.sin(self.theta2))

    @property
    def theta1(self):
        return self._theta1

    @property
    def theta2(self):
        # We want to represent angles in [-pi, pi]
        value = self._theta2
        while value > np.pi:
            value -= 2*np.pi
        return value

    def set_angles(self, theta1, theta2):
        while theta1 < -np.pi:
            theta1 += 2*np.pi
        while theta1 > np.pi:
            theta1 -= 2*np.pi
        while theta2 < -np.pi:
            theta2 += 2*np.pi
        while theta2 > np.pi:
            theta2 -= 2*np.pi
        self._theta1 = theta1
        if theta2 >= theta1:
            self._theta2 = theta2
        else:
            # If theta2 is smaller than theta1, it means we have done a full turn
            self._theta2 = theta2 + 2*np.pi
