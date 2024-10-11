# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

import numpy as np

import xobjects as xo

from ..c_init import GeomCInit
from ....general import _pkg_root


class CircularSegment(xo.Struct):
    """Circular arc segment, defined by a centre and radius, and the starting/end angles defined anti-clockwise"""
    R  = xo.Float64
    s  = xo.Float64  # s-coordinate of the centre
    x  = xo.Float64  # x-coordinate of the centre
    t1 = xo.Float64  # Starting angle
    t2 = xo.Float64  # Ending angle

    _depends_on = [GeomCInit]
    _extra_c_sources = [_pkg_root / 'scattering_routines' / 'geometry' / 'segments' / 'circular.h']

    def evaluate(self, t):
        R = self.R
        sC = self.s
        xC = self.x
        t1 = self.t1
        t2 = self.t2
        while t1 < -np.pi:
            t1 += 2*np.pi
        while t1 > np.pi:
            t1 -= 2*np.pi
        while t2 < -np.pi:
            t2 += 2*np.pi
        while t2 > np.pi:
            t2 -= 2*np.pi
        t = np.array(t)
        while np.any(t < -np.pi):
            t[t < -np.pi] += 2*np.pi
        while np.any(t > np.pi):
            t[t > np.pi] -= 2*np.pi
        if t1 < t2:
            mask = (t1 <= t) & (t <= t2)
        else:
            mask = (t1 <= t) | (t <= t2)
        return sC + R*np.cos(t[mask]), xC + R*np.sin(t[mask])
