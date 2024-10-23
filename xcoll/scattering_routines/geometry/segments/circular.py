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

    def __init__(self, *args, **kwargs):
        if 't1' in kwargs:
            while kwargs['t1'] < -np.pi:
                kwargs['t1'] += 2*np.pi
            while kwargs['t1'] > np.pi:
                kwargs['t1'] -= 2*np.pi
        if 't2' in kwargs:
            while kwargs['t2'] < -np.pi:
                kwargs['t2'] += 2*np.pi
            while kwargs['t2'] > np.pi:
                kwargs['t2'] -= 2*np.pi
        super().__init__(*args, **kwargs)

    def __repr__(self):
        p1, p2 = self.get_vertices()
        return f"CircularSegment(({p1[0]:.3}, {p1[1]:.3})-b-({np.rad2deg(self.t1):.0f}" + u'\xb0' \
             + f":{np.rad2deg(self.t2):.0f}" + u'\xb0' + f":{self.R:.3})-b-({p2[0]:.3}, {p2[1]:.3}))"

    def evaluate(self, t):
        R = self.R
        sC = self.s
        xC = self.x
        t1 = self.t1
        t2 = self.t2
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

    def get_vertices(self):
        s1 = self.round(self.s + self.R*np.cos(self.t1))
        x1 = self.round(self.x + self.R*np.sin(self.t1))
        s2 = self.round(self.s + self.R*np.cos(self.t2))
        x2 = self.round(self.x + self.R*np.sin(self.t2))
        return (s1, x1), (s2, x2)

