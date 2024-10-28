# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

import numpy as np

import xobjects as xo

from ...general import _pkg_root
from ..trajectories import all_trajectories, DriftTrajectory


class BezierSegment(xo.Struct):
    """Bézier segment, defined by a start and end point P1 and P2, and two control points that define the curve"""
    s1 = xo.Float64
    x1 = xo.Float64
    s2 = xo.Float64
    x2 = xo.Float64
    cs1 = xo.Float64
    cx1 = xo.Float64
    cs2 = xo.Float64
    cx2 = xo.Float64

    _depends_on = all_trajectories
    _extra_c_sources = [_pkg_root / 'geometry' / 'segments' / 'bezier.h']

    max_crossings = {DriftTrajectory: 3}

    def __repr__(self):
        return f"BezierSegment(({self.s1:.3}, {self.x1:.3})-c-({self.cs1:.3}, {self.cx1:.3}) -- " \
             + f"({self.cs2:.3}, {self.cx2:.3})-c-({self.s2:.3}, {self.x2:.3}))"

    def evaluate(self, t):
        s1  = self.s1
        x1  = self.x1
        s2  = self.s2
        x2  = self.x2
        cs1 = self.cs1
        cx1 = self.cx1
        cs2 = self.cs2
        cx2 = self.cx2
        t = np.array(t)
        mask = (t >= 0) & (t <= 1)
        return (1-t[mask])**3 * s1 + 3*t[mask]*(1-t[mask])**2 * cs1 + 3*(1-t[mask])*t[mask]**2 * cs2 + t[mask]**3 * s2, \
               (1-t[mask])**3 * x1 + 3*t[mask]*(1-t[mask])**2 * cx1 + 3*(1-t[mask])*t[mask]**2 * cx2 + t[mask]**3 * x2

    def get_vertices(self):
        return (self.s1, self.x1), (self.s2, self.x2)

