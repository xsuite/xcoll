# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

import numpy as np

import xobjects as xo

from ..c_init import GeomCInit
from ....general import _pkg_root


class HalfOpenLineSegment(xo.Struct):
    """Half-open line segment from a point (s, x) to +/-infinity along an angle t"""
    s = xo.Float64
    x = xo.Float64
    t = xo.Float64 # angle (wrt s-axis) towards inf

    _depends_on = [GeomCInit]
    _extra_c_sources = [_pkg_root / 'scattering_routines' / 'geometry' / 'segments' / 'halfopen_line.h']

    def evaluate(self, t):
        s1 = self.s
        x1 = self.x
        s2 = s1 + np.cos(self.t)
        x2 = x1 + np.sin(self.t)
        t = np.array(t)
        mask = t >= 0
        return s1*(1-t[mask]) + s2*t[mask], x1*(1-t[mask]) + x2*t[mask]
