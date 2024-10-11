# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

import numpy as np

import xobjects as xo

from ....general import _pkg_root
from ..c_init import GeomCInit


class LineSegment(xo.Struct):
    """Line segment between two points (s1, x1) -- (s2, x2)"""
    s1 = xo.Float64
    x1 = xo.Float64
    s2 = xo.Float64
    x2 = xo.Float64

    _depends_on = [GeomCInit]
    _extra_c_sources = [_pkg_root / 'scattering_routines' / 'geometry' / 'segments' / 'line.h']

    def evaluate(self, t):
        s1 = self.s1
        x1 = self.x1
        s2 = self.s2
        x2 = self.x2
        t = np.array(t)
        mask = (t >= 0) & (t <= 1)
        return s1*(1-t[mask]) + s2*t[mask], x1*(1-t[mask]) + x2*t[mask]
