# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

import numpy as np

import xobjects as xo

from ...general import _pkg_root
from ..trajectories import all_trajectories, DriftTrajectory


class LineSegment(xo.Struct):
    """Line segment between two points (s1, x1) -- (s2, x2)"""
    s1 = xo.Float64
    x1 = xo.Float64
    s2 = xo.Float64
    x2 = xo.Float64

    _depends_on = all_trajectories
    _extra_c_sources = [_pkg_root / 'geometry' / 'segments' / 'line.h']

    max_crossings = {DriftTrajectory: 2}

    def __repr__(self):
        return f"LineSegment(({self.s1:.3}, {self.x1:.3}) -- ({self.s2:.3}, {self.x2:.3}))"

    def evaluate(self, t):
        s1 = self.s1
        x1 = self.x1
        s2 = self.s2
        x2 = self.x2
        t = np.array(t)
        mask = (t >= 0) & (t <= 1)
        return s1*(1-t[mask]) + s2*t[mask], x1*(1-t[mask]) + x2*t[mask]

    def get_vertices(self):
        return (self.s1, self.x1), (self.s2, self.x2)
