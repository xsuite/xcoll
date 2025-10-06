# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

import numpy as np

import xobjects as xo

from ...general import _pkg_root
from ..c_init import BoundingBox


class DriftTrajectory(xo.Struct):
    """A trajectory describing a drift.

    In parametrised form, it is given by:
        s(𝜆) = s0 + 𝜆 cos(𝜃0)     𝜆  = 0..∞
        x(𝜆) = x0 + 𝜆 sin(𝜃0)     𝜃0 = -π/2..π/2
        𝜃(𝜆) = 𝜃0

    where (s0, x0) and 𝜃0 are the initial particle coordinates resp. angle, and
    𝜆 represents the distance travelled along the direction 𝜃0.

    In practice, we do not provide 𝜃0 but xp = tan(𝜃0).
    """

    s0 = xo.Float64
    x0 = xo.Float64
    sin_t0 = xo.Float64
    cos_t0 = xo.Float64
    tan_t0 = xo.Float64
    box = BoundingBox

    _extra_c_sources = [_pkg_root / 'geometry' / 'trajectories' / 'drift.h']

    _kernels = {'set_params': xo.Kernel(
                                c_name='DriftTrajectory_set_params',
                                args=[xo.Arg(xo.ThisClass, name="traj"),
                                      xo.Arg(xo.Float64, name="s0"),
                                      xo.Arg(xo.Float64, name="x0"),
                                      xo.Arg(xo.Float64, name="xp")],
                                ret=None),
                'update_box': xo.Kernel(
                                c_name='DriftTrajectory_update_box',
                                args=[xo.Arg(xo.ThisClass, name="traj"),
                                      xo.Arg(xo.Float64, name="l1"),
                                      xo.Arg(xo.Float64, name="l2")],
                                ret=None)}

    def __init__(self, *args, **kwargs):
        xp = kwargs.pop('xp', False)
        theta0 = kwargs.pop('theta0', False)
        super().__init__(*args, **kwargs)
        if xp is not False:
            self.set_params(s0=self.s0, x0=self.x0, xp=xp)
        elif theta0 is not False:
            self.set_params(s0=self.s0, x0=self.x0, xp=np.tan(theta0))
        self.box = BoundingBox()
        self.update_box(box=self.box, l1=0., l2=1.)

    def __str__(self):
        return f"DriftTrajectory(s0={self.s0}, x0={self.x0}, xp={self.xp})"

    @property
    def xp(self):
        return self.tan_t0

    @property
    def theta0(self):
        return self.round(np.arctan2(self.sin_t0, self.cos_t0))

    def update_box(self, l1, l2):
        if l1 >= l2:
            raise ValueError("l1 must be smaller than l2!")
        if l1 < 0 or l1 > 1:
            raise ValueError("l1 must be in [0, 1]!")
        if l2 < 0 or l2 > 1:
            raise ValueError("l2 must be in [0, 1]!")
        self.update_box(box=self.box, l1=l1, l2=l2)
