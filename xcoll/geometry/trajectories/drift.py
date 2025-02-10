# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

import numpy as np

import xobjects as xo

from ...general import _pkg_root
from ..c_init import GeomCInit


class DriftTrajectory(xo.Struct):
    """A trajectory describing a drift.

    In parametrised form, it is given by:
        s(𝜆) = s0 + 𝜆 cos(𝜃0)     𝜆  = 0..∞
        x(𝜆) = x0 + 𝜆 sin(𝜃0)     𝜃0 = -π..π
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

    _depends_on = [GeomCInit]
    _extra_c_sources = [_pkg_root / 'geometry' / 'trajectories' / 'drift.h']

    _kernels = {'set_params': xo.Kernel(
                                c_name='DriftTrajectory_set_params',
                                args=[xo.Arg(xo.ThisClass, name="traj"),
                                      xo.Arg(xo.Float64, name="s0"),
                                      xo.Arg(xo.Float64, name="x0"),
                                      xo.Arg(xo.Float64, name="xp")],
                                ret=None)}

    def __init__(self, *args, **kwargs):
        xp = kwargs.pop('xp', False)
        theta0 = kwargs.pop('theta0', False)
        super().__init__(*args, **kwargs)
        if xp is not False:
            self.set_params(s0=self.s0, x0=self.x0, xp=xp)
        elif theta0 is not False:
            self.sin_t0 = np.sin(theta0)
            self.cos_t0 = np.cos(theta0)
            self.tan_t0 = np.tan(theta0)

    def __str__(self):
        return f"DriftTrajectory(s0={self.s0}, x0={self.x0}, xp={self.xp})"

    @property
    def xp(self):
        return self.tan_t0

    @property
    def theta0(self):
        return self.round(np.arctan2(self.sin_t0, self.cos_t0))

#     args_hv = [
#             # The arguments that define the particle trajectory, common to both planes
#             xo.Arg(xo.Float64, pointer=False, name="s0")  # Particle s
#     ]
#     args_h = [
#             # The arguments that define the horizontal (after rotation) particle trajectory
#             xo.Arg(xo.Float64, pointer=False, name="x0"),  # Particle x
#             xo.Arg(xo.Float64, pointer=False, name="xm")   # Particle slope in the x direction (= xp = tan(theta_x))
#     ]
#     args_v = [
#             # The arguments that define the vertical (after rotation) particle trajectory
#             xo.Arg(xo.Float64, pointer=False, name="y0"),  # Particle y
#             xo.Arg(xo.Float64, pointer=False, name="ym")   # Particle slope in the y direction (= yp = tan(theta_y))
#     ]
