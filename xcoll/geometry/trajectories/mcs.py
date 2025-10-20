## copyright ############################### #
## This file is part of the Xcoll package.   #
## Copyright (c) CERN, 2025.                 #
## ######################################### #

import numpy as np

import xobjects as xo

from ...general import _pkg_root
from ..c_init.bounding_box import BoundingBox


class MultipleCoulombTrajectory(xo.Struct):
    """A trajectory describing a multiple coulomb scattering trajectory.

    In parametrised form, it is given by:
        s(𝜆) = s0 + 𝜆 cos(𝜃0) - (𝜉1/√12 + 𝜉2/2) 𝜆 𝛺(𝜆) sin(𝜃0)     𝜆  = 0..∞
        x(𝜆) = x0 + 𝜆 sin(𝜃0) + (𝜉1/√12 + 𝜉2/2) 𝜆 𝛺(𝜆) cos(𝜃0)     𝜃0 = -π..π
        𝜃(𝜆) = 𝜃0 + 𝜉2 𝛺(𝜆)

    where (s0, x0) and 𝜃0 are the initial particle coordinates resp. angle,
    𝜆 represents the travelled distance projected along the direction 𝜃0, and
    𝜉1 and 𝜉2 are two random normal variables. Finally, 𝛺(𝜆) represented the
    expected average scattering angle and is estimated by
        𝛺(𝜆) = (13.6 MeV) / (𝛽 pc) sqrt(𝜆 q^2 / (X0 𝛽^2)) (1 + 0.038 ln (𝜆 q^2 / (X0 𝛽^2)))

    where X0 is the material's radiation length, and q, 𝛽, and pc are the particle's charge,
    relativistic 𝛽, and momentum.

    In practice, we do not provide 𝜃0 but xp = tan(𝜃0).
    """

    s0 = xo.Float64
    x0 = xo.Float64
    sin_t0 = xo.Float64
    cos_t0 = xo.Float64
    tan_t0 = xo.Float64
    Xt0 = xo.Float64  #  X0 𝛽^2 / q^2
    A0 = xo.Float64   # (𝜉1/√12 + 𝜉2/2) (13.6 MeV) / (𝛽 pc)
    B0 = xo.Float64   # 𝜉2 (13.6 MeV) / (𝛽 pc)
    _l1 = xo.Float64  # start parameter along trajectory (default -5)
    _l2 = xo.Float64  # end parameter along trajectory (default 5)

    _depends_on = [BoundingBox]
    _extra_c_sources = [_pkg_root / 'geometry' / 'trajectories' / 'mcs.h']

    _kernels = {'set_params': xo.Kernel(
                                c_name='MultipleCoulombTrajectory_set_params',
                                args=[xo.Arg(xo.ThisClass, name="traj"),
                                      xo.Arg(xo.Float64, name="X0"),
                                      xo.Arg(xo.Float64, name="ran_1"),
                                      xo.Arg(xo.Float64, name="ran_2"),
                                      xo.Arg(xo.Float64, name="s0"),
                                      xo.Arg(xo.Float64, name="x0"),
                                      xo.Arg(xo.Float64, name="xp"),
                                      xo.Arg(xo.Float64, name="pc"),
                                      xo.Arg(xo.Float64, name="beta"),
                                      xo.Arg(xo.Float64, name="q")],
                                ret=None),
                'update_box': xo.Kernel(
                                c_name='MultipleCoulombTrajectory_update_box',
                                args=[xo.Arg(xo.ThisClass, name="traj"),
                                      xo.Arg(BoundingBox, name="box"),
                                      xo.Arg(xo.Float64, name="l1"),
                                      xo.Arg(xo.Float64, name="l2")],
                                ret=None)}

    def __init__(self, *args, **kwargs):
        X0 = kwargs.pop('X0', False)
        ran_1 = kwargs.pop('ran_1', False)
        ran_2 = kwargs.pop('ran_2', False)
        xp = kwargs.pop('xp', False)
        theta0 = kwargs.pop('theta0', False)
        pc = kwargs.pop('pc', False)
        beta = kwargs.pop('beta', False)
        q = kwargs.pop('q', False)
        super().__init__(*args, **kwargs)
        if pc is not False and beta is not False and q is not False and X0 is not False\
        and ran_1 is not False and ran_2 is not False:
            if xp is not False:
                self.set_params(X0=X0, ran_1=ran_1, ran_2=ran_2, s0=self.s0, x0=self.x0,
                                xp=xp, pc=pc, beta=beta, q=q)
            elif theta0 is not False:
                self.set_params(X0=X0, ran_1=ran_1, ran_2=ran_2, s0=self.s0, x0=self.x0,
                                xp=np.tan(theta0), pc=pc, beta=beta, q=q)

    def __str__(self):
        return f"MultipleCoulombTrajectory(s0={self.s0}, x0={self.x0}, xp={self.xp})"

    @property
    def xp(self):
        return self.tan_t0

    @property
    def theta0(self):
        return self.round(np.arctan2(self.sin_t0, self.cos_t0))
