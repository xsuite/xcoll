## copyright ############################### #
## This file is part of the Xcoll package.   #
## Copyright (c) CERN, 2025.                 #
## ######################################### #

import numpy as np

import xobjects as xo

from ...general import _pkg_root


class CircularTrajectory(xo.Struct):
    """A trajectory describing a particle along a circular path (as with bent crystals).

    In parametrised form, it is given by:
        s(𝜆) = sR + R cos(𝜆 + 𝜃I)     𝜆  = -π..π
        x(𝜆) = xR + R sin(𝜆 + 𝜃I)
        𝜃(𝜆) = 𝜃I + 𝜆 + chan. effects

    where (sR, xR) is the centre of the bend, R the distance of the particle to the centre
    R = sqrt((s0-sR)^2 + (x0-xR)^2), and 𝜃I the particle's polar angle wrt the centre
    𝜃I = tan-1((x0-xR)/(s0-sR)).

    In practice, we do not provide R nor 𝜃I, but they are calculated from (s0, x0).
    """

    R = xo.Float64
    sR = xo.Float64
    xR = xo.Float64
    sin_tI = xo.Float64
    cos_tI = xo.Float64
    tan_tI = xo.Float64

    _extra_c_sources = [_pkg_root / 'geometry' / 'trajectories' / 'circular.h']

    _kernels = {'set_params': xo.Kernel(
                                c_name='CircularTrajectory_set_params',
                                args=[xo.Arg(xo.ThisClass, name="traj"),
                                      xo.Arg(xo.Float64, name="sR"),
                                      xo.Arg(xo.Float64, name="xR"),
                                      xo.Arg(xo.Float64, name="s0"),
                                      xo.Arg(xo.Float64, name="x0")],
                                ret=xo.Float64)}

    def __init__(self, *args, **kwargs):
        s0 = False
        x0 = False
        if 's0' in kwargs and 'x0' in kwargs:
            s0 = kwargs.pop('s0')
            x0 = kwargs.pop('x0')
        super().__init__(*args, **kwargs)
        if s0 is not False and x0 is not False:
            self.set_initial_angle(s0, x0)

    def __str__(self):
        return f"CircularTrajectory(R={self.R}, sR={self.sR}, xR={self.xR}, tI={self.tI})"

    @property
    def tI(self):
        return self.round(np.arctan2(self.tan_tI))

    @tI.setter
    def tI(self, val):
        self.tan_tI = np.tan(val)
        self.sin_tI = np.sin(self.tI)
        self.cos_tI = np.cos(self.tI)

    def set_initial_angle(self, s0, x0):
        R = np.sqrt((s0-self.sR)**2 + (x0-self.xR)**2)
        self.R = R
        self.tan_tI = (x0-self.xR) / (s0-self.sR)
        self.sin_tI = (x0-self.xR) / R
        self.cos_tI = (s0-self.sR) / R

#     args_hv = [
#             # The arguments that define the particle trajectory, common to both planes
#             xo.Arg(xo.Float64, pointer=False, name="s0"),  # Particle s
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

