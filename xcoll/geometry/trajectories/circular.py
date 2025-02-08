## copyright ############################### #
## This file is part of the Xcoll package.   #
## Copyright (c) CERN, 2025.                 #
## ######################################### #

import numpy as np

import xobjects as xo

from ...general import _pkg_root
from ..c_init import GeomCInit


class CircularTrajectory(xo.Struct):
    """A trajectory describing a particle along a circular path (as with bent crystals).

    In parametrised form, it is given by:
        s(𝜆) = sR + R cos(𝜆 + 𝜃I)     𝜆  = -π..π
        x(𝜆) = xR + R sin(𝜆 + 𝜃I)
        𝜃(𝜆) = 𝜃I + 𝜆 + chan. effects

    where (sR, xR) is the centre of the bend, R the distance of the particle to the centre
    R = sqrt((s0-sR)^2 + (x0-xR)^2), and 𝜃I the particle's polar angle wrt the centre
    𝜃I = tan-1((x0-xR)/(s0-sR)).
    """

    R = xo.Float64
    sR = xo.Float64
    xR = xo.Float64
    sin_tI = xo.Float64
    cos_tI = xo.Float64
    tan_tI = xo.Float64

    _depends_on = [GeomCInit]
    _extra_c_sources = [_pkg_root / 'geometry' / 'trajectories' / 'circular.h']

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

