## copyright ############################### #
## This file is part of the Xcoll package.   #
## Copyright (c) CERN, 2025.                 #
## ######################################### #

import numpy as np

import xobjects as xo

from ...general import _pkg_root
from ..c_init import GeomCInit


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
        𝛺(𝜆) = (13.6 MeV) / (pc) sqrt(𝜆 q^2 / (X0 𝛽^2)) (1 + 0.038 ln (𝜆 q^2 / (X0 𝛽^2)))

    where X0 is the material's radiation length, and q, 𝛽, and pc are the particle's charge,
    relativistic 𝛽, and momentum.
    """

    s0 = xo.Float64
    x0 = xo.Float64
    sin_t0 = xo.Float64
    cos_t0 = xo.Float64
    tan_t0 = xo.Float64
    Xt0 = xo.Float64  #  X0 𝛽^2 / q^2
    A0 = xo.Float64   # (𝜉1/√12 + 𝜉2/2) (13.6 MeV) / (pc)
    B0 = xo.Float64   # 𝜉2 (13.6 MeV) / (pc)

    _depends_on = [GeomCInit]
    _extra_c_sources = [_pkg_root / 'geometry' / 'trajectories' / 'mcs.h']

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

