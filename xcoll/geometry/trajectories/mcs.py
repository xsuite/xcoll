## copyright ############################### #
## This file is part of the Xcoll package.   #
## Copyright (c) CERN, 2024.                 #
## ######################################### #

import numpy as np

import xobjects as xo
from ...general import _pkg_root
from ..c_init import GeomCInit

class MultipleCoulombTrajectory(xo.Struct):
    ''' Wrapper that holds C code and API for a multiple coulomb scattering trajectory'''

    _depends_on = [GeomCInit]
    _extra_c_sources = [_pkg_root / 'scattering_routines' / 'geometry' / 'trajectories' / 'mcs.h']

    args_hv = [
            # The arguments that define the particle trajectory, common to both planes
            xo.Arg(xo.Float64, pointer=False, name="s0"),  # Particle s
    ]
    args_h = [
            # The arguments that define the horizontal (after rotation) particle trajectory
            xo.Arg(xo.Float64, pointer=False, name="x0"),  # Particle x
            xo.Arg(xo.Float64, pointer=False, name="xm")   # Particle slope in the x direction (= xp = tan(theta_x))
    ]
    args_v = [
            # The arguments that define the vertical (after rotation) particle trajectory
            xo.Arg(xo.Float64, pointer=False, name="y0"),  # Particle y
            xo.Arg(xo.Float64, pointer=False, name="ym")   # Particle slope in the y direction (= yp = tan(theta_y))
    ]

