## copyright ############################### #
## This file is part of the Xcoll package.   #
## Copyright (c) CERN, 2025.                 #
## ######################################### #

import numpy as np

import xobjects as xo

from ...general import _pkg_root
from ..c_init import BoundingBox


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
    box = BoundingBox

    _extra_c_sources = [_pkg_root / 'geometry' / 'trajectories' / 'circular.h']

    _kernels = {'set_params': xo.Kernel(
                                c_name='CircularTrajectory_set_params',
                                args=[xo.Arg(xo.ThisClass, name="traj"),
                                      xo.Arg(xo.Float64, name="sR"),
                                      xo.Arg(xo.Float64, name="xR"),
                                      xo.Arg(xo.Float64, name="s0"),
                                      xo.Arg(xo.Float64, name="x0")],
                                ret=None),
                'init_bounding_box': xo.Kernel(
                                c_name='CircularTrajectory_init_bounding_box',
                                args=[xo.Arg(xo.ThisClass, name="traj"),
                                      xo.Arg(BoundingBox, name="box"), 
                                      xo.Arg(xo.Float64, name="l1"),
                                      xo.Arg(xo.Float64, name="l2")], # this is not parameters of mcs??
                                ret=None)}

    def __init__(self, *args, **kwargs):
        s0 = kwargs.pop('s0', False)
        x0 = kwargs.pop('x0', False)
        l1 = kwargs.pop('l1', 0)
        l2 = kwargs.pop('l2', np.pi)
        super().__init__(*args, **kwargs)
        if s0 is not False and x0 is not False:
            self.set_params(s0=s0, x0=x0, sR=self.sR, xR=self.xR)
        self.box = BoundingBox()
        self.init_bounding_box(box=self.box, l1=l1, l2=l2)

    def __str__(self):
        return f"CircularTrajectory(R={self.R}, sR={self.sR}, xR={self.xR}, tI={self.tI})"

    @property
    def tI(self):
        return self.round(np.arctan2(self.sin_tI, self.cos_tI))
