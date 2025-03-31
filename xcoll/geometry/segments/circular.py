# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

import numpy as np

import xobjects as xo

from ...general import _pkg_root
from ..c_init import BoundingBox


class CircularSegment(xo.Struct):
    """Circular arc segment, defined by a centre and radius, and the starting/end angles defined anti-clockwise"""
    R  = xo.Float64
    sR = xo.Float64  # s-coordinate of the centre
    xR = xo.Float64  # x-coordinate of the centre
    _theta1 = xo.Float64  # Starting angle
    _theta2 = xo.Float64  # Ending angle
    box = BoundingBox

    _extra_c_sources = [_pkg_root / 'geometry' / 'segments' / 'circular.h']
    _kernels = {'init_bounding_box': xo.Kernel(
                                        c_name='CircularSegment_init_bounding_box',
                                        args=[xo.Arg(xo.ThisClass, name="seg"),
                                              xo.Arg(xo.ThisClass, name="box"),
                                              xo.Arg(xo.Float64, name="t1"),
                                              xo.Arg(xo.Float64, name="t2")], # this is not parameters of mcs??
                                        ret=None)}
    def __init__(self, **kwargs):
        # Different ways to initialise a CircularSegment:
        # 1. Centre, radius, and angles: CircularSegment(R=..., sR=..., xR=..., theta1=..., theta2=...)
        # 2. Start point, radius, and angles: CircularSegment(R=..., s1=..., x1=..., theta1=..., theta2=...)
        # 3. End point, radius, and angles: CircularSegment(R=..., s2=..., x2=..., theta1=..., theta2=...)
        # 4. Start and end point, and the (possibly negative) curvature: CircularSegment(k=..., s1=..., x1=..., s2=..., x2=...)
        # 5. Centre, start point, and an angular shift: CircularSegment(sR=..., xR=..., s1=..., x1=..., delta_theta=...)
        # 6. Centre, end point, and an angular shift: CircularSegment(sR=..., xR=..., s2=..., x2=..., delta_theta=...)
        if 's1' in kwargs and 'x1' in kwargs and 's2' in kwargs and 'x2' in kwargs:
            if 'sR' in kwargs or 'xR' in kwargs:
                raise ValueError("Cannot provide start point, end point, and centre!")
            if 'theta1' in kwargs or 'theta2' in kwargs or 'delta_theta' in kwargs:
                raise ValueError("Cannot provide angles when providing start and end point!")
            if 'R' in kwargs:
                raise ValueError("Cannot provide radius when providing start and end point! "
                               + "Please use curvature 'k' = 1/R instead.")
            if 'k' not in kwargs:
                raise ValueError("Curvature 'k' must be provided when providing start and end point!")
            s1 = kwargs.pop('s1')
            x1 = kwargs.pop('x1')
            s2 = kwargs.pop('s2')
            x2 = kwargs.pop('x2')
            k = kwargs.pop('k')
            if np.isclose(k, 0.):
                raise ValueError("Curvature 'k' must be non-zero!")
            elif abs(k) > 1.:
                raise ValueError("Curvature 'k' must be in [-1, 1]!")
            m = np.sqrt((s1-s2)**2 + (x1-x2)**2)
            R = abs(1/k*m/2)
            ang_m = np.arctan2(x2-x1, s2-s1)
            gamma = np.arccos(1 - 2*k*k)
            if k > 0:
                theta1 = ang_m - np.pi/2 - gamma/2
                theta2 = ang_m - np.pi/2 + gamma/2
                kwargs['sR'] = s1 - R*np.cos(theta1)
                kwargs['xR'] = x1 - R*np.sin(theta1)
            else:
                theta1 = ang_m + np.pi/2 - gamma/2
                theta2 = ang_m + np.pi/2 + gamma/2
                kwargs['sR'] = s2 - R*np.cos(theta1)
                kwargs['xR'] = x2 - R*np.sin(theta1)
            kwargs['R'] = R
        elif 'delta_theta' in kwargs:
            if 'sR' not in kwargs or 'xR' not in kwargs:
                raise ValueError("Centre must be provided when providing angular shift!")
            sR = kwargs['sR']
            xR = kwargs['xR']
            if 's1' in kwargs or 'x1' in kwargs:
                if 's1' not in kwargs or 'x1' not in kwargs:
                    raise ValueError("If the start point is provided, both 's1' and 'x1' must be provided!")
                s1 = kwargs.pop('s1')
                x1 = kwargs.pop('x1')
                kwargs['R'] = np.sqrt((s1-sR)**2 + (x1-xR)**2)
                theta1 = np.arctan2(x1-xR, s1-sR)
                theta2 = theta1 + kwargs.pop('delta_theta')
            elif 's2' in kwargs or 'x2' in kwargs:
                if 's2' not in kwargs or 'x2' not in kwargs:
                    raise ValueError("If the end point is provided, both 's2' and 'x2' must be provided!")
                s2 = kwargs.pop('s2')
                x2 = kwargs.pop('x2')
                kwargs['R'] = np.sqrt((s2-sR)**2 + (x2-xR)**2)
                theta2 = np.arctan2(x2-xR, s2-sR)
                theta1 = theta2 - kwargs.pop('delta_theta')
        else:
            if 'R' not in kwargs:
                raise ValueError("Must provide radius, curvature, or 'delta_theta'!")
            if kwargs['R'] <= 0:
                raise ValueError("Radius must be strictly positive!")
            theta1 = kwargs.pop('theta1', -np.pi)
            theta2 = kwargs.pop('theta2', np.pi)
            if 'sR' in kwargs or 'xR' in kwargs:
                if 'sR' not in kwargs or 'xR' not in kwargs:
                    raise ValueError("If the centre is provided, both 'sR' and 'xR' must be provided!")
                if 's1' in kwargs or 'x1' in kwargs or 's2' in kwargs or 'x2' in kwargs:
                    raise ValueError("Centre and start/end points cannot be provided together!")
            elif 's1' in kwargs or 'x1' in kwargs:
                if 's1' not in kwargs or 'x1' not in kwargs:
                    raise ValueError("If the start point is provided, both 's1' and 'x1' must be provided!")
                if 's2' in kwargs or 'x2' in kwargs:
                    raise ValueError("Start and end points cannot be provided together when 'R' is provided!")
                s1 = kwargs.pop('s1')
                x1 = kwargs.pop('x1')
                kwargs['sR'] = s1 - kwargs['R']*np.cos(theta1)
                kwargs['xR'] = x1 - kwargs['R']*np.sin(theta1)
            elif 's2' in kwargs or 'x2' in kwargs:
                if 's2' not in kwargs or 'x2' not in kwargs:
                    raise ValueError("If the end point is provided, both 's2' and 'x2' must be provided!")
                s2 = kwargs.pop('s2')
                x2 = kwargs.pop('x2')
                kwargs['sR'] = s2 - kwargs['R']*np.cos(theta2)
                kwargs['xR'] = x2 - kwargs['R']*np.sin(theta2)
            else:
                raise ValueError("Must provide centre, start point, or end point when providing radius!")
        t1 = kwargs.pop('t1', 0.)
        t2 = kwargs.pop('t2', 2*np.pi)
        super().__init__(**kwargs)
        self.set_angles(theta1, theta2)
        self.box = BoundingBox()
        self.init_bounding_box(box=self.box, t1=t1, t2=t2)

    def __str__(self):
        p1, p2 = self.get_vertices()
        return f"CircularSegment(({p1[0]:.3}, {p1[1]:.3})-b-({np.rad2deg(self.theta1):.0f}" + u'\xb0' \
             + f":{np.rad2deg(self.theta2):.0f}" + u'\xb0' + f":{self.R:.3})-b-({p2[0]:.3}, {p2[1]:.3}))"

    def get_vertices(self):
        return (self.s1, self.x1), (self.s2, self.x2)

    def get_control_points(self):
        return (self.sR, self.xR),

    def _translate_inplace(self, ds, dx):
        self.sR += ds
        self.xR += dx

    def _rotate_inplace(self, ps, px, angle):
        c = np.cos(angle)
        s = np.sin(angle)
        self._translate_inplace(-ps, -px)
        new_sR = self.sR * c - self.xR * s
        new_xR = self.sR * s + self.xR * c
        self.sR = new_sR
        self.xR = new_xR
        self.set_angles(self.theta1 + angle, self.theta2 + angle)
        self._translate_inplace(ps, px)

    @property
    def s1(self):
        return self.round(self.sR + self.R*np.cos(self.theta1))

    @property
    def x1(self):
        return self.round(self.xR + self.R*np.sin(self.theta1))

    @property
    def s2(self):
        return self.round(self.sR + self.R*np.cos(self.theta2))

    @property
    def x2(self):
        return self.round(self.xR + self.R*np.sin(self.theta2))

    @property
    def theta1(self):
        return self._theta1

    @property
    def theta2(self):
        # We want to represent angles in [-pi, pi]
        value = self._theta2
        while value > np.pi:
            value -= 2*np.pi
        return value

    def set_angles(self, theta1, theta2):
        while theta1 < -np.pi:
            theta1 += 2*np.pi
        while theta1 > np.pi:
            theta1 -= 2*np.pi
        while theta2 < -np.pi:
            theta2 += 2*np.pi
        while theta2 > np.pi:
            theta2 -= 2*np.pi
        self._theta1 = theta1
        if theta2 >= theta1:
            self._theta2 = theta2
        else:
            # If theta2 is smaller than theta1, it means we have done a full turn
            self._theta2 = theta2 + 2*np.pi
