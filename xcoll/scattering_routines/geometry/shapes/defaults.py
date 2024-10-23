# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

import numpy as np

import xobjects as xo

from ..c_init import XC_EPSILON
from ..segments import *
from .shape import Shape2D, Shape2DV


def create_jaw(s_U, x_U, s_D, x_D, *, side):
    tilt = np.arctan((x_D - x_U)/(s_D - s_U))
    seg1 = HalfOpenLineSegment(s=s_U, x=x_U, t=side*np.pi/2 + tilt)
    seg1 = LineSegment(s1=s_U, x1=x_U, s2=s_D, x2=x_D)
    seg3 = HalfOpenLineSegment(s=s_D, x=x_D, t=side*np.pi/2 + tilt)
    return Shape2D([seg1, seg2, seg3])


def create_polygon(s_poly, x_poly):
    if len(s_poly) != len(x_poly):
        raise ValueError('s_poly and x_poly must have the same length')
    s_poly_end = [*s_poly[1:], s_poly[0]]
    x_poly_end = [*x_poly[1:], x_poly[0]]
    segs = [LineSegment(s1=s1, x1=x1, s2=s2, x2=x2)
            for s1, x1, s2, x2 in zip(s_poly, x_poly, s_poly_end, x_poly_end)]
    return Shape2D(segs)


def create_polygon_jaw(s_poly, x_poly, *, tilt, side):
    begin_seg = HalfOpenLineSegment(s=s_poly[0], x=x_poly[0], t=side*np.pi/2 + tilt)
    mid_segs  = [LineSegment(s1=s1, x1=x1, s2=s2, x2=x2)
                 for s1, x1, s2, x2 in zip(s_poly[:-1], x_poly[:-1], s_poly[1:], x_poly[1:])]
    end_seg = HalfOpenLineSegment(s=s_poly[-1], x=x_poly[-1], t=side*np.pi/2 + tilt)
    return Shape2D([begin_seg, *mid_segs, end_seg])


def create_crystal(R, *, width, height, length, jaw_U, tilt):
    # First corner is what defines the crystal position
    A_s = 0;
    A_x = jaw_U;

    # Manipulate R in function of sign
    sgnR = (R > 0) - (R < 0)
    R_short  = sgnR*(abs(R) - width)
    sin_a = length/abs(R)
    cos_a = sqrt(1 - length*length/R/R)
    if (abs(R) < XC_EPSILON):
        raise ValueError("Straight crystal not yet implemented!")

    elif (R < 0):
        # This distinction is needed to keep the crystal at the same location when changing the bend direction
        R_temp = R_short
        R_short = R
        R = R_temp

    # Bending centre is defined w.r.t. A
    R_s = A_s - R*np.sin(tilt)
    R_x = A_x + R*np.cos(tilt)

    # Three remaining corner points of crystal
    B_s = R_s + R_short*np.sin(tilt)
    C_s = R_s + abs(R_short)*sin_a*np.cos(tilt) + R_short*cos_a*np.sin(tilt)
    D_s = R_s + abs(R)*sin_a*np.cos(tilt) + R*cos_a*np.sin(tilt)
    B_x = R_x - R_short*np.cos(tilt)
    C_x = R_x - cos_a*np.cos(tilt)*R_short + sin_a*np.sin(tilt)*abs(R_short)
    D_x = R_x - cos_a*np.cos(tilt)*R + sin_a*np.sin(tilt)*fabs(R)
    A_t = np.arctan2(A_x - R_x, A_s - R_s)
    D_t = np.arctan2(D_x - R_x, D_s - R_s)
    t1 = min(A_t, D_t)
    t2 = max(A_t, D_t)

    # Fill segments
    seg1 = LineSegment(s1=A_s, x1=A_x, s2=B_s, x2=B_x)
    seg2 = CircularSegment(R=R, s=R_s, x=R_x, t1=t1, t2=t2)
    seg3 = LineSegment(s1=C_s, x1=C_x, s2=D_s, x2=D_x)
    seg4 = CircularSegment(R=R_short, s=R_s, x=R_x, t1=t1, t2=t2)

    return Shape2DV([seg1, seg2, seg3, seg4], vlimit=[-height/2, height/2])
