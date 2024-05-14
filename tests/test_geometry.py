# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import json

import xobjects as xo
import xcoll as xc

sys.path.insert(1, (xc._pkg_root.parent / 'tests' / 'data').as_posix())
from xcoll_geometry_test import XcollGeometryTest


def _init_kernels():
    geom = XcollGeometryTest()
    geom.compile_kernels(only_if_needed=True)#, save_source_as='geometry.c')
    return geom._context.kernels

kernels = _init_kernels()


def test_jaw():
    with open(xc._pkg_root.parent / "tests" / "data" / "geometry.json", "r") as fp:
        expected_s = json.load(fp)['expected_s_jaw']
    for tilt_L in expected_s.keys():
        tilt_tan_L = np.tan(np.deg2rad(int(tilt_L)))
        for tilt_R in expected_s[tilt_L].keys():
            tilt_tan_R = np.tan(np.deg2rad(int(tilt_R)))
            s_poly_L, x_poly_L, s_poly_R, x_poly_R = _generate_polygon_points(8, tilt_L=int(tilt_L), tilt_R=int(tilt_R))
            for part_ang in expected_s[tilt_L][tilt_R].keys():
                part_tan_x = np.tan(np.deg2rad(int(part_ang)))
                for part_x_cm in range(-90, 90):
                    part_x = part_x_cm/100
                    s_L = kernels.test_jaw(part_x=part_x, part_tan_x=part_tan_x, s_U=s_poly_L[1], x_U=x_poly_L[1],
                                           s_D=s_poly_L[2], x_D=x_poly_L[2], tilt_tan=tilt_tan_L, side=1)
                    s_R = kernels.test_jaw(part_x=part_x, part_tan_x=part_tan_x, s_U=s_poly_R[1], x_U=x_poly_R[1],
                                           s_D=s_poly_R[2], x_D=x_poly_R[2], tilt_tan=tilt_tan_R, side=-1)
                    s = min(s_L, s_R)
                    if s > 1.e10:
                        s = None
                        assert part_x_cm not in expected_s[tilt_L][tilt_R][part_ang]
                    else:
                        assert np.isclose(s, expected_s[tilt_L][tilt_R][part_ang][str(part_x_cm)])


def test_jaw_after_s():
    with open(xc._pkg_root.parent / "tests" / "data" / "geometry.json", "r") as fp:
        expected_s = json.load(fp)['expected_s_jaw_after_s']
    for tilt_L in expected_s.keys():
        tilt_tan_L = np.tan(np.deg2rad(int(tilt_L)))
        for tilt_R in expected_s[tilt_L].keys():
            tilt_tan_R = np.tan(np.deg2rad(int(tilt_R)))
            s_poly_L, x_poly_L, s_poly_R, x_poly_R = _generate_polygon_points(8, tilt_L=int(tilt_L), tilt_R=int(tilt_R))
            for part_ang in expected_s[tilt_L][tilt_R].keys():
                part_tan_x = np.tan(np.deg2rad(int(part_ang)))
                for part_x_cm in range(-90, 90):
                    part_x = part_x_cm/100
                    s_L = kernels.test_jaw_after_s(part_x=part_x, part_tan_x=part_tan_x, s_U=s_poly_L[1], x_U=x_poly_L[1],
                                           s_D=s_poly_L[2], x_D=x_poly_L[2], tilt_tan=tilt_tan_L, side=1, current_s=0.7)
                    s_R = kernels.test_jaw_after_s(part_x=part_x, part_tan_x=part_tan_x, s_U=s_poly_R[1], x_U=x_poly_R[1],
                                           s_D=s_poly_R[2], x_D=x_poly_R[2], tilt_tan=tilt_tan_R, side=-1, current_s=0.7)
                    s = min(s_L, s_R)
                    if s > 1.e10:
                        s = None
                        assert part_x_cm not in expected_s[tilt_L][tilt_R][part_ang]
                    else:
                        assert np.isclose(s, expected_s[tilt_L][tilt_R][part_ang][str(part_x_cm)])


def test_polygon():
    kernels = _init_kernels()
    with open(xc._pkg_root.parent / "tests" / "data" / "geometry.json", "r") as fp:
        expected_s = json.load(fp)['expected_s_polygon']
    for tilt in expected_s.keys():
        s_poly, x_poly, _, _ = _generate_polygon_points(8, tilt_L=int(tilt))
        for part_ang in expected_s[tilt].keys():
            part_tan_x = np.tan(np.deg2rad(int(part_ang)))
            for part_x_cm in range(-60, 80):
                part_x = part_x_cm/100
                s = kernels.test_polygon(part_x=part_x, part_tan_x=part_tan_x, s_poly=s_poly,
                                         x_poly=x_poly, num_polys=len(s_poly))
                if s > 1.e10:
                    s = None
                    assert part_x_cm not in expected_s[tilt][part_ang]
                else:
                    assert np.isclose(s, expected_s[tilt][part_ang][str(part_x_cm)])
                # _plot_poly(part_x, part_tan_x, s_poly, x_poly, s=s)   # Visual inspection of test validity


def test_polygon_after_s():
    kernels = _init_kernels()
    with open(xc._pkg_root.parent / "tests" / "data" / "geometry.json", "r") as fp:
        expected_s = json.load(fp)['expected_s_polygon_after_s']
    for tilt in expected_s.keys():
        s_poly, x_poly, _, _ = _generate_polygon_points(8, tilt_L=int(tilt))
        for part_ang in expected_s[tilt].keys():
            part_tan_x = np.tan(np.deg2rad(int(part_ang)))
            for part_x_cm in range(-60, 80):
                part_x = part_x_cm/100
                s = kernels.test_polygon_after_s(part_x=part_x, part_tan_x=part_tan_x, s_poly=s_poly,
                                                 x_poly=x_poly, num_polys=len(s_poly), current_s=0.7)
                if s > 1.e10:
                    s = None
                    assert part_x_cm not in expected_s[tilt][part_ang]
                else:
                    assert np.isclose(s, expected_s[tilt][part_ang][str(part_x_cm)])
                # _plot_poly(part_x, part_tan_x, s_poly, x_poly, s=s)   # Visual inspection of test validity


def test_open_polygon():
    kernels = _init_kernels()
    with open(xc._pkg_root.parent / "tests" / "data" / "geometry.json", "r") as fp:
        expected_s = json.load(fp)['expected_s_open_polygon']
    for tilt_L in expected_s.keys():
        tilt_tan_L = np.tan(np.deg2rad(int(tilt_L)))
        for tilt_R in expected_s[tilt_L].keys():
            tilt_tan_R = np.tan(np.deg2rad(int(tilt_R)))
            s_poly_L, x_poly_L, s_poly_R, x_poly_R = _generate_polygon_points(8, tilt_L=int(tilt_L), tilt_R=int(tilt_R))
            for part_ang in expected_s[tilt_L][tilt_R].keys():
                part_tan_x = np.tan(np.deg2rad(int(part_ang)))
                for part_x_cm in range(-90, 90):
                    part_x = part_x_cm/100
                    s_L = kernels.test_open_polygon(part_x=part_x, part_tan_x=part_tan_x, s_poly=s_poly_L, x_poly=x_poly_L,
                                                    num_polys=len(s_poly_L), tilt_tan=tilt_tan_L, side=1)
                    s_R = kernels.test_open_polygon(part_x=part_x, part_tan_x=part_tan_x, s_poly=s_poly_R, x_poly=x_poly_R,
                                                    num_polys=len(s_poly_R), tilt_tan=tilt_tan_R, side=-1)
                    s = min(s_L, s_R)
                    if s > 1.e10:
                        s = None
                        assert part_x_cm not in expected_s[tilt_L][tilt_R][part_ang]
                    else:
                        assert np.isclose(s, expected_s[tilt_L][tilt_R][part_ang][str(part_x_cm)])
                    # _plot_poly(part_x, part_tan_x, s_poly, x_poly, s=s)   # Visual inspection of test validity


def test_open_polygon_after_s():
    kernels = _init_kernels()
    with open(xc._pkg_root.parent / "tests" / "data" / "geometry.json", "r") as fp:
        expected_s = json.load(fp)['expected_s_open_polygon_after_s']
    for tilt_L in expected_s.keys():
        tilt_tan_L = np.tan(np.deg2rad(int(tilt_L)))
        for tilt_R in expected_s[tilt_L].keys():
            tilt_tan_R = np.tan(np.deg2rad(int(tilt_R)))
            s_poly_L, x_poly_L, s_poly_R, x_poly_R = _generate_polygon_points(8, tilt_L=int(tilt_L), tilt_R=int(tilt_R))
            for part_ang in expected_s[tilt_L][tilt_R].keys():
                part_tan_x = np.tan(np.deg2rad(int(part_ang)))
                for part_x_cm in range(-90, 90):
                    part_x = part_x_cm/100
                    s_L = kernels.test_open_polygon_after_s(part_x=part_x, part_tan_x=part_tan_x, s_poly=s_poly_L, x_poly=x_poly_L,
                                                    num_polys=len(s_poly_L), tilt_tan=tilt_tan_L, side=1, current_s=0.7)
                    s_R = kernels.test_open_polygon_after_s(part_x=part_x, part_tan_x=part_tan_x, s_poly=s_poly_R, x_poly=x_poly_R,
                                                    num_polys=len(s_poly_R), tilt_tan=tilt_tan_R, side=-1, current_s=0.7)
                    s = min(s_L, s_R)
                    if s > 1.e10:
                        s = None
                        assert part_x_cm not in expected_s[tilt_L][tilt_R][part_ang]
                    else:
                        assert np.isclose(s, expected_s[tilt_L][tilt_R][part_ang][str(part_x_cm)])
                    # _plot_poly(part_x, part_tan_x, s_poly, x_poly, s=s)   # Visual inspection of test validity


def test_crystal():
    kernels = _init_kernels()
    with open(xc._pkg_root.parent / "tests" / "data" / "geometry.json", "r") as fp:
        expected_s = json.load(fp)['expected_s_crystal']
    for tilt in expected_s.keys():
        tilt_sin = np.sin(np.deg2rad(int(tilt)))
        tilt_cos = np.cos(np.deg2rad(int(tilt)))
        for R in expected_s[tilt].keys():
            for part_ang in expected_s[tilt][R].keys():
                part_tan_x = np.tan(np.deg2rad(int(part_ang)))
                for part_x_cm in range(-60, 80):
                    part_x = part_x_cm/100
                    s = kernels.test_crystal(part_x=part_x, part_tan_x=part_tan_x, R=float(R), width=0.1476,
                                      length=1.6597, jaw_U=0.112, tilt_sin=tilt_sin, tilt_cos=tilt_cos)
                    if s > 1.e10:
                        s = None
                        assert part_x_cm not in expected_s[tilt][R][part_ang]
                    else:
                        assert np.isclose(s, expected_s[tilt][R][part_ang][str(part_x_cm)])


def test_crystal_after_s():
    kernels = _init_kernels()
    with open(xc._pkg_root.parent / "tests" / "data" / "geometry.json", "r") as fp:
        expected_s = json.load(fp)['expected_s_crystal_after_s']
    for tilt in expected_s.keys():
        tilt_sin = np.sin(np.deg2rad(int(tilt)))
        tilt_cos = np.cos(np.deg2rad(int(tilt)))
        for R in expected_s[tilt].keys():
            for part_ang in expected_s[tilt][R].keys():
                part_tan_x = np.tan(np.deg2rad(int(part_ang)))
                for part_x_cm in range(-60, 80):
                    part_x = part_x_cm/100
                    s = kernels.test_crystal_after_s(part_x=part_x, part_tan_x=part_tan_x, R=float(R), width=0.1476,
                                      length=1.6597, jaw_U=0.112, tilt_sin=tilt_sin, tilt_cos=tilt_cos, current_s=0.7)
                    if s > 1.e10:
                        s = None
                        assert part_x_cm not in expected_s[tilt][R][part_ang]
                    else:
                        assert np.isclose(s, expected_s[tilt][R][part_ang][str(part_x_cm)])


def _generate_polygon_points(num_poly, tilt_L=0, tilt_R=0):
    with open(xc._pkg_root.parent / "tests" / "data" / "geometry.json", "r") as fp:
        rans = json.load(fp)['rans']
    len_between = num_poly-4
    between = [[(i+1)/(len_between+1) + rans[i]*0.15-0.075, rans[i+len_between]*0.15+0.025]
               for i in range(len_between)]
    poly_L = [[0,0.4],  [0,0.1],  *between, [1,0.1],  [1,0.4]]
    between = [[(i+1)/(len_between+1) + rans[2*len_between+i]*0.15-0.075, -rans[3*len_between+i]*0.15-0.025]
               for i in range(len_between)]
    poly_R = [[0,-0.4], [0,-0.1], *between, [1,-0.1], [1,-0.4]]
    cos_L = np.cos(np.deg2rad(-tilt_L))
    sin_L = np.sin(np.deg2rad(-tilt_L))
    cos_R = np.cos(np.deg2rad(-tilt_R))
    sin_R = np.sin(np.deg2rad(-tilt_R))
    poly_s_L = np.array([(s-0.5)*cos_L  + (x-0.1)*sin_L + 0.5 for s,x in poly_L], dtype=np.float64)
    poly_x_L = np.array([-(s-0.5)*sin_L + (x-0.1)*cos_L + 0.1 for s,x in poly_L], dtype=np.float64)
    poly_s_R = np.array([(s-0.5)*cos_R  + (x+0.1)*sin_R + 0.5 for s,x in poly_R], dtype=np.float64)
    poly_x_R = np.array([-(s-0.5)*sin_R + (x+0.1)*cos_R - 0.1 for s,x in poly_R], dtype=np.float64)
    return poly_s_L, poly_x_L, poly_s_R, poly_x_R

def _plot_poly(part_x, part_tan, poly_s_L, poly_x_L, poly_s_R=None, poly_x_R=None, is_open=False, s=None):
    fig, ax = plt.subplots(1, 1, figsize=(8,5.6))
    if is_open:
        ax.plot(poly_s_L, poly_x_L, 'k-')
        if poly_s_R is not None and poly_x_R is not None:
            ax.plot(poly_s_R, poly_x_R, 'k-')
    else:
        ax.plot([*poly_s_L, poly_s_L[0]], [*poly_x_L, poly_x_L[0]], 'k-')
        if poly_s_R is not None and poly_x_R is not None:
            ax.plot([*poly_s_R, poly_s_R[0]], [*poly_x_R, poly_x_R[0]], 'k-')
    ax.plot([-0.5,1.5], [part_x-0.5*part_tan,part_x+1.5*part_tan], 'b-')
    ax.set_xlim((-0.5,1.5))
    ax.set_ylim((-0.7,0.7))
    if s is not None:
        ax.axvline(s, c='r', ls='--')
    plt.show()
