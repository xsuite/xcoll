# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import json

import xobjects as xo
import xcoll as xc


def test_closed_polygon():
    kernels = _create_geometry_kernel()
    with open(xc._pkg_root.parent / "tests" / "data" / "geometry.json", "r") as fp:
        expected_s = json.load(fp)['expected_s']
    for tilt in expected_s.keys():
        poly_s, poly_x, _, _ = _generate_polygon_points(8, tilt_L=int(tilt))
        for part_ang in expected_s[tilt].keys():
            part_tan = np.tan(np.deg2rad(int(part_ang)))
            for part_x_cm in range(-60, 80):
                part_x = part_x_cm/100
                s = kernels.check_poly(part_x=part_x, part_tan=part_tan, poly_s=poly_s,
                                       poly_x=poly_x, num_polys=len(poly_s), is_closed=True)
                if s > 1.e10:
                    s = None
                    assert part_x_cm not in expected_s[tilt][part_ang]
                else:
                    assert np.isclose(s, expected_s[tilt][part_ang][str(part_x_cm)])
                # _plot_poly(part_x, part_tan, poly_s, poly_x, s=s)   # Visual inspection of test validity

def test_open_polygon():
    kernels = _create_geometry_kernel()
    # with open(xc._pkg_root.parent / "tests" / "data" / "geometry.json", "r") as fp:
    #     expected_s = json.load(fp)['expected_s_open']
    expected_s = {}
    for tilt_L in expected_s.keys():
        for tilt_R in expected_s[tilt_L].keys():
            poly_s_L, poly_x_L, poly_s_R, poly_x_R = _generate_polygon_points(8, tilt_L=int(tilt_L), tilt_R=int(tilt_R))
            for part_ang in expected_s[tilt_L][tilt_R].keys():
                part_tan = np.tan(np.deg2rad(int(part_ang)))
                for part_x_cm in range(-90, 90):
                    part_x = part_x_cm/100
                    s_L = kernels.check_open_poly(part_x=part_x, part_tan=part_tan, poly_s=poly_s_L, poly_x=poly_x_L,
                                                  num_polys=len(poly_s_L), tan_tilt=np.tan(np.deg2rad(tilt_L)), side=1)
                    s_R = kernels.check_open_poly(part_x=part_x, part_tan=part_tan, poly_s=poly_s_R, poly_x=poly_x_R,
                                                  num_polys=len(poly_s_R), tan_tilt=np.tan(np.deg2rad(tilt_R)), side=-1)
                    s = min(s_L, s_R)
                    if s > 1.e10:
                        s = None
                        assert part_x_cm not in expected_s[tilt_L][tilt_R][part_ang]
                    else:
                        assert np.isclose(s, expected_s[tilt_L][tilt_R][part_ang][str(part_x_cm)])
                    # _plot_poly(part_x, part_tan, poly_s, poly_x, s=s)   # Visual inspection of test validity


def _create_geometry_kernel():
    src_poly = xc._pkg_root / 'scattering_routines' / 'geometry' / 'polygon.h'
    kernels_poly = {
        'check_poly': xo.Kernel(
                c_name='get_s_of_first_crossing_with_polygon',
                args=[
                    xo.Arg(xo.Float64, pointer=False, name='part_x'),
                    xo.Arg(xo.Float64, pointer=False, name='part_tan'),
                    xo.Arg(xo.Float64, pointer=True, name='poly_s'),
                    xo.Arg(xo.Float64, pointer=True, name='poly_x'),
                    xo.Arg(xo.Int8, name='num_polys'),
                    xo.Arg(xo.Int8, name='is_closed')
                ],
                ret=xo.Arg(xo.Float64, pointer=False, name='s')),
        'check_open_poly': xo.Kernel(
                c_name='get_s_of_first_crossing_with_open_polygon',
                args=[
                    xo.Arg(xo.Float64, pointer=False, name='part_x'),
                    xo.Arg(xo.Float64, pointer=False, name='part_tan'),
                    xo.Arg(xo.Float64, pointer=True, name='poly_s'),
                    xo.Arg(xo.Float64, pointer=True, name='poly_x'),
                    xo.Arg(xo.Int8, name='num_polys'),
                    xo.Arg(xo.Float64, name='tan_tilt'),
                    xo.Arg(xo.Int8, name='side')
                ],
                ret=xo.Arg(xo.Float64, pointer=False, name='s'))
    }
    context = xo.ContextCpu()
    context.add_kernels(sources=[src_poly], kernels=kernels_poly)
    return context.kernels

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
