# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import json

import xobjects as xo
import xcoll as xc

sys.path.insert(1, (xc._pkg_root.parent / 'tests').as_posix())
from test_geometry import _init_kernels, _generate_polygon_points


all_tilt_L = ['-89', '-63', '-37', '-24', '-13', '0', '11', '29', '44', '71', '89']
all_tilt_R = ['-88', '-68', '-32', '-28', '-12', '0', '9', '25', '50', '70', '89']
all_tilt_LR = ['[-89, 25]', '[-53, 0]', '[-24, -28]', '[-13, 9]', '[0, 0]', 
               '[11, -12]', '[29, 25]', '[0, 50]', '[71, -32]', '[89, -68]']
all_part_ang_x = ['-72', '-33', '-9', '0', '11', '42', '68']
all_part_ang_y = ['-70', '-38', '-8', '0', '14', '29', '69']
all_part_cm = np.array(['-90', '-89', '-88', '-87', '-86', '-85', '-84', '-83', '-82', '-81', '-80', '-79', '-78', '-77',
               '-76', '-75', '-74', '-73', '-72', '-71', '-70', '-69', '-68', '-67', '-66', '-65', '-64', '-63',
               '-62', '-61', '-60', '-59', '-58', '-57', '-56', '-55', '-54', '-53', '-52', '-51', '-50', '-49',
               '-48', '-47', '-46', '-45', '-44', '-43', '-42', '-41', '-40', '-39', '-38', '-37', '-36', '-35',
               '-34', '-33', '-32', '-31', '-30', '-29', '-28', '-27', '-26', '-25', '-24', '-23', '-22', '-21',
               '-20', '-19', '-18', '-17', '-16', '-15', '-14', '-13', '-12', '-11', '-10', '-9', '-8', '-7', '-6',
               '-5', '-4', '-3', '-2', '-1', '0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12',
               '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23', '24', '25', '26', '27', '28',
               '29', '30', '31', '32', '33', '34', '35', '36', '37', '38', '39', '40', '41', '42', '43', '44',
               '45', '46', '47', '48', '49', '50', '51', '52', '53', '54', '55', '56', '57', '58', '59', '60',
               '61', '62', '63', '64', '65', '66', '67', '68', '69', '70', '71', '72', '73', '74', '75', '76',
               '77', '78', '79', '80', '81', '82', '83', '84', '85', '86', '87', '88', '89', '90'])
all_R = ['-10', '-9', '-8', '-7', '-6', '-5', '-4', '-3', '-2', '-1', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10']


def _loop_1jaw_1partdim(name, func, num_polys):
    with open(xc._pkg_root.parent / "tests" / "data" / "geometry.json", "r") as fp:
        data = json.load(fp)
    expected_s = {}
    for tilt in all_tilt_L:
        expected_s[tilt] = {}
        tilt_tan = np.tan(np.deg2rad(int(tilt)))
        s_poly, x_poly, _, _ = _generate_polygon_points(num_polys, tilt_L=int(tilt))
        for part_ang in all_part_ang_x:
            expected_s[tilt][part_ang] = {}
            part_tan_x = np.tan(np.deg2rad(int(part_ang)))
            for part_x_cm in all_part_cm[35:]:
                part_x = int(part_x_cm)/100.
                s = func(part_x, part_tan_x, None, None, s_poly, x_poly, tilt_tan, 1)
                if s < 1.e10:
                    expected_s[tilt][part_ang][part_x_cm] = s
    data[name] = expected_s
    with open(xc._pkg_root.parent / "tests" / "data" / "geometry.json", "w") as fp:
        json.dump(data, fp)

def _loop_cry_1jaw_1partdim(name, func):
    with open(xc._pkg_root.parent / "tests" / "data" / "geometry.json", "r") as fp:
        data = json.load(fp)
    expected_s = {}
    for tilt in all_tilt_L:
        expected_s[tilt] = {}
        tilt_sin = np.sin(np.deg2rad(int(tilt)))
        tilt_cos = np.cos(np.deg2rad(int(tilt)))
        for R in all_R:
            expected_s[tilt][R] = {}
            for part_ang in all_part_ang_x:
                expected_s[tilt][R][part_ang] = {}
                part_tan_x = np.tan(np.deg2rad(int(part_ang)))
                for part_x_cm in all_part_cm[35:]:
                    part_x = int(part_x_cm)/100.
                    s = func(part_x, part_tan_x, None, None, float(R), tilt_sin, tilt_cos)
                    if s < 1.e10:
                        expected_s[tilt][R][part_ang][part_x_cm] = s
    data[name] = expected_s
    with open(xc._pkg_root.parent / "tests" / "data" / "geometry.json", "w") as fp:
        json.dump(data, fp)

def _loop_1jaw_2partdim(name, func, num_polys):
    with open(xc._pkg_root.parent / "tests" / "data" / "geometry.json", "r") as fp:
        data = json.load(fp)
    expected_s = {}
    for tilt in all_tilt_L:
        expected_s[tilt] = {}
        tilt_tan = np.tan(np.deg2rad(int(tilt)))
        s_poly, x_poly, _, _ = _generate_polygon_points(num_polys, tilt_L=int(tilt))
        for part_ang_x in all_part_ang_x:
            expected_s[tilt][part_ang_x] = {}
            part_tan_x = np.tan(np.deg2rad(int(part_ang_x)))
            for part_ang_y in all_part_ang_y:
                expected_s[tilt][part_ang_x][part_ang_y] = {}
                part_tan_y = np.tan(np.deg2rad(int(part_ang_y)))
                for part_x_cm in all_part_cm[35::10]:
                    expected_s[tilt][part_ang_x][part_ang_y][part_x_cm] = {}
                    part_x = int(part_x_cm)/100.
                    for part_y_cm in all_part_cm[35::10]:
                        part_y = int(part_y_cm)/100.
                        s = func(part_x, part_tan_x, part_y, part_tan_y, s_poly, x_poly, tilt_tan, 1)
                        if s < 1.e10:
                            expected_s[tilt][part_ang_x][part_ang_y][part_x_cm][part_y_cm] = s
    data[name] = expected_s
    with open(xc._pkg_root.parent / "tests" / "data" / "geometry.json", "w") as fp:
        json.dump(data, fp)

def _loop_cry_1jaw_2partdim(name, func):
    with open(xc._pkg_root.parent / "tests" / "data" / "geometry.json", "r") as fp:
        data = json.load(fp)
    expected_s = {}
    for tilt in all_tilt_L:
        expected_s[tilt] = {}
        tilt_sin = np.sin(np.deg2rad(int(tilt)))
        tilt_cos = np.cos(np.deg2rad(int(tilt)))
        for R in all_R:
            expected_s[tilt][R] = {}
            for part_ang_x in all_part_ang_x:
                expected_s[tilt][R][part_ang_x] = {}
                part_tan_x = np.tan(np.deg2rad(int(part_ang_x)))
                for part_ang_y in all_part_ang_y:
                    expected_s[tilt][R][part_ang_x][part_ang_y] = {}
                    part_tan_y = np.tan(np.deg2rad(int(part_ang_y)))
                    for part_x_cm in all_part_cm[35::10]:
                        expected_s[tilt][R][part_ang_x][part_ang_y][part_x_cm] = {}
                        part_x = int(part_x_cm)/100.
                        for part_y_cm in all_part_cm[35::10]:
                            part_y = int(part_y_cm)/100.
                            s = func(part_x, part_tan_x, part_y, part_tan_y, float(R), tilt_sin, tilt_cos)
                            if s < 1.e10:
                                expected_s[tilt][R][part_ang_x][part_ang_y][part_x_cm][part_y_cm] = s
    data[name] = expected_s
    with open(xc._pkg_root.parent / "tests" / "data" / "geometry.json", "w") as fp:
        json.dump(data, fp)

def _loop_2jaw_1partdim(name, func, num_polys):
    with open(xc._pkg_root.parent / "tests" / "data" / "geometry.json", "r") as fp:
        data = json.load(fp)
    expected_s = {}
    for tilt_L in all_tilt_L:
        expected_s[tilt_L] = {}
        tilt_tan_L = np.tan(np.deg2rad(int(tilt_L)))
        for tilt_R in all_tilt_R:
            expected_s[tilt_L][tilt_R] = {}
            tilt_tan_R = np.tan(np.deg2rad(int(tilt_R)))
            s_poly_L, x_poly_L, s_poly_R, x_poly_R = _generate_polygon_points(num_polys, tilt_L=int(tilt_L), tilt_R=int(tilt_R))
            for part_ang in all_part_ang_x:
                expected_s[tilt_L][tilt_R][part_ang] = {}
                part_tan_x = np.tan(np.deg2rad(int(part_ang)))
                for part_x_cm in all_part_cm:
                    part_x = int(part_x_cm)/100.
                    s_L = func(part_x, part_tan_x, None, None, s_poly_L, x_poly_L, tilt_tan_L, 1)
                    s_R = func(part_x, part_tan_x, None, None, s_poly_R, x_poly_R, tilt_tan_R, -1)
                    s = min(s_L, s_R)
                    if s < 1.e10:
                        expected_s[tilt_L][tilt_R][part_ang][part_x_cm] = s
    data[name] = expected_s
    with open(xc._pkg_root.parent / "tests" / "data" / "geometry.json", "w") as fp:
        json.dump(data, fp)

def _loop_doublejaw_2partdim(name, func, num_polys):
    with open(xc._pkg_root.parent / "tests" / "data" / "geometry.json", "r") as fp:
        data = json.load(fp)
    expected_s = {}
    for tilt_LR in all_tilt_LR:
        expected_s[tilt_LR] = {}
        tilt_L, tilt_R = tilt_LR.strip('][').split(', ')
        tilt_tan_L = np.tan(np.deg2rad(int(tilt_L)))
        tilt_tan_R = np.tan(np.deg2rad(int(tilt_R)))
        s_poly_L, x_poly_L, s_poly_R, x_poly_R = _generate_polygon_points(num_polys, tilt_L=int(tilt_L), tilt_R=int(tilt_R))
        for part_ang_x in all_part_ang_x:
            expected_s[tilt_LR][part_ang_x] = {}
            part_tan_x = np.tan(np.deg2rad(int(part_ang_x)))
            for part_ang_y in all_part_ang_y:
                expected_s[tilt_LR][part_ang_x][part_ang_y] = {}
                part_tan_y = np.tan(np.deg2rad(int(part_ang_y)))
                for part_x_cm in all_part_cm[::10]:
                    expected_s[tilt_LR][part_ang_x][part_ang_y][part_x_cm] = {}
                    part_x = int(part_x_cm)/100.
                    for part_y_cm in all_part_cm[::10]:
                        part_y = int(part_y_cm)/100.
                        s_L = func(part_x, part_tan_x, part_y, part_tan_y, s_poly_L, x_poly_L, tilt_tan_L, 1)
                        s_R = func(part_x, part_tan_x, part_y, part_tan_y, s_poly_R, x_poly_R, tilt_tan_R, -1)
                        s = min(s_L, s_R)
                        if s < 1.e10:
                            expected_s[tilt_LR][part_ang_x][part_ang_y][part_x_cm][part_y_cm] = s
    data[name] = expected_s
    with open(xc._pkg_root.parent / "tests" / "data" / "geometry.json", "w") as fp:
        json.dump(data, fp)


def create_jaw():
    def func(part_x, part_tan_x, part_y, part_tan_y, s_poly, x_poly, tilt_tan, side):
        return kernels.test_jaw(part_x=part_x, part_tan_x=part_tan_x, s_U=s_poly[1], x_U=x_poly[1],
                             s_D=s_poly[2], x_D=x_poly[2], tilt_tan=tilt_tan, side=side)
    _loop_2jaw_1partdim(name='expected_s_jaw', func=func, num_polys=4)

def create_jaw_after_s():
    def func(part_x, part_tan_x, part_y, part_tan_y, s_poly, x_poly, tilt_tan, side):
        return kernels.test_jaw_after_s(part_x=part_x, part_tan_x=part_tan_x, s_U=s_poly[1], x_U=x_poly[1],
                                        s_D=s_poly[2], x_D=x_poly[2], tilt_tan=tilt_tan, side=side, current_s=0.6)
    _loop_2jaw_1partdim(name='expected_s_jaw_after_s', func=func, num_polys=4)

def create_jaw_with_vlimit():
    def func(part_x, part_tan_x, part_y, part_tan_y, s_poly, x_poly, tilt_tan, side):
        return kernels.test_jaw_with_vlimit(part_x=part_x, part_tan_x=part_tan_x, part_y=part_y, part_tan_y=part_tan_y,
                                            s_U=s_poly[1], x_U=x_poly[1], s_D=s_poly[2], x_D=x_poly[2], tilt_tan=tilt_tan,
                                            side=side, y_min=-0.1, y_max=0.25)
    _loop_doublejaw_2partdim(name='expected_s_jaw_with_vlimit', func=func, num_polys=4)

def create_jaw_after_s_with_vlimit():
    def func(part_x, part_tan_x, part_y, part_tan_y, s_poly, x_poly, tilt_tan, side):
        return kernels.test_jaw_after_s_with_vlimit(part_x=part_x, part_tan_x=part_tan_x, part_y=part_y, part_tan_y=part_tan_y,
                                            s_U=s_poly[1], x_U=x_poly[1], s_D=s_poly[2], x_D=x_poly[2], tilt_tan=tilt_tan,
                                            side=side, y_min=-0.1, y_max=0.25, current_s=0.6)
    _loop_doublejaw_2partdim(name='expected_s_jaw_after_s_with_vlimit', func=func, num_polys=4)

def create_polygon():
    def func(part_x, part_tan_x, part_y, part_tan_y, s_poly, x_poly, tilt_tan, side):
        return kernels.test_polygon(part_x=part_x, part_tan_x=part_tan_x, s_poly=s_poly,
                                    x_poly=x_poly, num_polys=len(s_poly))
    _loop_1jaw_1partdim(name='expected_s_polygon', func=func, num_polys=8)

def create_polygon_after_s():
    def func(part_x, part_tan_x, part_y, part_tan_y, s_poly, x_poly, tilt_tan, side):
        return kernels.test_polygon_after_s(part_x=part_x, part_tan_x=part_tan_x, s_poly=s_poly,
                                            x_poly=x_poly, num_polys=len(s_poly), current_s=0.6)
    _loop_1jaw_1partdim(name='expected_s_polygon_after_s', func=func, num_polys=8)

def create_polygon_with_vlimit():
    def func(part_x, part_tan_x, part_y, part_tan_y, s_poly, x_poly, tilt_tan, side):
        return kernels.test_polygon_with_vlimit(part_x=part_x, part_tan_x=part_tan_x, part_y=part_y, part_tan_y=part_tan_y,
                                                s_poly=s_poly, x_poly=x_poly, num_polys=len(s_poly), y_min=-0.1, y_max=0.25)
    _loop_1jaw_2partdim(name='expected_s_polygon_with_vlimit', func=func, num_polys=8)

def create_polygon_after_s_with_vlimit():
    def func(part_x, part_tan_x, part_y, part_tan_y, s_poly, x_poly, tilt_tan, side):
        return kernels.test_polygon_after_s_with_vlimit(part_x=part_x, part_tan_x=part_tan_x, part_y=part_y, part_tan_y=part_tan_y,
                                                s_poly=s_poly, x_poly=x_poly, num_polys=len(s_poly), y_min=-0.1, y_max=0.25, current_s=0.6)
    _loop_1jaw_2partdim(name='expected_s_polygon_after_s_with_vlimit', func=func, num_polys=8)

def create_open_polygon():
    def func(part_x, part_tan_x, part_y, part_tan_y, s_poly, x_poly, tilt_tan, side):
        return kernels.test_open_polygon(part_x=part_x, part_tan_x=part_tan_x, s_poly=s_poly, x_poly=x_poly,
                                         num_polys=len(s_poly), tilt_tan=tilt_tan, side=side)
    _loop_2jaw_1partdim(name='expected_s_open_polygon', func=func, num_polys=8)

def create_open_polygon_after_s():
    def func(part_x, part_tan_x, part_y, part_tan_y, s_poly, x_poly, tilt_tan, side):
        return kernels.test_open_polygon_after_s(part_x=part_x, part_tan_x=part_tan_x, s_poly=s_poly, x_poly=x_poly,
                                                 num_polys=len(s_poly), tilt_tan=tilt_tan, side=side, current_s=0.6)
    _loop_2jaw_1partdim(name='expected_s_open_polygon_after_s', func=func, num_polys=8)

def create_open_polygon_with_vlimit():
    def func(part_x, part_tan_x, part_y, part_tan_y, s_poly, x_poly, tilt_tan, side):
        return kernels.test_open_polygon_with_vlimit(part_x=part_x, part_tan_x=part_tan_x, part_y=part_y, part_tan_y=part_tan_y,
                                                     s_poly=s_poly, x_poly=x_poly, num_polys=len(s_poly), tilt_tan=tilt_tan, side=side,
                                                     y_min=-0.1, y_max=0.25)
    _loop_doublejaw_2partdim(name='expected_s_open_polygon_with_vlimit', func=func, num_polys=8)

def create_open_polygon_after_s_with_vlimit():
    def func(part_x, part_tan_x, part_y, part_tan_y, s_poly, x_poly, tilt_tan, side):
        return kernels.test_open_polygon_after_s_with_vlimit(part_x=part_x, part_tan_x=part_tan_x, part_y=part_y, part_tan_y=part_tan_y,
                                                     s_poly=s_poly, x_poly=x_poly, num_polys=len(s_poly), tilt_tan=tilt_tan, side=side,
                                                     y_min=-0.1, y_max=0.25, current_s=0.6)
    _loop_doublejaw_2partdim(name='expected_s_open_polygon_after_s_with_vlimit', func=func, num_polys=8)

def create_crystal():
    def func(part_x, part_tan_x, part_y, part_tan_y, R, tilt_sin, tilt_cos):
        return kernels.test_crystal(part_x=part_x, part_tan_x=part_tan_x, R=R, width=0.15, length=0.27,
                                    jaw_U=0.11+1.e-12, tilt_sin=tilt_sin, tilt_cos=tilt_cos)
    _loop_cry_1jaw_1partdim(name='expected_s_crystal', func=func)

def create_crystal_after_s():
    def func(part_x, part_tan_x, part_y, part_tan_y, R, tilt_sin, tilt_cos):
        return kernels.test_crystal_after_s(part_x=part_x, part_tan_x=part_tan_x, R=R, width=0.15, length=0.27,
                                            jaw_U=0.11+1.e-12, tilt_sin=tilt_sin, tilt_cos=tilt_cos, current_s=0.6)
    _loop_cry_1jaw_1partdim(name='expected_s_crystal_after_s', func=func)

def create_crystal_with_vlimit():
    def func(part_x, part_tan_x, part_y, part_tan_y, R, tilt_sin, tilt_cos):
        return kernels.test_crystal_with_vlimit(part_x=part_x, part_tan_x=part_tan_x, part_y=part_y, part_tan_y=part_tan_y,
                                                R=R, width=0.15, length=0.27, jaw_U=0.11+1.e-12, tilt_sin=tilt_sin, tilt_cos=tilt_cos,
                                                y_min=-0.1, y_max=0.25)
    _loop_cry_1jaw_2partdim(name='expected_s_crystal_with_vlimit', func=func)

def create_crystal_after_s_with_vlimit():
    def func(part_x, part_tan_x, part_y, part_tan_y, R, tilt_sin, tilt_cos):
        return kernels.test_crystal_after_s_with_vlimit(part_x=part_x, part_tan_x=part_tan_x, part_y=part_y, part_tan_y=part_tan_y,
                                                R=R, width=0.15, length=0.27, jaw_U=0.11+1.e-12, tilt_sin=tilt_sin, tilt_cos=tilt_cos,
                                                y_min=-0.1, y_max=0.25, current_s=0.6)
    _loop_cry_1jaw_2partdim(name='expected_s_crystal_after_s_with_vlimit', func=func)


kernels = _init_kernels()

# create_jaw()
# create_jaw_after_s()
# create_jaw_with_vlimit()
# create_jaw_after_s_with_vlimit()

# create_polygon()
# create_polygon_after_s()
# create_polygon_with_vlimit()
# create_polygon_after_s_with_vlimit()

# create_open_polygon()
# create_open_polygon_after_s()
# create_open_polygon_with_vlimit()
# create_open_polygon_after_s_with_vlimit()

create_crystal()
create_crystal_after_s()
create_crystal_with_vlimit()
create_crystal_after_s_with_vlimit()
