# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

import sys
import json
import pytest
import numpy as np
from pathlib import Path
from functools import cache


import xcoll as xc
import xtrack as xt
from xobjects.test_helpers import skip_if_forbid_compile, for_all_test_contexts
from xobjects.context import get_context_from_string

sys.path.insert(1, (Path(__file__).parent / 'data').as_posix())
from xcoll_geometry_test import XcollGeometryTest


# Cache compilation to avoid duplication

@cache
def _get_geometry(context_name):
    context = get_context_from_string(context_name)
    geom = XcollGeometryTest(_context=context)
    geom.compile_kernels(only_if_needed=True)
    return geom

@pytest.fixture
def geometry(test_context):
    skip_if_forbid_compile()
    return _get_geometry(test_context)


# Arrays and particles are built on geometry.context below. This is important
# for GPU contexts, where all kernel arguments must belong to the exact same
# context instance as the cached test element.

@for_all_test_contexts
@pytest.mark.xcother
def test_jaw(geometry, test_context):
    context = geometry.context

    def run(particles, s_poly, x_poly, tilt_tan, side):
        geometry.test_jaw(
            particles, s_U=s_poly[1], x_U=x_poly[1],
            s_D=s_poly[2], x_D=x_poly[2], tilt_tan=tilt_tan, side=side)

    _check_2jaw_1partdim(
        context, name='expected_s_jaw', run=run, num_polys=4)


@for_all_test_contexts
@pytest.mark.xcother
def test_jaw_after_s(geometry, test_context):
    context = geometry.context

    def run(particles, s_poly, x_poly, tilt_tan, side):
        geometry.test_jaw_after_s(
            particles, s_U=s_poly[1], x_U=x_poly[1],
            s_D=s_poly[2], x_D=x_poly[2], tilt_tan=tilt_tan, side=side,
            current_s=0.6)

    _check_2jaw_1partdim(
        context, name='expected_s_jaw_after_s', run=run, num_polys=4)


@for_all_test_contexts
@pytest.mark.xcother
def test_jaw_with_vlimit(geometry, test_context):
    context = geometry.context

    def run(particles, s_poly, x_poly, tilt_tan, side):
        geometry.test_jaw_with_vlimit(
            particles, s_U=s_poly[1], x_U=x_poly[1],
            s_D=s_poly[2], x_D=x_poly[2], tilt_tan=tilt_tan, side=side,
            y_min=-0.1, y_max=0.25)

    _check_doublejaw_2partdim(
        context, name='expected_s_jaw_with_vlimit', run=run,
        num_polys=4)


@for_all_test_contexts
@pytest.mark.xcother
def test_jaw_after_s_with_vlimit(geometry, test_context):
    context = geometry.context

    def run(particles, s_poly, x_poly, tilt_tan, side):
        geometry.test_jaw_after_s_with_vlimit(
            particles, s_U=s_poly[1], x_U=x_poly[1],
            s_D=s_poly[2], x_D=x_poly[2], tilt_tan=tilt_tan, side=side,
            y_min=-0.1, y_max=0.25, current_s=0.6)

    _check_doublejaw_2partdim(
        context, name='expected_s_jaw_after_s_with_vlimit', run=run,
        num_polys=4)


@for_all_test_contexts
@pytest.mark.xcother
def test_polygon(geometry, test_context):
    context = geometry.context

    def run(particles, s_poly, x_poly, tilt_tan, side):
        num_polys = len(s_poly)
        s_poly, x_poly = _polygon_to_context(context, s_poly, x_poly)
        geometry.test_polygon(
            particles, s_poly=s_poly, x_poly=x_poly, num_polys=num_polys)

    _check_1jaw_1partdim(
        context, name='expected_s_polygon', run=run, num_polys=8)


@for_all_test_contexts
@pytest.mark.xcother
def test_polygon_after_s(geometry, test_context):
    context = geometry.context

    def run(particles, s_poly, x_poly, tilt_tan, side):
        num_polys = len(s_poly)
        s_poly, x_poly = _polygon_to_context(context, s_poly, x_poly)
        geometry.test_polygon_after_s(
            particles, s_poly=s_poly, x_poly=x_poly, num_polys=num_polys,
            current_s=0.6)

    _check_1jaw_1partdim(
        context, name='expected_s_polygon_after_s', run=run, num_polys=8)


@for_all_test_contexts
@pytest.mark.xcother
def test_polygon_with_vlimit(geometry, test_context):
    context = geometry.context

    def run(particles, s_poly, x_poly, tilt_tan, side):
        num_polys = len(s_poly)
        s_poly, x_poly = _polygon_to_context(context, s_poly, x_poly)
        geometry.test_polygon_with_vlimit(
            particles, s_poly=s_poly, x_poly=x_poly, num_polys=num_polys,
            y_min=-0.1, y_max=0.25)

    _check_1jaw_2partdim(
        context, name='expected_s_polygon_with_vlimit', run=run,
        num_polys=8)


@for_all_test_contexts
@pytest.mark.xcother
def test_polygon_after_s_with_vlimit(geometry, test_context):
    context = geometry.context

    def run(particles, s_poly, x_poly, tilt_tan, side):
        num_polys = len(s_poly)
        s_poly, x_poly = _polygon_to_context(context, s_poly, x_poly)
        geometry.test_polygon_after_s_with_vlimit(
            particles, s_poly=s_poly, x_poly=x_poly, num_polys=num_polys,
            y_min=-0.1, y_max=0.25, current_s=0.6)

    _check_1jaw_2partdim(
        context, name='expected_s_polygon_after_s_with_vlimit', run=run,
        num_polys=8)


@for_all_test_contexts
@pytest.mark.xcother
def test_open_polygon(geometry, test_context):
    context = geometry.context

    def run(particles, s_poly, x_poly, tilt_tan, side):
        num_polys = len(s_poly)
        s_poly, x_poly = _polygon_to_context(context, s_poly, x_poly)
        geometry.test_open_polygon(
            particles, s_poly=s_poly, x_poly=x_poly, num_polys=num_polys,
            tilt_tan=tilt_tan, side=side)

    _check_2jaw_1partdim(
        context, name='expected_s_open_polygon', run=run, num_polys=8)


@for_all_test_contexts
@pytest.mark.xcother
def test_open_polygon_after_s(geometry, test_context):
    context = geometry.context

    def run(particles, s_poly, x_poly, tilt_tan, side):
        num_polys = len(s_poly)
        s_poly, x_poly = _polygon_to_context(context, s_poly, x_poly)
        geometry.test_open_polygon_after_s(
            particles, s_poly=s_poly, x_poly=x_poly, num_polys=num_polys,
            tilt_tan=tilt_tan, side=side, current_s=0.6)

    _check_2jaw_1partdim(
        context, name='expected_s_open_polygon_after_s', run=run,
        num_polys=8)


@for_all_test_contexts
@pytest.mark.xcother
def test_open_polygon_with_vlimit(geometry, test_context):
    context = geometry.context

    def run(particles, s_poly, x_poly, tilt_tan, side):
        num_polys = len(s_poly)
        s_poly, x_poly = _polygon_to_context(context, s_poly, x_poly)
        geometry.test_open_polygon_with_vlimit(
            particles, s_poly=s_poly, x_poly=x_poly, num_polys=num_polys,
            tilt_tan=tilt_tan, side=side, y_min=-0.1, y_max=0.25)

    _check_doublejaw_2partdim(
        context, name='expected_s_open_polygon_with_vlimit', run=run,
        num_polys=8)


@for_all_test_contexts
@pytest.mark.xcother
def test_open_polygon_after_s_with_vlimit(geometry, test_context):
    context = geometry.context

    def run(particles, s_poly, x_poly, tilt_tan, side):
        num_polys = len(s_poly)
        s_poly, x_poly = _polygon_to_context(context, s_poly, x_poly)
        geometry.test_open_polygon_after_s_with_vlimit(
            particles, s_poly=s_poly, x_poly=x_poly, num_polys=num_polys,
            tilt_tan=tilt_tan, side=side, y_min=-0.1, y_max=0.25,
            current_s=0.6)

    _check_doublejaw_2partdim(
        context, name='expected_s_open_polygon_after_s_with_vlimit',
        run=run, num_polys=8)


@for_all_test_contexts
@pytest.mark.xcother
def test_crystal(geometry, test_context):
    context = geometry.context

    def run(particles, R, tilt_sin, tilt_cos):
        geometry.test_crystal(
            particles, R=R, width=0.15, length=0.27, jaw_U=0.11 + 1.e-12,
            tilt_sin=tilt_sin, tilt_cos=tilt_cos)

    _check_cry_1jaw_1partdim(
        context, name='expected_s_crystal', run=run)


@for_all_test_contexts
@pytest.mark.xcother
def test_crystal_after_s(geometry, test_context):
    context = geometry.context

    def run(particles, R, tilt_sin, tilt_cos):
        geometry.test_crystal_after_s(
            particles, R=R, width=0.15, length=0.27, jaw_U=0.11 + 1.e-12,
            tilt_sin=tilt_sin, tilt_cos=tilt_cos, current_s=0.6)

    _check_cry_1jaw_1partdim(
        context, name='expected_s_crystal_after_s', run=run)


@for_all_test_contexts
@pytest.mark.xcother
def test_crystal_with_vlimit(geometry, test_context):
    context = geometry.context

    def run(particles, R, tilt_sin, tilt_cos):
        geometry.test_crystal_with_vlimit(
            particles, R=R, width=0.15, length=0.27, jaw_U=0.11 + 1.e-12,
            tilt_sin=tilt_sin, tilt_cos=tilt_cos,
            y_min=-0.1, y_max=0.25)

    _check_cry_1jaw_2partdim(
        context, name='expected_s_crystal_with_vlimit', run=run)


@for_all_test_contexts
@pytest.mark.xcother
def test_crystal_after_s_with_vlimit(geometry, test_context):
    context = geometry.context

    def run(particles, R, tilt_sin, tilt_cos):
        geometry.test_crystal_after_s_with_vlimit(
            particles, R=R, width=0.15, length=0.27, jaw_U=0.11 + 1.e-12,
            tilt_sin=tilt_sin, tilt_cos=tilt_cos,
            y_min=-0.1, y_max=0.25, current_s=0.6)

    _check_cry_1jaw_2partdim(
        context, name='expected_s_crystal_after_s_with_vlimit', run=run)


@cache
def _load_geometry_data():
    with open(xc._pkg_root.parent / "tests" / "data" / "geometry.json", "r") as fp:
        return json.load(fp)

def _generate_polygon_points(num_poly, tilt_L=0, tilt_R=0):
    rans = _load_geometry_data()['rans']
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
    import matplotlib.pyplot as plt
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


def _particles_from_1d_points(context, points):
    part_x, part_tan_x, expected_s = (
        np.asarray(values, dtype=np.float64) for values in zip(*points))
    particles = xt.Particles(
        _context=context, p0c=1.e9, x=part_x, px=part_tan_x)
    return particles, expected_s

def _polygon_to_context(context, s_poly, x_poly):
    return (
        context.nparray_to_context_array(s_poly),
        context.nparray_to_context_array(x_poly),
    )

def _particles_from_2d_points(context, points):
    part_x, part_tan_x, part_y, part_tan_y, expected_s = (
        np.asarray(values, dtype=np.float64) for values in zip(*points))
    particles = xt.Particles(
        _context=context, p0c=1.e9,
        x=part_x, px=part_tan_x, y=part_y, py=part_tan_y)
    return particles, expected_s

def _particle_s(context, particles):
    return context.nparray_from_context_array(particles.s)

def _assert_particle_s(context, particles, expected_s):
    np.testing.assert_allclose(
        _particle_s(context, particles), expected_s, rtol=1.e-5, atol=1.e-8)


def _check_1jaw_1partdim(context, name, run, num_polys):
    expected_s = _load_geometry_data()[name]
    for tilt, points_by_angle in expected_s.items():
        points = []
        for part_ang, points_by_x in points_by_angle.items():
            part_tan_x = np.tan(np.deg2rad(int(part_ang)))
            points.extend(
                (int(part_x_cm) / 100., part_tan_x, s)
                for part_x_cm, s in points_by_x.items())
        if not points:
            continue
        particles, reference_s = _particles_from_1d_points(context, points)
        tilt_tan = np.tan(np.deg2rad(int(tilt)))
        s_poly, x_poly, _, _ = _generate_polygon_points(
            num_polys, tilt_L=int(tilt))
        run(particles, s_poly, x_poly, tilt_tan, 1)
        _assert_particle_s(context, particles, reference_s)

def _check_cry_1jaw_1partdim(context, name, run):
    expected_s = _load_geometry_data()[name]
    for tilt, points_by_radius in expected_s.items():
        tilt_sin = np.sin(np.deg2rad(int(tilt)))
        tilt_cos = np.cos(np.deg2rad(int(tilt)))
        for radius, points_by_angle in points_by_radius.items():
            points = []
            for part_ang, points_by_x in points_by_angle.items():
                part_tan_x = np.tan(np.deg2rad(int(part_ang)))
                points.extend(
                    (int(part_x_cm) / 100., part_tan_x, s)
                    for part_x_cm, s in points_by_x.items())
            if not points:
                continue
            particles, reference_s = _particles_from_1d_points(context, points)
            run(particles, float(radius), tilt_sin, tilt_cos)
            _assert_particle_s(context, particles, reference_s)

def _check_1jaw_2partdim(context, name, run, num_polys):
    expected_s = _load_geometry_data()[name]
    for tilt, points_by_angle_x in expected_s.items():
        points = []
        for part_ang_x, points_by_angle_y in points_by_angle_x.items():
            part_tan_x = np.tan(np.deg2rad(int(part_ang_x)))
            for part_ang_y, points_by_x in points_by_angle_y.items():
                part_tan_y = np.tan(np.deg2rad(int(part_ang_y)))
                for part_x_cm, points_by_y in points_by_x.items():
                    points.extend(
                        (int(part_x_cm) / 100., part_tan_x,
                         int(part_y_cm) / 100., part_tan_y, s)
                        for part_y_cm, s in points_by_y.items())
        if not points:
            continue
        particles, reference_s = _particles_from_2d_points(context, points)
        tilt_tan = np.tan(np.deg2rad(int(tilt)))
        s_poly, x_poly, _, _ = _generate_polygon_points(
            num_polys, tilt_L=int(tilt))
        run(particles, s_poly, x_poly, tilt_tan, 1)
        _assert_particle_s(context, particles, reference_s)

def _check_cry_1jaw_2partdim(context, name, run):
    expected_s = _load_geometry_data()[name]
    for tilt, points_by_radius in expected_s.items():
        tilt_sin = np.sin(np.deg2rad(int(tilt)))
        tilt_cos = np.cos(np.deg2rad(int(tilt)))
        for radius, points_by_angle_x in points_by_radius.items():
            points = []
            for part_ang_x, points_by_angle_y in points_by_angle_x.items():
                part_tan_x = np.tan(np.deg2rad(int(part_ang_x)))
                for part_ang_y, points_by_x in points_by_angle_y.items():
                    part_tan_y = np.tan(np.deg2rad(int(part_ang_y)))
                    for part_x_cm, points_by_y in points_by_x.items():
                        points.extend(
                            (int(part_x_cm) / 100., part_tan_x,
                             int(part_y_cm) / 100., part_tan_y, s)
                            for part_y_cm, s in points_by_y.items())
            if not points:
                continue
            particles, reference_s = _particles_from_2d_points(context, points)
            run(particles, float(radius), tilt_sin, tilt_cos)
            _assert_particle_s(context, particles, reference_s)

def _check_2jaw_1partdim(context, name, run, num_polys):
    expected_s = _load_geometry_data()[name]
    for tilt_L, points_by_tilt_R in expected_s.items():
        tilt_tan_L = np.tan(np.deg2rad(int(tilt_L)))
        for tilt_R, points_by_angle in points_by_tilt_R.items():
            points = []
            for part_ang, points_by_x in points_by_angle.items():
                part_tan_x = np.tan(np.deg2rad(int(part_ang)))
                points.extend(
                    (int(part_x_cm) / 100., part_tan_x, s)
                    for part_x_cm, s in points_by_x.items())
            if not points:
                continue
            particles_L, reference_s = _particles_from_1d_points(
                context, points)
            particles_R = particles_L.copy()
            tilt_tan_R = np.tan(np.deg2rad(int(tilt_R)))
            polygons = _generate_polygon_points(
                num_polys, tilt_L=int(tilt_L), tilt_R=int(tilt_R))
            s_poly_L, x_poly_L, s_poly_R, x_poly_R = polygons
            run(particles_L, s_poly_L, x_poly_L, tilt_tan_L, 1)
            run(particles_R, s_poly_R, x_poly_R, tilt_tan_R, -1)
            actual_s = np.minimum(
                _particle_s(context, particles_L),
                _particle_s(context, particles_R))
            np.testing.assert_allclose(
                actual_s, reference_s, rtol=1.e-5, atol=1.e-8)

def _check_doublejaw_2partdim(context, name, run, num_polys):
    expected_s = _load_geometry_data()[name]
    for tilt_LR, points_by_angle_x in expected_s.items():
        points = []
        for part_ang_x, points_by_angle_y in points_by_angle_x.items():
            part_tan_x = np.tan(np.deg2rad(int(part_ang_x)))
            for part_ang_y, points_by_x in points_by_angle_y.items():
                part_tan_y = np.tan(np.deg2rad(int(part_ang_y)))
                for part_x_cm, points_by_y in points_by_x.items():
                    points.extend(
                        (int(part_x_cm) / 100., part_tan_x,
                         int(part_y_cm) / 100., part_tan_y, s)
                        for part_y_cm, s in points_by_y.items())
        if not points:
            continue
        particles_L, reference_s = _particles_from_2d_points(context, points)
        particles_R = particles_L.copy()
        tilt_L, tilt_R = tilt_LR.strip('][').split(', ')
        tilt_tan_L = np.tan(np.deg2rad(int(tilt_L)))
        tilt_tan_R = np.tan(np.deg2rad(int(tilt_R)))
        polygons = _generate_polygon_points(
            num_polys, tilt_L=int(tilt_L), tilt_R=int(tilt_R))
        s_poly_L, x_poly_L, s_poly_R, x_poly_R = polygons
        run(particles_L, s_poly_L, x_poly_L, tilt_tan_L, 1)
        run(particles_R, s_poly_R, x_poly_R, tilt_tan_R, -1)
        actual_s = np.minimum(
            _particle_s(context, particles_L),
            _particle_s(context, particles_R))
        np.testing.assert_allclose(
            actual_s, reference_s, rtol=1.e-5, atol=1.e-8)
