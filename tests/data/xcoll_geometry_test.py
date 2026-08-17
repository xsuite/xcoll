# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

import xobjects as xo
import xtrack as xt
from xcoll import InteractionRecord

# Example jaw code:
# -----------------
#
# void test_jaw_local_particle(XcollGeometryTestData el, LocalParticle* part0, \
#                 double s_U, double x_U, double s_D, double x_D, \
#                 double tilt_tan, int8_t side){
#     START_PER_PARTICLE_BLOCK(part0, part);
#     Segment segments[3];
#     create_jaw(segments, s_U, x_U, s_D, x_D, tilt_tan, side);
#     double part_x = LocalParticle_get_x(part);
#     double part_tan = LocalParticle_get_xp(part);
#     double s = get_s_of_first_crossing(part_x, part_tan, segments, 3);
#     LocalParticle_set_s(part, s);
#     END_PER_PARTICLE_BLOCK;
# }


# Auto-generate code for all objects and methods
# ----------------------------------------------

def _create_c_vars(vars, after_s, vlimit):
    vars_end = ''
    if vlimit != '':
        vars_end += ', double y_min, double y_max'
    if after_s != '':
        vars_end += ', double current_s'
    return f'XcollGeometryTestData el, LocalParticle* part0, {vars}{vars_end}'


def _create_c_particle_vars(vlimit):
    lines = [
        'double part_x = LocalParticle_get_x(part);',
        'double part_tan_x = LocalParticle_get_xp(part);',
    ]
    if vlimit != '':
        lines.extend([
            'double part_y = LocalParticle_get_y(part);',
            'double part_tan_y = LocalParticle_get_yp(part);',
        ])
    return '\n    '.join(lines)

def _create_c_crossing_func(num_segments, after_s, vlimit):
    if after_s == '' and vlimit == '':
        return f'double s = get_s_of_first_crossing(part_x, part_tan_x, segments, {num_segments});'
    elif vlimit == '':
        return f'double s = get_s_of_crossing_after_s(part_x, part_tan_x, segments, {num_segments}, current_s);'
    elif after_s == '':
        return f'double s = get_s_of_first_crossing_with_vlimit(part_x, part_tan_x, part_y, part_tan_y, segments, {num_segments}, y_min, y_max);'
    else:
        return f'double s = get_s_of_crossing_after_s_with_vlimit(part_x, part_tan_x, part_y, part_tan_y, segments, {num_segments}, y_min, y_max, current_s);'


# Production tracking creates segments once per element. These test kernels
# recreate them for each particle because Segments cannot be returned to
# Python, keeping geometry construction and crossing detection in one function.

src_geomtest = [
    '#include "xobjects/headers/common.h"',
    '#include "xtrack/headers/track.h"',
    '#include "xcoll/scattering_routines/geometry/segments.h"',
    '#include "xcoll/scattering_routines/geometry/objects.h"',
    '#include "xcoll/scattering_routines/geometry/rotation.h"',
    '#include "xcoll/scattering_routines/geometry/methods.h"',
    '#include "xcoll/scattering_routines/geometry/get_s.h"',
    '#include "xcoll/interaction_record/interaction_record_src/interaction_record.h"',
    '#include "xtrack/beam_elements/elements_src/track_drift.h"',
    '#include "xtrack/beam_elements/elements_src/track_srotation.h"',
    '#include "xtrack/beam_elements/elements_src/track_xyshift.h"',
]

# Jaw
for vlimit in ['', '_with_vlimit']:
    for after_s in ['', '_after_s']:
        vars = _create_c_vars("double s_U, double x_U, double s_D, double x_D, double tilt_tan, int8_t side", after_s, vlimit)
        src_geomtest.append(f"GPUFUN")
        src_geomtest.append(f"void test_jaw_local_particle{after_s}{vlimit}({vars}){{")
        src_geomtest.append(f"    START_PER_PARTICLE_BLOCK(part0, part);")
        src_geomtest.append(f"    {_create_c_particle_vars(vlimit)}")
        src_geomtest.append(f"    Segment segments[3];")
        src_geomtest.append(f"    create_jaw(segments, s_U, x_U, s_D, x_D, tilt_tan, side);")
        src_geomtest.append(f"    {_create_c_crossing_func(3, after_s, vlimit)}")
        src_geomtest.append(f"    LocalParticle_set_s(part, s);")
        src_geomtest.append(f"    END_PER_PARTICLE_BLOCK;")
        src_geomtest.append(f"}}")
        src_geomtest.append(f"")
# Polygon
for vlimit in ['', '_with_vlimit']:
    for after_s in ['', '_after_s']:
        vars = _create_c_vars("GPUGLMEM double* s_poly, GPUGLMEM double* x_poly, int8_t num_polys", after_s, vlimit)
        src_geomtest.append(f"GPUFUN")
        src_geomtest.append(f"void test_polygon_local_particle{after_s}{vlimit}({vars}){{")
        src_geomtest.append(f"    START_PER_PARTICLE_BLOCK(part0, part);")
        src_geomtest.append(f"    {_create_c_particle_vars(vlimit)}")
        src_geomtest.append(f"    double local_s_poly[XC_MAX_SEGMENTS];")
        src_geomtest.append(f"    double local_x_poly[XC_MAX_SEGMENTS];")
        src_geomtest.append(f"    for (int8_t i = 0; i < num_polys; i++) {{")
        src_geomtest.append(f"        local_s_poly[i] = s_poly[i];")
        src_geomtest.append(f"        local_x_poly[i] = x_poly[i];")
        src_geomtest.append(f"    }}")
        src_geomtest.append(f"    Segment segments[XC_MAX_SEGMENTS];")
        src_geomtest.append(f"    create_polygon(segments, local_s_poly, local_x_poly, num_polys);")
        src_geomtest.append(f"    {_create_c_crossing_func('num_polys', after_s, vlimit)}")
        src_geomtest.append(f"    LocalParticle_set_s(part, s);")
        src_geomtest.append(f"    END_PER_PARTICLE_BLOCK;")
        src_geomtest.append(f"}}")
        src_geomtest.append(f"")
# Open polygon
for vlimit in ['', '_with_vlimit']:
    for after_s in ['', '_after_s']:
        vars = _create_c_vars("GPUGLMEM double* s_poly, GPUGLMEM double* x_poly, int8_t num_polys, double tilt_tan, int8_t side", after_s, vlimit)
        src_geomtest.append(f"GPUFUN")
        src_geomtest.append(f"void test_open_polygon_local_particle{after_s}{vlimit}({vars}){{")
        src_geomtest.append(f"    START_PER_PARTICLE_BLOCK(part0, part);")
        src_geomtest.append(f"    {_create_c_particle_vars(vlimit)}")
        src_geomtest.append(f"    double local_s_poly[XC_MAX_SEGMENTS];")
        src_geomtest.append(f"    double local_x_poly[XC_MAX_SEGMENTS];")
        src_geomtest.append(f"    for (int8_t i = 0; i < num_polys; i++) {{")
        src_geomtest.append(f"        local_s_poly[i] = s_poly[i];")
        src_geomtest.append(f"        local_x_poly[i] = x_poly[i];")
        src_geomtest.append(f"    }}")
        src_geomtest.append(f"    Segment segments[XC_MAX_SEGMENTS];")
        src_geomtest.append(f"    create_open_polygon(segments, local_s_poly, local_x_poly, num_polys, tilt_tan, side);")
        src_geomtest.append(f"    {_create_c_crossing_func('num_polys+1', after_s, vlimit)}")
        src_geomtest.append(f"    LocalParticle_set_s(part, s);")
        src_geomtest.append(f"    END_PER_PARTICLE_BLOCK;")
        src_geomtest.append(f"}}")
        src_geomtest.append(f"")
# Crystal
for vlimit in ['', '_with_vlimit']:
    for after_s in ['', '_after_s']:
        vars = _create_c_vars("double R, double width, double length, double jaw_U, double tilt_sin, double tilt_cos", after_s, vlimit)
        src_geomtest.append(f"GPUFUN")
        src_geomtest.append(f"void test_crystal_local_particle{after_s}{vlimit}({vars}){{")
        src_geomtest.append(f"    START_PER_PARTICLE_BLOCK(part0, part);")
        src_geomtest.append(f"    {_create_c_particle_vars(vlimit)}")
        src_geomtest.append(f"    Segment segments[XC_MAX_SEGMENTS];  // by value, no allocation")
        src_geomtest.append(f"    create_crystal(segments, R, width, length, jaw_U, tilt_sin, tilt_cos);")
        src_geomtest.append(f"    {_create_c_crossing_func(4, after_s, vlimit)}")
        src_geomtest.append(f"    LocalParticle_set_s(part, s);")
        src_geomtest.append(f"    END_PER_PARTICLE_BLOCK;")
        src_geomtest.append(f"}}")
        src_geomtest.append(f"")
src_geomtest = '\n'.join(src_geomtest)


# Multiply object kernels for all methods
# ---------------------------------------

def mult_kernels(kernel_dct):
    new_kernel_dct = {}
    for name, ker in kernel_dct.items():
        new_kernel_dct[name] = ker
        # Add after_s kernels
        new_name = f'{name}_after_s'
        new_args = [*ker.args, xo.Arg(xo.Float64, name='current_s')]
        new_c_name = f'{ker.c_name}_after_s'
        new_kernel_dct[new_name] = xo.Kernel(c_name=new_c_name, args=new_args)
        # Add with_vlimit kernels
        new_name = f'{name}_with_vlimit'
        new_args = [*ker.args, xo.Arg(xo.Float64, name='y_min'), xo.Arg(xo.Float64, name='y_max')]
        new_c_name = f'{ker.c_name}_with_vlimit'
        new_kernel_dct[new_name] = xo.Kernel(c_name=new_c_name, args=new_args)
        # Add both
        new_name = f'{name}_after_s_with_vlimit'
        new_args = [*new_args, xo.Arg(xo.Float64, name='current_s')]
        new_c_name = f'{ker.c_name}_after_s_with_vlimit'
        new_kernel_dct[new_name] = xo.Kernel(c_name=new_c_name, args=new_args)
    return new_kernel_dct

class XcollGeometryTest(xt.BeamElement):
    _xofields = {}

    allow_track = False
    allow_no_prebuilt_kernel = True

    _depends_on = [InteractionRecord]
    _extra_c_sources = [src_geomtest]

    _per_particle_kernels = mult_kernels({
        'test_jaw': xo.Kernel(
                c_name='test_jaw_local_particle',
                args=[
                    xo.Arg(xo.Float64, name='s_U'),
                    xo.Arg(xo.Float64, name='x_U'),
                    xo.Arg(xo.Float64, name='s_D'),
                    xo.Arg(xo.Float64, name='x_D'),
                    xo.Arg(xo.Float64, name='tilt_tan'),
                    xo.Arg(xo.Int8, name='side')
                ]),
        'test_polygon': xo.Kernel(
                c_name='test_polygon_local_particle',
                args=[
                    xo.Arg(xo.Float64, pointer=True, name='s_poly'),
                    xo.Arg(xo.Float64, pointer=True, name='x_poly'),
                    xo.Arg(xo.Int8, name='num_polys')
                ]),
        'test_open_polygon': xo.Kernel(
                c_name='test_open_polygon_local_particle',
                args=[
                    xo.Arg(xo.Float64, pointer=True, name='s_poly'),
                    xo.Arg(xo.Float64, pointer=True, name='x_poly'),
                    xo.Arg(xo.Int8, name='num_polys'),
                    xo.Arg(xo.Float64, name='tilt_tan'),
                    xo.Arg(xo.Int8, name='side')
                ]),
        'test_crystal': xo.Kernel(
                c_name='test_crystal_local_particle',
                args=[
                    xo.Arg(xo.Float64, name='R'),
                    xo.Arg(xo.Float64, name='width'),
                    xo.Arg(xo.Float64, name='length'),
                    xo.Arg(xo.Float64, name='jaw_U'),
                    xo.Arg(xo.Float64, name='tilt_sin'),
                    xo.Arg(xo.Float64, name='tilt_cos')
                ])
    })
