# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

import xobjects as xo
import xtrack as xt
from xcoll.geometry.old_geometry import XcollGeometry


# Example jaw code:
# -----------------
#
# double test_jaw(double part_x, double part_tan, double s_U, double x_U, double s_D, \
#                 double x_D, double tilt_tan, int8_t side){
#     Segment* segments = create_jaw(s_U, x_U, s_D, x_D, tilt_tan, side);
#     double s = get_s_of_first_crossing(part_x, part_tan, segments, 3);
#     destroy_jaw(segments);  // Important to free memory!!
#     return s;
# }


# Auto-generate code for all objects and methods
# ----------------------------------------------

def _create_c_vars(vars, after_s, vlimit):
    vars_part = 'double part_x, double part_tan_x, '
    if vlimit != '':
        vars_part += 'double part_y, double part_tan_y, '
    vars_end = ''
    if vlimit != '':
        vars_end += ', double y_min, double y_max'
    if after_s != '':
        vars_end += ', double current_s'
    return f'{vars_part}{vars}{vars_end}'

def _create_c_crossing_func(num_segments, after_s, vlimit):
    if after_s == '' and vlimit == '':
        return f'double s = get_s_of_first_crossing(part_x, part_tan_x, segments, {num_segments});'
    elif vlimit == '':
        return f'double s = get_s_of_crossing_after_s(part_x, part_tan_x, segments, {num_segments}, current_s);'
    elif after_s == '':
        return f'double s = get_s_of_first_crossing_with_vlimit(part_x, part_tan_x, part_y, part_tan_y, segments, {num_segments}, y_min, y_max);'
    else:
        return f'double s = get_s_of_crossing_after_s_with_vlimit(part_x, part_tan_x, part_y, part_tan_y, segments, {num_segments}, y_min, y_max, current_s);'


# These tests would be very inefficient in tracking, as the segments should be
# created before tracking each particle. For ease of testing (and because we
# cannot return Segments to Python), we combine both into one function.

src_geomtest = []
# Jaw
for vlimit in ['', '_with_vlimit']:
    for after_s in ['', '_after_s']:
        vars = _create_c_vars("double s_U, double x_U, double s_D, double x_D, double tilt_tan, int8_t side", after_s, vlimit)
        src_geomtest.append(f"/*gpufun*/")
        src_geomtest.append(f"double test_jaw{after_s}{vlimit}({vars}){{")
        src_geomtest.append(f"    Segment* segments = create_jaw(s_U, x_U, s_D, x_D, tilt_tan, side);")
        src_geomtest.append(f"    {_create_c_crossing_func(3, after_s, vlimit)}")
        src_geomtest.append(f"    destroy_jaw(segments);  // Important to free memory!!")
        src_geomtest.append(f"    return s;")
        src_geomtest.append(f"}}")
        src_geomtest.append(f"")
# Polygon
for vlimit in ['', '_with_vlimit']:
    for after_s in ['', '_after_s']:
        vars = _create_c_vars("double* s_poly, double* x_poly, int8_t num_polys", after_s, vlimit)
        src_geomtest.append(f"/*gpufun*/")
        src_geomtest.append(f"double test_polygon{after_s}{vlimit}({vars}){{")
        src_geomtest.append(f"    Segment* segments = create_polygon(s_poly, x_poly, num_polys);")
        src_geomtest.append(f"    {_create_c_crossing_func('num_polys', after_s, vlimit)}")
        src_geomtest.append(f"    destroy_polygon(segments, num_polys);  // Important to free memory!!")
        src_geomtest.append(f"    return s;")
        src_geomtest.append(f"}}")
        src_geomtest.append(f"")
# Open polygon
for vlimit in ['', '_with_vlimit']:
    for after_s in ['', '_after_s']:
        vars = _create_c_vars("double* s_poly, double* x_poly, int8_t num_polys, double tilt_tan, int8_t side", after_s, vlimit)
        src_geomtest.append(f"/*gpufun*/")
        src_geomtest.append(f"double test_open_polygon{after_s}{vlimit}({vars}){{")
        src_geomtest.append(f"    Segment* segments = create_open_polygon(s_poly, x_poly, num_polys, tilt_tan, side);")
        src_geomtest.append(f"    {_create_c_crossing_func('num_polys+1', after_s, vlimit)}")
        src_geomtest.append(f"    destroy_open_polygon(segments, num_polys);  // Important to free memory!!")
        src_geomtest.append(f"    return s;")
        src_geomtest.append(f"}}")
        src_geomtest.append(f"")
# Crystal
for vlimit in ['', '_with_vlimit']:
    for after_s in ['', '_after_s']:
        vars = _create_c_vars("double R, double width, double length, double jaw_U, double tilt_sin, double tilt_cos", after_s, vlimit)
        src_geomtest.append(f"/*gpufun*/")
        src_geomtest.append(f"double test_crystal{after_s}{vlimit}({vars}){{")
        src_geomtest.append(f"    Segment* segments = create_crystal(R, width, length, jaw_U, tilt_sin, tilt_cos);")
        src_geomtest.append(f"    {_create_c_crossing_func(4, after_s, vlimit)}")
        src_geomtest.append(f"    destroy_crystal(segments);  // Important to free memory!!")
        src_geomtest.append(f"    return s;")
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
        new_kernel_dct[new_name] = xo.Kernel(c_name=new_name, args=new_args, ret=ker.ret)
        # Add with_vlimit kernels
        new_name = f'{name}_with_vlimit'
        new_args = [*ker.args, xo.Arg(xo.Float64, name='y_min'), xo.Arg(xo.Float64, name='y_max')]
        new_args.insert(2, xo.Arg(xo.Float64, name='part_tan_y'))
        new_args.insert(2, xo.Arg(xo.Float64, name='part_y'))
        new_kernel_dct[new_name] = xo.Kernel(c_name=new_name, args=new_args, ret=ker.ret)
        # Add both
        new_name = f'{name}_after_s_with_vlimit'
        new_args = [*new_args, xo.Arg(xo.Float64, name='current_s')]
        new_kernel_dct[new_name] = xo.Kernel(c_name=new_name, args=new_args, ret=ker.ret)
    return new_kernel_dct

class XcollGeometryTest(xt.BeamElement):
    _xofields = {}

    allow_track = False

    _depends_on = [XcollGeometry]

    _extra_c_sources = [src_geomtest]

    _kernels = mult_kernels({
        'test_jaw': xo.Kernel(
                c_name='test_jaw',
                args=[
                    xo.Arg(xo.Float64, name='part_x'),
                    xo.Arg(xo.Float64, name='part_tan_x'),
                    xo.Arg(xo.Float64, name='s_U'),
                    xo.Arg(xo.Float64, name='x_U'),
                    xo.Arg(xo.Float64, name='s_D'),
                    xo.Arg(xo.Float64, name='x_D'),
                    xo.Arg(xo.Float64, name='tilt_tan'),
                    xo.Arg(xo.Int8, name='side')
                ],
                ret=xo.Arg(xo.Float64, pointer=False, name='s')),
        'test_polygon': xo.Kernel(
                c_name='test_polygon',
                args=[
                    xo.Arg(xo.Float64, pointer=False, name='part_x'),
                    xo.Arg(xo.Float64, pointer=False, name='part_tan_x'),
                    xo.Arg(xo.Float64, pointer=True, name='s_poly'),
                    xo.Arg(xo.Float64, pointer=True, name='x_poly'),
                    xo.Arg(xo.Int8, name='num_polys')
                ],
                ret=xo.Arg(xo.Float64, pointer=False, name='s')),
        'test_open_polygon': xo.Kernel(
                c_name='test_open_polygon',
                args=[
                    xo.Arg(xo.Float64, pointer=False, name='part_x'),
                    xo.Arg(xo.Float64, pointer=False, name='part_tan_x'),
                    xo.Arg(xo.Float64, pointer=True, name='s_poly'),
                    xo.Arg(xo.Float64, pointer=True, name='x_poly'),
                    xo.Arg(xo.Int8, name='num_polys'),
                    xo.Arg(xo.Float64, name='tilt_tan'),
                    xo.Arg(xo.Int8, name='side')
                ],
                ret=xo.Arg(xo.Float64, pointer=False, name='s')),
        'test_crystal': xo.Kernel(
                c_name='test_crystal',
                args=[
                    xo.Arg(xo.Float64, name='part_x'),
                    xo.Arg(xo.Float64, name='part_tan_x'),
                    xo.Arg(xo.Float64, name='R'),
                    xo.Arg(xo.Float64, name='width'),
                    xo.Arg(xo.Float64, name='length'),
                    xo.Arg(xo.Float64, name='jaw_U'),
                    xo.Arg(xo.Float64, name='tilt_sin'),
                    xo.Arg(xo.Float64, name='tilt_cos')
                ],
                ret=xo.Arg(xo.Float64, pointer=False, name='s'))
    })
