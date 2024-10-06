# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

import xobjects as xo
import xtrack as xt
from xcoll.scattering_routines.geometry import XcollGeometry


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

def _create_c_vars(vars, vlimit, after_s):
    vars_part = 'double part_s, double part_x, double part_tan_x, '
    if vlimit == '_vlimit':
        vars_part += 'double part_y, double part_tan_y, '
    vars_end = ''
    if vlimit == '_vlimit':
        vars_end += ', double y_min, double y_max'
    if after_s == '_after_s':
        vars_end += ', double after_s'
    return f'{vars_part}{vars}{vars_end}'

def _create_c_crossing_func(num_segments, vlimit, after_s):
    if after_s == '_first' and vlimit == '':
        return f'double s = crossing_drift_first(segments, {num_segments}, part_s, part_x, part_tan_x);'
    elif vlimit == '':
        return f'double s = crossing_drift_after_s(segments, {num_segments}, part_s, part_x, part_tan_x, after_s);'
    elif after_s == '_first':
        return f'double s = crossing_drift_vlimit_first(segments, {num_segments}, part_s, part_x, part_tan_x, part_y, part_tan_y, y_min, y_max);'
    else:
        return f'double s = crossing_drift_vlimit_after_s(segments, {num_segments}, part_s, part_x, part_tan_x, part_y, part_tan_y, y_min, y_max, after_s);'


# These tests would be very inefficient in tracking, as the segments should be
# created before tracking each particle. For ease of testing (and because we
# cannot return Segments to Python), we combine both into one function.

src_geomtest = []
# Jaw
for vlimit in ['', '_vlimit']:
    for after_s in ['_first', '_after_s']:
        vars = _create_c_vars("double s_U, double x_U, double s_D, double x_D, double tilt_tan, int8_t side", vlimit, after_s)
        src_geomtest.append(f"/*gpufun*/")
        src_geomtest.append(f"double test_jaw{vlimit}{after_s}({vars}){{")
        src_geomtest.append(f"    Segment* segments = create_jaw(s_U, x_U, s_D, x_D, tilt_tan, side);")
        src_geomtest.append(f"    {_create_c_crossing_func(3, vlimit, after_s)}")
        src_geomtest.append(f"    destroy_jaw(segments);  // Important to free memory!!")
        src_geomtest.append(f"    return s;")
        src_geomtest.append(f"}}")
        src_geomtest.append(f"")
# Polygon
for vlimit in ['', '_vlimit']:
    for after_s in ['_first', '_after_s']:
        vars = _create_c_vars("double* s_poly, double* x_poly, int8_t num_polys", vlimit, after_s)
        src_geomtest.append(f"/*gpufun*/")
        src_geomtest.append(f"double test_polygon{vlimit}{after_s}({vars}){{")
        src_geomtest.append(f"    Segment* segments = create_polygon(s_poly, x_poly, num_polys);")
        src_geomtest.append(f"    {_create_c_crossing_func('num_polys', vlimit, after_s)}")
        src_geomtest.append(f"    destroy_polygon(segments, num_polys);  // Important to free memory!!")
        src_geomtest.append(f"    return s;")
        src_geomtest.append(f"}}")
        src_geomtest.append(f"")
# Open polygon
for vlimit in ['', '_vlimit']:
    for after_s in ['_first', '_after_s']:
        vars = _create_c_vars("double* s_poly, double* x_poly, int8_t num_polys, double tilt_tan, int8_t side", vlimit, after_s)
        src_geomtest.append(f"/*gpufun*/")
        src_geomtest.append(f"double test_open_polygon{vlimit}{after_s}({vars}){{")
        src_geomtest.append(f"    Segment* segments = create_open_polygon(s_poly, x_poly, num_polys, tilt_tan, side);")
        src_geomtest.append(f"    {_create_c_crossing_func('num_polys+1', vlimit, after_s)}")
        src_geomtest.append(f"    destroy_open_polygon(segments, num_polys);  // Important to free memory!!")
        src_geomtest.append(f"    return s;")
        src_geomtest.append(f"}}")
        src_geomtest.append(f"")
# Crystal
for vlimit in ['', '_vlimit']:
    for after_s in ['_first', '_after_s']:
        vars = _create_c_vars("double R, double width, double length, double jaw_U, double tilt_sin, double tilt_cos", vlimit, after_s)
        src_geomtest.append(f"/*gpufun*/")
        src_geomtest.append(f"double test_crystal{vlimit}{after_s}({vars}){{")
        src_geomtest.append(f"    Segment* segments = create_crystal(R, width, length, jaw_U, tilt_sin, tilt_cos);")
        src_geomtest.append(f"    {_create_c_crossing_func(4, vlimit, after_s)}")
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
        new_name = f'{name}_first'
        new_kernel_dct[new_name] = xo.Kernel(c_name=new_name, args=ker.args, ret=ker.ret)
        # Add after_s kernels
        new_name = f'{name}_after_s'
        new_args = [*ker.args, xo.Arg(xo.Float64, name='after_s')]
        new_kernel_dct[new_name] = xo.Kernel(c_name=new_name, args=new_args, ret=ker.ret)
        # Add with_vlimit kernels
        new_name = f'{name}_vlimit_first'
        new_args = [*ker.args, xo.Arg(xo.Float64, name='y_min'), xo.Arg(xo.Float64, name='y_max')]
        new_args.insert(3, xo.Arg(xo.Float64, name='part_tan_y'))
        new_args.insert(3, xo.Arg(xo.Float64, name='part_y'))
        new_kernel_dct[new_name] = xo.Kernel(c_name=new_name, args=new_args, ret=ker.ret)
        # Add both
        new_name = f'{name}_vlimit_after_s'
        new_args = [*new_args, xo.Arg(xo.Float64, name='after_s')]
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
                    xo.Arg(xo.Float64, name='part_s'),
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
                    xo.Arg(xo.Float64, pointer=False, name='part_s'),
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
                    xo.Arg(xo.Float64, pointer=False, name='part_s'),
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
                    xo.Arg(xo.Float64, name='part_s'),
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
