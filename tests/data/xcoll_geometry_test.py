# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

import xobjects as xo
import xtrack as xt
from xcoll.scattering_routines.geometry import XcollGeometry

# These tests would be very inefficient in tracking, as the segments should be
# created before tracking each particle. For ease of testing (and because we
# cannot return Segments to Python), we combine both into one function.



# Example jaw code:
# -----------------
#
# double test_jaw_first(double part_s, double part_x, double part_tan, double s_U, double x_U, \
#                       double s_D, double x_D, double tilt_tan, int8_t side){
#     Segment* segments = create_jaw(s_U, x_U, s_D, x_D, tilt_tan, side);
#     double s = get_s_of_first_crossing(part_x, part_tan, segments, 3);
#     destroy_jaw(segments);  // Important to free memory!!
#     return s;
# }


# Auto-generate code and kernels for all objects and methods
# ----------------------------------------------------------

def _create_source_per_method_for_object(name, object_vars, object_constructor, object_destructor, num_segments):
    return f'''
/*gpufun*/
double test_{name}_first(double part_s, double part_x, double part_tan_x, {object_vars}){{
    {object_constructor}
    double s = crossing_drift_first(segments, {num_segments}, part_s, part_x, part_tan_x);
    {object_destructor}  // Important to free memory!!
    return s;
}}

/*gpufun*/
double test_{name}_after_s(double part_s, double part_x, double part_tan_x, {object_vars}, double after_s){{
    {object_constructor}
    double s = crossing_drift_after_s(segments, {num_segments}, part_s, part_x, part_tan_x, after_s);
    {object_destructor}  // Important to free memory!!
    return s;
}}

/*gpufun*/
double test_{name}_vlimit_first(double part_s, double part_x, double part_tan_x, double part_y, double part_tan_y, {object_vars}, double y_min, double y_max){{
    {object_constructor}
    double s = crossing_drift_vlimit_first(segments, {num_segments}, part_s, part_x, part_tan_x, part_y, part_tan_y, y_min, y_max);
    {object_destructor}  // Important to free memory!!
    return s;
}}

/*gpufun*/
double test_{name}_vlimit_after_s(double part_s, double part_x, double part_tan_x, double part_y, double part_tan_y, {object_vars}, double y_min, double y_max, double after_s){{
    {object_constructor}
    double s = crossing_drift_vlimit_after_s(segments, {num_segments}, part_s, part_x, part_tan_x, part_y, part_tan_y, y_min, y_max, after_s);
    {object_destructor}  // Important to free memory!!
    return s;
}}
'''

def _create_kernels_per_method_for_object(name, object_args):
    return {
        f'test_{name}_first': xo.Kernel(
                c_name=f'test_{name}_first',
                args=[
                    xo.Arg(xo.Float64, name='part_s'),     # particle args
                    xo.Arg(xo.Float64, name='part_x'),     # particle args
                    xo.Arg(xo.Float64, name='part_tan_x'), # particle args
                    *object_args
                ],
                ret=xo.Arg(xo.Float64, pointer=False, name='s')),
        f'test_{name}_after_s': xo.Kernel(
                c_name=f'test_{name}_after_s',
                args=[
                    xo.Arg(xo.Float64, name='part_s'),     # particle args
                    xo.Arg(xo.Float64, name='part_x'),     # particle args
                    xo.Arg(xo.Float64, name='part_tan_x'), # particle args
                    *object_args,
                    xo.Arg(xo.Float64, name='after_s')     # method args
                ],
                ret=xo.Arg(xo.Float64, pointer=False, name='s')),
        f'test_{name}_vlimit_first': xo.Kernel(
                c_name=f'test_{name}_vlimit_first',
                args=[
                    xo.Arg(xo.Float64, name='part_s'),     # particle args
                    xo.Arg(xo.Float64, name='part_x'),     # particle args
                    xo.Arg(xo.Float64, name='part_tan_x'), # particle args
                    xo.Arg(xo.Float64, name='part_y'),     # particle args
                    xo.Arg(xo.Float64, name='part_tan_y'), # particle args
                    *object_args,
                    xo.Arg(xo.Float64, name='y_min'),      # method args
                    xo.Arg(xo.Float64, name='y_max')       # method args
                ],
                ret=xo.Arg(xo.Float64, pointer=False, name='s')),
        f'test_{name}_vlimit_after_s': xo.Kernel(
                c_name=f'test_{name}_vlimit_after_s',
                args=[
                    xo.Arg(xo.Float64, name='part_s'),     # particle args
                    xo.Arg(xo.Float64, name='part_x'),     # particle args
                    xo.Arg(xo.Float64, name='part_tan_x'), # particle args
                    xo.Arg(xo.Float64, name='part_y'),     # particle args
                    xo.Arg(xo.Float64, name='part_tan_y'), # particle args
                    *object_args,
                    xo.Arg(xo.Float64, name='y_min'),      # method args
                    xo.Arg(xo.Float64, name='y_max'),      # method args
                    xo.Arg(xo.Float64, name='after_s')     # method args
                ],
                ret=xo.Arg(xo.Float64, pointer=False, name='s'))
    }



class XcollGeometryTest(xt.BeamElement):
    _xofields = {}

    allow_track = False

    _depends_on = [XcollGeometry]

    _extra_c_sources = [
        _create_source_per_method_for_object(
            name="jaw", num_segments=3,
            object_vars="double s_U, double x_U, double s_D, double x_D, double tilt_tan, int8_t side",
            object_constructor="Segment* segments = create_jaw(s_U, x_U, s_D, x_D, tilt_tan, side);",
            object_destructor="destroy_jaw(segments);") \
        + _create_source_per_method_for_object(
            name="polygon", num_segments='num_polys',
            object_vars="double* s_poly, double* x_poly, int8_t num_polys",
            object_constructor="Segment* segments = create_polygon(s_poly, x_poly, num_polys);",
            object_destructor="destroy_polygon(segments, num_polys);") \
        + _create_source_per_method_for_object(
            name="open_polygon", num_segments='num_polys',
            object_vars="double* s_poly, double* x_poly, int8_t num_polys, double tilt_tan, int8_t side",
            object_constructor="Segment* segments = create_open_polygon(s_poly, x_poly, num_polys, tilt_tan, side);",
            object_destructor="destroy_open_polygon(segments, num_polys);") \
        + _create_source_per_method_for_object(
            name="crystal", num_segments=4,
            object_vars="double R, double width, double length, double jaw_U, double tilt_sin, double tilt_cos",
            object_constructor="Segment* segments = create_crystal(R, width, length, jaw_U, tilt_sin, tilt_cos);",
            object_destructor="destroy_crystal(segments);")
    ]

    _kernels = {
        **_create_kernels_per_method_for_object(name="jaw", object_args=[
                    xo.Arg(xo.Float64, name='s_U'),
                    xo.Arg(xo.Float64, name='x_U'),
                    xo.Arg(xo.Float64, name='s_D'),
                    xo.Arg(xo.Float64, name='x_D'),
                    xo.Arg(xo.Float64, name='tilt_tan'),
                    xo.Arg(xo.Int8,    name='side'),
            ]),
        **_create_kernels_per_method_for_object(name="polygon", object_args=[
                    xo.Arg(xo.Float64, pointer=True, name='s_poly'),
                    xo.Arg(xo.Float64, pointer=True, name='x_poly'),
                    xo.Arg(xo.Int8,    name='num_polys')
            ]),
        **_create_kernels_per_method_for_object(name="open_polygon", object_args=[
                    xo.Arg(xo.Float64, pointer=True, name='s_poly'),
                    xo.Arg(xo.Float64, pointer=True, name='x_poly'),
                    xo.Arg(xo.Int8,    name='num_polys'),
                    xo.Arg(xo.Float64, name='tilt_tan'),
                    xo.Arg(xo.Int8,    name='side')
            ]),
        **_create_kernels_per_method_for_object(name="crystal", object_args=[
                    xo.Arg(xo.Float64, name='R'),
                    xo.Arg(xo.Float64, name='width'),
                    xo.Arg(xo.Float64, name='length'),
                    xo.Arg(xo.Float64, name='jaw_U'),
                    xo.Arg(xo.Float64, name='tilt_sin'),
                    xo.Arg(xo.Float64, name='tilt_cos')
            ]),
    }
