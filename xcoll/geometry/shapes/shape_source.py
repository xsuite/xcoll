# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

import numpy as np

import xobjects as xo

from ..c_init import xo_to_ctypes, xo_to_cnames
from ..segments import get_max_crossings
from ..trajectories import all_trajectories, DriftTrajectory, args_cross_h, args_cross_v, args_vlimit


# =================
# == s positions ==
# =================

# Different types of s selections (first crossing, last one before a given s, first one after a given s, last crossing)
all_s_positions = {
    'first': {
        'args': [],
        'code': f"return Trajectory_get_first({xo_to_cnames(args_cross_h)});"
        },
    'before_s': {
        'args': [xo.Arg(xo.Float64, name="before_s")],
        'code': f"return Trajectory_get_before_s({xo_to_cnames(args_cross_h)}, before_s);"
        },
    'after_s': {
        'args': [xo.Arg(xo.Float64, name="after_s")],
        'code': f"return Trajectory_get_after_s({xo_to_cnames(args_cross_h)}, after_s);"
        },
    'last': {
        'args': [],
        'code': f"return Trajectory_get_last({xo_to_cnames(args_cross_h)});"
        }
}

# # Sanity check to assert the drift trajectory code has all s position functions
# def assert_s_pos_sources(tra):
#     for s_pos, s_args in all_s_positions.items():
#         cnames = f"{args_cross_h[0].get_c_type()[:-1]} {args_cross_h[0].name}, {xo_to_ctypes(args_cross_h[1:])}"
#         if s_args['args']:
#             cnames += f", {xo_to_ctypes(s_args['args'])}"
#         header = f"/*gpufun*/\ndouble Trajectory_get_{s_pos}({cnames})"
#         header_found = False
#         for src in tra._extra_c_sources:
#             if isinstance(src, str):
#                 if header in src:
#                     header_found = True
#                     break
#             else:
#                 with open(src) as f:
#                     if header in f.read():
#                         header_found = True
#                         break
#         if not header_found:
#             raise SystemError(f"Missing or corrupt C Trajectory_get_s function for {s_pos}.")

# assert_s_pos_sources(DriftTrajectory)


# ==========================
# == C source for Shape2D ==
# ==========================

shape_source = []
# for tra in all_trajectories:
#     c_types = f"{xo_to_ctypes(tra.args_hv)}, {xo_to_ctypes(tra.args_h)}"
#     c_names = f"{xo_to_cnames(tra.args_hv)}, {xo_to_cnames(tra.args_h)}"
#     # Function to get all crossings with all segments (loop and sort)
#     shape_source.append(f"""
# /*gpufun*/
# void Shape2D_crossing_{tra.name}(Shape2D shape, {xo_to_ctypes(args_cross_h)}, {c_types}){{
#     int64_t n_segments = Shape2D_len_segments(shape);
#     for (int8_t i=0; i<n_segments;i++) {{
#         LocalSegment seg = Shape2D_getp1_segments(shape, i);
#         LocalSegment_crossing_{tra.name}(seg, {xo_to_cnames(args_cross_h)}, {c_names});
#     }}
#     sort_array_of_double(s, (int64_t) *n_hit);
# }}""")

#     # Functions to get one crossing: (first crossing, last one before a given s, first one after a given s, last crossing)
#     # This is just a template; the actual code will be generated at object instantiations and input inside the switch cases
#     for s_pos, s_vals in all_s_positions.items():
#         s_pos_types = ''
#         if s_vals['args']:
#             s_pos_types = ", " + ", ".join([f"{arg.get_c_type()} {arg.name}" for arg in s_vals['args']])
#         shape_source.append(f"""
# /*gpufun*/
# double Shape2D_crossing_{tra.name}_{s_pos}(Shape2D shape, {c_types}{s_pos_types}){{
#     int8_t n_hit = 0;
#     int64_t seg_id = Shape2D_get__seg_id(shape);
#     switch (seg_id){{
#         /*START_SEG_ID_CASES_Shape2D_{tra.__name__}_{s_pos}*/
#         /*END_SEG_ID_CASES_Shape2D_{tra.__name__}_{s_pos}*/
#         default:
#             printf("Unknown seg_id for Shape2D: %ld\\nPlease recompile.", seg_id); fflush(stdout);  //only_for_context cpu_serial
#             (void) n_hit;  // Avoid warning
#             return XC_S_MAX;
#     }}
# }}""")


# ===========================
# == C source for Shape2DV ==
# ===========================

shape_v_source = []
# for tra in all_trajectories:
#     c_types = f"{xo_to_ctypes(tra.args_hv)}, {xo_to_ctypes(tra.args_h)}, {xo_to_ctypes(tra.args_v)}"
#     c_names_h = f"{xo_to_cnames(tra.args_hv)}, {xo_to_cnames(tra.args_h)}"
#     c_names_v = f"{xo_to_cnames(tra.args_hv)}, {xo_to_cnames(tra.args_v)}"

#     # Function to get all crossings with all segments, considering vertical limits
#     shape_v_source.append(f"""
# /*gpufun*/
# void Shape2DV_crossing_{tra.name}(Shape2DV shape, {xo_to_ctypes(args_cross_h)}, {c_types}){{
#     double {xo_to_cnames(args_cross_v)}[2];
#     {xo_to_ctypes(args_vlimit[0])} = Shape2DV_get_vlimit(shape, 0);
#     {xo_to_ctypes(args_vlimit[1])} = Shape2DV_get_vlimit(shape, 1);
#     int8_t v_result = {tra.__name__}_vlimit({xo_to_cnames(args_cross_v)}, {c_names_v}, {xo_to_cnames(args_vlimit)});
#     if (v_result == 0){{
#         return;  // Completely outside - no crossing possible
#     }} else {{
#         int64_t n_segments = Shape2DV_len_segments(shape);
#         for (int8_t i=0; i<n_segments;i++) {{
#             LocalSegment seg = Shape2DV_getp1_segments(shape, i);
#             LocalSegment_crossing_{tra.name}(seg, {xo_to_cnames(args_cross_h)}, {c_names_h});
#         }}
#         sort_array_of_double(s, (int64_t) *n_hit);
#         if (v_result == 1){{
#             calculate_overlap_array_interval(s, n_hit, restrict_s);
#         }}
#     }}
# }}""")

#     # Functions to get one crossing: (first crossing, last one before a given s, first one after a given s, last crossing)
#     # This is just a template; the actual code will be generated at object instantiations and input inside the switch cases
#     for s_pos, s_vals in all_s_positions.items():
#         s_pos_types = ''
#         if s_vals['args']:
#             s_pos_types = ", " + ", ".join([f"{arg.get_c_type()} {arg.name}" for arg in s_vals['args']])
#         shape_v_source.append(f"""
# /*gpufun*/
# double Shape2DV_crossing_{tra.name}_{s_pos}(Shape2DV shape, {c_types}{s_pos_types}){{
#     int8_t n_hit = 0;
#     int64_t seg_id = Shape2DV_get__seg_id(shape);
#     switch (seg_id){{
#         /*START_SEG_ID_CASES_Shape2DV_{tra.__name__}_{s_pos}*/
#         /*END_SEG_ID_CASES_Shape2DV_{tra.__name__}_{s_pos}*/
#         default:
#             printf("Unknown seg_id for Shape2DV: %ld\\nPlease recompile.", seg_id); fflush(stdout);  //only_for_context cpu_serial
#             (void) n_hit;  // Avoid warning
#             return XC_S_MAX;
#     }}
# }}""")


# ==================
# == C source API ==
# ==================

# Function to get the seg_id and max_crossings for each existing object type
def get_seg_ids(shape):
    sources = shape.__class__._extra_c_sources
    seg_ids = {}
    for src in sources:
        if f'SEG_ID_CASES_{shape.__class__.__name__}_DriftTrajectory' in src:
            cases = src.split('SEG_ID_CASES')[1].split(' case ')
            if len(cases) > 1:
                for case in cases[1:]:
                    seg_id = int(case.split(':')[0])
                    max_crossings = int(case.split(':')[2].split('\n')[0].strip())
                    if max_crossings not in seg_ids:
                        seg_ids[max_crossings] = seg_id
                    elif seg_id != seg_ids[max_crossings]:
                        raise SystemError(f"C code for Xcoll geometry is corrupted! Please inspect.")
    return seg_ids


# Function to inject the switch cases code for the new object type into the Shape2D source
def create_cases_in_source(shape, trajectory):
    cls = shape.__class__
    max_crossings = get_max_crossings(shape, trajectory)
    sources_new = []
    # for src in cls._extra_c_sources:
    #     for s_pos, s_code in all_s_positions.items():
    #         c_names = f"{xo_to_cnames(trajectory.args_hv)}, {xo_to_cnames(trajectory.args_h)}"
    #         if cls.__name__.endswith('V'):
    #             c_names += f", {xo_to_cnames(trajectory.args_v)}"
    #         src = src.replace(f"/*END_SEG_ID_CASES_{cls.__name__}_{trajectory.__name__}_{s_pos}*/",
    # f"""case {shape._seg_id}: {{  // SIZE: {max_crossings}
    #         double {xo_to_cnames(args_cross_h[1])}[{max_crossings}];
    #         {cls.__name__}_crossing_{trajectory.name}(shape, &{xo_to_cnames(args_cross_h)}, {c_names});
    #         {s_code['code']}
    #     }}
    #     /*END_SEG_ID_CASES_{cls.__name__}_{trajectory.__name__}_{s_pos}*/""")
    #     sources_new.append(src)

    cls._extra_c_sources = sources_new
