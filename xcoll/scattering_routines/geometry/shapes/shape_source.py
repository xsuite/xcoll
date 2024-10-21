# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

import numpy as np

import xobjects as xo

from ..trajectories import trajectories, trajectories_vlimit_sources, get_max_crossings


# =================
# == s positions ==
# =================

# Different types of s selections (first crossing, last one before a given s, first one after a given s, last crossing)
all_s_positions = {
    'first': {
        'args': [],
        'code': """if (n_hit>0){
                return s[0];
            }
            return XC_S_MAX;"""},
    'before_s': {
        'args': [xo.Arg(xo.Float64, name="before_s")],
        'code': """for (int8_t i=n_hit-1; i>=0; i--){
                if (s[i] <= before_s){
                    return s[i];
                }
            }
            return XC_S_MAX;"""},
    'after_s': {
        'args': [xo.Arg(xo.Float64, name="after_s")],
        'code': """for (int8_t i=0; i<n_hit; i++){
                if (s[i] >= after_s){
                    return s[i];
                }
            }
            return XC_S_MAX;"""},
    'last': {
        'args': [],
        'code': """if (n_hit>0){
                return s[n_hit-1];
            }
            return XC_S_MAX;"""}
}


# ==========================
# == C source for Shape2D ==
# ==========================

shape_source = []
for trajectory, vals in trajectories.items():
    c_types = ", ".join([f"{arg.get_c_type()} {arg.name}" for arg in vals['args']])
    c_names = ", ".join([f"{arg.name}" for arg in vals['args']])

    # Function to get all crossings with all segments (loop and sort)
    shape_source.append(f"""
/*gpufun*/
void Shape2D_crossing_{trajectory}(Shape2D shape, int8_t* n_hit, double* s, {c_types}){{
    int64_t n_segments = Shape2D_len_segments(shape);
    for (int8_t i=0; i<n_segments;i++) {{
        LocalSegment seg = Shape2D_getp1_segments(shape, i);
        LocalSegment_crossing_{trajectory}(seg, n_hit, s, {c_names});
    }}
    sort_array_of_double(s, (int64_t) *n_hit);
}}""")

    # Functions to get one crossing: (first crossing, last one before a given s, first one after a given s, last crossing)
    # This is just a template; the actual code will be generated at object instantiations and input inside the switch cases
    for s_pos, s_vals in all_s_positions.items():
        s_pos_types = ''
        if s_vals['args']:
            s_pos_types = ", " + ", ".join([f"{arg.get_c_type()} {arg.name}" for arg in s_vals['args']])
        shape_source.append(f"""
/*gpufun*/
double Shape2D_crossing_{trajectory}_{s_pos}(Shape2D shape, {c_types}{s_pos_types}){{
    int8_t n_hit = 0;
    int64_t seg_id = Shape2D_get__seg_id(shape);
    switch (seg_id){{
        /*START_SEG_ID_CASES_Shape2D_{trajectory}_{s_pos}*/
        /*END_SEG_ID_CASES_Shape2D_{trajectory}_{s_pos}*/
        default:
            printf("Unknown seg_id for Shape2D: %ld\\nPlease recompile.", seg_id); fflush(stdout);  //only_for_context cpu_serial
            (void) n_hit;  // Avoid warning
            return XC_S_MAX;
    }}
}}""")


# ===========================
# == C source for Shape2DV ==
# ===========================

shape_v_source = []
for trajectory, vals in trajectories.items():
    c_types_all    = ", ".join([f"{arg.get_c_type()} {arg.name}" for arg in vals['args'] + vals['args_vlimit_extra']])
    c_types_vlimit = ", ".join([f"{arg.get_c_type()} {arg.name}" for arg in vals['args_vlimit_common'] + vals['args_vlimit_extra']])
    c_names        = ", ".join([f"{arg.name}" for arg in vals['args']])
    c_names_vlimit = ", ".join([f"{arg.name}" for arg in vals['args_vlimit_common'] + vals['args_vlimit_extra']])

    # Assert that all trajectories have a vlimit function and add it to the source code
    header = f"/*gpufun*/\nint8_t vlimit_{trajectory}(double* restrict_s, {c_types_vlimit}, double ymin, double ymax)"
    if header not in trajectories_vlimit_sources[trajectory]:
        raise SystemError(f"Missing or corrupt C vlimit function for {trajectory}.")
    shape_v_source.append(trajectories_vlimit_sources[trajectory])

    # Function to get all crossings with all segments, considering vertical limits
    shape_v_source.append(f"""
/*gpufun*/
void Shape2DV_crossing_{trajectory}(Shape2DV shape, int8_t* n_hit, double* s, {c_types_all}){{
    double restrict_s[2];
    double ymin = Shape2DV_get_vlimit(shape, 0);
    double ymax = Shape2DV_get_vlimit(shape, 1);
    int8_t v_result = vlimit_{trajectory}(restrict_s, {c_names_vlimit}, ymin, ymax);
    if (v_result == 0){{
        return;  // Completely outside - no crossing possible
    }} else {{
        int64_t n_segments = Shape2DV_len_segments(shape);
        for (int8_t i=0; i<n_segments;i++) {{
            LocalSegment seg = Shape2DV_getp1_segments(shape, i);
            LocalSegment_crossing_{trajectory}(seg, n_hit, s, {c_names});
        }}
        sort_array_of_double(s, (int64_t) *n_hit);
        if (v_result == 1){{
            calculate_overlap_array_interval(s, n_hit, restrict_s);
        }}
    }}
}}""")

    # Functions to get one crossing: (first crossing, last one before a given s, first one after a given s, last crossing)
    # This is just a template; the actual code will be generated at object instantiations and input inside the switch cases
    for s_pos, s_vals in all_s_positions.items():
        s_pos_types = ''
        if s_vals['args']:
            s_pos_types = ", " + ", ".join([f"{arg.get_c_type()} {arg.name}" for arg in s_vals['args']])
        shape_v_source.append(f"""
/*gpufun*/
double Shape2DV_crossing_{trajectory}_{s_pos}(Shape2DV shape, {c_types_all}{s_pos_types}){{
    int8_t n_hit = 0;
    int64_t seg_id = Shape2DV_get__seg_id(shape);
    switch (seg_id){{
        /*START_SEG_ID_CASES_Shape2DV_{trajectory}_{s_pos}*/
        /*END_SEG_ID_CASES_Shape2DV_{trajectory}_{s_pos}*/
        default:
            printf("Unknown seg_id for Shape2DV: %ld\\nPlease recompile.", seg_id); fflush(stdout);  //only_for_context cpu_serial
            (void) n_hit;  // Avoid warning
            return XC_S_MAX;
    }}
}}""")


# ==================
# == C source API ==
# ==================

# Function to get the seg_id and max_crossings for each existing object type
def get_seg_ids(shape):
    sources = shape.__class__._extra_c_sources
    seg_ids = {}
    for src in sources:
        if f'SEG_ID_CASES_{shape.__class__.__name__}_drift' in src:
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
    for src in cls._extra_c_sources:
        for s_pos, s_code in all_s_positions.items():
            c_args = trajectories[trajectory]['args'].copy()
            if cls.__name__.endswith('V'):
                c_args += trajectories[trajectory]['args_vlimit_extra']
            c_names = ", ".join([f"{arg.name}" for arg in c_args])
            src = src.replace(f"/*END_SEG_ID_CASES_{cls.__name__}_{trajectory}_{s_pos}*/",
    f"""case {shape._seg_id}: {{  // SIZE: {max_crossings}
            double s[{max_crossings}];
            {cls.__name__}_crossing_{trajectory}(shape, &n_hit, s, {c_names});
            {s_code['code']}
        }}
        /*END_SEG_ID_CASES_{cls.__name__}_{trajectory}_{s_pos}*/""")
        sources_new.append(src)

    cls._extra_c_sources = sources_new
