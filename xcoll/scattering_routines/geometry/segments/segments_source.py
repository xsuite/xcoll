# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

import numpy as np

import xobjects as xo

from ..trajectories import trajectories_c_args, get_max_crossings


segments_source = []
for trajectory, c_args in trajectories_c_args.items():
    # Function to get all crossings with all segments
    segments_source.append(f"""
/*gpufun*/
void Segments_crossing_{trajectory}(Segments segs, int8_t* n_hit, double* s, {c_args[0]}){{
    int64_t n_segments = Segments_len_data(segs);
    for (int8_t i=0; i<n_segments;i++) {{
        Segment seg = Segments_getp1_data(segs, i);
        Segment_crossing_{trajectory}(seg, n_hit, s, {c_args[1]});
    }}
    sort_array_of_double(s, (int64_t) *n_hit);
}}
""")

    # Function to get all crossings with all segments, considering vertical limits - TODO TODO
    segments_source.append(f"""
/*gpufun*/
void Segments_crossing_{trajectory}_vlimit(Segments segs, int8_t* n_hit, double* s, {c_args[2]}){{
    int64_t n_segments = Segments_len_data(segs);
    for (int8_t i=0; i<n_segments;i++) {{
        Segment seg = Segments_getp1_data(segs, i);
        // Segment_crossing_{trajectory}(seg, n_hit, s, {c_args[3]});
    }}
    sort_array_of_double(s, (int64_t) *n_hit);
}}
""")

    for s_pos in ['first', 'after_s']:
        for vlimit in ['', '_vlimit']:
            this_c_args = c_args[0] if vlimit == '' else c_args[2]
            this_c_args = f"{this_c_args}, double after_s" if s_pos == 'after_s' else this_c_args
            segments_source.append(f"""
/*gpufun*/
double Segments_crossing_{trajectory}{vlimit}_{s_pos}(Segments segs, {this_c_args}){{
    int8_t n_hit = 0;
    int64_t seg_id = Segments_get__seg_id(segs);
    switch (seg_id){{
        /*START_SEG_ID_CASES_{trajectory}{vlimit}_{s_pos}*/
        /*END_SEG_ID_CASES_{trajectory}{vlimit}_{s_pos}*/
        default:
            printf("Unknown seg_id for Segment: %ld\\nPlease recompile.", seg_id);
            return XC_S_MAX;
    }}
}}
""")

def get_seg_ids(sources):
    seg_ids = {}
    for src in sources:
        if 'SEG_ID_CASES_drift' in src:
            cases = src.split('SEG_ID_CASES')[1].split(' case ')
            if len(cases) > 1:
                for case in cases[1:]:
                    seg_id = int(case.split(':')[0])
                    max_crossings = int(case.split(':')[2].split('\n')[0].strip())
                    if max_crossings not in seg_ids:
                        seg_ids[max_crossings] = seg_id
                    elif seg_id != seg_ids[max_crossings]:
                        raise ValueError(f"C code for Xcoll geometry iscorrupted! Please inspect.")
    return seg_ids


_source_s = {
    'first': """if (n_hit>0){
                return s[0];
            }
            return XC_S_MAX;""",
    'after_s': """for (int8_t i=0; i<n_hit; i++){
                if (s[i] >= after_s){
                    return s[i];
                }
            }
            return XC_S_MAX;"""
}

def create_cases_in_source(segments, trajectory):
    max_crossings = get_max_crossings(segments, trajectory)
    sources_new = []
    for src in segments._extra_c_sources:
        for s_pos in ['first', 'after_s']:
            for vlimit in ['', '_vlimit']:
                if vlimit == '':
                    c_args = trajectories_c_args[trajectory][1]
                else:
                    c_args = trajectories_c_args[trajectory][3]
                src = src.replace(f"/*END_SEG_ID_CASES_{trajectory}{vlimit}_{s_pos}*/",
    f"""case {segments._seg_id}: {{  // SIZE: {max_crossings}
            double s[{max_crossings}];
            Segments_crossing_{trajectory}{vlimit}(segs, &n_hit, s, {c_args});
            {_source_s[s_pos]}
        }}
        /*END_SEG_ID_CASES_{trajectory}{vlimit}_{s_pos}*/""")
        sources_new.append(src)

    return sources_new
