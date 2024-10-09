# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

import numpy as np

import xobjects as xo

from .trajectories import trajectories, trajectories_c_args, get_max_crossings
from ..c_init import GeomCInit
from .line import LineSegment
from .halfopen_line import HalfOpenLineSegment
from .circular import CircularSegment
from .bezier import BezierSegment


all_segments = (LineSegment, HalfOpenLineSegment, CircularSegment, BezierSegment)


# # Sanity check to assert Segment crossing functions are correctly defined for all trajectories
for trajectory, c_args in trajectories_c_args.items():
    for seg in all_segments:
        header = f"/*gpufun*/\nvoid {seg.__name__}_crossing_{trajectory}({seg.__name__} seg, int8_t* n_hit, double* s, {c_args[0]})"
        if not np.any([header in src for src in seg._extra_c_sources]):
            raise ValueError(f"Missing or corrupt C crossing function for {trajectory} in {seg.__name__}.")


class Segment(xo.UnionRef):
    _reftypes = all_segments
    _methods = [xo.Method(
                    c_name=f"crossing_{trajectory}",
                    args=[
                        xo.Arg(xo.Int8,    pointer=True,  name="n_hit"),
                        xo.Arg(xo.Float64, pointer=True,  name="s"),
                        *args["crossing_args"]
                    ],
                    ret=None)
                for trajectory, args in trajectories.items()]


class Segments(xo.Struct):
    data = Segment[:]
    _extra_c_sources = [f"""
/*gpufun*/
void Segments_crossing_{trajectory}(Segments segs, int8_t* n_hit, double* s, {c_args[0]}){{
    int64_t n_segments = Segments_len_data(segs);
    for (int8_t i=0; i<n_segments;i++) {{
        Segment seg = Segments_getp1_data(segs, i);
        Segment_crossing_{trajectory}(seg, n_hit, s, {c_args[1]});
    }}
    sort_array_of_double(s, (int64_t) *n_hit);
}}

/*gpufun*/
void Segments_crossing_{trajectory}_vlimit(Segments segs, int8_t* n_hit, double* s, {c_args[0]}){{
    int64_t n_segments = Segments_len_data(segs);
    for (int8_t i=0; i<n_segments;i++) {{
        Segment seg = Segments_getp1_data(segs, i);
        Segment_crossing_{trajectory}(seg, n_hit, s, {c_args[1]});
    }}
    sort_array_of_double(s, (int64_t) *n_hit);
}}
"""
            for trajectory, args in trajectories.items()]

    def __getitem__(self, i):
        return self.data[i]

    def __iter__(self):
        return iter(self.data)

    def evaluate(self, t):
        s = []
        x = []
        for seg in self.data:
            this_s, this_x = seg.evaluate(t)
            s.append(this_s)
            x.append(this_x)
        return np.concatenate(s), np.concatenate(x)

class GeomObject(xo.Struct):
    segments = Segments

    def __init__(self, segments):
        self.segments = Segments(data=segments)
        self._extra_sources = [
            """
/*gpufun*/
void GeomObject_crossing_drift(GeomObject obj, int8_t* n_hit, double* s, double s0, double x0, double m){
    Segments_crossing_drift(GeomObject_get_segments(obj), n_hit, s, s0, x0, m);
}
"""
        ]

# class Jaw
#     segments = Segments

# get_c_type
