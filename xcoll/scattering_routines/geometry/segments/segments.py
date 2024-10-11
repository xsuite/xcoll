# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

import numpy as np

import xobjects as xo

from ..trajectories import trajectories, trajectories_c_args
from ..c_init import GeomCInit
from .line import LineSegment
from .halfopen_line import HalfOpenLineSegment
from .circular import CircularSegment
from .bezier import BezierSegment
from .segments_source import segments_source


all_segments = (LineSegment, HalfOpenLineSegment, CircularSegment, BezierSegment)


# Sanity check to assert Segment crossing functions are correctly defined for all trajectories
for trajectory, c_args in trajectories_c_args.items():
    for seg in all_segments:
        header = f"/*gpufun*/\nvoid {seg.__name__}_crossing_{trajectory}({seg.__name__} seg, int8_t* n_hit, double* s, {c_args[0]})"
        header_found = False
        for src in seg._extra_c_sources:
            if isinstance(src, str) and header in src:
                header_found = True
                break
            with open(src) as f:
                if header in f.read():
                    header_found = True
                    break
        if not header_found:
            raise ValueError(f"Missing or corrupt C crossing function for {trajectory} in {seg.__name__}.")


class Segment(xo.UnionRef):
    """General segment, acting as a xobject-style parent class for all segment types"""
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
    """Array of segments, representing an object in the geometry"""
    data    = Segment[:]
    _seg_id = xo.Int64  # This links the object to the correct array size for the crossings s

    _extra_c_sources = segments_source

    def __init__(self, segments=None, **kwargs):
        if segments is not None:
            if 'data' in kwargs:
                raise ValueError("Cannot provide 'segments' and 'data' at the same time")
            kwargs['data'] = segments
        elif 'data' in kwargs:
            segments = kwargs['data']
        kwargs['_seg_id'] = 0
        super().__init__(**kwargs)
        if segments is not None:
            print(len(segments)) # Need to get seg_id, make cases in C, recompile/assign compilation

    def __repr__(self):
        return f"Segments([{', '.join([seg.__class__.__name__ + '(...)' for seg in self])}])"

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
