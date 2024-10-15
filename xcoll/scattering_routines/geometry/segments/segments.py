# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

import numpy as np

import xobjects as xo

from ..trajectories import trajectories, trajectories_c_args, get_max_crossings
from ..c_init import GeomCInit
from .line import LineSegment
from .halfopen_line import HalfOpenLineSegment
from .circular import CircularSegment
from .bezier import BezierSegment
from .segments_source import segments_source, get_seg_ids, create_cases_in_source

all_segments = (LineSegment, HalfOpenLineSegment, CircularSegment, BezierSegment)


# Sanity check to assert Segment crossing functions are correctly defined for each segment type and for all trajectories
for trajectory, c_args in trajectories_c_args.items():
    for seg in all_segments:
        header = f"/*gpufun*/\nvoid {seg.__name__}_crossing_{trajectory}({seg.__name__} seg, int8_t* n_hit, double* s, {c_args['c_types']})"
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


# ===========================
# == General segment class ==
# ===========================

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


# TODO TODO Need to recompile/assign compilation
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
        # Each different object type will get its own seg_id, by inspecting the source code
        # First we check if code for this object type already exists
        seg_ids = get_seg_ids(Segments._extra_c_sources)
        max_crossings = get_max_crossings(segments, 'drift')
        add_code = False
        if max_crossings in seg_ids:
            kwargs['_seg_id'] = seg_ids[max_crossings]
        else:
            add_code = True
            kwargs['_seg_id'] = max(seg_ids.values()) + 1 if len(seg_ids) > 0 else 0
        super().__init__(**kwargs)
        if segments is not None and add_code:
            for trajectory in trajectories.keys():
                Segments._extra_c_sources = create_cases_in_source(self, trajectory)

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
