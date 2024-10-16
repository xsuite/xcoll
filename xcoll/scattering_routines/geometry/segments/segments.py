# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

import numpy as np

import xobjects as xo

from ..trajectories import trajectories, get_max_crossings
from .line import LineSegment
from .halfopen_line import HalfOpenLineSegment
from .circular import CircularSegment
from .bezier import BezierSegment
from .segments_source import segments_source, segments_vlimit_source, get_seg_ids, \
                             create_cases_in_source, assert_localsegment_sources

all_segments = (LineSegment, HalfOpenLineSegment, CircularSegment, BezierSegment)


# Sanity check to assert all segment types have crossing functions for all trajectories
for seg in all_segments:
    assert_localsegment_sources(seg)


# ===========================
# == General segment class ==
# ===========================

class LocalSegment(xo.UnionRef):
    """General segment, acting as a xobject-style parent class for all segment types"""
    _reftypes = all_segments
    _methods = [xo.Method(
                    c_name=f"crossing_{trajectory}",
                    args=[
                        xo.Arg(xo.Int8,    pointer=True,  name="n_hit"),
                        xo.Arg(xo.Float64, pointer=True,  name="s"),
                        *vals["args"]
                    ],
                    ret=None)
                for trajectory, vals in trajectories.items()]


# TODO TODO Need to recompile/assign compilation
class Segments(xo.Struct):
    """Array of segments, representing an object in the geometry"""
    segments = LocalSegment[:]
    _seg_id  = xo.Int64  # This links the object to the correct array size for the crossings s

    _extra_c_sources = segments_source

    def __init__(self, segments=None, **kwargs):
        if not segments:
            raise ValueError("Need to provide `segments`.")
        kwargs['segments'] = segments
        _init_segments_class(self, **kwargs)

    def __repr__(self):
        return f"Segments([{', '.join([seg.__class__.__name__ + '(...)' for seg in self])}])"

    def __getitem__(self, i):
        return self.segments[i]

    def __iter__(self):
        return iter(self.segments)

    def evaluate(self, t):
        s = []
        x = []
        for seg in self.segments:
            this_s, this_x = seg.evaluate(t)
            s.append(this_s)
            x.append(this_x)
        return np.concatenate(s), np.concatenate(x)


class SegmentsVLimit(xo.Struct):
    """Array of segments, representing an object in the geometry, with vertical limits"""
    segments = LocalSegment[:]
    vlimit   = xo.Float64[2]
    _seg_id  = xo.Int64  # This links the object to the correct array size for the crossings s

    _depends_on = [Segments]
    _extra_c_sources = segments_vlimit_source

    def __init__(self, segments=None, vlimit=None, **kwargs):
        if not segments:
            raise ValueError("Need to provide `segments`.")
        kwargs['segments'] = segments
        if not vlimit or isinstance(vlimit, str) \
        or not hasattr(vlimit, '__iter__') or len(vlimit) != 2:
            raise ValueError("Need to provide `vlimit` as [ymin, ymax].")
        kwargs['vlimit'] = vlimit
        _init_segments_class(self, **kwargs)

    def __repr__(self):
        return f"SegmentsVLimit([{', '.join([seg.__class__.__name__ + '(...)' for seg in self])}], vlimit=[{self.vlimit[0], self.vlimit[1]}])"

    def __getitem__(self, i):
        return self.segments[i]

    def __iter__(self):
        return iter(self.segments)

    def evaluate(self, t):
        s = []
        x = []
        for seg in self.segments:
            this_s, this_x = seg.evaluate(t)
            s.append(this_s)
            x.append(this_x)
        return np.concatenate(s), np.concatenate(x)


def _init_segments_class(seg, **kwargs):
    # Each different object type will get its own seg_id, by inspecting the source code
    # First we check if code for this object type already exists
    seg_ids = get_seg_ids(seg)
    max_crossings = get_max_crossings(kwargs['segments'], 'drift')  # test with drift to check if seg_id already has source
    add_code = False
    if max_crossings in seg_ids:
        kwargs['_seg_id'] = seg_ids[max_crossings]
    else:
        add_code = True
        kwargs['_seg_id'] = max(seg_ids.values()) + 1 if len(seg_ids) > 0 else 0
    super(seg.__class__, seg).__init__(**kwargs)
    if add_code:
        for trajectory in trajectories.keys():
            create_cases_in_source(seg, trajectory)
