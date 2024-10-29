# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

import numpy as np

import xobjects as xo

from ..c_init import  XC_EPSILON, xo_to_ctypes
from ..trajectories import all_trajectories, DriftTrajectory, args_cross_h, args_cross_v, args_vlimit

from .line import LineSegment
from .halfopen_line import HalfOpenLineSegment
from .circular import CircularSegment
from .bezier import BezierSegment


all_segments = (LineSegment, HalfOpenLineSegment, CircularSegment, BezierSegment)


class LocalSegment(xo.UnionRef):
    """General segment, acting as a xobject-style parent class for all segment types"""
    _reftypes = all_segments
    _methods = [xo.Method(
                    c_name=f"crossing_{tra.name}",
                    args=[*args_cross_h, *tra.args_hv, *tra.args_h],
                    ret=None)
                for tra in all_trajectories]

    def __init__(self, *args, **kwargs):
        raise ValueError("LocalSegment is an abstract class and should not be instantiated")

    @classmethod
    def from_dict(cls, dct, **kwargs):
        """Returns the correct segment object from a dictionary in the same style as a HybridClass"""
        this_dct = dct.copy()
        this_cls = this_dct.pop('__class__')
        class_found = False
        for cls in all_segments:
            if this_cls == cls.__name__:
                class_found = True
                break
        if not class_found:
            raise ValueError(f"Not a segment class: {this_cls}")
        return cls(**this_dct, **kwargs)


# Sanity check to assert all segment types have crossing functions for all trajectories
def assert_localsegment_sources(seg):
    for tra in all_trajectories:
        header = f"/*gpufun*/\nvoid {seg.__name__}_crossing_{tra.name}({seg.__name__} seg, {xo_to_ctypes(args_cross_h)}, " \
               + f"{xo_to_ctypes(tra.args_hv)}, {xo_to_ctypes(tra.args_h)})"
        header_found = False
        for src in seg._extra_c_sources:
            if isinstance(src, str):
                if header in src:
                    header_found = True
                    break
            else:
                with open(src) as f:
                    if header in f.read():
                        header_found = True
                        break
        if not header_found:
            raise SystemError(f"Missing or corrupt C crossing function for {tra.__name__} in {seg.__name__}.")


for seg in all_segments:
    assert_localsegment_sources(seg)
    assert hasattr(seg, 'evaluate')
    assert hasattr(seg, 'get_vertices')
    assert hasattr(seg, 'max_crossings')


# Define common methods for all segments
def seg__eq__(self, other):
    """Check if two segments are equal"""
    return self.to_dict() == other.to_dict()

def to_dict(self):
    """Returns a dictionary in the same style as a HybridClass"""
    return {'__class__': self.__class__.__name__, **self._to_json()}

@classmethod
def from_dict(cls, dct, **kwargs):
    """Returns the object from a dictionary in the same style as a HybridClass"""
    this_dct = dct.copy()
    this_cls = this_dct.pop('__class__')
    if this_cls != cls.__name__:
        raise ValueError(f"Expected class {cls.__name__}, got {this_cls}")
    return cls(**this_dct, **kwargs)

def seg_round(self, val):
    """Built-in to provide rounding to Xcoll precision"""
    return round(val, -int(np.log10(XC_EPSILON)))

def is_open(self):
    """Check if the segment is an open segment"""
    return len(self.get_vertices()) == 1

def connection_to(self, other):
    """Get the point(s) at which the segment is connected to another segment"""
    overlap = [vert1 for vert1 in self.get_vertices()
                if np.any([np.allclose(vert1, vert2, atol=XC_EPSILON)
                            for vert2 in other.get_vertices()])]
    return overlap

def is_connected_to(self, other):
    """Check if this segment is connected to another segment"""
    return len(self.connection_to(other)) > 0

for seg in all_segments:
    seg.__eq__ = seg__eq__
    seg.to_dict = to_dict
    seg.from_dict = from_dict
    seg.round = seg_round
    seg.is_open = is_open
    seg.connection_to = connection_to
    seg.is_connected_to = is_connected_to


# Add some missing docstrings
for seg in all_segments:
    seg.evaluate.__doc__ = """Evaluate the segment over t using the parametric equation"""
    seg.get_vertices.__doc__ = """Get the vertices of the segment"""


# Function to get the maximum number of crossings for a given object type
def get_max_crossings(segments, trajectory=DriftTrajectory):
    if hasattr(segments, '__iter__') and all(isinstance(seg, all_segments) for seg in segments):
        max_crossings = 0
        for seg in segments:
            max_crossings += seg.max_crossings[trajectory]
        return max_crossings
    elif isinstance(segments, all_segments):
        return segments.max_crossings[trajectory]
    else:
        raise ValueError(f"Unexpected type for segments: {type(segments)}")
