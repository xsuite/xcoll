# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

import numpy as np

import xobjects as xo

from ..c_init import  XC_EPSILON
from ..trajectories import trajectories
from .line import LineSegment
from .halfopen_line import HalfOpenLineSegment
from .circular import CircularSegment
from .bezier import BezierSegment


all_segments = (LineSegment, HalfOpenLineSegment, CircularSegment, BezierSegment)

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



# Sanity check to assert all segment types have crossing functions for all trajectories
def assert_localsegment_sources(seg):
    for trajectory, vals in trajectories.items():
        c_types = ", ".join([f"{arg.get_c_type()} {arg.name}" for arg in vals['args']])
        header = f"/*gpufun*/\nvoid {seg.__name__}_crossing_{trajectory}({seg.__name__} seg, int8_t* n_hit, double* s, {c_types})"
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
            raise SystemError(f"Missing or corrupt C crossing function for {trajectory} in {seg.__name__}.")


for seg in all_segments:
    assert_localsegment_sources(seg)
    assert hasattr(seg, 'evaluate')
    assert hasattr(seg, 'get_vertices')


# Define common methods for all segments
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
    seg.round = seg_round
    seg.is_open = is_open
    seg.connection_to = connection_to
    seg.is_connected_to = is_connected_to


# Add some missing docstrings
for seg in all_segments:
    seg.evaluate.__doc__ = """Evaluate the segment over t using the parametric equation"""
    seg.get_vertices.__doc__ = """Get the vertices of the segment"""
