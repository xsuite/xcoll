# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

import numpy as np

import xobjects as xo

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
