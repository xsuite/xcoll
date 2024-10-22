# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

import numpy as np

import xobjects as xo

from ..segments import LocalSegment
from ..trajectories import trajectories
from .shape import Shape2D, _init_shape
from .shape_source import all_s_positions, shape_v_source


# TODO: need to rename kernels (remove Shape2DV_ prefix) when xobjects PR #109 is merged
class Shape2DV(xo.Struct):
    """Array of segments, representing an object in the geometry, with vertical limits"""
    segments = LocalSegment[:]
    vlimit   = xo.Float64[2]
    _seg_id  = xo.Int64  # This links the object to the correct array size for the crossings s

    _needs_compilation = True

    _depends_on = [Shape2D]
    _extra_c_sources = shape_v_source

    _kernels = {**{
                f"Shape2DV_crossing_{trajectory}": xo.Kernel(
                    c_name=f"Shape2DV_crossing_{trajectory}",
                    args=[
                        xo.Arg(xo.ThisClass, name="shape"),
                        xo.Arg(xo.Int8,    pointer=True, name="n_hit"),
                        xo.Arg(xo.Float64, pointer=True, name="s"),
                        *vals["args"],
                        *vals["args_vlimit_extra"]
                    ],
                    ret=None)
                for trajectory, vals in trajectories.items()},
                **{
                f"Shape2DV_crossing_{trajectory}_{s_pos}": xo.Kernel(
                    c_name=f"Shape2DV_crossing_{trajectory}_{s_pos}",
                    args=[
                        xo.Arg(xo.ThisClass, name="shape"),
                        *vals["args"],
                        *vals["args_vlimit_extra"],
                        *s_vals["args"]
                    ],
                    ret=xo.Arg(xo.Float64, pointer=False, name='s'))
                for s_pos, s_vals in all_s_positions.items() for trajectory, vals in trajectories.items()}
                }

    def __init__(self, segments=None, vlimit=None, **kwargs):
        if not segments:
            raise ValueError("Need to provide `segments`.")
        if not hasattr(segments, '__iter__'):
            raise ValueError("The variable `segments` should be a list of segments.")
        kwargs['segments'] = segments
        if not vlimit or isinstance(vlimit, str) \
        or not hasattr(vlimit, '__iter__') or len(vlimit) != 2:
            raise ValueError("Need to provide `vlimit` as [vmin, vmax].")
        kwargs['vlimit'] = vlimit
        _init_shape(self, **kwargs)

    def __repr__(self):
        shape2d = Shape2D.__repr__(self)
        return shape2d[:-1] + f", vlimit=[{self.vlimit[0]}, {self.vlimit[1]}])"
        # segs = ',\n          '.join([seg.__repr__() for seg in self])
        # return "Shape2DV([" + segs + "\n         ])"

    def __getitem__(self, i):
        return Shape2D.__getitem__(self, i)
        # return self.segments[i]

    def __iter__(self):
        return Shape2D.__iter__(self)
        # return iter(self.segments)

    def __len__(self):
        return Shape2D.__len__(self)
        # return len(self.segments)

    def __getattr__(self, attr):
        return Shape2D.__getattr__(self, attr)

    def is_composite(self):
        return Shape2D.is_composite(self)

    def is_open(self):
        return Shape2D.is_open(self)

    def get_shapes(self):
        for shape in self._shapes:
            yield Shape2DV([self.segments[i] for i in shape], vlimit=self.vlimit)

    def plot(self, axes=None):
        return Shape2D.plot(self, axes=axes)

    def plot3d(self, axes=None, num_points=100):
        return Shape2D.plot3d(self, axes=axes, num_points=num_points)
