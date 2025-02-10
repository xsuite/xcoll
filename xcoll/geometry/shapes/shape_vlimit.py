# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

import numpy as np
import pandas as pd

import xobjects as xo

from ..segments import LocalSegment
from ..trajectories import all_trajectories, args_cross_h
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

    # _kernels = {**{
    #             f"Shape2DV_crossing_{tra.name}": xo.Kernel(
    #                 c_name=f"Shape2DV_crossing_{tra.name}",
    #                 args=[
    #                     xo.Arg(xo.ThisClass, name="shape"),
    #                     *args_cross_h, *tra.args_hv, *tra.args_h, *tra.args_v
    #                 ],
    #                 ret=None)
    #             for tra in all_trajectories},
    #             **{
    #             f"Shape2DV_crossing_{tra.name}_{s_pos}": xo.Kernel(
    #                 c_name=f"Shape2DV_crossing_{tra.name}_{s_pos}",
    #                 args=[
    #                     xo.Arg(xo.ThisClass, name="shape"),
    #                     *tra.args_hv, *tra.args_h, *tra.args_v, *s_vals["args"]
    #                 ],
    #                 ret=xo.Arg(xo.Float64, pointer=False, name='s'))
    #             for s_pos, s_vals in all_s_positions.items() for tra in all_trajectories}
    #             }

    def __init__(self, segments=None, vlimit=None, **kwargs):
        if not segments:
            raise ValueError("Need to provide `segments`.")
        if not isinstance(segments, (list, set, tuple, np.ndarray, pd.Series)):
            raise ValueError("The variable `segments` should be a list of segments.")
        kwargs['segments'] = segments
        if not vlimit or isinstance(vlimit, str) \
        or not hasattr(vlimit, '__iter__') or len(vlimit) != 2:
            raise ValueError("Need to provide `vlimit` as [vmin, vmax].")
        kwargs['vlimit'] = vlimit
        _init_shape(self, **kwargs)

    def __repr__(self):
        shape2d = Shape2D.__repr__(self)
        sp_last = "\n         ], " if len(self._shapes) == 1 else "], "
        return shape2d[:-2] + sp_last + f"vlimit=[{self.vlimit[0]}, {self.vlimit[1]}])"

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

    def to_dict(self, attr):
        return Shape2D.to_dict(self, attr)

    @classmethod
    def from_dict(cls, dct):
        return Shape2D.from_dict(cls, dct)

    def is_composite(self):
        return Shape2D.is_composite(self)

    def is_open(self):
        return Shape2D.is_open(self)

    def get_shapes(self):
        if self.is_composite():
            for shape in self._shapes:
                yield Shape2DV([self[i] for i in shape], vlimit=list(self.vlimit))
        else:
            yield self

    def get_vertex_tree(self):
        return Shape2D.get_vertex_tree(self)

    def get_vertices(self, **kwargs):
        return Shape2D.get_vertices(self, **kwargs)

    def plot(self, **kwargs):
        return Shape2D.plot(self, **kwargs)

    def plot3d(self, **kwargs):
        return Shape2D.plot3d(self, **kwargs)


Shape2DV.is_composite.__doc__ = Shape2D.is_composite.__doc__
Shape2DV.is_open.__doc__ = Shape2D.is_open.__doc__
Shape2DV.get_shapes.__doc__ = Shape2D.get_shapes.__doc__
Shape2DV.get_vertex_tree.__doc__ = Shape2D.get_vertex_tree.__doc__
Shape2DV.get_vertices.__doc__ = Shape2D.get_vertices.__doc__
Shape2DV.plot.__doc__ = Shape2D.plot.__doc__
Shape2DV.plot3d.__doc__ = Shape2D.plot3d.__doc__
