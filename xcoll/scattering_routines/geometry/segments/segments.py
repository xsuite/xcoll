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
from .segments_source import all_s_positions, segments_source, segments_vlimit_source, \
                        get_seg_ids, create_cases_in_source, assert_localsegment_sources

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


# =========================================
# == Segments class to create 2D objects ==
# =========================================

# TODO: need to rename kernels (remove Segments_ prefix) when xobjects PR #142 is merged
class Segments(xo.Struct):
    """Array of segments, representing an object in the geometry"""
    segments = LocalSegment[:]
    _seg_id  = xo.Int64  # This links the object to the correct array size for the crossings s
    _needs_compilation = True

    _extra_c_sources = segments_source

    _kernels = {**{
                f"Segments_crossing_{trajectory}": xo.Kernel(
                    c_name=f"Segments_crossing_{trajectory}",
                    args=[
                        xo.Arg(xo.ThisClass, name="segs"),
                        xo.Arg(xo.Int8,    pointer=True, name="n_hit"),
                        xo.Arg(xo.Float64, pointer=True, name="s"),
                        *vals["args"]
                    ],
                    ret=None)
                for trajectory, vals in trajectories.items()},
                **{
                f"Segments_crossing_{trajectory}_{s_pos}": xo.Kernel(
                    c_name=f"Segments_crossing_{trajectory}_{s_pos}",
                    args=[
                        xo.Arg(xo.ThisClass, name="segs"),
                        *vals["args"],
                        *s_vals["args"]
                    ],
                    ret=xo.Arg(xo.Float64, pointer=False, name='s'))
                for s_pos, s_vals in all_s_positions.items() for trajectory, vals in trajectories.items()}
                }

    def __init__(self, segments=None, **kwargs):
        if not segments:
            raise ValueError("Need to provide `segments`.")
        if not hasattr(segments, '__iter__'):
            raise ValueError("The variable `segments` should be a list of segments.")
        kwargs['segments'] = segments
        _init_segments_class(self, **kwargs)

    def __repr__(self):
        return f"Segments([{', '.join([seg.__class__.__name__ + '(...)' for seg in self])}])"

    def __getitem__(self, i):
        return self.segments[i]

    def __iter__(self):
        return iter(self.segments)

    def __getattr__(self, attr):
        # TODO: to be removed when xobjects PR #142 is merged
        if not attr.startswith(self.__class__.__name__):
            attr = f"{self.__class__.__name__}_{attr}"
        if attr in self._kernels:
            return PyMethod(kernel_name=attr, element=self)
        raise ValueError(f"Attribute {attr} not found in {self.__class__.__name__}")

    def plot(self, axes=None):
        if axes is None:
            import matplotlib.pyplot as plt
            _, axes = plt.subplots()
        s = []
        x = []
        if np.any([seg.__class__ == HalfOpenLineSegment or seg.__class__ == CircularSegment
                   for seg in self.segments]):
            t = np.linspace(-4, 4, 4000)
        else:
            t = np.linspace(0, 1, 4000)
        for seg in self.segments:
            this_s, this_x = seg.evaluate(t)
            s.append(this_s)
            x.append(this_x)
        axes.scatter(np.concatenate(s), np.concatenate(x), s=1)
        return axes.figure, axes


# ==============================================================
# == Segments class to create 2D objects with vertical limits ==
# ==============================================================

# TODO: need to rename kernels (remove SegmentsVLimit_ prefix) when xobjects PR #109 is merged
class SegmentsVLimit(xo.Struct):
    """Array of segments, representing an object in the geometry, with vertical limits"""
    segments = LocalSegment[:]
    vlimit   = xo.Float64[2]
    _seg_id  = xo.Int64  # This links the object to the correct array size for the crossings s
    _needs_compilation = True

    _depends_on = [Segments]
    _extra_c_sources = segments_vlimit_source

    _kernels = {**{
                f"SegmentsVLimit_crossing_{trajectory}": xo.Kernel(
                    c_name=f"SegmentsVLimit_crossing_{trajectory}",
                    args=[
                        xo.Arg(xo.ThisClass, name="segs"),
                        xo.Arg(xo.Int8,    pointer=True, name="n_hit"),
                        xo.Arg(xo.Float64, pointer=True, name="s"),
                        *vals["args"]
                    ],
                    ret=None)
                for trajectory, vals in trajectories.items()},
                **{
                f"SegmentsVLimit_crossing_{trajectory}_{s_pos}": xo.Kernel(
                    c_name=f"Segments_crossing_{trajectory}_{s_pos}",
                    args=[
                        xo.Arg(xo.ThisClass, name="segs"),
                        *vals["args"],
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
            raise ValueError("Need to provide `vlimit` as [ymin, ymax].")
        kwargs['vlimit'] = vlimit
        _init_segments_class(self, **kwargs)

    def __repr__(self):
        return f"SegmentsVLimit([{', '.join([seg.__class__.__name__ + '(...)' for seg in self])}], vlimit=[{self.vlimit[0], self.vlimit[1]}])"

    def __getitem__(self, i):
        return self.segments[i]

    def __iter__(self):
        return iter(self.segments)

    def __getattr__(self, attr):
        return Segments.__getattr__(self, attr)

    def plot(self, axes=None):
        return Segments.plot(self, axes=axes)


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
            seg.__class__._needs_compilation = True


class PyMethod:
    # Similar class as for the xt.BeamElement, but without the Metaclass magic
    # (and hence no need for PyMethodDescriptor)
    def __init__(self, kernel_name, element):
        self.kernel_name = kernel_name
        self.element = element

    def __call__(self, **kwargs):
        instance = self.element
        context = instance._context
        # import pdb; pdb.set_trace()

        if instance.__class__._needs_compilation:
            # We don't have the HybridClass metaclass magic, so we need to manually replace ThisClass
            for ker in instance.__class__._kernels.values():
                for arg in ker.args:
                    if arg.atype == xo.ThisClass:
                        arg.atype = instance.__class__
            instance.__class__.compile_kernels(instance, save_source_as="temp2.c")
            instance.__class__._needs_compilation = False
        kernel = context.kernels[self.kernel_name]
        kwargs['segs'] = instance

        return kernel( **kwargs)