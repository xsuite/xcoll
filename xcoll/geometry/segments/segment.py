# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

import numpy as np

import xobjects as xo

from ..c_init import GeomCInit, PyMethod, XC_GEOM_EPSILON

from .line import LineSegment
from .halfopen_line import HalfOpenLineSegment
from .circular import CircularSegment
from .bezier import BezierSegment


all_segments = (LineSegment, HalfOpenLineSegment, CircularSegment, BezierSegment)


segment_methods = {
    'func_s': xo.Method(
        c_name=f"func_s",
        args=[xo.Arg(xo.Float64, name="t")],
        ret=xo.Arg(xo.Float64, name="s")),
    'func_x': xo.Method(
        c_name=f"func_x",
        args=[xo.Arg(xo.Float64, name="t")],
        ret=xo.Arg(xo.Float64, name="x")),
    'deriv_s': xo.Method(
        c_name=f"deriv_s",
        args=[xo.Arg(xo.Float64, name="t")],
        ret=xo.Arg(xo.Float64, name="s")),
    'deriv_x': xo.Method(
        c_name=f"deriv_x",
        args=[xo.Arg(xo.Float64, name="t")],
        ret=xo.Arg(xo.Float64, name="x"))
}


class LocalSegment(xo.UnionRef):
    """General segment, acting as a xobject-style parent class for all segment types"""
    _reftypes = all_segments
    _methods = list(segment_methods.values())

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


# Add kernels for func_ and deriv_ functions to all segments
def __getattr(self, attr):
    # Prepend the segment name to the kernel names to avoid duplication conflicts
    kernel_name = f"{self.__class__.__name__}_{attr}"
    if kernel_name in self._kernels:
        return PyMethod(kernel_name=kernel_name, element=self, element_name='seg')
    raise ValueError(f"Attribute {attr} not found in {self.__class__.__name__}")

for seg in all_segments:
    this_kernels = getattr(seg, '_kernels', {})
    _kernels = {key: xo.Kernel(c_name=f"{seg.__name__}_{val.c_name}",
                               ret=val.ret, args=[xo.Arg(xo.ThisClass, name="seg"), *val.args])
                for key, val in segment_methods.items()}
    this_kernels.update(_kernels)
    # Prepend the segment name to the kernel names to avoid duplication conflicts
    this_kernels = {f"{seg.__name__}_{key}": val for key, val in this_kernels.items()}
    seg._kernels = this_kernels
    seg.__getattr__ = __getattr
    seg._needs_compilation = True


# Define common methods for all segments
from ..trajectories.trajectory import __eq, __repr, to_dict, from_dict, __copy, __round

def is_open(self):
    """Check if the segment is an open segment"""
    return len(self.get_vertices()) == 1

def connection_to(self, other):
    """Get the point(s) at which the segment is connected to another segment"""
    overlap = [vert1 for vert1 in self.get_vertices()
                if np.any([np.allclose(vert1, vert2, atol=XC_GEOM_EPSILON)
                            for vert2 in other.get_vertices()])]
    return overlap

def is_connected_to(self, other):
    """Check if this segment is connected to another segment"""
    return len(self.connection_to(other)) > 0

def translate(self, ds, dx, *, inplace=False):
    """Translate the segment by (ds, dx). If `inplace` is False, the method returns a
    new segment, otherwise the segment is modified in place and nothing is returned.
    """
    if inplace:
        self._translate_inplace(ds, dx)
    else:
        new_seg = self.copy()
        new_seg._translate_inplace(ds, dx)
        return new_seg

def rotate(self, ps, px, angle, *, inplace=False):
    """Rotates the segment over an angle at a pivot point (ps, px). If `inplace` is False,
    the method returns a new segment, otherwise the segment is modified in place and
    nothing is returned.
    """
    if inplace:
        self._rotate_inplace(ps, px, angle)
    else:
        new_seg = self.copy()
        new_seg._rotate_inplace(ps, px, angle)
        return new_seg

for seg in all_segments:
    seg.name = seg.__name__.lower()[:-7]
    seg.__eq__ = __eq
    if not '__repr__' in seg.__dict__:
        seg.__repr__ = __repr
    if not '__str__' in seg.__dict__:
        seg.__str__ = __repr
    seg.to_dict = to_dict
    seg.from_dict = from_dict
    seg.copy = __copy
    seg.round = __round
    seg.is_open = is_open
    seg.connection_to = connection_to
    seg.is_connected_to = is_connected_to
    seg.translate = translate
    seg.rotate = rotate


# Sanity check to assert all segments have C code for func_ and deriv_ functions
def assert_segment_sources(tra):
    assert seg in all_segments
    name = seg.__name__
    for func in ['func_s', 'func_x', 'deriv_s', 'deriv_x']:
        header = f"/*gpufun*/\ndouble {name}_{func}({name} seg, double t)"
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
            raise SystemError(f"Missing or corrupt C function:  double {name}_{func}"
                            + f"({name} seg, double t).")

for seg in all_segments:
    assert_segment_sources(seg)
    assert hasattr(seg, 'get_vertices')
    assert hasattr(seg, '_translate_inplace')
    assert hasattr(seg, '_rotate_inplace')

# Add some missing docstrings
for seg in all_segments:
    seg.get_vertices.__doc__ = """Get the vertices of the segment"""


# Function to get the maximum number of crossings for a given object type
def get_max_crossings(segments, trajectory):
    if hasattr(segments, '__iter__') and all(isinstance(seg, all_segments) for seg in segments):
        max_crossings = 0
        for seg in segments:
            if hasattr(seg, '_max_crossings') and trajectory in seg._max_crossings:
                max_crossings += seg._max_crossings[trajectory]
            else:
                max_crossings += 2
        return max_crossings
    elif isinstance(segments, all_segments):
        return segments.max_crossings[trajectory]
    else:
        raise ValueError(f"Unexpected type for segments: {type(segments)}")
