# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

import numpy as np
import matplotlib.pyplot as plt

import xobjects as xo

from ..segments import LocalSegment, CircularSegment
from ..trajectories import trajectories, get_max_crossings
from .shape_source import all_s_positions, shape_source, get_seg_ids, create_cases_in_source


# TODO: need to rename kernels (remove Shape2D_ prefix) when xobjects PR #142 is merged
class Shape2D(xo.Struct):
    """Array of segments, representing an object in the geometry"""
    segments = LocalSegment[:]
    _seg_id  = xo.Int64  # This links the object to the correct array size for the crossings s

    _needs_compilation = True
    _extra_c_sources = shape_source

    _kernels = {**{
                f"Shape2D_crossing_{trajectory}": xo.Kernel(
                    c_name=f"Shape2D_crossing_{trajectory}",
                    args=[
                        xo.Arg(xo.ThisClass, name="shape"),
                        xo.Arg(xo.Int8,    pointer=True, name="n_hit"),
                        xo.Arg(xo.Float64, pointer=True, name="s"),
                        *vals["args"]
                    ],
                    ret=None)
                for trajectory, vals in trajectories.items()},
                **{
                f"Shape2D_crossing_{trajectory}_{s_pos}": xo.Kernel(
                    c_name=f"Shape2D_crossing_{trajectory}_{s_pos}",
                    args=[
                        xo.Arg(xo.ThisClass, name="shape"),
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
        _init_shape(self, **kwargs)

    def __repr__(self):
        cls = self.__class__.__name__
        sp = (1 + len(cls) )* ' '
        segs = ('\n' + sp + '], [\n ' + sp).join([
                    (',\n '+ sp).join([self[i].__repr__() for i in shape])
                for shape in self._shapes])
        return cls + "([" + segs + "\n" + sp + "])"

    def __getitem__(self, i):
        return self.segments[i]

    def __iter__(self):
        return iter(self.segments)

    def __len__(self):
        return len(self.segments)

    def __getattr__(self, attr):
        if not attr.startswith(self.__class__.__name__):
            # TODO: to be removed when xobjects PR #142 is merged
            attr = f"{self.__class__.__name__}_{attr}"
        if attr in self._kernels:
            return PyMethod(kernel_name=attr, element=self)
        raise ValueError(f"Attribute {attr} not found in {self.__class__.__name__}")

    def is_composite(self):
        return len(self._shapes) > 1

    def is_open(self):
        return [self[shape[0]].is_open() for shape in self._shapes]

    def get_shapes(self):
        for shape in self._shapes:
            yield Shape2D([self[i] for i in shape])

    def plot(self, axes=None):
        if axes is None:
            _, axes = plt.subplots()
        colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
        for i, shape in enumerate(self._shapes):
            color = colors[i % len(colors)]
            for j in shape:
                seg = self.segments[j]
                if seg.__class__ == CircularSegment:
                    t = np.linspace(-np.pi, np.pi, 1000)
                elif seg.is_open():
                    t = np.linspace(0, 0.6, 1000)
                else:
                    t = np.linspace(0, 1, 1000)
                s, x = seg.evaluate(t)
                axes.plot(s, x, c=color)
                if seg.is_open():
                    t = np.linspace(0.6, 1, 1000)
                    s, x = seg.evaluate(t)
                    axes.plot(s, x, c=color, ls='--')
        axes.set_aspect('equal')
        return axes.figure, axes

    def plot3d(self, axes=None, num_points=100):
        if axes is None:
            import matplotlib.pyplot as plt
            fig = plt.figure(clear=True)
            axes = fig.add_subplot(projection='3d')


def _init_shape(shape, **kwargs):
    # Each different object type will get its own seg_id, by inspecting the source code
    # First we check if code for this object type already exists
    seg_ids = get_seg_ids(shape)
    max_crossings = get_max_crossings(kwargs['segments'], 'drift')  # test with drift to check if seg_id already has source
    add_code = False
    if max_crossings in seg_ids:
        kwargs['_seg_id'] = seg_ids[max_crossings]
    else:
        add_code = True
        kwargs['_seg_id'] = max(seg_ids.values()) + 1 if len(seg_ids) > 0 else 0
    allow_single = kwargs.pop('allow_single', False)
    super(shape.__class__, shape).__init__(**kwargs)
    shape._shapes = _get_shapes(shape, allow_single=allow_single)
    if add_code:
        for trajectory in trajectories.keys():
            create_cases_in_source(shape, trajectory)
            shape.__class__._needs_compilation = True

# Recursively find all connected segments
# Note that this is order-dependent, i.e. we cannot find the cycle [0,1,2,0] when starting from 1 or 2
def _get_cycle(i, segs, group=None, *, allow_single=False):
    if group is None:
        # Start of loop
        group = [i]
    elif segs[i].is_open():
        # We have finished an open cycle!
        # Sanity check: if the last segment of the cycle is HalfOpen, the first one must be too
        if not segs[group[0]].is_open() or i == group[0]:
            raise ValueError(f"Shape malformed: HalfOpen segment in the middle of cycle (Segment IDs {group})!")
        return group
    found_next = False
    for next_i in range(i+1, len(segs)):
        if segs[i].is_connected_to(segs[next_i]):
            found_next = True
            break
    if found_next:
        group.append(next_i)
        return _get_cycle(next_i, segs, group, allow_single=allow_single)
    else:
        # Check for errors
        if len(group) == 1 and not allow_single:
            raise ValueError(f"Shape malformed: Segment {i} is not connected to any other segment! " \
                          + "If this is intended, set allow_single=True.")
        elif len(group) == 2:
            if len(segs[i].connection_to(segs[group[0]])) != 2:
                raise ValueError("Shape malformed: A shape of two segments is only allowed if both vertices " \
                              + f"connect (Segment IDs {group})!")
        elif not segs[i].is_connected_to(segs[group[0]]):
            raise ValueError(f"Shape malformed: Cycle (Segment IDs {group}) is not closed!")
        # We finished a closed cycle!
        return group

def _get_shapes(shape, allow_single=False):
    segs = {i: seg for i, seg in enumerate(shape)}
    shapes = []
    visited = set()
    for i, seg in segs.items():
        if i not in visited:
            cycle = _get_cycle(i, shape, allow_single=allow_single)
            visited.update(cycle)
            shapes.append(cycle)
    return shapes


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
            instance.__class__.compile_kernels(instance) #, save_source_as="temp.c")
            instance.__class__._needs_compilation = False
        kernel = context.kernels[self.kernel_name]
        kwargs['shape'] = instance

        return kernel( **kwargs)
