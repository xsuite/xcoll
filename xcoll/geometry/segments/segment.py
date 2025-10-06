# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

import numpy as np

import xobjects as xo

from ..c_init import PyMethod, XC_GEOM_EPSILON

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
        ret=xo.Arg(xo.Float64, name="x")),
    'update_box': xo.Method(
        c_name=f"update_box",
        args=[xo.Arg(xo.Float64, name="t1"),
              xo.Arg(xo.Float64, name="t2")],
        ret=None),
    'bounded_below': xo.Method(
        c_name=f"bounded_below",
        args=[],
        ret=xo.Arg(xo.Int8, name="bounded_below")),
    'bounded_above': xo.Method(
        c_name=f"bounded_above",
        args=[],
        ret=xo.Arg(xo.Int8, name="bounded_above"))
}


class LocalSegment(xo.UnionRef):
    """General segment, acting as a xobject-style parent class for all segment types"""
    _reftypes = all_segments
    _methods = list(segment_methods.values())

    #def __init__(self, *args, **kwargs):
    #    raise ValueError("LocalSegment is an abstract class and should not be instantiated")

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


def is_open(self):
    """Check if the segment is an open segment"""
    return not bool(self.bounded_above()) or not bool(self.bounded_below())

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

def plot(self, t1=0, t2=1, ax=None, plot_bounding_box=True, plot_control_points=True):
    """Plot the segment and its bounding box"""
    import matplotlib.pyplot as plt
    if ax is None:
        fig, ax  = plt.subplots(figsize=(8, 8))
    else:
        fig = ax.figure

    # Plot the seg
    t_values = np.linspace(0, 1, 100)
    s_values = np.array([self.func_s(t=t) for t in t_values])
    x_values = np.array([self.func_x(t=t) for t in t_values])
    ax.plot(s_values, x_values, 'b-', label=f"{self.name} segment")

    # Get and plot the bounding box
    t1 = max(t1, 0)
    if not self.is_open():
        t2 = min(t2, 1)
    if plot_bounding_box:
        vertices = np.vstack([np.array(self.box.vertices), 
                              np.array(self.box.vertices)[0]])
        ax.plot(vertices.T[0], vertices.T[1], 'k--', label='Bounding Box')

    # Get vertices and control points
    s_start, x_start = self.func_s(t=0), self.func_x(t=0)
    s_t1, x_t1 = self.func_s(t=t1), self.func_x(t=t1)
    if not self.is_open():
        s_end, x_end = self.func_s(t=1), self.func_x(t=1)
    s_t2, x_t2 = self.func_s(t=t2), self.func_x(t=t2)
    cp = self.get_control_points()

    # Plot the control lines
    if cp and plot_control_points:
        ax.plot([s_start, cp[0][0]], [x_start, cp[0][1]], c='lightgray', lw=1)
        if not self.is_open():
            ax.plot([cp[-1][0], s_end], [cp[-1][1], x_end], c='lightgray', lw=1)

    # Plot the vertices and control points
    ax.plot([s_t1, s_t2], [x_t1, x_t2], 'bo')
    ax.plot([s_start], [x_start], 'go', label='Endpoints')
    if not self.is_open():
        ax.plot([s_end], [x_end], 'go', label='Endpoints')
    if plot_control_points:
        for s, x in cp:
            ax.plot([s], [x], 'ro', label='Control Points')

    ax.set_xlabel('s')
    ax.set_ylabel('x')
    ax.set_aspect('equal', 'box')
    return fig, ax

@classmethod
def _inspect(cls, plot_bounding_box=True, plot_control_points=True, **kwargs):
    # Quick method to plot the segment and its bounding box in an interactive way, for testing
    # kwargs needs to have all arguments as keys, and val should be [min, max, initial_value]
    import matplotlib.pyplot as plt
    from matplotlib.widgets import Slider
    plt.ion()  # Enable interactive mode
    fig, ax = cls(**{kk: vv[-1] for kk, vv in kwargs.items()}).plot(
                        plot_bounding_box=plot_bounding_box,
                        plot_control_points=plot_control_points)
    plt.subplots_adjust(left=0.1, bottom=0.1+0.025*len(kwargs))  # Space for the sliders

    state = {'ax': ax} # Use a mutable object to store the axis (ax)

    def update_plot(val):
        this_kwargs = {arg: sliders[arg].val for arg in kwargs}
        # Clear the existing plot and redraw
        state['ax'].clear()
        fig, state['ax'] = cls(**this_kwargs).plot(t1=sliders['t1'].val, t2=sliders['t2'].val, ax=state['ax'],
                        plot_bounding_box=plot_bounding_box,
                        plot_control_points=plot_control_points)
        plt.draw()

    all_kwargs = kwargs
    all_kwargs['t1'] = [0, 1, 0]
    all_kwargs['t2'] = [0, 1, 1]
    ax_sliders = {arg: plt.axes([0.1, 0.025*(len(all_kwargs)-i-1), 0.8, 0.03], facecolor='lightgrey')
                  for i, arg in enumerate(all_kwargs)}
    sliders = {arg: Slider(ax_sliders[arg], arg, val[0], val[1], valinit=val[2])
               for arg, val in all_kwargs.items()}
    for slider in sliders.values():
        slider.on_changed(update_plot)
    plt.show(block=False)  # Make plt.show() non-blocking
    while plt.fignum_exists(fig.number):
        plt.pause(0.1) # Keep the plot alive and responsive

for seg in all_segments:
    # old_init = seg.__dict__.get('__init__', None)
    # def __init(self, *args, **kwargs):
    #     if old_init:
    #         old_init(self, *args, **kwargs)
    #     else:
    #         super().__init__(*args, **kwargs)
    #     if self.is_open():
    #         assert len(self.get_vertices()) == 1
    #     else:
    #         assert len(self.get_vertices()) == 2
    # seg.__init__ = __init
    seg.is_open = is_open
    seg.connection_to = connection_to
    seg.is_connected_to = is_connected_to
    seg.translate = translate
    seg.rotate = rotate
    seg.plot = plot
    seg._inspect = _inspect

# Add some missing docstrings
for seg in all_segments:
    seg.get_vertices.__doc__ = """Get the vertices of the segment"""


# Sanity check to assert all segments have C code for func_ and deriv_ functions
def assert_segment_sources(tra):
    assert seg in all_segments
    def _check_source(header, seg):
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
            raise SystemError(f"Missing or corrupt C function in segment {seg.__name__}:  {header}.")
    name = seg.__name__
    for func in ['bounded_above', 'bounded_below']:
        header = f"/*gpufun*/\nint8_t {name}_{func}({name} seg)"
        _check_source(header, seg)
    for func in ['func_s', 'func_x', 'deriv_s', 'deriv_x']:
        header = f"/*gpufun*/\ndouble {name}_{func}({name} seg, double t)"
        _check_source(header, seg)

def get_control_points(self):
    return ()

for seg in all_segments:
    assert_segment_sources(seg)
    assert hasattr(seg, 'get_vertices')
    assert hasattr(seg, '_translate_inplace')
    assert hasattr(seg, '_rotate_inplace')
    if not hasattr(seg, 'get_control_points'):
        seg.get_control_points = get_control_points


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
