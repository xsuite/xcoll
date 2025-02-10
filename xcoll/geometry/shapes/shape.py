# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

import xobjects as xo

from ..c_init import XC_EPSILON, PyMethod
from ..segments import LocalSegment, LineSegment, HalfOpenLineSegment, CircularSegment, get_max_crossings
from ..trajectories import all_trajectories, DriftTrajectory, args_cross_h
from .shape_source import all_s_positions, shape_source, get_seg_ids, create_cases_in_source


class Shape2D(xo.Struct):
    """Array of segments, representing objects in the geometry.
       The array can contain multiple shapes, as long as each shape is either a closed loop
       (in which case the segments don't need to be ordered) or an open loop which starts
       and ends with a half-open segment (in which the order of the segments matters)."""
    segments = LocalSegment[:]
    _seg_id  = xo.Int64  # This links the object to the correct array size for the crossings s

    _needs_compilation = True
    _extra_c_sources = shape_source

    # _kernels = {**{
    #             f"crossing_{tra.name}": xo.Kernel(
    #                 c_name=f"Shape2D_crossing_{tra.name}",
    #                 args=[
    #                     xo.Arg(xo.ThisClass, name="shape"),
    #                     *args_cross_h, *tra.args_hv, *tra.args_h
    #                 ],
    #                 ret=None)
    #             for tra in all_trajectories},
    #             **{
    #             f"crossing_{tra.name}_{s_pos}": xo.Kernel(
    #                 c_name=f"Shape2D_crossing_{tra.name}_{s_pos}",
    #                 args=[
    #                     xo.Arg(xo.ThisClass, name="shape"),
    #                     *tra.args_hv, *tra.args_h, *s_vals["args"]
    #                 ],
    #                 ret=xo.Arg(xo.Float64, pointer=False, name='s'))
    #             for s_pos, s_vals in all_s_positions.items() for tra in all_trajectories}
    #             }

    def __init__(self, segments=None, **kwargs):
        if not segments:
            raise ValueError("Need to provide `segments`.")
        if not isinstance(segments, (list, set, tuple, np.ndarray, pd.Series)):
            raise ValueError("The variable `segments` should be a list of segments.")
        kwargs['segments'] = segments
        _init_shape(self, **kwargs)

    def __repr__(self):
        cls = self.__class__.__name__
        sp = (1 + len(cls) )* ' '
        segs = ('\n' + sp + '], [\n ' + sp).join([
                    (',\n '+ sp).join([self[i].__repr__() for i in shape])
                for shape in self._shapes])
        sp_last = "\n" + sp if len(self._shapes) > 1 else ""
        return cls + "([" + segs + sp_last +"])"

    def __getitem__(self, i):
        return self.segments[i]

    def __iter__(self):
        return iter(self.segments)

    def __len__(self):
        return len(self.segments)

    def __getattr__(self, attr):
        if attr in self._kernels:
            return PyMethod(kernel_name=attr, element=self, element_name='shape')
        raise ValueError(f"Attribute {attr} not found in {self.__class__.__name__}")

    def __eq__(self, other):
        """Check if two shapes are equal"""
        this_dct = self.to_dict()
        other_dct = other.to_dict()
        # Order of segments is not important
        this_dct['segments'] = set(this_dct['segments'])
        other_dct['segments'] = set(other_dct['segments'])
        return this_dct == other_dct

    def to_dict(self):
        """Returns a dictionary in the same style as a HybridClass"""
        this_dict = {'__class__': self.__class__.__name__}
        this_dict.update(self._to_json())
        this_dict['segments'] = [{'__class__': kk, **vv} for kk, vv in this_dict['segments']]
        return this_dict

    @classmethod
    def from_dict(cls, dct, **kwargs):
        """Returns the object from a dictionary in the same style as a HybridClass"""
        this_dct = dct.copy()
        this_cls = this_dct.pop('__class__')
        if this_cls != cls.__name__:
            raise ValueError(f"Expected class {cls.__name__}, got {this_cls}")
        this_dct['segments'] = [LocalSegment.from_dict(seg) for seg in this_dct['segments']]
        return cls(**this_dct, **kwargs)

    def is_composite(self):
        """Returns True if the shape is a composite of multiple shapes"""
        return len(self._shapes) > 1

    def is_open(self):
        """Returns for each sub-shape whether or not it is open"""
        return [self[shape[0]].is_open() for shape in self._shapes]

    def get_shapes(self):
        """Returns the shape's sub-shapes"""
        if self.is_composite():
            for shape in self._shapes:
                yield Shape2D([self[i] for i in shape])
        else:
            yield self

    def get_vertex_tree(self):
        """Returns for each vertex the indices of the segments that are connected to it"""
        segs = {i: seg for i, seg in enumerate(self)}
        vertices = {}
        for i, seg in segs.items():
            for j in range(i+1, len(segs)):
                for vert in seg.connection_to(segs[j]):
                    # Small hack to find the vertex up to numerical precision
                    vert_found = False
                    for vert_key in vertices.keys():
                        if np.allclose(vert, vert_key, atol=XC_EPSILON):
                            vert_found = True
                            break
                    if vert_found:
                        vertices[vert_key].add(i)
                        vertices[vert_key].add(j)
                    else:
                        vertices[vert] = {i, j}
        return vertices

    def get_vertices(self, *, smooth_points=None, scale_open_points=None):
        """Returns the coordinates of the vertices of the shapes.
           If scale_open_points is provided, the open segments get an extra vertex
           at this distance, in the direction of the segment.
           If smooth_points is provided, curved segments are interpolated."""
        all_coords = []
        if smooth_points is True:
            smooth_points = 20
        elif smooth_points is None or smooth_points is False:
            smooth_points = 0
        elif smooth_points > 0 and smooth_points < 2:
            smooth_points += 2
        smooth_points = int(smooth_points)
        if scale_open_points is True:
            scale_open_points = 1
        elif scale_open_points is None or scale_open_points is False:
            scale_open_points = 0
        for i, (shape, is_open) in enumerate(zip(self._shapes, self.is_open())):
            coords = []
            if is_open:
                for seg1, seg2 in zip(shape[:-1], shape[1:]):
                    coords.append(*self[seg1].connection_to(self[seg2]))
                    _interpolate(self[seg2], coords, smooth_points)
                if scale_open_points:
                    assert scale_open_points > 0
                    coords_initial = tuple(coord[0] for coord in self[shape[0]].evaluate(scale_open_points))
                    coords.insert(0, coords_initial)
                    coords_final = tuple(coord[0] for coord in self[shape[-1]].evaluate(scale_open_points))
                    coords.append(coords_final)
            elif len(shape) == 1:
                seg = self[shape[0]]
                if seg.is_open():
                    coords.append(seg.get_vertices()[0])
                    if scale_open_points:
                        assert scale_open_points > 0
                        coords.append(tuple(coord[0] for coord in self[shape[0]].evaluate(scale_open_points)))
                else:
                    coords.append(seg.get_vertices()[0])
                    _interpolate(seg, coords, smooth_points)
                    coords.append(seg.get_vertices()[-1])
            elif len(shape) == 2:
                for i, vert in enumerate(self[shape[0]].connection_to(self[shape[1]])):
                    coords.append(vert)
                    _interpolate(self[shape[i]], coords, smooth_points)
                coords.append(coords[0])
            else:
                for seg1, seg2 in zip([shape[-1], *shape[:-1]], shape):
                    coords.append(*self[seg1].connection_to(self[seg2]))
                    _interpolate(self[seg2], coords, smooth_points)
                coords.append(coords[0])
            all_coords.append(coords)
        return all_coords

    def plot(self, *, axes=None, smooth_points=100, scale_open_points=0.6):
        """Returns (fig, axes) with a plot of the shape. Axes can optionally be provided."""
        if axes is None:
            _, axes = plt.subplots(clear=True)
        colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
        vertices = self.get_vertices(smooth_points=smooth_points, scale_open_points=scale_open_points)
        for i, (verts, shape) in enumerate(zip(vertices, self._shapes)):
            c = colors[i % len(colors)]
            axes.plot(*np.array(verts).T, c=c)
            for seg in shape:
                if self[seg].is_open():
                    s, x = self[seg].evaluate([scale_open_points, 1])
                    axes.plot(s, x, c=c, ls='--')
        axes.set_aspect('equal')
        return axes.figure, axes

    def plot3d(self, *, axes=None, smooth_points=100, scale_open_points=0.6):
        """Returns (fig, axes) with a 3d plot of the shape. Axes can optionally be provided."""
        if axes is None:
            fig = plt.figure(clear=True)
            axes = fig.add_subplot(projection='3d')
            axes.view_init(elev=25, azim=-110)
            axes.set_xlabel('s')
            axes.set_ylabel('x')
        vlimit = 'vlimit' in [field.name for field in self._fields]
        min_s = 0
        max_s = 0
        min_x = 0
        max_x = 0
        vmin = self.vlimit[0] if vlimit else -0.2
        vmax = self.vlimit[1] if vlimit else 0.2
        colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
        vertices = self.get_vertices(smooth_points=smooth_points, scale_open_points=scale_open_points)
        for i, verts in enumerate(vertices):
            c = colors[i % len(colors)]
            s = np.array(verts).T[0]
            x = np.array(verts).T[1]
            shape_B = np.array([s, x, [vmin for _ in range(len(s))]]).T
            shape_T = np.array([s, x, [vmax for _ in range(len(s))]]).T
            axes.plot_surface(*np.array([shape_B, shape_T]).T, alpha=0.35, color=c)
            if vlimit:
                face_B = Poly3DCollection([shape_B], alpha=0.35, color=c)
                face_T = Poly3DCollection([shape_T], alpha=0.35, color=c)
                axes.add_collection3d(face_B)
                axes.add_collection3d(face_T)
            min_s = min(min_s, min(s))
            max_s = max(max_s, max(s))
            min_x = min(min_x, min(x))
            max_x = max(max_x, max(x))
        axes.set_box_aspect([max_s-min_s, max_x-min_x, vmax-vmin])
        return axes.figure, axes


def _init_shape(shape, **kwargs):
    # Each different object type will get its own seg_id, by inspecting the source code
    # First we check if code for this object type already exists
    seg_ids = get_seg_ids(shape)
    max_crossings = get_max_crossings(kwargs['segments'])
    add_code = False
    if max_crossings in seg_ids:
        kwargs['_seg_id'] = seg_ids[max_crossings]
    else:
        add_code = True
        kwargs['_seg_id'] = max(seg_ids.values()) + 1 if len(seg_ids) > 0 else 0
    super(shape.__class__, shape).__init__(**kwargs)
    # Sort the segments into shapes and verify the structure
    shape._shapes = _get_shapes(shape)
    # Shapes are defined anti-clockwise
    if len(shape._shapes) > 1:
        shape._shapes = [list(reversed(shape)) if is_clockwise(verts) else shape
                        for shape, verts in zip(shape._shapes, shape.get_vertices())]
    if add_code:
        for tra in all_trajectories:
            create_cases_in_source(shape, tra)
            shape.__class__._needs_compilation = True


def _get_shapes(shape):
    # We create the vertex tree to verify the shapes are well-formed
    verts = shape.get_vertex_tree()
    for vert, vert_segs in verts.items():
        # Verify, for each vertex, how many segments are connected to it
        if len(vert_segs) == 1 :
            raise ShapeMalformedError(shape, f"Fatal: {vert} connected to one segment. " \
                                            + "This cannot happen (as it should not count " \
                                            + "as a connection). Check the code.", vert_segs)
        elif len(vert_segs) > 2:
            raise ShapeMalformedError(shape, f"Vertex {vert} connected to more than 2 segments: " \
                                           + f"{vert_segs}!", vert_segs)
    for seg_id in range(len(shape)):
        # Verify, for each segment, that all of its vertices are connected to another segment
        n_vert = len([ss for _, ss in verts.items() if seg_id in ss])
        if n_vert == 0 and len(verts) > 0:
            # A disconnected segment is only allowed if it is alone (for testing)
            raise ShapeMalformedError(shape, f"Segment {seg_id} is completely disconnected!", seg_id)
        if n_vert == 1 and not shape[seg_id].is_open():
            raise ShapeMalformedError(shape, f"Cycle is not closed: Segment {seg_id} is " \
                                            + "connected to only one other segment!", seg_id)
        if n_vert > 3 :
            raise ShapeMalformedError(shape, f"Fatal: Segment {seg_id} has more than 2 " \
                                            + "vertices. This is not supported. Check the code.", seg_id)
    # Next we loop over all cycles of segments
    if len(verts) == 0:
        # Single segment
        shapes = [[0]]
    else:
        shapes = []
        # We start with open segments, as the open cycles depend on the correct order
        open_segs = [seg_id for seg_id, seg in enumerate(shape) if seg.is_open()]
        closed_segs = [seg_id for seg_id, seg in enumerate(shape) if not seg.is_open()]
        for start_id in open_segs + closed_segs:
            if start_id in [seg for group in shapes for seg in group]:
                # Already in a group
                continue
            group = []
            visited = set()
            next_id = {start_id}
            while len(next_id) > 0:
                assert len(next_id) == 1
                next_id = next_id.pop()
                group.append(next_id)
                if next_id != start_id and shape[next_id].is_open():
                    # End of an open cycle
                    assert shape[start_id].is_open()
                    break
                # We create a dictionary of the unvisited vertices that connect to the segment next_id
                # If there are multiple (like for the initial choice) we take the one with the smallest index
                dd = {vv: sum(ss) for vv, ss in verts.items() if vv not in visited and next_id in ss}
                next_vert = min(dd, key=dd.get)
                # Get the id of the other segment that connects to next_vert
                next_id = verts[next_vert] - set(group)
                visited.add(next_vert)
            shapes.append(group)
    return shapes


def is_clockwise(vertices):
    """Returns True if the vertices are defined in a clockwise order"""
    if vertices[-1] != vertices[0]:
        vertices.append(vertices[0])
    area = 0
    for i in range(len(vertices)-1):
        x1, y1 = vertices[i]
        x2, y2 = vertices[i + 1]
        area += (x2 - x1) * (y2 + y1)
    # If area is positive, it's clockwise, otherwise it's counter-clockwise
    return area > 0


def _interpolate(segment, coords, smooth_points):
    if smooth_points:
        if isinstance(segment, CircularSegment):
            t1 = segment.t1
            t2 = segment.t2
            if t2 < t1:
                t2 += 2*np.pi
            t = np.linspace(t1, t2, smooth_points)
        elif isinstance(segment, (LineSegment, HalfOpenLineSegment)):
            return
        else:
            t = np.linspace(0, 1, smooth_points)
        s, x = segment.evaluate(t)
        interp = [(ss,xx) for ss,xx in zip(s,x)]
        if np.allclose(interp[0], coords[-1], atol=XC_EPSILON):
            coords.extend(interp[1:-1])
        else:
            coords.extend(list(reversed(interp[1:-1])))


class ShapeMalformedError(ValueError):
    def __init__(self, shape=None, message='', bad_seg=[]):
        if shape:
            segs = {i: seg for i, seg in enumerate(shape)}
            message += "\n" + "\n".join([f"{i}: {seg}" for i, seg in segs.items()])
            if not hasattr(bad_seg, '__iter__'):
                bad_seg = [bad_seg]
            plt.close()
            fig, ax = plt.subplots()
            for i, seg in enumerate(shape):
                t = np.linspace(-4, 4, 1000)
                s, x = seg.evaluate(t)
                c = 'tab:blue'
                if i in bad_seg:
                    c = 'tab:red'
                ax.plot(s, x, c=c)
            fig.show()
        super().__init__(message)
