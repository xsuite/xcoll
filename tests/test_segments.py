# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

import numpy as np

import xobjects as xo
from xobjects.test_helpers import for_all_test_contexts
import xcoll as xc
from xcoll.scattering_routines.geometry.c_init import XC_EPSILON


num_seg = 1000

@for_all_test_contexts
def test_line_segment(test_context):
    # Test the python API of LineSegment
    t = np.linspace(-4, 4, 4001)
    for s1, x1, s2, x2 in zip(*np.random.uniform(-8, 8, [4, num_seg])):
        # Test instantiation
        seg = xc.LineSegment(s1=s1, x1=x1, s2=s2, x2=x2, _context=test_context)
        print(seg)  # test __repr__
        assert np.allclose((seg.s1, seg.x1, seg.s2, seg.x2), (s1, x1, s2, x2), atol=XC_EPSILON)
        # Test vertices
        verts = seg.get_vertices()
        assert len(verts) == 2
        assert np.allclose(verts[0], (s1, x1), atol=XC_EPSILON)
        assert np.allclose(verts[1], (s2, x2), atol=XC_EPSILON)
        # Test evaluation
        s, x = seg.evaluate(t)
        # First and last evaluated points have to equal the vertices
        assert np.allclose((s[0], x[0], s[-1], x[-1]), (s1, x1, s2, x2), atol=XC_EPSILON)
        # The evaluated points' slope has to be correct
        slope = (x2 - x1) / (s2 - s1)
        assert np.allclose((x[1:] - x1) / (s[1:] - s1), slope, atol=XC_EPSILON)
        # The evaluated points have to be within the segment
        assert np.all(min(s1, s2) <= s)
        assert np.all(s <= max(s1, s2))
        assert np.all(min(x1, x2) <= x)
        assert np.all(x <= max(x1, x2))



@for_all_test_contexts
def test_halfopen_line_segment(test_context):
    # Test the python API of HalfOpenLineSegment
    t = np.linspace(-4, 4, 4001)
    for s1, x1, t1 in zip(*np.random.uniform(-8, 8, [3, num_seg])):
        # Test instantiation
        seg = xc.HalfOpenLineSegment(s=s1, x=x1, t=t1, _context=test_context)
        print(seg)  # test __repr__
        assert np.allclose((seg.s, seg.x, seg.t), (s1, x1, t1), atol=XC_EPSILON)
        # Test vertices
        verts = seg.get_vertices()
        assert len(verts) == 1
        assert np.allclose(verts[0], (s1, x1), atol=XC_EPSILON)
        # Test evaluation
        s, x = seg.evaluate(t)
        # First evaluated point has to equal the vertices
        assert np.allclose((s[0], x[0]), (s1, x1), atol=XC_EPSILON)
        # The evaluated points' slope has to be correct
        assert np.allclose((x[1:] - x1) / (s[1:] - s1), np.tan(t1), atol=XC_EPSILON)
        # The evaluated points have to be within the segment
        if np.cos(t1) >= 0:
            assert np.all(s1 <= s)
        else:
            assert np.all(s <= s1)
        if np.sin(t1) >= 0:
            assert np.all(x1 <= x)
        else:
            assert np.all(x <= x1)


@for_all_test_contexts
def test_circular_segment(test_context):
    # Test the python API of CircularSegment
    t = np.linspace(-4, 4, 4001)
    for R, s0, x0 , t1, t2 in zip(*np.random.uniform(-8, 8, [5, num_seg])):
        # Test instantiation
        R=abs(R)
        seg = xc.CircularSegment(R=R, s=s0, x=x0, t1=t1, t2=t2, _context=test_context)
        print(seg)  # test __repr__
        t1 = np.mod(t1 + np.pi, 2*np.pi) - np.pi  # Move to [-pi, pi]
        t2 = np.mod(t2 + np.pi, 2*np.pi) - np.pi  # Move to [-pi, pi]
        assert np.allclose((seg.R, seg.s, seg.x), (R, s0, x0), atol=XC_EPSILON)
        assert np.allclose((seg.t1, seg.t2), (t1, t2), atol=XC_EPSILON)
        # Test vertices
        verts = seg.get_vertices()
        assert len(verts) == 2
        assert np.allclose(verts[0], (s0 + R*np.cos(t1), x0 + R*np.sin(t1)), atol=XC_EPSILON)
        assert np.allclose(verts[1], (s0 + R*np.cos(t2), x0 + R*np.sin(t2)), atol=XC_EPSILON)
        # Test evaluation
        s, x = seg.evaluate(t)
        # First and last evaluated points are not necessarily close to the vertices, because
        # of cyclicity (as t was sampled from -pi to pi, but if e.g. -pi < t2 < t1, the t array
        # will start at -pi instead of t1).
        # The evaluated points have to lie on the circle
        assert np.allclose((s-s0)**2 + (x-x0)**2, R**2, atol=XC_EPSILON)
        # The evaluated points have to be within the segment
        angs = np.arctan2(x - x0, s - s0)
        if t1 < t2:
            assert np.all((t1 <= angs) & (angs <= t2))
        else:
            assert np.all((t1 <= angs) | (angs <= t2))


@for_all_test_contexts
def test_bezier_segment(test_context):
    t = np.linspace(-4, 4, 4001)
    for s1, x1, s2, x2, cs1, cx1, cs2, cx2 in zip(*np.random.uniform(-8, 8, [8, num_seg])):
        seg = xc.BezierSegment(s1=s1, x1=x1, s2=s2, x2=x2, cs1=cs1, cx1=cx1, cs2=cs2, cx2=cx2, _context=test_context)
        print(seg)  # test __repr__
        assert np.allclose((seg.s1, seg.x1, seg.s2, seg.x2, seg.cs1, seg.cx1, seg.cs2, seg.cx2),
                            (s1, x1, s2, x2, cs1, cx1, cs2, cx2), atol=XC_EPSILON)
        # Test vertices
        verts = seg.get_vertices()
        assert len(verts) == 2
        assert np.allclose(verts[0], (s1, x1), atol=XC_EPSILON)
        assert np.allclose(verts[1], (s2, x2), atol=XC_EPSILON)
        # Test evaluation
        s, x = seg.evaluate(t)
        # First and last evaluated points have to equal the vertices
        assert np.allclose((s[0], x[0], s[-1], x[-1]), (s1, x1, s2, x2), atol=XC_EPSILON)


# from xcoll.scattering_routines.geometry.trajectories import get_max_crossings
# from xobjects.test_helpers import for_all_test_contexts

# path = Path(__file__).parent / 'data_test_segments'
# # add kernels 
# context = xo.ContextCpu()
# context.add_kernels(
#     kernels={
#         "Segments_crossing_drift": xo.Kernel(
#             c_name="Segments_crossing_drift",
#             args=[
#                 xo.Arg(xcSegments, name="segs"),
#                 xo.Arg(xo.Int8,    pointer=True, name="n_hit"),
#                 xo.Arg(xo.Float64, pointer=True,  name="s"),
#                 xo.Arg(xo.Float64, pointer=False, name="s0"),
#                 xo.Arg(xo.Float64, pointer=False, name="x0"),
#                 xo.Arg(xo.Float64, pointer=False, name="m")
#             ],
#             ret=None,
#         )
#     }, save_source_as="temp.c"
# )


# # Segment(s) to be tested #############################
# line_segment      = xc.LineSegment(s1=0,x1=0,s2=0.2,x2=1)
# half_open_segment = xc.HalfOpenLineSegment(s=2, x=1, t=np.pi/4)
# bezier_segment    = xc.BezierSegment(s1=0.2, x1=1, s2=1+0.5*np.cos(3*np.pi/4), x2=0.5+0.5*np.sin(3*np.pi/4), cs1=0.5, cx1=2.5, cs2=1+0.5*np.cos(3*np.pi/4)-0.6, cx2=0.5+0.5*np.sin(3*np.pi/4)-0.6)

# half_circle       = xc.CircularSegment(R=0.5, s=1, x=0.5, t1=-np.pi/2, t2=np.pi/2)
# half_circle_2     = xc.CircularSegment(R=0.5, s=1, x=0.5, t1=np.pi/2, t2=3*np.pi/2)
# half_circle_3     = xc.CircularSegment(R=0.5, s=1, x=0.5, t1=0, t2=np.pi)
# half_circle_4     = xc.CircularSegment(R=0.5, s=1, x=0.5, t1=np.pi, t2=2*np.pi)

# seg        = xc.Segments([xc.LineSegment(s1=0, x1=0, s2=0.2, x2=1),
#                           xc.BezierSegment(s1=0.2, x1=1, s2=1+0.5*np.cos(3*np.pi/4), x2=0.5+0.5*np.sin(3*np.pi/4), cs1=0.5, cx1=2.5, cs2=1+0.5*np.cos(3*np.pi/4)-0.6, cx2=0.5+0.5*np.sin(3*np.pi/4)-0.6),
#                           xc.CircularSegment(R=0.5, s=1, x=0.5, t1=-np.pi/2, t2=3*np.pi/4),
#                           xc.LineSegment(s1=1, x1=0, s2=0, x2=0)])

# seg_bezier = xc.Segments([xc.BezierSegment(s1=0.2, x1=1, s2=1+0.5*np.cos(3*np.pi/4), x2=0.5+0.5*np.sin(3*np.pi/4), cs1=0.5, cx1=2.5, cs2=1+0.5*np.cos(3*np.pi/4)-0.6, cx2=0.5+0.5*np.sin(3*np.pi/4)-0.6),
#                           xc.BezierSegment(s1=1+0.5*np.cos(3*np.pi/4), x1=0.5+0.5*np.sin(3*np.pi/4), s2=1.2, x2=0.5, cs1=1.4, cx1=1.6, cs2=1.1, cx2=0.9),
#                           xc.BezierSegment(s1=1.2, x1=0.5, s2=0.5, x2=0.2, cs1=1.2, cx1=0.1, cs2=0.4, cx2=0.1),
#                           xc.BezierSegment(s1=0.5, x1=0.2, s2=0.2, x2=1, cs1=0.7, cx1=0.5, cs2=0.1, cx2=0.2)])

# # Test functions #####################################
# def test_line_segment(): 
#     _check_ref('line',line_segment)

# def test_half_open_segment():
#     _check_ref('half_open', half_open_segment)

# def test_bezier_segment():
#     _check_ref('bezier', bezier_segment)

# def test_half_circle():
#     _check_ref('half_circle', half_circle)
#     _check_ref('half_circle_2', half_circle_2)   
#     _check_ref('half_circle_3', half_circle_3)
#     _check_ref('half_circle_4', half_circle_4)

# def test_segments(): 
#     _check_ref('segments', seg)

# def test_segments_bezier():
#     _check_ref('bezier_segments', seg_bezier)

# def _check_ref(name, seg, _context=None):
#     with open(Path(path, 'ref_'+name+'.json'), 'r') as fid:
#         ref = json.load(fid)
#     ref_s     = np.array(ref['hit_s'])
#     ref_x     = np.array(ref['hit_x'])
#     ref_m     = np.array(ref['m'])
#     ref_point = np.array(ref['point'])

#     n_hit = np.zeros(1, dtype=np.int8)
#     s = np.zeros(get_max_crossings(seg, "drift"), dtype=np.float64)

#     s0, x0 = ref_point
#     n_hit[0] = 0
#     hits_s = np.array([], dtype=np.float64)
#     hits_x = np.array([], dtype=np.float64)
#     context.kernels.Segments_crossing_drift(segs=seg, n_hit=n_hit, s=s, s0=s0, x0=x0, m=ref_m[0])
#     print('yes')
#     hits_s = np.concatenate((hits_s, s[:n_hit[0]]))
#     hits_x = np.concatenate((hits_x, x0 + ref_m[0] * (s[:n_hit[0]] - s0)))
#     print('yes')
#     assert np.allclose(hits_s, ref_s, atol=1e-12, rtol=0), "s does not match."
#     assert np.allclose(hits_x, ref_x, atol=1e-12, rtol=0), "x does not match"
#     print(f"Test for {name} passed.")
