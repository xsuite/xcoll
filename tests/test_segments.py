# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

import sys
import numpy as np
import pandas as pd
import json
from pathlib import Path

import xobjects as xo
import xcoll as xc
from xcoll.scattering_routines.geometry.segments import *
from xcoll.scattering_routines.geometry.c_init import *
from xcoll.scattering_routines.geometry.segments.trajectories import get_max_crossings
from xobjects.test_helpers import for_all_test_contexts

path = Path(__file__).parent / 'data_test_segments'
# add kernels 
context = xo.ContextCpu()
context.add_kernels(
    kernels={
        "Segments_crossing_drift": xo.Kernel(
            c_name="Segments_crossing_drift",
            args=[
                xo.Arg(Segments,   name="segs"),
                xo.Arg(xo.Int8,    pointer=True, name="n_hit"),
                xo.Arg(xo.Float64, pointer=True,  name="s"),
                xo.Arg(xo.Float64, pointer=False, name="s0"),
                xo.Arg(xo.Float64, pointer=False, name="x0"),
                xo.Arg(xo.Float64, pointer=False, name="m")
            ],
            ret=None,
        )
    }, save_source_as="temp.c"
)


# Segment(s) to be tested #############################
line_segment      = LineSegment(s1=0,x1=0,s2=0.2,x2=1)
half_open_segment = HalfOpenLineSegment(s=2, x=1, t=np.pi/4)
bezier_segment    = BezierSegment(s1=0.2, x1=1, s2=1+0.5*np.cos(3*np.pi/4), x2=0.5+0.5*np.sin(3*np.pi/4), cs1=0.5, cx1=2.5, cs2=1+0.5*np.cos(3*np.pi/4)-0.6, cx2=0.5+0.5*np.sin(3*np.pi/4)-0.6)

half_circle       = CircularSegment(R=0.5, s=1, x=0.5, t1=-np.pi/2, t2=np.pi/2)
half_circle_2     = CircularSegment(R=0.5, s=1, x=0.5, t1=np.pi/2, t2=3*np.pi/2)
half_circle_3     = CircularSegment(R=0.5, s=1, x=0.5, t1=0, t2=np.pi)
half_circle_4     = CircularSegment(R=0.5, s=1, x=0.5, t1=np.pi, t2=2*np.pi)

seg        = Segments(data=[LineSegment(s1=0, x1=0, s2=0.2, x2=1),
                            BezierSegment(s1=0.2, x1=1, s2=1+0.5*np.cos(3*np.pi/4), x2=0.5+0.5*np.sin(3*np.pi/4), cs1=0.5, cx1=2.5, cs2=1+0.5*np.cos(3*np.pi/4)-0.6, cx2=0.5+0.5*np.sin(3*np.pi/4)-0.6),
                            CircularSegment(R=0.5, s=1, x=0.5, t1=-np.pi/2, t2=3*np.pi/4),
                            LineSegment(s1=1, x1=0, s2=0, x2=0)])

seg_bezier = Segments(data=[BezierSegment(s1=0.2, x1=1, s2=1+0.5*np.cos(3*np.pi/4), x2=0.5+0.5*np.sin(3*np.pi/4), cs1=0.5, cx1=2.5, cs2=1+0.5*np.cos(3*np.pi/4)-0.6, cx2=0.5+0.5*np.sin(3*np.pi/4)-0.6),
                            BezierSegment(s1=1+0.5*np.cos(3*np.pi/4), x1=0.5+0.5*np.sin(3*np.pi/4), s2=1.2, x2=0.5, cs1=1.4, cx1=1.6, cs2=1.1, cx2=0.9),
                            BezierSegment(s1=1.2, x1=0.5, s2=0.5, x2=0.2, cs1=1.2, cx1=0.1, cs2=0.4, cx2=0.1),
                            BezierSegment(s1=0.5, x1=0.2, s2=0.2, x2=1, cs1=0.7, cx1=0.5, cs2=0.1, cx2=0.2)])

# Test functions #####################################
def test_line_segment(): 
    _check_ref('line',line_segment)

def test_half_open_segment():
    _check_ref('half_open', half_open_segment)

def test_bezier_segment():
    _check_ref('bezier', bezier_segment)

def test_half_circle():
    _check_ref('half_circle', half_circle)
    _check_ref('half_circle_2', half_circle_2)   
    _check_ref('half_circle_3', half_circle_3)
    _check_ref('half_circle_4', half_circle_4)

def test_segments(): 
    _check_ref('segments', seg)

def test_segments_bezier():
    _check_ref('bezier_segments', seg_bezier)

def _check_ref(name, seg, _context=None):
    print(f"Testing {name}")
    # if _context is None:
    #     _context = xo.ContextCpu()
    with open(Path(path, 'ref_'+name+'.json'), 'r') as fid:
        ref = json.load(fid)
    print('yes')
    ref_s     = np.array(ref['hit_s'])
    ref_x     = np.array(ref['hit_x'])
    ref_m     = np.array(ref['m'])
    ref_point = np.array(ref['point'])
    print('yes')
    n_hit = np.zeros(1, dtype=np.int8)
    s = np.zeros(get_max_crossings(seg, "drift"), dtype=np.float64)
    print('yes')
    s0, x0 = ref_point
    n_hit[0] = 0
    hits_s = np.array([], dtype=np.float64)
    hits_x = np.array([], dtype=np.float64)
    context.kernels.Segments_crossing_drift(segs=seg, n_hit=n_hit, s=s, s0=s0, x0=x0, m=ref_m[0])
    print('yes')
    hits_s = np.concatenate((hits_s, s[:n_hit[0]]))
    hits_x = np.concatenate((hits_x, x0 + ref_m[0] * (s[:n_hit[0]] - s0)))
    print('yes')
    assert np.allclose(hits_s, ref_s, atol=1e-12, rtol=0), "s does not match."
    assert np.allclose(hits_x, ref_x, atol=1e-12, rtol=0), "x does not match"
    print(f"Test for {name} passed.")

test_line_segment()
test_half_open_segment()
test_bezier_segment()
test_half_circle()
test_segments()
test_segments_bezier()