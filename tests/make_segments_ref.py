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
path = Path(__file__).parent / 'data_test_segments'

def create_ref(segs, name, point, part_theta):
    n_hit = np.zeros(1, dtype=np.int8)
    m = np.tan(np.radians(part_theta))
    s0, x0 = point
    s = np.zeros(get_max_crossings(segs, "drift"), dtype=np.float64)

    n_hit[0] = 0
    hits_s = []
    hits_x = []
    n_hit[0] = 0
    context.kernels.Segments_crossing_drift(segs=segs, n_hit=n_hit, s=s, s0=s0, x0=x0, m=m)
    hits_s += list(s[:n_hit[0]])
    hits_x += list(x0 + m*(s[:n_hit[0]] - s0))


    data = {
        "point": point,
        "m": [m],
        "hit_s": hits_s,
        "hit_x": hits_x,
    } 
    with open(Path(path,f"ref_{name}.json"), "w") as file:
        json.dump(data, file)
    print(f"Created reference for {name}.")

# single segments to create references for 
line_segment      = LineSegment(s1=0,x1=0,s2=0.2,x2=1)
half_open_segment = HalfOpenLineSegment(s=0.8, x=1, t=np.pi/2)
bezier_segment    = BezierSegment(s1=0.2, x1=1, s2=1+0.5*np.cos(3*np.pi/4), x2=0.5+0.5*np.sin(3*np.pi/4), cs1=0.5, cx1=2.5, cs2=1+0.5*np.cos(3*np.pi/4)-0.6, cx2=0.5+0.5*np.sin(3*np.pi/4)-0.6)
half_circle       = CircularSegment(R=0.5, s=1, x=0.5, t1=-np.pi/2, t2=np.pi/2)
half_circle_2     = CircularSegment(R=0.5, s=1, x=0.5, t1=np.pi/2, t2=3*np.pi/2)
half_circle_3     = CircularSegment(R=0.5, s=1, x=0.5, t1=0, t2=np.pi)
half_circle_4     = CircularSegment(R=0.5, s=1, x=0.5, t1=np.pi, t2=2*np.pi)

create_ref(Segments(data=[line_segment]), name="line", point=[0.67,0.53], part_theta=30)
create_ref(Segments(data=[half_open_segment]), name="half_open", point=[1.5,0.5], part_theta=120)
create_ref(Segments(data=[bezier_segment]), name="bezier", point=[0.5,0.8], part_theta=120)
create_ref(Segments(data=[half_circle]), name="half_circle", point=[1.131,0.151], part_theta=77)
create_ref(Segments(data=[half_circle_2]), name="half_circle_2", point=[1.4,0.67], part_theta=43)
create_ref(Segments(data=[half_circle_3]), name="half_circle_3", point=[1.0,0.8], part_theta=120)
create_ref(Segments(data=[half_circle_4]), name="half_circle_4", point=[0.75,0.8], part_theta=170)

# Create reference for Segments
seg        = Segments(data=[LineSegment(s1=0, x1=0, s2=0.2, x2=1),
                            BezierSegment(s1=0.2, x1=1, s2=1+0.5*np.cos(3*np.pi/4), x2=0.5+0.5*np.sin(3*np.pi/4), cs1=0.5, cx1=2.5, cs2=1+0.5*np.cos(3*np.pi/4)-0.6, cx2=0.5+0.5*np.sin(3*np.pi/4)-0.6),
                            CircularSegment(R=0.5, s=1, x=0.5, t1=-np.pi/2, t2=3*np.pi/4),
                            LineSegment(s1=1, x1=0, s2=0, x2=0)])

seg_bezier = Segments(data=[BezierSegment(s1=0.2, x1=1, s2=1+0.5*np.cos(3*np.pi/4), x2=0.5+0.5*np.sin(3*np.pi/4), cs1=0.5, cx1=2.5, cs2=1+0.5*np.cos(3*np.pi/4)-0.6, cx2=0.5+0.5*np.sin(3*np.pi/4)-0.6),
                            BezierSegment(s1=1+0.5*np.cos(3*np.pi/4), x1=0.5+0.5*np.sin(3*np.pi/4), s2=1.2, x2=0.5, cs1=1.4, cx1=1.6, cs2=1.1, cx2=0.9),
                            BezierSegment(s1=1.2, x1=0.5, s2=0.5, x2=0.2, cs1=1.2, cx1=0.1, cs2=0.4, cx2=0.1),
                            BezierSegment(s1=0.5, x1=0.2, s2=0.2, x2=1, cs1=0.7, cx1=0.5, cs2=0.1, cx2=0.2)])

create_ref(seg, name="segments", point=[0.249,0.934], part_theta=160)
create_ref(seg_bezier, name="bezier_segments", point=[0.931,0.83], part_theta=7)
