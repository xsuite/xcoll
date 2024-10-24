# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

import xobjects as xo

from ..c_init import xo_to_ctypes
from .drift import DriftTrajectory


all_trajectories = [DriftTrajectory]


args_cross_h = [
    xo.Arg(xo.Int8,    pointer=True, name="n_hit"),
    xo.Arg(xo.Float64, pointer=True, name="s"),
]
args_cross_v = [
    xo.Arg(xo.Float64, pointer=True, name="restrict_s"),
]
args_vlimit = [
    xo.Arg(xo.Float64, pointer=False, name="ymin"),
    xo.Arg(xo.Float64, pointer=False, name="ymax")
]
args_length = [
    xo.Arg(xo.Float64, pointer=False, name="s1"),
    xo.Arg(xo.Float64, pointer=False, name="s2")
]


# Sanity check to assert all trajectory types have a vlimit function and a length function
def assert_trajectory_sources(seg):
    for tra in all_trajectories:
        header = f"/*gpufun*/\nint8_t {tra.__name__}_vlimit({xo_to_ctypes(args_cross_v)}, {xo_to_ctypes(tra.args_hv)}, " \
               + f"{xo_to_ctypes(tra.args_v)}, {xo_to_ctypes(args_vlimit)})"
        header_found = False
        for src in tra._extra_c_sources:
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
            raise SystemError(f"Missing or corrupt C vlimit function for {tra.__name__}.")

        header = f"/*gpufun*/\ndouble {tra.__name__}_length({xo_to_ctypes(tra.args_hv)}, {xo_to_ctypes(tra.args_h)}, " \
                 f"{xo_to_ctypes(tra.args_v)}, {xo_to_ctypes(args_length)})"
        header_found = False
        for src in tra._extra_c_sources:
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
            raise SystemError(f"Missing or corrupt C length function for {tra.__name__}.")


for tra in all_trajectories:
    assert_trajectory_sources(tra)
    assert hasattr(tra, 'args_hv')
    assert hasattr(tra, 'args_h')
    assert hasattr(tra, 'args_v')
    tra.name = tra.__name__.lower()[:-10]
