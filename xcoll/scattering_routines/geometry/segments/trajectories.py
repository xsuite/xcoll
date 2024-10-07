# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

import xobjects as xo


trajectories = {
    "drift": {
        "crossing_args": [
            xo.Arg(xo.Float64, pointer=False, name="s0"),
            xo.Arg(xo.Float64, pointer=False, name="x0"),
            xo.Arg(xo.Float64, pointer=False, name="m")
        ],
        "max_crossings": {
            "line": 2,          # 2 crossings in case of parallel trajectory
            "halfopenline": 2,
            "circular": 2,
            "bezier": 3
        }
    }
}


trajectories_c_args = {trajectory: [", ".join([f"{arg.get_c_type()} {arg.name}" for arg in args['crossing_args']]),
                                    ", ".join([f"{arg.name}" for arg in args['crossing_args']])]
                       for trajectory, args in trajectories.items()}

def get_max_crossings(segments, trajectory):
    from xcoll.scattering_routines.geometry.segments import all_segments, Segments
    if isinstance(segments, Segments):
        max_crossings = 0
        for seg in segments:
            max_crossings += trajectories[trajectory]["max_crossings"][seg.__class__.__name__.lower()[:-7]]
        return max_crossings
    elif isinstance(segments, all_segments):
        return trajectories[trajectory]["max_crossings"][segments.__class__.__name__.lower()[:-7]]
    else:
        raise ValueError(f"Unexpected type for segments: {type(segments)}")

