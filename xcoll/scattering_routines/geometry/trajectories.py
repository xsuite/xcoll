# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

import xobjects as xo


# Trajectory definitions
# ----------------------

trajectories = {
    "drift": {
        "args": [
            # The arguments that define the horizontal (after rotation) particle trajectory
            xo.Arg(xo.Float64, pointer=False, name="s0"),  # Particle s
            xo.Arg(xo.Float64, pointer=False, name="x0"),  # Particle x
            xo.Arg(xo.Float64, pointer=False, name="xm")   # Particle slope in the x direction (= xp = tan(theta_x))
        ],
        "args_vlimit_common": [
            # The arguments that define the vertical (after rotation) particle trajectory, shared with the above
            xo.Arg(xo.Float64, pointer=False, name="s0")   # Particle s
        ],
        "args_vlimit_extra": [
            # The arguments that define the vertical (after rotation) particle trajectory, new
            xo.Arg(xo.Float64, pointer=False, name="y0"),  # Particle y
            xo.Arg(xo.Float64, pointer=False, name="ym")   # Particle slope in the y direction (= yp = tan(theta_y))
        ],
        "max_crossings": {
            "line": 2,          # 2 crossings in case of parallel trajectory
            "halfopenline": 2,
            "circular": 2,
            "bezier": 3
        }
    ## Add more trajectories here ##
    }
}


trajectories_vlimit_sources = {
    'drift': """
/*gpufun*/
int8_t vlimit_drift(double* restrict_s, double s0, double y0, double ym, double ymin, double ymax){
    if (fabs(ym) < XC_EPSILON){
        // Trajectory parallel to s axis
        if (y0 < ymin || y0 > ymax){
            return 0;  // Completely outside - no crossing possible
        } else {
            return -1; // Completely inside - no vertical check needed
        }
    } else {
        restrict_s[0] = (ymin - y0)/ym + s0;
        restrict_s[1] = (ymax - y0)/ym + s0;
        SWAP(restrict_s, 0, 1);   // To make sure these are sorted
        return 1;  // Default behavior: check overlap with horizontal crossings
    }
}
"""
    ## Add vlimit functions for each trajectory here ##
}


# Function to get the maximum number of crossings for a given object type
def get_max_crossings(segments, trajectory):
    from xcoll.scattering_routines.geometry.segments import all_segments
    if hasattr(segments, '__iter__') and all(isinstance(seg, all_segments) for seg in segments):
        max_crossings = 0
        for seg in segments:
            max_crossings += trajectories[trajectory]["max_crossings"][seg.__class__.__name__.lower()[:-7]]
        return max_crossings
    elif isinstance(segments, all_segments):
        return trajectories[trajectory]["max_crossings"][segments.__class__.__name__.lower()[:-7]]
    else:
        raise ValueError(f"Unexpected type for segments: {type(segments)}")
