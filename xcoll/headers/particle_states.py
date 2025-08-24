# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

from ..general import _pkg_root


with open(_pkg_root / 'headers' / 'particle_states.h', 'r') as fid:
    particle_states_source = fid.read()

particle_states = {
    line.split()[1][3:]: int(line.split()[2])
    for line in particle_states_source.split('\n')
    if len(line.split()) > 1 and line.split()[1][:3] == 'XC_' # select the source lines with the definitions
}

def __getattr__(name):
    if name in particle_states:
        return particle_states[name]
    else:
        raise AttributeError(f"module '{__name__}' has no attribute '{name}'")
