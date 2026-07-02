# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2026.                 #
# ######################################### #

import time
import numpy as np
import pandas as pd
from pathlib import Path

import xtrack as xt
import xcoll as xc
import xpart as xp


beam      = 1
plane     = "H"
output    = "output"
num_bunch = 40  # One of the main bunches contributing to losses.

num_particles = 5000
momentum      = 6.8e12

start_time = time.time()

path_in = Path(__file__).parent
path_out = Path.cwd()


def async_beam_dump_per_turn(sync=False, oprefire=False, aprefire=False,
                             bunch=0, turn=0):
    phase_data = (df_angles['phase_data'].values)[~np.isnan(df_angles['phase_data'])]
    kicker_magnets = {'b1': [f"mkd.{i}5l6.b1" for i in "abcdefghijklmno"],
                      'b2': [f"mkd.{i}5r6.b2" for i in "abcdefghijklmno"]}
    sign = -1 if beam == 1 else 1
    if sync:
        if turn == 0:
            for mlp in kicker_magnets[f'b{beam}']:
                line[mlp].knl[0] = sign*phase_data[bunch]*1e-3
        else:
            for mlp in kicker_magnets[f'b{beam}']:
                line[mlp].knl[0] = sign*phase_data[-1]*1e-3
    elif aprefire:
        if turn == 0:
            for index, mlp in enumerate(kicker_magnets[f'b{beam}']):
                line[mlp].knl[0] = sign*df_angles[f'a_data{index+1}'][bunch]*1e-3
        else:
            for index, mlp in enumerate(kicker_magnets[f'b{beam}']):
                line[mlp].knl[0] = sign*df_angles[f'a_data{index+1}'].values[-1]*1e-3
    elif oprefire:
        n_mkd = len(kicker_magnets[f'b{beam}'])
        if turn == 0:
            for index, mlp in enumerate(kicker_magnets[f'b{beam}']):
                line[mlp].knl[0] = sign*df_angles[f'a_data{n_mkd-index}'][bunch]*1e-3
        else:
            for index, mlp in enumerate(kicker_magnets[f'b{beam}']):
                line[mlp].knl[0] = sign*df_angles[f'a_data{n_mkd-index}'].values[-1]*1e-3


# Main script
# ===========

# Load from json
env = xt.load(path_in / 'machines' / f'lhc_run3_b{beam}.json')
line = env[f'lhcb{beam}']


# Initialise colldb
colldb = xc.CollimatorDatabase.from_yaml(path_in / 'colldbs' / f'lhc_run3.yaml', beam=beam)


# Open angles file
df_angles = pd.read_csv(path_in / 'extras' / f'data_asd_phases.csv')


# Install collimators into line
colldb.install_fluka_collimators(line=line, verbose=True)


# Aperture model check
print('\nAperture model check after introducing collimators:')
df_with_coll = line.check_aperture()
print(df_with_coll)
assert not np.any(df_with_coll.has_aperture_problem)


# Prepare engine
line.particle_ref = xt.Particles('proton', p0c=momentum)
xc.fluka.engine.particles_ref = line.particle_ref
xc.fluka.engine.capacity = 5*num_particles
xc.fluka.engine.relative_capacity = 2


# Build the tracker
line.build_tracker()


# Assign the optics to deduce the gap settings
line.xcoll.collimators.assign_optics()
xc.fluka.engine.start(line=line, cwd='run_fluka_temp', clean=False, verbose=True, touches=True)


# Create initial particles
x_norm, px_norm = xp.generate_2D_gaussian(num_particles=num_particles)
y_norm, py_norm = xp.generate_2D_gaussian(num_particles=num_particles)
sigma_z = 2.5e-9
tw =  line.twiss()
zeta_co = tw['zeta', 'ip4']
delta_co = tw['delta','ip4']
zeta, delta = xp.generate_longitudinal_coordinates(
    line=line,
    num_particles=num_particles,
    distribution='gaussian',
    sigma_z=sigma_z,
    particle_ref=line.particle_ref)
part = line.build_particles(x_norm=x_norm, y_norm=y_norm,
                            _capacity=3*num_particles,
                            px_norm=px_norm,
                            py_norm=py_norm,
                            zeta=zeta + zeta_co,
                            delta = delta + delta_co,
                            nemitt_x = 2.5e-6,
                            nemitt_y = 2.5e-6,
                            at_element='ip4')
part.start_tracking_at_element = -1


# Track
# =====

line.xcoll.scattering.enable()

sync=False
oprefire=True
aprefire=False

# First turn
async_beam_dump_per_turn(sync=sync, oprefire=oprefire, aprefire=aprefire, bunch=num_bunch, turn=0)
line.track(part, num_turns=1, time=True, with_progress=1, ele_start='ip4', ele_stop='ip4')

# Second turn
async_beam_dump_per_turn(sync=sync, oprefire=oprefire, aprefire=aprefire, bunch=num_bunch, turn=1)
line.track(part, num_turns=1, time=True, with_progress=1, ele_start='ip4', ele_stop='ip4')

line.xcoll.scattering.disable()

# Stop the FLUKA connection (and return to the previous directory)
xc.fluka.engine.stop(clean=True)

line_is_reversed = True if f'{beam}' == '2' else False
ThisLM = xc.LossMap(line, line_is_reversed=line_is_reversed, part=part)
ThisLM.to_json(file=path_out / 'results' / f'lossmap_fluka_ASD_B{beam}{plane}.json')

# Save a summary of the collimator losses to a text file
ThisLM.save_summary(file=path_out / 'results' / f'coll_summary_fluka_ASD_B{beam}{plane}.out')
print(ThisLM.summary)

print(f"Total calculation time {time.time()-start_time}s")

# Plot loss map
ThisLM.plot(savefig=path_out / 'plots' / 'lossmaps' / f'lossmap_fluka_ASD_B{beam}{plane}.pdf', zoom='betatron')
