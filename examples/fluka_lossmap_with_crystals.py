# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

import numpy as np
from pathlib import Path
import time
start_time = time.time()
import matplotlib.pyplot as plt

import xobjects as xo
import xpart as xp
import xtrack as xt
import xcoll as xc


beam              = 1
plane             = 'H'
num_turns         = 10
num_particles     = 50
capacity          = 500*num_particles
relative_capacity = 500
particle_ref      = xt.Particles('Pb208', p0c=6.8e12*82)

path_in = Path(__file__).parent
path_out = Path.cwd()


# Load from json
env = xt.load(path_in / 'machines' / f'lhc_run3_b{beam}.json')
line = env[f'lhcb{beam}']


# Initialise colldb
colldb = xc.CollimatorDatabase.from_yaml(path_in / 'colldbs' / f'lhc_run3_crystals.yaml',
                                               beam=beam, ignore_crystals=False)


# Install collimators into line
colldb.install_fluka_collimators(line=line, verbose=True)

# Aperture model check
print('\nAperture model check after introducing collimators:')
df_with_coll = line.check_aperture()
assert not np.any(df_with_coll.has_aperture_problem)


# Primary crystal
tcpc = f"tcpc{plane.lower()}.a{6 if plane=='V' else 4 if f'{beam}'=='1' else 5}{'l' if f'{beam}'=='1' else 'r'}7.b{beam}"


# Assign the optics to deduce the gap settings
line.collimators.assign_optics()
line.collimators.align_to_beam_divergence()


# Connect to FLUKA
xc.fluka.engine.particle_ref = particle_ref
xc.fluka.engine.capacity = capacity
xc.fluka.engine.relative_capacity = relative_capacity
xc.fluka.engine.seed = 5656565
xc.fluka.engine.start(line=line, capacity=xc.fluka.engine.capacity, cwd='run_fluka_temp', clean=False, verbose=True, return_ions=True)


# # Generate initial pencil distribution on crystal
# part = line[tcp].generate_pencil(num_particles)
# Generate initial halo
x_norm, px_norm, _, _ = xp.generate_2D_uniform_circular_sector(r_range=(5, 5.04), num_particles=num_particles)
y_norm  = np.random.normal(scale=0.01, size=num_particles)
py_norm = np.random.normal(scale=0.01, size=num_particles)
part = line.build_particles(
            x_norm=x_norm, px_norm=px_norm, y_norm=y_norm, py_norm=py_norm,
            nemitt_x=line[tcpc].nemitt_x, nemitt_y=line[tcpc].nemitt_y,
            at_element=tcpc, particle_ref=xc.fluka.engine.particle_ref,
            _capacity=xc.fluka.engine.capacity)


# Move the line to an OpenMP context to be able to use all cores
line.discard_tracker()
line.build_tracker(_context=xo.ContextCpu(omp_num_threads='auto'))
# Should move iobuffer as well in case of impacts


# Track!
line.scattering.enable()
line.track(part, num_turns=num_turns, time=True, with_progress=1)
line.scattering.disable()
print(f"Done tracking in {line.time_last_track:.1f}s.")


# Move the line back to the default context to be able to use all prebuilt kernels for the aperture interpolation
line.discard_tracker()
line.build_tracker(_context=xo.ContextCpu())

# Save loss map to json
lm_time = time.time()
line_is_reversed = True if f'{beam}' == '2' else False
ThisLM = xc.LossMap(line, line_is_reversed=line_is_reversed, part=part)
print(f"Loss map created in {time.time()-lm_time:.1f}s.")
ThisLM.to_json(file=path_out / 'results' / f'lossmap_fluka_crystals_B{beam}{plane}.json')

# Save a summary of the collimator losses to a text file
ThisLM.save_summary(file=path_out / 'results' / f'coll_summary_fluka_crystals_B{beam}{plane}.out')
print(ThisLM.summary)


# Stop the FLUKA connection  (this should be after creating the loss map to get the correct particle_ref)
xc.fluka.engine.stop(clean=False)


print(f"Total calculation time {time.time()-start_time}s")

# Plot loss map
ThisLM.plot(savefig=path_out / 'plots' / 'lossmaps' / f'lossmap_fluka_crystals_B{beam}{plane}.pdf', zoom='betatron')
plt.show()

