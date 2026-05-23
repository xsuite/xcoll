# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

import numpy as np
from pathlib import Path
import time
start_time = time.time()
import matplotlib.pyplot as plt

import xobjects as xo
import xtrack as xt
import xcoll as xc


beam          = 1
plane         = 'H'
num_turns     = 10
num_particles = 5000

path_in = Path(__file__).parent
path_out = Path.cwd()


# Load from json
env = xt.load(path_in / 'machines' / f'lhc_run3_b{beam}.json')
line = env[f'lhcb{beam}']


# Initialise colldb
colldb = xc.CollimatorDatabase.from_yaml(path_in / 'colldbs' / f'lhc_run3.yaml', beam=beam)


# Install collimators into line
colldb.install_geant4_collimators(line=line, verbose=True)


# Aperture model check
print('\nAperture model check after introducing collimators:')
df_with_coll = line.check_aperture()
assert not np.any(df_with_coll.has_aperture_problem)


# Assign the optics to deduce the gap settings
line.collimators.assign_optics()


# Connect to Geant4
xc.geant4.engine.start(line=line, cwd='run_geant4_temp', clean=True, verbose=True)


# Generate initial pencil distribution on horizontal collimator
tcp  = f"tcp.{'c' if plane=='H' else 'd'}6{'l' if f'{beam}'=='1' else 'r'}7.b{beam}"
part = line[tcp].generate_pencil(num_particles)


# # Treat warnings as errors to debug
# np.set_printoptions(threshold=np.inf)
# import warnings
# warnings.filterwarnings("error")


# Track!
line.scattering.enable()
line.track(part, num_turns=num_turns, time=True, with_progress=1)
line.scattering.disable()
print(f"Done tracking in {line.time_last_track:.1f}s.")


# Save loss map to json
lm_time = time.time()
line_is_reversed = True if f'{beam}' == '2' else False
ThisLM = xc.LossMap(line, line_is_reversed=line_is_reversed, part=part)
print(f"Loss map created in {time.time()-lm_time:.1f}s.")
ThisLM.to_json(file=path_out / 'results' / f'lossmap_geant4_B{beam}{plane}.json')

# Save a summary of the collimator losses to a text file
ThisLM.save_summary(file=path_out / 'results' / f'coll_summary_geant4_B{beam}{plane}.out')
print(ThisLM.summary)


# Stop the Geant4 connection (this should be after creating the loss map to get the correct particle_ref)
xc.geant4.engine.stop(clean=True)


print(f"Total calculation time {time.time()-start_time}s")

# Plot loss map
ThisLM.plot(savefig=path_out / 'plots' / 'lossmaps' / f'lossmap_geant4_B{beam}{plane}.pdf', zoom='betatron')
plt.show()

