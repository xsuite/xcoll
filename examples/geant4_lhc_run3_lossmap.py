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
path_out = Path.cwd() / 'plots'


# Load from json
line = xt.Line.from_json(path_in / 'machines' / f'lhc_run3_b{beam}.json')


# Initialise colldb
colldb = xc.CollimatorDatabase.from_yaml(path_in / 'colldb' / f'lhc_run3.yaml', beam=beam)


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


# Stop the Geant4 connection (and return to the previous directory)
xc.geant4.engine.stop(clean=True)


# Save lossmap to json, which can be loaded, combined (for more statistics),
# and plotted with the 'lossmaps' package
line_is_reversed = True if f'{beam}' == '2' else False
ThisLM = xc.LossMap(line, line_is_reversed=line_is_reversed, part=part)
ThisLM.to_json(file=Path(path_out, f'lossmap_B{beam}{plane}.json'))

# Save a summary of the collimator losses to a text file
ThisLM.save_summary(file=Path(path_out, f'coll_summary_B{beam}{plane}.out'))
print(ThisLM.summary)

print(f"Total calculation time {time.time()-start_time}s")

ThisLM.plot(savefig=Path(path_out, f'lossmap_geant4_B{beam}{plane}.pdf'))
plt.show()

