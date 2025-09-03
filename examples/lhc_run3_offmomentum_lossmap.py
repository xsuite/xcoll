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
import xtrack as xt
import xpart as xp
import xobjects as xo
import xcoll as xc


# We do the majority of the script on the default context to be able to use prebuilt kernels
context = xo.ContextCpu()


# This script takes around 8 minutes on a modern CPU (80s preparation+interpolation, 400s tracking)
beam = 1
plane = 'DPpos'

num_particles  = 500
sweep          = 300
sweep          = -abs(sweep) if plane == 'DPpos' else abs(sweep)
num_turns      = int(20*abs(sweep))

path_in  = xc._pkg_root.parent / 'examples'
path_out = Path.cwd()


# Load from json
line = xt.Line.from_json(path_in / 'machines' / f'lhc_run3_b{beam}.json')


# Initialise colldb
colldb = xc.CollimatorDatabase.from_yaml(path_in / 'colldb' / f'lhc_run3.yaml', beam=beam)


# Install collimators into line
colldb.install_everest_collimators(line=line, verbose=True)


# Aperture model check
print('\nAperture model check after introducing collimators:')
df_with_coll = line.check_aperture()
assert not np.any(df_with_coll.has_aperture_problem)


# Build the tracker
line.build_tracker()


# Assign the optics to deduce the gap settings
line.collimators.assign_optics()


# Generate initial matched bunch
part = xp.generate_matched_gaussian_bunch(nemitt_x=colldb.nemitt_x,
                                          nemitt_y=colldb.nemitt_y,
                                          sigma_z=7.55e-2, num_particles=num_particles, line=line)


# Move the line to an OpenMP context to be able to use all cores
line.discard_tracker()
line.build_tracker(_context=xo.ContextCpu(omp_num_threads='auto'))
# Should move iobuffer as well in case of impacts


# Print some info of the RF sweep
rf_sweep = xc.RFSweep(line)
rf_sweep.prepare(sweep_per_turn=sweep/num_turns)
rf_sweep.info()


# Track during RF sweep:
line.scattering.enable()
line.track(particles=part, num_turns=num_turns, time=True, with_progress=5)
line.scattering.disable()
print(f"Done sweeping RF in {line.time_last_track:.1f}s.")


# Move the line back to the default context to be able to use all prebuilt kernels for the aperture interpolation
line.discard_tracker()
line.build_tracker(_context=xo.ContextCpu())


# Let's visualise how the losses move from IR7 to IR3 during the sweep
tcp3 = line.element_names.index('tcp.6l3.b1')
turns_ip3 = part.at_turn[part.at_element == tcp3]
tcp7v = line.element_names.index('tcp.d6l7.b1')
turns_ip7v = part.at_turn[part.at_element == tcp7v]
tcp7h = line.element_names.index('tcp.c6l7.b1')
turns_ip7h = part.at_turn[part.at_element == tcp7h]
tcp7s = line.element_names.index('tcp.b6l7.b1')
turns_ip7s = part.at_turn[part.at_element == tcp7s]

plt.figure(figsize=(8,5))
plt.hist([turns_ip7v, turns_ip7h, turns_ip7s, turns_ip3], 100, stacked=True,
         label=['Lost on TCP.D7', 'Lost on TCP.C7', 'Lost on TCP.B7', 'Lost on TCP.3'])
plt.xlabel('turn')
plt.ylabel('number of lost particles (stacked)')
plt.legend()
plt.tight_layout()
plt.savefig('OffMomentumPosition.png', dpi=300)
plt.show()


# In this example, the IR3 region gets populated very fast. When this is not the case, e.g. for injection,
# one should keep only the last part of the sweep (as this is what happens in reality), as follows:


# In reality the loss map only shows the last second (1.3s BLM integration time)
# A typical RF sweep in reality shifts around 25Hz per second, i.e. around 2.2mHz/turn.
# We sweep 50mHz/turn, so we are ~22.5 times faster than in reality.
# We keep the equivalent of the last 3 seconds (slightly higher integration time, in order not to lose too much statistics).
# We want to keep the last 3/22.5 seconds = 1500 turns.
# This is equal to:
seconds_to_keep = 3
turns_to_keep = seconds_to_keep*25*num_turns/abs(sweep)
last_turn = part.at_turn.max()
part2 = part.filter(part.at_turn > last_turn - turns_to_keep)
print(f"Keeping last {int(turns_to_keep)} turns (equivalent integration time of {seconds_to_keep} seconds).")
print(f"This means we use {len(part2.x)} particles (of which {len(part2.x[part2.state > 0])} survived).")


# Save lossmap to json, which can be loaded, combined (for more statistics),
# and plotted with the 'lossmaps' package
line_is_reversed = True if f'{beam}' == '2' else False
ThisLM = xc.LossMap(line, line_is_reversed=line_is_reversed, part=part)
ThisLM.to_json(file=Path(path_out, f'lossmap_B{beam}{plane}.json'))

# Save a summary of the collimator losses to a text file
ThisLM.save_summary(file=Path(path_out, f'coll_summary_B{beam}{plane}.out'))
print(ThisLM.summary)

print(f"Total calculation time {time.time()-start_time}s")

ThisLM.plot()
plt.show()
