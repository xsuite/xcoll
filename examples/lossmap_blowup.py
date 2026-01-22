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
import xcoll as xc


beam          = 1
plane         = 'V'
num_turns     = 1500
num_particles = 10_000

path_in = Path(__file__).parent
path_out = Path.cwd()


# Load from json
env = xt.load(path_in / 'machines' / f'lhc_run3_b{beam}.json')
line = env[f'lhcb{beam}']


# Initialise colldb
colldb = xc.CollimatorDatabase.from_yaml(path_in / 'colldbs' / f'lhc_run3.yaml', beam=beam)


# Install collimators into line
colldb.install_everest_collimators(line=line, verbose=True)


# Install ADT into line
# ADT kickers in LHC are named adtk[hv].[abcd]5[lr]4.b1 (with the position 5l4 (B1H or B2V) or 5r4 (B1V or B2H)
# These are not in the line, but their tank names are: adtk[hv].[abcd]5[lr]4.[abcd].b1  (32 markers)
pos = 'b5l4' if f'{beam}' == '1' and plane == 'H' else 'b5r4'
pos = 'b5l4' if f'{beam}' == '2' and plane == 'V' else pos
name = f'adtk{plane.lower()}.{pos}.b{beam}'
tank_start = f'adtk{plane.lower()}.{pos}.a.b{beam}'
tank_end   = f'adtk{plane.lower()}.{pos}.d.b{beam}'
adt_pos = 0.5*line.get_s_position(tank_start) + 0.5*line.get_s_position(tank_end)
adt = xc.BlowUp.install(line, name=f'{name}_blowup', at_s=adt_pos, plane=plane, stop_at_turn=num_turns,
                        amplitude=0.5, use_individual_kicks=True)


# Aperture model check
print('\nAperture model check after introducing elements:')
df_with_coll = line.check_aperture()
assert not np.any(df_with_coll.has_aperture_problem)


# Assign the optics to deduce the gap settings, and calibrate the ADT
tw = line.twiss(reverse=False)
line.collimators.assign_optics(twiss=tw)
if plane == 'H':
    adt.calibrate_by_emittance(nemitt=colldb.nemitt_x, twiss=tw)
else:
    adt.calibrate_by_emittance(nemitt=colldb.nemitt_y, twiss=tw)


# Optimise the line
line.optimize_for_tracking()


# Generate a matched Gaussian bunch
part = xp.generate_matched_gaussian_bunch(num_particles=num_particles, total_intensity_particles=1.6e11,
                                          nemitt_x=colldb.nemitt_x, nemitt_y=colldb.nemitt_y, sigma_z=7.55e-2, line=line)


# Move the line to an OpenMP context to be able to use all cores
line.discard_tracker()
line.build_tracker(_context=xo.ContextCpu(omp_num_threads='auto'))
# Should move iobuffer as well in case of impacts


# Track!
line.scattering.enable()
adt.activate()
line.track(part, num_turns=num_turns, time=True, with_progress=1)
adt.deactivate()
line.scattering.disable()
print(f"Done tracking in {line.time_last_track:.1f}s.")


# Move the line back to the default context to be able to use all prebuilt kernels for the aperture interpolation
line.discard_tracker()
line.build_tracker(_context=xo.ContextCpu())


# Save loss map to json
line_is_reversed = True if f'{beam}' == '2' else False
ThisLM = xc.LossMap(line, line_is_reversed=line_is_reversed, part=part)
ThisLM.to_json(file=path_out / 'results' / f'lossmap_blowup_B{beam}{plane}.json')

# Save a summary of the collimator losses to a text file
ThisLM.save_summary(file=path_out / 'results' / f'coll_summary_blowup_B{beam}{plane}.out')
print(ThisLM.summary)

print(f"Total calculation time {time.time()-start_time}s")

# Plot loss map
ThisLM.plot(savefig=path_out / 'plots' / 'lossmaps' / f'lossmap_blowup_B{beam}{plane}.pdf', zoom='betatron')
plt.show()


# Plot
collimators_with_most_losses = ThisLM.summary.sort_values('n', ascending=False).name.values[:5]

bin_centers = {}
losses = {}
for coll in collimators_with_most_losses:
    idx = line.element_names.index(coll)
    mask = part.at_element == idx
    n, x, _ = plt.hist(part.at_turn[mask], 500, cumulative=False)
    bin_centers[coll] = 0.5*(x[1:]+x[:-1])
    losses[coll] = n
    plt.close()

_, ax = plt.subplots(figsize=(10,5))
for coll in collimators_with_most_losses:
    ax.plot(bin_centers[coll], losses[coll], label=coll)
ax.set_ylabel('number of lost particles')
ax.set_xlabel("turn")
ax.set_yscale('log')
ax.legend()
ax.set_title(f"Loss map B{beam}{plane} with blow-up ")
plt.savefig(path_out / 'plots' / 'lossmaps' / f'lossmap_blowup_over_turns_B{beam}{plane}.png', dpi=300)
plt.show()
