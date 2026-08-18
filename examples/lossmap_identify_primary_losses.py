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


context = xo.ContextCpu(omp_num_threads='auto')  # For CPU
# context = xo.ContextCupy()                     # For CUDA GPUs
# context = xo.ContextPyopencl()                 # For OpenCL GPUs

beam          = 1
plane         = 'V'
num_turns     = 500
num_particles = 25_000

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
adt = xc.BlowUp.install(line, name=f'{name}_blowup', at=adt_pos, plane=plane, stop_at_turn=num_turns,
                        amplitude=1.5, use_individual_kicks=False)


# Aperture model check
print('\nAperture model check after introducing elements:')
df_with_coll = line.check_aperture()
assert not np.any(df_with_coll.has_aperture_problem)


# Assign the optics to deduce the gap settings, and calibrate the ADT
tw = line.twiss()
line.xcoll.collimators.assign_optics(twiss=tw)
if plane == 'H':
    adt.calibrate_by_emittance(nemitt=colldb.nemitt_x, twiss=tw)
else:
    adt.calibrate_by_emittance(nemitt=colldb.nemitt_y, twiss=tw)


# Mark secondary particles during tracking, to be able to distinguish them
# from primaries in the loss map later
line.xcoll.scattering.identify_primary_losses()


# Bring one secondary very close to the primaries to induce hierarchy breaking
line['tcsg.d4l7.b1'].gap = 5.2


# Optimise the line
line.optimize_for_tracking()


# Move to a more efficient context for tracking
line.discard_tracker()
line.build_tracker(_context=context)


# Generate a matched Gaussian bunch
part = xp.generate_matched_gaussian_bunch(num_particles=num_particles,
                                          total_intensity_particles=1.6e11,
                                          nemitt_x=colldb.nemitt_x,
                                          nemitt_y=colldb.nemitt_y,
                                          sigma_z=7.55e-2,
                                          line=line,
                                          _context=context)


# Track!
line.xcoll.scattering.enable()
adt.activate()
line.track(part, num_turns=num_turns, time=True, with_progress=1)
adt.deactivate()
line.xcoll.scattering.disable()
print(f"Done tracking in {line.time_last_track:.1f}s.")


# Move the line back to CPU to be able to use all prebuilt kernels for the aperture interpolation
if not isinstance(context, xo.ContextCpu):
    line.discard_tracker()
    line.build_tracker(_context=xo.ContextCpu())


# Save loss map to json
lm_time = time.time()
line_is_reversed = True if f'{beam}' == '2' else False
ThisLM = xc.LossMap(line, line_is_reversed=line_is_reversed, part=part)
print(f"Loss map created in {time.time()-lm_time:.1f}s.")
ThisLM.to_json(file=path_out / 'results' / f'lossmap_primary_B{beam}{plane}.json')

# Save a summary of the collimator losses to a text file
ThisLM.save_summary(file=path_out / 'results' / f'coll_summary_primary_B{beam}{plane}.out')
print(ThisLM.summary)

print(f"Total calculation time {time.time()-start_time}s")

# Plot loss map
ThisLM.plot(savefig=path_out / 'plots' / 'lossmaps' / f'lossmap_primary_B{beam}{plane}.pdf', zoom='betatron')
plt.show()

