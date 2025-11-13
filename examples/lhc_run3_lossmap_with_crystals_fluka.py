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


# We do the majority of the script on the default context to be able to use prebuilt kernels
context = xo.ContextCpu()


# This script takes around 3 minutes on a modern CPU (90s preparation+interpolation, 90s tracking)
beam = 1
plane = 'H'
tilt = 0  # rad !!

num_turns = 2 # 00
num_particles = 50 #000

path_in  = xc._pkg_root.parent / 'examples'
path_out = Path.cwd()


# Load from json
line = xt.Line.from_json(path_in / 'machines' / f'lhc_run3_b{beam}.json')


# Initialise colldb
colldb = xc.CollimatorDatabase.from_yaml(path_in / 'colldb' / f'lhc_run3_crystals.yaml',
                                               beam=beam, ignore_crystals=False)


# Install collimators into line
colldb.install_fluka_collimators(line=line, verbose=True)

# Aperture model check
print('\nAperture model check after introducing collimators:')
df_with_coll = line.check_aperture()
assert not np.any(df_with_coll.has_aperture_problem)


# Primary crystal
tcpc = f"tcpc{plane.lower()}.a{6 if plane=='V' else 4 if f'{beam}'=='1' else 5}{'l' if f'{beam}'=='1' else 'r'}7.b{beam}"

tw = line.twiss4d()

alfa_x = tw.rows[tcpc].alfx[0]
alfa_y = tw.rows[tcpc].alfy[0]
beta_x  = tw.rows[tcpc].betx[0]
beta_y  = tw.rows[tcpc].bety[0]
gap_tcpc = colldb[tcpc]['gap']
gamma0 = line.particle_ref.gamma0
beta0 = line.particle_ref.beta0

# import pdb; pdb.set_trace()

geo_emitt_x = colldb.nemitt_x/(gamma0*beta0)
geo_emitt_y = colldb.nemitt_y/(gamma0*beta0)

ch_x = -gap_tcpc * alfa_x * np.sqrt(geo_emitt_x/beta_x)
ch_y = -gap_tcpc * alfa_y * np.sqrt(geo_emitt_y/beta_y)

# Apply settings
line[tcpc].bending_angle = ch_x # 40.e-6
line[tcpc].width         = 0.002
line[tcpc].height        = 0.05
line[tcpc].length         = 0.004
# line[tcpc].align_to_beam_divergence()

# Build the tracker
line.build_tracker()


# Assign the optics to deduce the gap settings
line.collimators.assign_optics()

# Connect to FLUKA
xc.fluka.engine.particle_ref = xt.Particles.reference_from_pdg_id(pdg_id='Pb208', p0c=6.8e12*82)
xc.fluka.engine.capacity = 20*num_particles
xc.fluka.engine.seed = 5656565
xc.fluka.engine.start(line=line, capacity=xc.fluka.engine.capacity, cwd='run_fluka_temp', clean=False, verbose=True, return_ions=True)


# Optimise the line
# line.optimize_for_tracking()


# # Generate initial pencil distribution on crystal
# part = line[tcp].generate_pencil(num_particles)
# Generate initial halo
x_norm, px_norm, _, _ = xp.generate_2D_uniform_circular_sector(r_range=(5, 5.04), num_particles=num_particles)
y_norm  = np.random.normal(scale=0.01, size=num_particles)
py_norm = np.random.normal(scale=0.01, size=num_particles)
part = line.build_particles(
            x_norm=x_norm, px_norm=px_norm, y_norm=y_norm, py_norm=py_norm,
            nemitt_x=line[tcpc].nemitt_x, nemitt_y=line[tcpc].nemitt_y,
            at_element=tcpc, particle_ref=xc.fluka.engine.particle_ref, _context=context, _capacity=xc.fluka.engine.capacity)


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

# Stop the FLUKA connection (and return to the previous directory)
xc.fluka.engine.stop(clean=False)

# Save lossmap to json, which can be loaded, combined (for more statistics),
# and plotted with the 'lossmaps' package
line_is_reversed = True if f'{beam}' == '2' else False
ThisLM = xc.LossMap(line, line_is_reversed=line_is_reversed, part=part)
ThisLM.to_json(file=Path(path_out, f'lossmap_B{beam}{plane}.json'))

# Save a summary of the collimator losses to a text file
ThisLM.save_summary(file=Path(path_out, f'coll_summary_B{beam}{plane}.out'))
print(ThisLM.summary)


# Impacts
# =======
#nabs = summary.loc[summary.collname==tcpc, 'nabs'].values[0]
#imp = colldb.impacts.to_pandas()
#outfile = Path(path_out, f'cry_{round(tilt*1e6)}urad_impacts_B{beam}{plane}.json')
#imp.to_json(outfile)

# Only keep the first impact
#imp.drop_duplicates(subset=['parent_id'], inplace=True)
#assert np.unique(imp.interaction_type) == ['Enter Jaw']

# Only keep those first impacts that are on the crystal
#first_impacts = imp[imp.collimator==tcpc].parent_id
#num_first_impacts = len(first_impacts)
#assert num_first_impacts==len(np.unique(first_impacts))

#ineff = nabs/num_first_impacts
#outfile = Path(path_out, f'cry_{round(tilt*1e6)}urad_ineff_B{beam}{plane}.json')
#with outfile.open('w') as fid:
#    json.dump({'first_impacts': num_first_impacts, 'nabs': nabs, 'ineff': ineff}, fid)

#print(f"Out of {num_first_impacts} particles hitting the crystal {tcpc} (angle {round(tilt*1e6)}urad) first, {round(nabs)} are absorbed in the crystal (for an inefficiency of {ineff:5}.")
#print(f"Critical angle is {round(line[tcpc].critical_angle*1.e6, 1)}urad.")

print(f"Total calculation time {time.time()-start_time}s")

ThisLM.plot()
plt.show()
