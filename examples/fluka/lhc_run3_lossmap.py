# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

import numpy as np
from pathlib import Path
import time
start_time = time.time()
import os, contextlib

import xobjects as xo
import xtrack as xt
import xpart as xp
import xcoll as xc


context = xo.ContextCpu()

beam          = 1
plane         = 'H'

num_turns     = 2
num_particles = 5000
engine        = 'fluka'

path_in  = xc._pkg_root.parent / 'examples'
path_out = Path.cwd()


# Start FLUKA server
xc.FlukaEngine.start_server("lhc_run3_30cm.inp", n_alloc=2*num_particles)


# Load from json
with open(os.devnull, 'w') as fid:
    with contextlib.redirect_stdout(fid):
        line = xt.Line.from_json(path_in / 'machines' / f'lhc_run3_b{beam}.json')


# Aperture model check
print('\nAperture model check on imported model:')
with open(os.devnull, 'w') as fid:
    with contextlib.redirect_stdout(fid):
        df_imported = line.check_aperture()
assert not np.any(df_imported.has_aperture_problem)


# Initialise collmanager
coll_manager = xc.CollimatorManager.from_yaml(path_in / 'colldb' / f'lhc_run3.yaml', line=line, beam=beam, _context=context)


# Install collimators into line
if engine == 'fluka':
    coll_manager.install_fluka_collimators(verbose=True)
else:
    raise ValueError(f"Unknown scattering engine {engine}!")


# Aperture model check
print('\nAperture model check after introducing collimators:')
with open(os.devnull, 'w') as fid:
    with contextlib.redirect_stdout(fid):
        df_with_coll = line.check_aperture()
assert not np.any(df_with_coll.has_aperture_problem)


# Build the tracker
coll_manager.build_tracker()


# Set FLUKA reference particle
particle_ref = xp.Particles.reference_from_pdg_id(pdg_id='proton', p0c=6.8e12)
xc.FlukaEngine().set_particle_ref(particle_ref)


coll_manager.set_openings()
# Assign optics manually
tw = line.twiss()
for coll in coll_manager.collimator_names:
    idx = line.element_names.index(coll)
    line[coll].ref_x = tw.x[idx]
    line[coll].ref_y = tw.y[idx]


# Generate initial pencil distribution on horizontal collimator
tcp  = f"tcp.{'c' if plane=='H' else 'd'}6{'l' if f'{beam}'=='1' else 'r'}7.b{beam}"
part = coll_manager.generate_pencil_on_collimator(tcp, num_particles=num_particles)


# Optimise the line
line.optimize_for_tracking()
idx = line.element_names.index(tcp)
part.at_element = idx
part.start_tracking_at_element = idx

print(part.at_element)
print(f"{np.mean(part.x)} +- {np.std(part.x)}")

# Track
coll_manager.enable_scattering()
line.track(part, num_turns=num_turns, time=True)
coll_manager.disable_scattering()
print(f"Done tracking in {line.time_last_track:.1f}s.")


# Save lossmap to json, which can be loaded, combined (for more statistics),
# and plotted with the 'lossmaps' package
_ = coll_manager.lossmap(part, file=Path(path_out,f'lossmap_B{beam}{plane}.json'))


# Save a summary of the collimator losses to a text file
summary = coll_manager.summary(part, file=Path(path_out,f'coll_summary_B{beam}{plane}.out'))
print(summary)
print(f"Total calculation time {time.time()-start_time}s")

exit()
