import json
import numpy as np
from pathlib import Path
import time
start_time = time.time()
import sys, os, contextlib

import xpart as xp
import xtrack as xt
import xcoll as xc

if len(sys.argv) < 4:
    raise ValueError("Need three arguments: beam, plane, and tilt.")
beam  = int(sys.argv[1])
plane = str(sys.argv[2])
tilt  = int(sys.argv[3])*1.e-6


# On a modern CPU, we get ~5000 particle*turns/s
# So this script should take a bit less than 40 minutes
# beam          =  1
# plane         = 'H'
# tilt          = 0         # rad !!

num_turns     = 200
num_particles = 5000
engine        = 'everest'

path_in  = xc._pkg_root.parent / 'examples'
path_out = Path.cwd()


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
coll_manager = xc.CollimatorManager.from_yaml(path_in / 'colldb' / f'lhc_run3_crystals.yaml', line=line,
                                              beam=beam, ignore_crystals=False, record_impacts=True)


# Install collimators into line
if engine == 'everest':
    coll_manager.install_everest_collimators(verbose=True)
else:
    raise ValueError(f"Unknown scattering engine {engine}!")


# Aperture model check
print('\nAperture model check after introducing collimators:')
with open(os.devnull, 'w') as fid:
    with contextlib.redirect_stdout(fid):
        df_with_coll = line.check_aperture()
assert not np.any(df_with_coll.has_aperture_problem)


# Primary crystal
tcpc = f"tcpc{plane.lower()}.a{6 if plane=='V' else 4 if beam==1 else 5}{'l' if beam==1 else 'r'}7.b{beam}"


# Build the tracker
coll_manager.build_tracker()


# Set the collimator openings based on the colldb,
# or manually override with the option gaps={collname: gap}
coll_manager.colldb.active        = {tcpc: True}
coll_manager.colldb.gap           = {tcpc: 5}
coll_manager.set_openings()


# # Generate initial pencil distribution on crystal
# part = coll_manager.generate_pencil_on_collimator(tcpc, num_particles=num_particles)
# Generate initial halo
x_norm, px_norm, _, _ = xp.generate_2D_uniform_circular_sector(r_range=(5, 5.04), num_particles=num_particles)
y_norm  = np.random.normal(scale=0.01, size=num_particles)
py_norm = np.random.normal(scale=0.01, size=num_particles)
part = line.build_particles(
            x_norm=x_norm, px_norm=px_norm, y_norm=y_norm, py_norm=py_norm,
            nemitt_x=coll_manager.colldb.emittance[0], nemitt_y=coll_manager.colldb.emittance[1],
            at_element=tcpc, match_at_s=coll_manager.s_active_front[tcpc],
            _buffer=coll_manager._part_buffer
)

# Optimise the line
line.optimize_for_tracking()
idx = line.element_names.index(tcpc)
part.at_element = idx
part.start_tracking_at_element = idx


# Apply settings
line[tcpc].tilt_L        = tilt
line[tcpc].bending_angle = 40.e-6
line[tcpc].xdim          = 0.002
line[tcpc].ydim          = 0.05


# Track
coll_manager.enable_scattering()
line.track(part, num_turns=num_turns, time=True)
coll_manager.disable_scattering()
print(f"Done tracking in {line.time_last_track:.1f}s.")


# Save lossmap to json, which can be loaded, combined (for more statistics),
# and plotted with the 'lossmaps' package
# _ = coll_manager.lossmap(part, file=Path(path_out,f'lossmap_B{beam}{plane}.json'))


# Save a summary of the collimator losses to a text file
summary = coll_manager.summary(part, file=Path(path_out, f'cry_{round(tilt*1e6)}urad_summary_B{beam}{plane}.out'))
nabs = summary.loc[summary.collname==tcpc, 'nabs'].values[0]
print(summary)


# Impacts
# =======
imp = coll_manager.impacts.to_pandas()
outfile = Path(path_out, f'cry_{round(tilt*1e6)}urad_impacts_B{beam}{plane}.json')
imp.to_json(outfile)

# Only keep the first impact
imp.drop_duplicates(subset=['parent_id'], inplace=True)
assert np.unique(imp.interaction_type) == ['Enter Jaw']

# Only keep those first impacts that are on the crystal
first_impacts = imp[imp.collimator==tcpc].parent_id
num_first_impacts = len(first_impacts)
assert num_first_impacts==len(np.unique(first_impacts))

ineff = nabs/num_first_impacts
outfile = Path(path_out, f'cry_{round(tilt*1e6)}urad_ineff_B{beam}{plane}.json')
with outfile.open('w') as fid:
    json.dump({'first_impacts': num_first_impacts, 'nabs': nabs, 'ineff': ineff}, fid)

print(f"Out of {num_first_impacts} particles hitting the crystal {tcpc} (angle {round(tilt*1e6)}urad) first, {round(nabs)} are absorbed in the crystal (for an inefficiency of {ineff:5}.")
print(f"Critical angle is {round(line[tcpc].critical_angle*1.e6, 1)}urad.")
print(f"Total calculation time {time.time()-start_time}s")

exit()
