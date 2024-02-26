import numpy as np
from pathlib import Path
import time
start_time = time.time()
import sys, os, contextlib

import xobjects as xo
import xpart as xp
import xtrack as xt
import xcoll as xc


context = xo.ContextCpu(omp_num_threads='auto')

# On a modern CPU, we get ~5000 particle*turns/s
# So this script should take around half an hour
beam          = 2
plane         = 'H'
num_turns     = 200
num_particles = 50000
engine        = 'everest'

path_in  = xc._pkg_root.parent / 'examples'
path_out = Path.cwd()


# Load from json
with open(os.devnull, 'w') as fid:
    with contextlib.redirect_stdout(fid):
        line = xt.Line.from_json(path_in / 'machines' / f'lhc_run3_b{beam}.json')


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
# part = xc.generate_pencil_on_collimator(line, tcpc, num_particles=num_particles,
#                                         nemitt_x=3.5e-6, nemitt_y=3.5e-6)
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
line_is_reversed = True if beam == 2 else False
ThisLM = xc.LossMap(line, line_is_reversed=line_is_reversed, part=part)
ThisLM.to_json(file=Path(path_out, f'lossmap_B{beam}{plane}.json'))

# Save a summary of the collimator losses to a text file
ThisLM.save_summary(file=Path(path_out, f'coll_summary_B{beam}{plane}.out'))
print(ThisLM.summary)



# Impacts
# =======
#nabs = summary.loc[summary.collname==tcpc, 'nabs'].values[0]
#imp = coll_manager.impacts.to_pandas()
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

exit()
