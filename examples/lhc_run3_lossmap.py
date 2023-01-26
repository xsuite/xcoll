import json
import numpy as np
from pathlib import Path
import xtrack as xt
import xcoll as xc

path_in  = Path.cwd()
path_out = Path.cwd()

beam   = '1'
plane  = 'h'
engine = 'everest'



# Load from json
with open(Path(path_in,'machines',f'lhc_run3_b{beam}.json'), 'r') as fid:
    loaded_dct = json.load(fid)
line = xt.Line.from_dict(loaded_dct)


# Aperture model check
print('\nAperture model check on imported model:')
df_imported = line.check_aperture()
assert not np.any(df_imported.has_aperture_problem)


# Initialise collmanager
line_is_reversed = True if beam=='2' else False
coll_manager = xc.CollimatorManager(
    line=line, line_is_reversed=line_is_reversed,
    colldb=xc.load_SixTrack_colldb(Path(path_in,'colldb',f'lhc_run3_b{beam}.dat'), emit=3.5e-6)
    )


# Install collimators in line as everest collimators
if engine == 'everest':
    coll_manager.install_everest_collimators(verbose=True)
elif engine == 'pyeverest':
    coll_manager.install_pyeverest_collimators(verbose=True)
elif engine == 'k2':
    coll_manager.install_k2_collimators(verbose=True)
else:
    raise ValueError()

    
# Build the tracker
coll_manager.build_tracker()


# Set the collimator openings based on the colldb,
# or manually override with the option gaps={collname: gap}
coll_manager.set_openings()


# Aperture model check
print('\nAperture model check after introducing collimators:')
df_with_coll = line.check_aperture()
assert not np.any(df_with_coll.has_aperture_problem)


# Generate initial pencil distribution on horizontal collimator
tcp  = f"tcp.{'c' if plane=='h' else 'd'}6{'l' if beam=='1' else 'r'}7.b{beam}"
part = coll_manager.generate_pencil_on_collimator(tcp, num_particles=50000)


# Track
coll_manager.track(part, num_turns=200)


# Get losses for lossmap
coll_manager.create_lossmap(part)


# Save to json
# These files can be loaded, combined (for more statistics), and plotted with the 'lossmaps' package
with open(Path(path_out,f'lossmap_B{beam+plane}.json'), 'w') as fid:
    json.dump(coll_manager.lossmap, fid, indent=True)


summary = coll_manager.coll_summary(part)
with open(Path(path_out,f'coll_summary_B{beam+plane}.json'), 'w') as fid:
    json.dump(summary.to_json(), fid, indent=True)

print(summary)

exit()
