import json
import numpy    as np
import matplotlib.pyplot as plt

import xobjects as xo
import xtrack   as xt
import xpart    as xp
import xcoll    as xc


# Load from json
with open('machines/lhc_run3_b2.json', 'r') as fid:
    loaded_dct = json.load(fid)
line = xt.Line.from_dict(loaded_dct)

# Aperture model check
print('\nAperture model check on imported model:')
df_imported = line.check_aperture()
assert not np.any(df_imported.has_aperture_problem)

# Initialise collmanager
coll_manager = xc.CollimatorManager(
    line=line,
    line_is_reversed=True,
    colldb=xc.load_SixTrack_colldb('colldb/lhc_run3_b2.dat', emit=3.5e-6)
    )

# Install collimators in line as everest collimators
coll_manager.install_everest_collimators(verbose=True)

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
part = coll_manager.generate_pencil_on_collimator('tcp.d6r7.b2', num_particles=50000)

# Track
coll_manager.track(part, num_turns=200)

# Get losses for lossmap
coll_manager.create_lossmap(part)

# Save to json
# These files can be loaded, combined (for more statistics), and plotted with the 'lossmaps' package
with open('lossmap_B2V.json', 'w') as fid:
    json.dump(coll_manager.lossmap, fid, indent=True)

print(coll_manager.coll_summary(part))

exit()