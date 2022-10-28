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

# Initialise collmanager,on the specified buffer
coll_manager = xc.CollimatorManager(
    line=line,
    line_is_reversed=True,
    colldb=xc.load_SixTrack_colldb('colldb/lhc_run3_b2.dat', emit=3.5e-6)
    )

# Install collimators in line as black absorbers
coll_manager.install_everest_collimators(verbose=True)

# Build the tracker
coll_manager.build_tracker()

# Align the collimators
coll_manager.align_collimators_to('front')

# Set the collimator openings based on the colldb,
# or manually override with the option gaps={collname: gap}
coll_manager.set_openings()

# Aperture model check
print('\nAperture model check after introducing collimators:')
df_with_coll = line.check_aperture()
assert not np.any(df_with_coll.has_aperture_problem)


# Horizontal loss map
num_particles = 50000
coll = 'tcp.c6r7.b2'

# Collimator plane: generate pencil distribution in normalized coordinates
x_norm, px_norm, _, _ = xp.generate_2D_pencil(
                             num_particles=num_particles,
                             pos_cut_sigmas=coll_manager.colldb.gap[coll],
                             dr_sigmas=0.002,
                             side='+-')

# Other plane: generate gaussian distribution in normalized coordinates
y_norm = np.random.normal(scale=0.01, size=num_particles)
py_norm = np.random.normal(scale=0.01, size=num_particles)

part = xp.build_particles(
            tracker=coll_manager.line.tracker,
            x_norm=x_norm, px_norm=px_norm,
            y_norm=y_norm, py_norm=py_norm,
            scale_with_transverse_norm_emitt=coll_manager.colldb.emittance,
            at_element=coll,
            match_at_s=coll_manager.s_match[coll])

# Track
coll_manager.track(part, num_turns=1)

# Get losses for lossmap
coll_manager.create_lossmap(part)

# Save to json
# These files can be loaded, combined (for more statistics), and plotted with the 'lossmaps' package
with open('lossmap_B2H.json', 'w') as fid:
    json.dump(coll_manager.lossmap, fid, indent=True)

print(coll_manager.coll_summary(part))

exit()