import json
import numpy as np
from pathlib import Path
import xtrack as xt
import xpart as xp
import xcoll as xc


# On a modern CPU, we get ~5000 particle*turns/s
# So this script should take around 45 minutes
beam      = 1
lmtype    = 'DPpos'

num_particles  = 2000
sweep          = 300
sweep          = -sweep if lmtype == 'DPpos' else sweep
pretrack_turns = 500
num_turns      = 20*abs(sweep)
at_element     = 'ip3'
engine         = 'everest'

path_in  = xc._pkg_root.parent / 'examples'
path_out = Path.cwd()


# Load from json
line = xt.Line.from_json(path_in / 'machines' / f'lhc_run3_b{beam}.json')


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


# Install collimators into line
if engine == 'everest':
    coll_manager.install_everest_collimators(verbose=True)
else:
    raise ValueError(f"Unknown scattering engine {engine}!")


# Aperture model check
print('\nAperture model check after introducing collimators:')
df_with_coll = line.check_aperture()
assert not np.any(df_with_coll.has_aperture_problem)

    
# Build the tracker
coll_manager.build_tracker()


# Set the collimator openings based on the colldb,
# or manually override with the option gaps={collname: gap}
coll_manager.set_openings()


# Generate initial halo distribution in IP3  (better matching)
# This creates amplitudes on a circle with R = [4, 6], e.g.:
# plt.scatter(np.sqrt(x_norm**2 + px_norm**2), np.sqrt(y_norm**2 + py_norm**2))
x_norm, px_norm, _, _ = xp.generate_2D_uniform_circular_sector(num_particles, r_range=(4, 6))
y_norm, py_norm, _, _ = xp.generate_2D_uniform_circular_sector(num_particles, r_range=(4, 6))
zeta, delta = xp.generate_longitudinal_coordinates(
                    num_particles=num_particles, distribution='gaussian', sigma_z=7.55e-2, line=line
            )
theta = np.random.uniform(0, np.pi/2, num_particles)
x_norm  *= np.cos(theta)
px_norm *= np.cos(theta)
y_norm  *= np.sin(theta)
py_norm *= np.sin(theta)
part = xp.build_particles(
        x_norm=x_norm, px_norm=px_norm, y_norm=y_norm, py_norm=py_norm, zeta=zeta, delta=delta,
        nemitt_x=3.5e-6, nemitt_y=3.5e-6, line=line, at_element=at_element
)


# Optimise the line
line.optimize_for_tracking(keep_markers=[at_element])
idx = line.element_names.index(at_element)
part.at_element = idx
part.start_tracking_at_element = idx


# Do a pre-tracking to empty the halo
coll_manager.enable_scattering()
line.track(part, num_turns=pretrack_turns, time=True)
coll_manager.disable_scattering()
print(f"Done pre-tracking to empty halo in {line.time_last_track:.1f}s.")

# We remove the lost particles, as we don't want them in the loss map
part = part.filter(part.state == 1)

# Track during RF sweep:
coll_manager.enable_scattering()
coll_manager.rf_sweep(sweep=sweep, num_turns=num_turns, particles=part, time=True)
coll_manager.disable_scattering()
print(f"Done sweeping RF in {line.time_last_track:.1f}s.")


# Save lossmap to json, which can be loaded, combined (for more statistics),
# and plotted with the 'lossmaps' package
_ = coll_manager.lossmap(part, file=Path(path_out,f'lossmap_B{beam+plane}.json'))

# Save a summary of the collimator losses to a text file
summary = coll_manager.summary(part, file=Path(path_out,f'coll_summary_B{beam+plane}.out'))
print(summary)

exit()
