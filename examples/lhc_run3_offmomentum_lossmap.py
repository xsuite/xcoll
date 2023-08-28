import json
import numpy as np
from pathlib import Path
import xtrack as xt
import xpart as xp
import xcoll as xc


# On a modern CPU, we get ~5000 particle*turns/s
# This script typically takes around 1 hour
beam      = 1
plane    = 'DPpos'

num_particles  = 5000
sweep          = 300
sweep          = -abs(sweep) if plane == 'DPpos' else abs(sweep)
num_turns      = int(20*abs(sweep))
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
coll_manager = xc.CollimatorManager.from_yaml(path_in / 'colldb' / f'lhc_run3.yaml', line=line, beam=beam)


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


# Optimise the line
line.optimize_for_tracking()


# Generate initial matched bunch
part = xp.generate_matched_gaussian_bunch(nemitt_x=coll_manager.colldb.emittance[0],
                                          nemitt_y=coll_manager.colldb.emittance[1],
                                          sigma_z=7.55e-2, num_particles=num_particles, line=line)


# Track during RF sweep:
coll_manager.enable_scattering()
coll_manager.rf_sweep(sweep=sweep, num_turns=num_turns, particles=part, time=True)
coll_manager.disable_scattering()
print(f"Done sweeping RF in {line.time_last_track:.1f}s.")

# Save lossmap to json, which can be loaded, combined (for more statistics),
# and plotted with the 'lossmaps' package
coll_manager.lossmap(part, file=Path(path_out,f'lossmap_B{beam}{plane}.json'))


# Save a summary of the collimator losses to a text file
summary = coll_manager.summary(part, file=Path(path_out,f'coll_summary_B{beam}{plane}.out'))
print(summary)


# Let's visualise how the losses move from IR7 to IR3 during the sweep
tcp3 = line.element_names.index('tcp.6l3.b1')
turns_ip3 = part.at_turn[part.at_element == tcp3]
tcp7v = line.element_names.index('tcp.d6l7.b1')
turns_ip7v = part.at_turn[part.at_element == tcp7v]
tcp7h = line.element_names.index('tcp.c6l7.b1')
turns_ip7h = part.at_turn[part.at_element == tcp7h]
tcp7s = line.element_names.index('tcp.b6l7.b1')
turns_ip7s = part.at_turn[part.at_element == tcp7s]

plt.figure(figsize=(8,5))
plt.hist([turns_ip7v, turns_ip7h, turns_ip7s, turns_ip3], 100, stacked=True,
         label=['Lost on TCP.D7', 'Lost on TCP.C7', 'Lost on TCP.B7', 'Lost on TCP.3'])
plt.xlabel('turn')
plt.ylabel('number of lost particles (stacked)')
plt.legend()
plt.tight_layout()
plt.savefig('OffMomentumPosition.png', dpi=300)
plt.show()


# In this example, the IR3 region gets populated very fast. When this is not the case, e.g. for injection,
# one should keep only the last part of the sweep (as this is what happens in reality), as follows:


# In reality the loss map only shows the last second (1.3s BLM integration time)
# A typical RF sweep in reality shifts around 25Hz per second, i.e. around 2.2mHz/turn.
# We sweep 50mHz/turn, so we are ~22.5 times faster than in reality.
# We keep the equivalent of the last 3 seconds (slightly higher integration time, in order not to lose too much statistics).
# We want to keep the last 3/22.5 seconds = 1500 turns.
# This is equal to:
seconds_to_keep = 3
turns_to_keep = seconds_to_keep*25*num_turns/abs(sweep)
last_turn = part.at_turn.max()
part2 = part.filter(part.at_turn > last_turn - turns_to_keep)
print(f"Keeping last {int(turns_to_keep)} turns (equivalent integration time of {seconds_to_keep} seconds).")
print(f"This means we use {len(part2.x)} particles (of which {len(part2.x[part2.state > 0])} survived).")


# Save lossmap to json, which can be loaded, combined (for more statistics),
# and plotted with the 'lossmaps' package
coll_manager.lossmap(part2, file=Path(path_out,f'lossmap_B{beam}{plane}_end.json'))


# Save a summary of the collimator losses to a text file
summary = coll_manager.summary(part2, file=Path(path_out,f'coll_summary_B{beam}{plane}_end.out'))
print(summary)

