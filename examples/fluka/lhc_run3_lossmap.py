import json
import numpy as np
from pathlib import Path
import xpart as xp
import xtrack as xt
import xcoll as xc


# On a modern CPU, we get ~5000 particle*turns/s
# So this script should take around half an hour
beam          = '1'
plane         = 'H'

num_turns     = 200
num_particles = 5000
engine        = 'fluka'

path_in  = xc._pkg_root.parent / 'examples'
path_out = Path.cwd()


# Start FLUKA server
xc.FlukaEngine.start_server("lhc_run3_30cm.inp")


# Load from json
line = xt.Line.from_json(path_in / 'machines' / f'lhc_run3_b{beam}.json')


# Aperture model check
print('\nAperture model check on imported model:')
df_imported = line.check_aperture()
assert not np.any(df_imported.has_aperture_problem)


# Initialise collmanager
coll_manager = xc.CollimatorManager.from_yaml(path_in / 'colldb' / f'lhc_run3.yaml', line=line, beam=beam)


# Hack: because of how FLUKA extends the collimators, it might overlap with some
#       elements, so we manually remove them
to_remove = []
ss = line.get_s_position()
for coll, val in xc.FlukaEngine().collimators.items():
    l = val['length']
    idx = line.element_names.index(coll)
    s_coll = ss[idx]
    i = idx - 1
    # We remove everything between the beginning and end of the collimator except drifts
    while ss[i] >= s_coll - l/2:
        el = line[i]
        if el.__class__.__name__ == 'Drift':
            i -= 1
            continue
        nn = line.element_names[i]
        if nn==f'{coll}_aper' or nn==f'{coll}_aper_patch':
            i -= 1
            continue
        if hasattr(el, 'length') and el.length > 0:
            raise ValueError()
        to_remove.append(nn)
        i -= 1
    i = idx + 1
    while ss[i] <= s_coll + l/2:
        el = line[i]
        if el.__class__.__name__ == 'Drift':
            i += 1
            continue
        nn = line.element_names[i]
        if nn==f'{coll}_aper' or nn==f'{coll}_aper_patch':
            i += 1
            continue
        if hasattr(el, 'length') and el.length > 0:
            raise ValueError()
        to_remove.append(nn)
        i += 1
new_line_names = []
for name in line.element_names:
    if name not in to_remove:
        new_line_names.append(name)
line.element_names = new_line_names


# Install collimators into line
if engine == 'fluka':
    coll_manager.install_fluka_collimators(verbose=True)
else:
    raise ValueError(f"Unknown scattering engine {engine}!")


# Aperture model check
print('\nAperture model check after introducing collimators:')
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
tcp  = f"tcp.{'c' if plane=='H' else 'd'}6{'l' if beam=='1' else 'r'}7.b{beam}"
part = coll_manager.generate_pencil_on_collimator(tcp, num_particles=num_particles)
# x_norm, px_norm, _, _ = xp.generate_2D_uniform_circular_sector(r_range=(5, 5.04), num_particles=num_particles)
# y_norm  = np.random.normal(scale=0.01, size=num_particles)
# py_norm = np.random.normal(scale=0.01, size=num_particles)
# part = line.build_particles(
#             x_norm=x_norm, px_norm=px_norm, y_norm=y_norm, py_norm=py_norm,
#             nemitt_x=coll_manager.colldb.emittance[0], nemitt_y=coll_manager.colldb.emittance[1],
#             at_element=tcp, match_at_s=coll_manager.s_active_front[tcp],
#             _context=coll_manager._buffer.context
# )


# Optimise the line
# line.optimize_for_tracking()
# idx = line.element_names.index(tcp)
# part.at_element = idx
# part.start_tracking_at_element = idx

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


exit()
