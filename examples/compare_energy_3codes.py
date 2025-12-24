# This is a script to test the FLUKA integration into xtrack.
# Before running the script, the source code needs to be compiled. In the root
# package folder, run ./compile_fluka.sh to do so.

import numpy as np
import xpart as xp
import xtrack as xt
import xcoll as xc
import time

import matplotlib.pyplot as plt


num_part = 1_000_000
capacity = 2*num_part
particle_ref = xt.Particles('proton', p0c=6.8e12)
nbins = 200
E_low = 1.e12
E_high = particle_ref.energy0[0]
unit = 1.e12

# num_part = 1_000_000
# capacity = 2*num_part
# particle_ref = xt.Particles('proton', p0c=450e9)
# nbins = 200
# E_low = 100.e9
# E_high = particle_ref.energy0[0]
# unit = 1.e9


# Stop any running engines
if xc.fluka.engine.is_running():
    xc.fluka.engine.stop(clean=True)
if xc.geant4.engine.is_running():
    xc.geant4.engine.stop(clean=True)

coll1 = xc.EverestCollimator(length=0.6, material=xc.materials.MolybdenumGraphite, jaw=0.001)
coll2 = xc.FlukaCollimator(length=0.6, material='mogr', jaw=0.001)
coll3 = xc.Geant4Collimator(length=0.6, material=xc.materials.MolybdenumGraphite, jaw=0.001)

# Connect to FLUKA
xc.fluka.engine.particle_ref = particle_ref
xc.fluka.engine.capacity = capacity
xc.fluka.engine.seed = 5656565
xc.fluka.engine.start(elements=coll2, clean=True, verbose=False)
xc.fluka.engine.init_tracking(max_particle_id=num_part)  # Prepares FLUKA for tracking, needed because we will track in chunks

# Connect to Geant4
xc.geant4.engine.particle_ref = particle_ref
xc.geant4.engine.seed = 5656565
xc.geant4.engine.start(elements=coll3, clean=True, verbose=False)


# Create an initial distribution of particles, random in 4D, on the left jaw (with the
# longitudinal coordinates set to zero)
np.random.seed(seed=23823842)
x_init   = np.random.normal(loc=0.001, scale=0.2e-3, size=num_part)
px_init  = np.random.normal(loc=0., scale=5.e-6, size=num_part)
y_init   = np.random.normal(loc=0., scale=1e-3, size=num_part)
py_init  = np.random.normal(loc=0., scale=5.e-6, size=num_part)
part_init = xp.build_particles(x=x_init, px=px_init, y=y_init, py=py_init,
                               particle_ref=xc.fluka.engine.particle_ref, _capacity=capacity)


# Do the tracking in Everest
print(f"Tracking {num_part} particles (Everest)...   ", end='', flush=True)
start = time.time()
part = part_init.copy()
coll1.track(part)
energy_everest = part.energy[part.state > 0]/unit
print(f"Done in {round(1000*(time.time()-start), 3)}ms")


# Do the tracking in FLUKA
start = time.time()
energy_fluka = np.array([], dtype=np.float64)
# We split it in small chunks to be more efficient
num_step = 10_000
steps = num_part // num_step
for i in range(steps):
    print(f"Tracking {num_part} particles (FLUKA)...     Step {i}/{steps}", flush=True, end='\r')
    mask = np.zeros(capacity, dtype=np.bool)
    mask[i*num_step:(i+1)*num_step] = True
    mask[num_part + i*num_step:num_part + (i+1)*num_step] = True  # For secondaries
    part = part_init.filter(mask)
    coll2.track(part)
    mask = (part.state > 0) & (part.energy >= E_low)
    parents, inv = np.unique(part.parent_particle_id[mask], return_inverse=True)
    this_energy_fluka = np.zeros_like(parents, dtype=np.float64)
    np.add.at(this_energy_fluka, inv, part.energy[mask]/unit)
    energy_fluka = np.concatenate((energy_fluka, this_energy_fluka))
if num_part > steps*num_step:  # Some leftover particles
    print(f"Tracking {num_part} particles (FLUKA)...     Step {steps}/{steps}", flush=True, end='\r')
    mask = np.zeros(capacity, dtype=np.bool)
    mask[steps*num_step:num_part] = True
    mask[num_part + steps*num_step:num_part + num_part] = True  # For secondaries
    part = part_init.filter(mask)
    coll2.track(part)
    mask = (part.state > 0) & (part.energy >= E_low)
    parents, inv = np.unique(part.parent_particle_id[mask], return_inverse=True)
    this_energy = np.zeros_like(parents, dtype=np.float64)
    np.add.at(this_energy, inv, part.energy[mask]/unit)
    energy_fluka = np.concatenate((energy_fluka, this_energy))
print(f"Tracking {num_part} particles (FLUKA)...     Done in {round(time.time()-start, 3)}s.")


# Do the tracking in Geant4
start = time.time()
energy_geant4 = np.array([], dtype=np.float64)
# We split it in small chunks to be more efficient
num_step = 10_000
steps = num_part // num_step
for i in range(steps):
    print(f"Tracking {num_part} particles (Geant4)...     Step {i}/{steps}", flush=True, end='\r')
    mask = np.zeros(capacity, dtype=np.bool)
    mask[i*num_step:(i+1)*num_step] = True
    mask[num_part + i*num_step:num_part + (i+1)*num_step] = True  # For secondaries
    part = part_init.filter(mask)
    coll3.track(part)
    mask = (part.state > 0) & (part.energy >= E_low)
    parents, inv = np.unique(part.parent_particle_id[mask], return_inverse=True)
    this_energy = np.zeros_like(parents, dtype=np.float64)
    np.add.at(this_energy, inv, part.energy[mask]/unit)
    energy_geant4 = np.concatenate((energy_geant4, this_energy))
if num_part > steps*num_step:  # Some leftover particles
    print(f"Tracking {num_part} particles (Geant4)...     Step {steps}/{steps}", flush=True, end='\r')
    mask = np.zeros(capacity, dtype=np.bool)
    mask[steps*num_step:num_part] = True
    mask[num_part + steps*num_step:num_part + num_part] = True  # For secondaries
    part = part_init.filter(mask)
    coll3.track(part)
    mask = (part.state > 0) & (part.energy >= E_low)
    parents, inv = np.unique(part.parent_particle_id[mask], return_inverse=True)
    this_energy = np.zeros_like(parents, dtype=np.float64)
    np.add.at(this_energy, inv, part.energy[mask]/unit)
    energy_geant4 = np.concatenate((energy_geant4, this_energy))
print(f"Tracking {num_part} particles (Geant4)...     Done in {round(time.time()-start, 3)}s.")


# Create bins based on the min and max energy
bins = np.logspace(np.log10(E_low/unit), np.log10(E_high/unit), nbins + 1)
bin_centres = np.sqrt(bins[:-1] * bins[1:])
log_edges = np.log10(bins)
dlog = np.diff(log_edges)

fig, ax = plt.subplots(figsize=(14, 4))
for ee, label in zip([energy_everest, energy_fluka, energy_geant4],
                     ['Everest', 'FLUKA', 'Geant4']):
    counts, edges = np.histogram(ee, bins=bins)
    bin_widths = np.diff(edges)
    dNdlogE = counts / (len(ee) * dlog)
    ax.step(bin_centres, dNdlogE, where='mid', label=label)

# Set horizontal axis to logarithmic scale
ax.set_xscale('log')
ax.set_yscale('log')

# Labels and legend
ax.set_xlabel(fr'Energy [$10^{{{int(np.log10(unit))}}}$eV]')
ax.set_ylabel('Normalised frequency ' + r'$\frac{dN}{d\log E}$')
ax.legend()
ax.grid(True, which='both', linestyle='--', linewidth=0.5)
plt.tight_layout()
plt.savefig('plots/compare_energy_3codes.png', dpi=300)
plt.show()


# Stop the engines
xc.fluka.engine.stop(clean=True)
xc.geant4.engine.stop(clean=True)
