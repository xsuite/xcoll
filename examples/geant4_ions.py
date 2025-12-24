import numpy as np
import xpart as xp
import xtrack as xt
import xtrack.particles.pdg as pdg
import xcoll as xc
import time

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm


if xc.geant4.engine.is_running():
    xc.geant4.engine.stop()

num_part = 10_000
capacity = 25*num_part
particle_ref = xt.Particles('Pb208', p0c=6.8e12*82)

# Create a Geant4 collimator
coll = xc.Geant4Collimator(length=0.05, material='mogr', jaw=0.001)

# Connect to Geant4
xc.geant4.engine.particle_ref = particle_ref
xc.geant4.engine.return_none = True
xc.geant4.engine.return_ions = True
xc.geant4.engine.start(elements=coll, relative_energy_cut=1e-3, clean=True, verbose=False)

# Create an initial distribution of particles, random in 4D, on the left jaw (with the
# longitudinal coordinates set to zero)
x_init   = np.random.normal(loc=0.011, scale=0.2e-3, size=num_part)
px_init  = np.random.normal(loc=0, scale=5.e-6, size=num_part)
y_init   = np.random.normal(loc=0., scale=1e-3, size=num_part)
py_init  = np.random.normal(loc=0., scale=5.e-6, size=num_part)
part_init = xp.build_particles(x=x_init, px=px_init, y=y_init, py=py_init,
                                particle_ref=xc.geant4.engine.particle_ref,
                                _capacity=capacity)
part = part_init.copy()

# Do the tracking in Geant4
print(f"Tracking {num_part} {pdg.get_name_from_pdg_id(particle_ref.pdg_id[0])}s...     ", end='', flush=True)
start = time.time()
coll.track(part)
print(f"Done in {round(time.time()-start, 3)}s.", flush=True)

# Stop the Geant4 server
xc.geant4.engine.stop(clean=True)

_, A, Z, _ = pdg.get_properties_from_pdg_id(part.pdg_id[part.particle_id >= num_part])
N = A - Z

Zmin, Zmax = Z.min(), Z.max()
Nmin, Nmax = N.min(), N.max()

# +1 because bins define edges
H, Nedges, Zedges = np.histogram2d(
    N, Z,
    bins=(np.arange(Nmin, Nmax+2), np.arange(Zmin, Zmax+2))
)

fig, ax = plt.subplots(figsize=(10, 7))
pcm = ax.pcolormesh(
    Nedges, Zedges, H.T,
    norm=LogNorm(vmin=1, vmax=H.max()),
    cmap='viridis'
)
cbar = fig.colorbar(pcm, ax=ax)
cbar.set_label("Counts (log scale)")
ax.set_xticks(np.arange(Nedges[0], Nedges[-1] + 3, 5))
ax.set_yticks(np.arange(Zedges[0], Zedges[-1] + 3, 5))
ax.grid(which='both', linestyle=':', linewidth=0.5, color='grey')
ax.set_aspect('equal', adjustable='box')
ax.set_xlabel("Neutron number N")
ax.set_ylabel("Proton number Z")
ax.set_title("Isotope map")
plt.tight_layout()
plt.show()

