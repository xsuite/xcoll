import numpy as np
import xpart as xp
import xtrack as xt
import xcoll as xc
import time


if xc.geant4.engine.is_running():
    xc.geant4.engine.stop()


# Make some Geant4 and Everest collimators to preload materials
g4_colls = []
ev_colls = []
for material in ['Air', 'Water', 'Beryllium', 'MoGR', 'Inermet180']:
    for length in [1.e-3, 2.e-2, 0.1, 0.4, 1]:
        g4_colls.append(xc.Geant4Collimator(length=length, material=material, jaw=0.001))
        ev_colls.append(xc.EverestCollimator(length=length, material=material, jaw=0.001))


# Start the engine
xc.geant4.engine.particle_ref = xt.Particles('proton', p0c=6.8e12)
xc.geant4.engine.start(elements=g4_colls, clean=True, verbose=False)


# Create an initial distribution of particles
num_part   = 5000
x_init     = np.random.normal(loc=0.05, scale=0.2e-3, size=num_part)
px_init    = np.random.normal(loc=0.,   scale=5.e-6,  size=num_part)
y_init     = np.random.normal(loc=0.,   scale=1e-3,   size=num_part)
py_init    = np.random.normal(loc=0.,   scale=5.e-6,  size=num_part)
zeta_init  = np.random.normal(loc=0.,   scale=1.e-2,  size=num_part)
delta_init = np.random.normal(loc=0.,   scale=5.e-4,  size=num_part)
part_init = xp.build_particles(x=x_init, px=px_init, y=y_init, py=py_init,
                               zeta=zeta_init, delta=delta_init,
                               particle_ref=xc.geant4.engine.particle_ref,
                               _capacity=num_part*2)

for ev_coll, g4_coll in zip(ev_colls, g4_colls):
    print(f"\nBenchmarking collimator: length={g4_coll.length} m, material={g4_coll.material.name}")
    # Copy initial particles for each collimator
    part = part_init.copy()
    print(f"    Tracking {num_part} protons with Everest...  ", end='', flush=True)
    start = time.time()
    ev_coll.track(part)
    print(f"Done in {(time.time()-start)*1000:.3}ms. Dead: {np.sum((part.state < 1) & (part.state > -999999))}. "
          f"Alive: {np.sum(part.state == 1)}", flush=True)
    part = part_init.copy()
    print(f"    Tracking {num_part} protons with Geant4...   ", end='', flush=True)
    start = time.time()
    g4_coll.track(part)
    print(f"Done in {time.time()-start:.3}s. Dead: {np.sum((part.state < 1) & (part.state > -999999))}. "
          f"Alive: {np.sum(part.state == 1)}", flush=True)
