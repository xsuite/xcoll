import numpy as np
import xpart as xp
import xpart.pdg as pdg
import xcoll as xc
import time

import matplotlib.pyplot as plt

num_part = int(100)
_capacity = num_part*1000

if xc.FlukaEngine.is_running():
    xc.FlukaEngine.stop()


def run_many_particles(particle_ref):
    # Create a FLUKA collimator
    coll = xc.FlukaCollimator(length=0.4, assembly='fcc_tcp')
    coll.jaw = 0.001


    # Connect to FLUKA
    xc.FlukaEngine.particle_ref = particle_ref
    xc.FlukaEngine.capacity = _capacity
    xc.FlukaEngine.start(elements=coll, clean=True, verbose=False, return_all=True,
                        return_neutral=True, electron_lower_momentum_cut=1.e9)


    # Create an initial distribution of particles, random in 4D, on the left jaw (with the
    # longitudinal coordinates set to zero)
    x_init   = np.random.normal(loc=0.0012, scale=0.2e-3, size=num_part)
    px_init  = np.random.normal(loc=-1.e-5, scale=5.e-6, size=num_part)
    y_init   = np.random.normal(loc=0., scale=1e-3, size=num_part)
    py_init  = np.random.normal(loc=0., scale=5.e-6, size=num_part)
    part_init = xp.build_particles(x=x_init, px=px_init, y=y_init, py=py_init,
                                particle_ref=xc.FlukaEngine.particle_ref,
                                _capacity=xc.FlukaEngine.capacity)
    part = part_init.copy()

    # Do the tracking in FLUKA
    print(f"Tracking {num_part} {pdg.get_name_from_pdg_id(particle_ref.pdg_id[0])}s  (FLUKA)...     ", end='')
    start = time.time()
    coll.track(part)
    print(f"Done in {round(time.time()-start, 3)}s.")
    # print(f"Incoming particles survived: {len(part.state[(part.state>0) & (part.particle_id<10)])}/{num_part}")
    pdg_ids = np.unique(part.pdg_id[part.particle_id>=10], return_counts=True)
    print("Returned particles:")
    for pdg_id, num in zip(*pdg_ids):
        try:
            name = pdg.get_name_from_pdg_id(pdg_id)
        except ValueError:
            name = 'unknown'
        print(f"  {num:6} {name:14}  (PDG ID: {pdg_id}, mass: {part.mass[part.pdg_id==pdg_id][0]} eV)")
    print()

    # Stop the FLUKA server
    xc.FlukaEngine.stop(clean=True)


run_many_particles(xp.Particles.reference_from_pdg_id(pdg_id='proton', p0c=6.8e12))
run_many_particles(xp.Particles.reference_from_pdg_id(pdg_id='Pb208', p0c=6.8e12*82))
run_many_particles(xp.Particles.reference_from_pdg_id(pdg_id='electron', p0c=200e9))
