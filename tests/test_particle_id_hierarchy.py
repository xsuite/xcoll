# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

import time
import pytest
import numpy as np

import xpart as xp
import xtrack as xt
import xcoll as xc

from xobjects.test_helpers import for_all_test_contexts
try:
    import rpyc
except ImportError as e:
    rpyc = None


path = xc._pkg_root.parent / 'tests' / 'data'
particle_ref = xt.Particles('proton', p0c=6.8e12)


@for_all_test_contexts(
    excluding=('ContextCupy', 'ContextPyopencl')  # Geant4 only on CPU
)
@pytest.mark.skipif(not xc.geant4.environment.compiled, reason="BDSIM+Geant4 installation not found")
@pytest.mark.skipif(rpyc is None, reason="rpyc not installed")
def test_geant4(test_context):
    length = 1.2
    num_slices = 2000
    coll = xc.Geant4Collimator(length=length/num_slices, jaw=0.001, material='Ti', _context=test_context)
    xc.geant4.engine.particle_ref = particle_ref
    xc.geant4.engine.start(elements=coll, relative_energy_cut=0.1,
                           bdsim_config_file=path / f'geant4_protons.gmad')

    step = 2e-5
    x = np.arange(0.002, 0.003+step, step)
    y = np.arange(-0.001, 0.001+step, step)
    X, Y = np.meshgrid(x,y)
    coords = np.vstack([X.ravel(), Y.ravel()]).T
    _capacity = 2*len(coords[:,0])
    part_init = xp.build_particles(x=coords[:,0], y=coords[:,1], particle_ref=particle_ref,
                                   _capacity=_capacity, _context=test_context)
    parents = part_init.particle_id[part_init.state==1]


    part = part_init.copy()
    t_start = time.time()
    for i in range(num_slices):
        if i % 100 == 0:
            print(f"Slice {i}/{num_slices}")
        coll.track(part)
    print(f"Time per track: {(time.time()-t_start)/num_slices*1e3:.2f}ms for "
        + f"{len(part_init.x)} protons through {coll.length/1000:.2f}mm")
    xc.geant4.engine.stop(clean=True)

    mask_child = part.particle_id > parents.max()
    child_id   = part.particle_id[mask_child]
    assert not np.any(part.parent_particle_id[mask_child] == child_id)
    assert np.all(part.parent_particle_id[~mask_child] == part.particle_id[~mask_child])
    print(f"{len(child_id)} children created from {len(parents)} primaries")

    ppart = xc.ParticlesTree(part)
    child_idx  = ppart.indices_of(child_id)     # vectorised lookup
    primary_id = ppart.root_ids[child_idx]      # roots per child, fast
    assert np.all(primary_id <= parents.max())
    assert np.all(np.isin(primary_id, parents))

    id_to_index = {pid_val: idx for idx, pid_val in enumerate(part.particle_id)}
    mask = np.array([id_to_index[p] for p in primary_id])  # Masks children to their primaries
    assert np.all(part_init.parent_particle_id[primary_id] == part.parent_particle_id[mask])

    dx = part_init.x[primary_id] - part.x[mask]
    dy = part_init.y[primary_id] - part.y[mask]
    distance = np.sqrt(dx*dx + dy*dy)  # Distance from primary to its children

    grid_sep = x[1]-x[0]
    assert np.all(distance <= 5*grid_sep)  # All children should be close
    print(f"Max distance between child and primary: {distance.max():.2e} m")
    print(f"RMS distance between child and primary: {distance.std():.2e} m")
    print(f"Average distance between child and primary: {distance.mean():.2e} m")
    print(f"Percentage of children within {grid_sep/2:.2e} m (stay within 1 "
          f"grid cell): { (distance <= grid_sep/2).sum() / len(distance) * 100:.2f}%")
    assert (distance <= grid_sep/2).sum() / len(distance) > 0.95  # Allow only 5% of children to leave the grid cell
