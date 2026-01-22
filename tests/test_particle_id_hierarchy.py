# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

import time
import pytest
import numpy as np
from pathlib import Path

import xpart as xp
import xtrack as xt
import xcoll as xc
from xpart.test_helpers import flaky_assertions, retry

from _common_api import engine_params


path = Path(__file__).parent / 'data'
particle_ref = xt.Particles('proton', p0c=6.8e12)


@pytest.mark.parametrize("engine", engine_params)
@retry()
def test_hierarchy(engine):
    if engine == "fluka":
        if xc.fluka.engine.is_running():
            xc.fluka.engine.stop(clean=True)
    elif engine == "geant4":
        if xc.geant4.engine.is_running():
            xc.geant4.engine.stop(clean=True)
    length = 1.2
    num_slices = 2000
    jaw = 0.001
    material = 'Ti'
    seed = 1336
    step = 2e-5

    x = np.arange(jaw + 0.001, jaw + 0.003 + step, step)
    y = np.arange(-0.001, 0.001+step, step)
    X, Y = np.meshgrid(x,y)
    coords = np.vstack([X.ravel(), Y.ravel()]).T
    capacity = 50*len(coords[:,0])

    if engine == "fluka":
        coll = xc.FlukaCollimator(length=length/num_slices, jaw=jaw, material=material)
        xc.fluka.engine.particle_ref = particle_ref
        xc.fluka.engine.start(elements=coll, seed=seed)
        xc.fluka.engine.return_all = True
        xc.fluka.engine.capacity = capacity
        part_init = xp.build_particles(x=coords[:,0], y=coords[:,1],
                                       particle_ref=xc.fluka.engine.particle_ref,
                                       _capacity=xc.fluka.engine.capacity)

    elif engine == "geant4":
        coll = xc.Geant4Collimator(length=length/num_slices, jaw=jaw, material=material)
        xc.geant4.engine.particle_ref = particle_ref
        xc.geant4.engine.start(elements=coll, seed=seed)
        xc.geant4.engine.return_all = True
        part_init = xp.build_particles(x=coords[:,0], y=coords[:,1],
                                       particle_ref=xc.geant4.engine.particle_ref,
                                       _capacity=capacity)

    parents = part_init.particle_id[part_init.state==1]
    part = part_init.copy()

    t_start = time.time()
    for i in range(num_slices):
        if i % 100 == 0:
            print(f"Slice {i}/{num_slices}")
        coll.track(part)
    print(f"Time per track: {(time.time()-t_start)/num_slices*1e3:.2f}ms for "
        + f"{len(part_init.x)} protons through {coll.length/1000:.2f}mm")

    mask_child = part.particle_id > parents.max()
    child_id   = part.particle_id[mask_child]

    with flaky_assertions():
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
        print(f"Percentage of children within {grid_sep/2:.2e} m (1 grid cell): "
            f"{(distance <= grid_sep/2).sum() / len(distance) * 100:.2f}%")
        print(' '*23 + f"within {grid_sep:.2e} m  (2 grid cells): "
            f"{(distance <= grid_sep).sum() / len(distance) * 100:.2f}%")
        print(' '*23 + f"within {1.5*grid_sep:.2e} m  (3 grid cells): "
            f"{(distance <= 1.5*grid_sep).sum() / len(distance) * 100:.2f}%")
        print(' '*23 + f"within {2*grid_sep:.2e} m  (4 grid cells): "
            f"{(distance <= 2*grid_sep).sum() / len(distance) * 100:.2f}%")
        print(' '*23 + f"within {2.5*grid_sep:.2e} m  (5 grid cells): "
            f"{(distance <= 2.5*grid_sep).sum() / len(distance) * 100:.2f}%")
        assert (distance <= grid_sep/2).sum() / len(distance) > 0.90    # Allow only 10% of children to leave the grid cell
        assert (distance <= grid_sep).sum() / len(distance) > 0.95      # Allow only 5% of children to leave two grid cells
        assert (distance <= 1.5*grid_sep).sum() / len(distance) > 0.98  # Allow only 2% of children to leave three grid cells
        assert (distance <= 2*grid_sep).sum() / len(distance) > 0.99  # Allow only 1% of children to leave four grid cells

    if engine == "fluka":
        xc.fluka.engine.stop(clean=True)
    elif engine == "geant4":
        xc.geant4.engine.stop(clean=True)

# TODO
# @pytest.mark.fluka
# def test_hierarchy_donadon():
#     if xc.fluka.engine.is_running():
#         xc.fluka.engine.stop(clean=True)
#     raise NotImplementedError("Need to implement test for child particle ids with DONADON collimator")
#     coll = xc.FlukaCollimator(length=0.0001, assembly='donadon', jaw=0)
#     xc.fluka.engine.particle_ref = xt.Particles.reference_from_pdg_id(pdg_id='electron', p0c=200e9)
#     xc.fluka.engine.capacity = 100_000
#     xc.fluka.engine.seed = 7856231
#     xc.fluka.engine.start(elements=coll, clean=False, verbose=True)
#     x_init, y_init = np.array(np.meshgrid(np.linspace(-0.01, 0.01, 21), np.linspace(-0.01, 0.01, 21))).reshape(2,-1)
#     px_init = 0
#     py_init = 0
#     part_init = xp.build_particles(x=x_init, px=px_init, y=y_init, py=py_init,
#                                    particle_ref=xc.fluka.engine.particle_ref,
#                                    _capacity=xc.fluka.engine.capacity)

#     part = part_init.copy()
#     xc.fluka.engine.stop(clean=True)
