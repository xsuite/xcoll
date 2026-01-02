# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

import time
import pytest
import numpy as np
from pathlib import Path

import xtrack as xt
import xpart as xp
import xcoll as xc
import xtrack.particles.pdg as pdg
from  xcoll import constants as xcc

from xobjects.test_helpers import for_all_test_contexts
try:
    import rpyc
except ImportError as e:
    rpyc = None


path = Path(__file__).parent / 'data'
particle_ref = xt.Particles('Pu-239', p0c=94*7.0e12)


@for_all_test_contexts(
    excluding=('ContextCupy', 'ContextPyopencl')  # Geant4 only on CPU
)
@pytest.mark.skipif(rpyc is None, reason="rpyc not installed")
@pytest.mark.skipif(not xc.geant4.environment.ready, reason="BDSIM+Geant4 installation not found")
def test_geant4(test_context):
    num_part = 500
    coll = xc.Geant4Collimator(length=0.01, jaw=0.001, material='Ti', _context=test_context)
    xc.geant4.engine.particle_ref = particle_ref
    part = xp.build_particles(x=np.ones(num_part)*2*coll.jaw_L,
                                   particle_ref=particle_ref, _capacity=10*num_part)
    xc.geant4.engine.start(elements=coll, seed=1336)

    t_start = time.time()
    coll.track(part)
    print(f"Time per track: {(time.time()-t_start)*1e3:.2f}ms for "
        + f"{num_part} Pu-239 ions through {coll.length:.2f}m")
    xc.geant4.engine.stop(clean=True)

    # Get only the initial particles that survived and all new particles (even if dead, as neutral particles will be flagged dead)
    mask = (part.state > 0) | (part.particle_id >= num_part)
    pdg_ids, counts = np.unique(part.pdg_id[mask], return_counts=True)
    for pdg_id, num in zip(pdg_ids, counts):
        try:
            name = pdg.get_name_from_pdg_id(pdg_id, long_name=False)
        except ValueError:
            name = 'unknown'
        if part.state[part.pdg_id==pdg_id][0] == xcc.MASSLESS_OR_NEUTRAL:
            mass = 0
        else:
            mass = part.mass[part.pdg_id==pdg_id][0]
        E = part.energy[mask & (part.pdg_id==pdg_id)]
        en = f"{E[~np.isnan(E)].mean():.1e} ± {E[~np.isnan(E)].std():.1e} eV"
        print(f"  {num:6} {name:12}{en:21}  (PDG ID: {pdg_id}, mass: {mass} eV)")

    # Check that we have some fission products (PDG ID > 1 billion
    assert np.sum((pdg_ids > 1000000000) & (pdg_ids < 1000942390)) > 0


@pytest.mark.skipif(not xc.fluka.environment.ready, reason="FLUKA installation not found")
def test_fluka():
    num_part = 500
    coll = xc.FlukaCollimator(length=0.01, jaw=0.001, material='Ti')
    xc.fluka.engine.particle_ref = particle_ref
    part = xp.build_particles(x=np.ones(num_part)*2*coll.jaw_L,
                                   particle_ref=particle_ref, _capacity=10*num_part)
    xc.fluka.engine.start(elements=coll, seed=1336)

    t_start = time.time()
    coll.track(part)
    print(f"Time per track: {(time.time()-t_start)*1e3:.2f}ms for "
        + f"{num_part} Pu-239 ions through {coll.length:.2f}m")
    xc.fluka.engine.stop(clean=True)

    # Get only the initial particles that survived and all new particles (even if dead, as neutral particles will be flagged dead)
    mask = (part.state > 0) | (part.particle_id >= num_part)
    pdg_ids, counts = np.unique(part.pdg_id[mask], return_counts=True)
    for pdg_id, num in zip(pdg_ids, counts):
        try:
            name = pdg.get_name_from_pdg_id(pdg_id, long_name=False)
        except ValueError:
            name = 'unknown'
        if part.state[part.pdg_id==pdg_id][0] == xcc.MASSLESS_OR_NEUTRAL:
            mass = 0
        else:
            mass = part.mass[part.pdg_id==pdg_id][0]
        E = part.energy[mask & (part.pdg_id==pdg_id)]
        en = f"{E[~np.isnan(E)].mean():.1e} ± {E[~np.isnan(E)].std():.1e} eV"
        print(f"  {num:6} {name:12}{en:21}  (PDG ID: {pdg_id}, mass: {mass} eV)")

    # Check that we have some fission products (PDG ID > 1 billion
    assert np.sum((pdg_ids > 1000000000) & (pdg_ids < 1000942390)) > 0
