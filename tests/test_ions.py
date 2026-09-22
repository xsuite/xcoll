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
import xtrack.particles.masses as xpm
from  xcoll import constants as xcc

from xobjects.test_helpers import for_all_test_contexts

from _common_api import engine_params


path = Path(__file__).parent / 'data'
particle_ref = xt.Particles('Pu-239', p0c=94*7.0e12)


@pytest.mark.parametrize("engine", engine_params)
def test_ions(engine):
    num_part = 2000
    capacity = 100*num_part
    length = 0.05
    jaw = 0.001
    material = 'Ti'
    seed = 1336

    if engine == "fluka":
        xc_engine = xc.fluka.engine
        xc_engine.relative_length_fortran_array = 100
        coll = xc.FlukaCollimator(length=length, jaw=jaw, material=material)

    elif engine == "geant4":
        xc_engine = xc.geant4.engine
        coll = xc.Geant4Collimator(length=length, jaw=jaw, material=material)

    if xc_engine.is_running():
        xc_engine.stop(clean=True)

    xc_engine.particle_ref = particle_ref
    xc_engine.return_none = True
    xc_engine.return_ions = True
    xc_engine.start(elements=coll, seed=seed)
    part = xp.build_particles(x=np.ones(num_part)*2*coll.jaw_L,
                                particle_ref=xc_engine.particle_ref,
                                _capacity=capacity)

    t_start = time.time()
    coll.track(part)
    xc_engine.stop(clean=True)

    print(f"Time per track: {(time.time()-t_start)*1e3:.2f}ms for "
        + f"{num_part} Pu-239 ions through {coll.length:.2f}m")

    mask_children = part.particle_id >= num_part
    pdg_ids = part.pdg_id[mask_children]
    print(f"Generated {len(pdg_ids)} children.")

    # We explicitly requested ions only
    assert len(pdg_ids) > 0
    assert np.all(pdg_ids > 1_000_000_000)
    _, A, Z, _ = pdg.get_properties_from_pdg_id(pdg_ids)
    assert np.all(A >= 2)
    assert np.all(Z >= 1)
    assert np.all(A >= Z)

    # Check masses
    A_from_mass = np.rint(part.mass[mask_children] / xpm.U_MASS_EV).astype(int)
    assert np.array_equal(A_from_mass, A)

    # Require a broad fragmentation spectrum
    assert len(np.unique(A)) > 100
    assert len(np.unique(Z)) > 40

    A_core = np.unique(A[(A >= 20) & (A < 100)])
    coverage = len(A_core) / 80
    print(f"A coverage in [20, 100): {coverage:.1%}")
    assert coverage > 0.9  # FLUKA 100% Geant4 95%
    assert np.max(np.diff(A_core)) <= 4
    assert set(A_core % 10) == set(range(10)) # Test that the last digit of A is not only 0 (was a bug)

    # Check iostopes coverage in [100, 240[
    A_tail = np.unique(A[(A >= 100) & (A < 240)])
    coverage = len(A_tail) / (np.max(A_tail) - 100 + 1) if len(A_tail) > 0 else 0
    print(f"A coverage in [100, {np.max(A_tail) if len(A_tail) > 0 else 100}): {coverage:.1%}")
    assert coverage > 0.5  # FLUKA 63% Geant4 64%

    # print info
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
