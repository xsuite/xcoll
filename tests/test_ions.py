# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

import numpy as np
import pytest
import time

import xtrack as xt
import xpart as xp
import xcoll as xc

from xobjects.test_helpers import for_all_test_contexts

# try the import here and skip tests if missing
# also need the import here in case of pytest --forked
try:
    import collimasim as cs
except ImportError:
    cs = None

path = xc._pkg_root.parent / 'tests' / 'data'
particle_ref = xt.Particles.reference_from_pdg_id(pdg_id='Pu-239', p0c=7.0e12)

@for_all_test_contexts(
    excluding=('ContextCupy', 'ContextPyopencl')  # Geant4 only on CPU
)
@pytest.mark.skipif(cs is None, reason="Geant4 tests need collimasim installed")
def test_bdsim_ions(test_context):
    num_part = 100
    coll = xc.Geant4Collimator(length=0.01, jaw=0.001, material='Ti', _context=test_context)
    xc.geant4.engine.particle_ref = particle_ref
    part = xp.build_particles(x=np.ones(num_part)*2*coll.jaw_L,
                                   particle_ref=particle_ref, _capacity=10*num_part)
    xc.geant4.engine.start(elements=coll, seed=1336,
                           bdsim_config_file=path / 'geant4_osmium.gmad')
    coll.track(part)
    xc.geant4.engine.stop()
    pdg_ids = np.unique(part.pdg_id)
    assert np.sum((pdg_ids > 1000000000) & (pdg_ids < 1000942390)) > 0
