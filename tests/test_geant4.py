import json
from pathlib import Path
import numpy as np
import pytest

import xobjects as xo
import xpart as xp
import xcoll as xc
from xobjects.test_helpers import for_all_test_contexts

# try the import here and skip tests if missing
# also need the import here in case of pytest --forked
try:
    import collimasim as cs
except ImportError:
    cs = None

path = Path.cwd() / 'data_test_geant4'

@for_all_test_contexts(
    excluding=('ContextCupy', 'ContextPyopencl')  # Geant4 only on CPU
)
@pytest.mark.skipif(cs is None, reason="Geant4 tests need collimasim installed")
def test_multiple_tracking(test_context):
    xc.Geant4Engine(random_generator_seed=1993, 
                    reference_pdg_id=2212, 
                    reference_kinetic_energy=7e12, 
                    relative_energy_cut=0.15, 
                    bdsim_config_file=str(path / f'settings_black_absorber_protons.gmad'))

    g4_collimators = _make_geant4_collimators(_context=test_context)
    ba_collimators = _make_black_absorbers(_context=test_context)

    part = _generate_particles(_context=test_context)

    part_ba = part.copy()

    for coll in g4_collimators:
        coll.track(part)

    for coll in ba_collimators:
        coll.track(part_ba)

    part.sort(interleave_lost_particles=True)
    part_ba.sort(interleave_lost_particles=True)

    assert np.all(part.filter(part.state==1).particle_id 
                  == part_ba.filter(part_ba.state==1).particle_id)


def _make_geant4_collimators(angles=[0,45,90], tilts=[0,0], _context=None):
    if _context is None:
        _context = xo.ContextCpu()
    jaws = [0.03, -0.02]
    co = [-0.01, 0.01]
    L = 0.873

    collimators = []
    for ii, angle in enumerate(angles):
        shift = co[0]*np.cos(angle) + co[1]*np.sin(angle)
        g4coll = xc.Geant4Collimator(length=L, angle=angle,
                                     jaw=jaws+shift, tilt=tilts,
                                     _context=_context, material='cu',
                                     collimator_id=f'g4coll_{ii}')
        collimators.append(g4coll)

    return collimators

def _make_black_absorbers(angles=[0,45,90], tilts=[0,0], _context=None):
    if _context is None:
        _context = xo.ContextCpu()
    jaws = [0.03, -0.02]
    co = [-0.01, 0.01]
    L = 0.873

    collimators = []
    for ii, angle in enumerate(angles):
        shift = co[0]*np.cos(angle) + co[1]*np.sin(angle)
        bacoll = xc.BlackAbsorber(length=L, angle=angle, 
                                     jaw=jaws+shift, tilt=tilts,
                                     _context=_context)
        collimators.append(bacoll)

    return collimators


def _generate_particles(four_dim=True, angle=0, _context=None):
    if _context is None:
        _context = xo.ContextCpu()
    # Make particles
    n_part = 10000
    x = np.random.uniform(-0.1, 0.1, n_part)
    y = np.random.uniform(-0.1, 0.1, n_part)
    if four_dim:
        px = np.random.uniform(-0.1, 0.1, n_part)
        py = np.random.uniform(-0.1, 0.1, n_part)
    else:
        px = 0
        py = 0

    ref = xp.Particles(mass0=xp.PROTON_MASS_EV, q0=1, p0c=7e12, _context=_context)
    part = xp.build_particles(x=x, y=y, px=px, py=py, particle_ref=ref, _context=_context)

    return part
