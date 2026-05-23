import time
import numpy as np
import pytest

import xtrack as xt
import xtrack.particles.pdg as pdg
from xtrack.particles import LAST_INVALID_STATE
from xtrack.particles import masses as xpm
from xpart.test_helpers import flaky_assertions, retry
import xcoll as xc
from  xcoll import constants as xcc

from _common_api import check_skip_old_bdsim, engine_params, engine_params


# # Treat warnings as errors to debug
# np.set_printoptions(threshold=np.inf)
# import warnings
# warnings.filterwarnings("error")


@pytest.mark.parametrize("engine", engine_params)
@retry()
def test_return_photons(engine):
    check_skip_old_bdsim(engine)
    print("Testing return_none=True, return_photons=True")
    num_part = 5000
    capacity = 100_000
    particle_ref = xt.Particles('proton', p0c=6.8e12)
    part = _run(engine, num_part, capacity, particle_ref, True, do_assert=False, return_type='photons')
    pdg = part.pdg_id[part.particle_id >= num_part]
    assert np.all(pdg == 22)
    with flaky_assertions():
       assert pdg.size > 0  # Some photons should be created to make this test meaningful


# TODO: does not work for FLUKA
@pytest.mark.geant4
@retry()
def test_return_electrons():
    engine = "geant4"
    check_skip_old_bdsim(engine)
    print("Testing return_none=True, return_electrons=True")
    num_part = 2500
    capacity = 250_000
    particle_ref = xt.Particles('proton', p0c=6.8e12)
    part = _run(engine, num_part, capacity, particle_ref, True, do_assert=False, return_type='electrons')
    pdg = part.pdg_id[part.particle_id >= num_part]
    assert np.all((pdg == 11) | (pdg == -11))
    with flaky_assertions():
        assert pdg.size > 0  # Some electrons should be created to make this test meaningful


@pytest.mark.parametrize("engine", engine_params)
@retry()
def test_return_protons(engine):
    check_skip_old_bdsim(engine)
    print("Testing return_none=True, return_protons=True")
    num_part = 2500
    capacity = 50_000
    particle_ref = xt.Particles('proton', p0c=6.8e12)
    part = _run(engine, num_part, capacity, particle_ref, True, do_assert=False, return_type='protons')
    pdg = part.pdg_id[part.particle_id >= num_part]
    assert np.all((pdg == 2212) | (pdg == -2212))
    with flaky_assertions():
        assert pdg.size > 0  # Some protons should be created to make this test meaningful


@pytest.mark.parametrize("engine", engine_params)
@retry()
def test_return_neutrons(engine):
    check_skip_old_bdsim(engine)
    print("Testing return_none=True, return_neutrons=True")
    num_part = 2500
    capacity = 50_000
    particle_ref = xt.Particles('proton', p0c=6.8e12)
    part = _run(engine, num_part, capacity, particle_ref, True, do_assert=False, return_type='neutrons')
    pdg = part.pdg_id[part.particle_id >= num_part]
    assert np.all((pdg == 2112) | (pdg == -2112))
    with flaky_assertions():
        assert pdg.size > 0  # Some neutrons should be created to make this test meaningful


@pytest.mark.parametrize("engine", engine_params)
@retry()
def test_return_mesons(engine):
    check_skip_old_bdsim(engine)
    print("Testing return_none=True, return_mesons=True")
    num_part = 2500
    capacity = 50_000
    particle_ref = xt.Particles('proton', p0c=6.8e12)
    part = _run(engine, num_part, capacity, particle_ref, True, do_assert=False, return_type='mesons')
    pdg = part.pdg_id[part.particle_id >= num_part]
    assert not np.any(pdg == 22)
    assert not np.any((pdg == 11) | (pdg == -11))
    assert not np.any((pdg == 2212) | (pdg == -2212))
    assert not np.any((pdg == 2112) | (pdg == -2112))
    assert np.all(np.abs(pdg) // 10 % 10 != 0)
    assert np.all(np.abs(pdg) // 100 % 10 != 0)
    assert np.all(np.abs(pdg) // 1000 % 10 == 0)
    with flaky_assertions():
        assert pdg.size > 0  # Some mesons should be created to make this test meaningful


# TODO: does not work for Geant4
@pytest.mark.fluka
@retry()
def test_return_baryons():
    engine = "fluka"
    check_skip_old_bdsim(engine)
    print("Testing return_none=True, return_other_baryons=True")
    num_part = 100
    capacity = 50_000
    particle_ref = xt.Particles('He4', p0c=450e9*4)
    part = _run(engine, num_part, capacity, particle_ref, True, do_assert=False, return_type='other_baryons')
    pdg = part.pdg_id[part.particle_id >= num_part]
    assert not np.any(pdg == 22)
    assert not np.any((pdg == 11) | (pdg == -11))
    assert not np.any((pdg == 2212) | (pdg == -2212))
    assert not np.any((pdg == 2112) | (pdg == -2112))
    assert np.all(np.abs(pdg) // 10 % 10 != 0)
    assert np.all(np.abs(pdg) // 100 % 10 != 0)
    assert np.all(np.abs(pdg) // 1000 % 10 != 0)
    with flaky_assertions():
        assert pdg.size > 0  # Some baryons should be created to make this test meaningful


@pytest.mark.parametrize("engine", engine_params)
@pytest.mark.parametrize('hit', [True, False], ids=['hit', 'miss'])
@retry()
def test_protons(engine, hit):
    print(f"Testing protons in {engine.capitalize()} with hit={hit}.")
    _run(engine, 500, 2500, xt.Particles('proton', p0c=6.8e12), hit)


@pytest.mark.parametrize("engine", engine_params)
@pytest.mark.parametrize('hit', [True, False], ids=['hit', 'miss'])
@retry()
def test_lead(engine, hit):
    print(f"Testing lead in {engine.capitalize()} with hit={hit}.")
    _run(engine, 100, 100_000, xt.Particles('Pb208', p0c=6.8e12*82), hit)


@pytest.mark.parametrize("engine", engine_params)
@pytest.mark.parametrize('hit', [True, False], ids=['hit', 'miss'])
@pytest.mark.parametrize('proton_ref', [True, False], ids=['proton_ref', 'antiproton_ref'])
@retry()
def test_antiprotons(engine, proton_ref, hit):
    print(f"Testing antiprotons in {engine.capitalize()} with hit={hit} and proton_ref={proton_ref}")
    p0c = 6.8e12
    if proton_ref:
        _run(engine, 100, 500, xt.Particles('proton', p0c=p0c), hit, tol=3e-11,
             mass_ratio=1, charge_ratio=-1, pdg_id=-2212)
    else:
        _run(engine, 100, 500, xt.Particles('antiproton', p0c=p0c), hit)


@pytest.mark.parametrize("engine", engine_params)
@pytest.mark.parametrize('hit', [True, False], ids=['hit', 'miss'])
@pytest.mark.parametrize('proton_ref', [True, False], ids=['proton_ref', 'electron_ref'])
@retry()
def test_electrons(engine, proton_ref, hit):
    check_skip_old_bdsim(engine, check_old_bdsim=proton_ref)
    print(f"Testing electrons in {engine.capitalize()} with hit={hit} and proton_ref={proton_ref}")
    p0c = 200e9
    if proton_ref:
        pdg_id = 11
        if engine == 'fluka':
            ref_mass = xc.fluka.particle_masses[pdg_id] or xpm.ELECTRON_MASS_EV
        elif engine == 'geant4':
            ref_mass = xc.geant4.particle_masses[pdg_id] or xpm.ELECTRON_MASS_EV
        _run(engine, 500, 50_000, xt.Particles('proton', p0c=p0c), hit, tol=3e-11,
             ref_mass=ref_mass, charge_ratio=-1, pdg_id=pdg_id)
    else:
        _run(engine, 500, 50_000, xt.Particles('electron', p0c=p0c), hit)


@pytest.mark.parametrize("engine", engine_params)
@pytest.mark.parametrize('hit', [True, False], ids=['hit', 'miss'])
@pytest.mark.parametrize('proton_ref', [True, False], ids=['proton_ref', 'positron_ref'])
@retry()
def test_positrons(engine, proton_ref, hit):
    check_skip_old_bdsim(engine, check_old_bdsim=proton_ref)
    print(f"Testing positrons in {engine.capitalize()} with hit={hit} and proton_ref={proton_ref}")
    p0c = 200e9
    if proton_ref:
        pdg_id = -11
        if engine == 'fluka':
            ref_mass = xc.fluka.particle_masses[pdg_id] or xpm.ELECTRON_MASS_EV
        elif engine == 'geant4':
            ref_mass = xc.geant4.particle_masses[pdg_id] or xpm.ELECTRON_MASS_EV
        _run(engine, 500, 50_000, xt.Particles('proton', p0c=p0c), hit, tol=3e-11,
             ref_mass=ref_mass, charge_ratio=1, pdg_id=pdg_id)
    else:
        _run(engine, 500, 50_000, xt.Particles('positron', p0c=p0c), hit)


@pytest.mark.parametrize("engine", engine_params)
@pytest.mark.parametrize('hit', [True, False], ids=['hit', 'miss'])
@pytest.mark.parametrize('proton_ref', [True, False], ids=['proton_ref', 'muon_ref'])
@retry()
def test_muons(engine, proton_ref, hit):
    check_skip_old_bdsim(engine)
    print(f"Testing muons in {engine.capitalize()} with hit={hit} and proton_ref={proton_ref}")
    p0c = 200e9
    if proton_ref:
        pdg_id = 13
        if engine == 'fluka':
            ref_mass = xc.fluka.particle_masses[pdg_id] or xpm.MUON_MASS_EV
        elif engine == 'geant4':
            ref_mass = xc.geant4.particle_masses[pdg_id] or xpm.MUON_MASS_EV
        _run(engine, 500, 50_000, xt.Particles('proton', p0c=p0c), hit, tol=3e-11,
             ref_mass=ref_mass, charge_ratio=-1, pdg_id=pdg_id)
    else:
        _run(engine, 500, 50_000, xt.Particles('muon', p0c=p0c), hit)


@pytest.mark.parametrize("engine", engine_params)
@pytest.mark.parametrize('hit', [True, False], ids=['hit', 'miss'])
@pytest.mark.parametrize('proton_ref', [True, False], ids=['proton_ref', 'antimuon_ref'])
@retry()
def test_antimuons(engine, proton_ref, hit):
    check_skip_old_bdsim(engine)
    print(f"Testing antimuons in {engine.capitalize()} with hit={hit} and proton_ref={proton_ref}")
    p0c = 200e9
    if proton_ref:
        pdg_id = -13
        if engine == 'fluka':
            ref_mass = xc.fluka.particle_masses[pdg_id] or xpm.MUON_MASS_EV
        elif engine == 'geant4':
            ref_mass = xc.geant4.particle_masses[pdg_id] or xpm.MUON_MASS_EV
        _run(engine, 500, 50_000, xt.Particles('proton', p0c=p0c), hit, tol=3e-11,
             ref_mass=ref_mass, charge_ratio=1, pdg_id=pdg_id)
    else:
        _run(engine, 500, 50_000, xt.Particles('antimuon', p0c=p0c), hit)


@pytest.mark.parametrize("engine", engine_params)
@pytest.mark.parametrize('hit', [True, False], ids=['hit', 'miss'])
@pytest.mark.parametrize('proton_ref', [True, False], ids=['proton_ref', 'pion_ref'])
@retry()
def test_positive_pions(engine, proton_ref, hit):
    check_skip_old_bdsim(engine)
    print(f"Testing positive pions in {engine.capitalize()} with hit={hit} and proton_ref={proton_ref}")
    p0c = 450e9
    if proton_ref:
        pdg_id = 211
        if engine == 'fluka':
            ref_mass = xc.fluka.particle_masses[pdg_id] or xpm.PION_MASS_EV
        elif engine == 'geant4':
            ref_mass = xc.geant4.particle_masses[pdg_id] or xpm.PION_MASS_EV
        _run(engine, 500, 50_000, xt.Particles('proton', p0c=p0c), hit, tol=3e-11,
             ref_mass=ref_mass, charge_ratio=1, pdg_id=pdg_id)
    else:
        _run(engine, 500, 50_000, xt.Particles('pi+', p0c=p0c), hit)


@pytest.mark.parametrize("engine", engine_params)
@pytest.mark.parametrize('hit', [True, False], ids=['hit', 'miss'])
@pytest.mark.parametrize('proton_ref', [True, False], ids=['proton_ref', 'pion_ref'])
@retry()
def test_negative_pions(engine, proton_ref, hit):
    check_skip_old_bdsim(engine)
    print(f"Testing negative pions in {engine.capitalize()} with hit={hit} and proton_ref={proton_ref}")
    p0c = 450e9
    if proton_ref:
        pdg_id = -211
        if engine == 'fluka':
            ref_mass = xc.fluka.particle_masses[pdg_id] or xpm.PION_MASS_EV
        elif engine == 'geant4':
            ref_mass = xc.geant4.particle_masses[pdg_id] or xpm.PION_MASS_EV
        _run(engine, 500, 50_000, xt.Particles('proton', p0c=p0c), hit, tol=3e-11,
             ref_mass=ref_mass, charge_ratio=-1, pdg_id=pdg_id)
    else:
        _run(engine, 500, 50_000, xt.Particles('pi-', p0c=p0c), hit)


def _run(engine, num_part, capacity, particle_ref, hit, tol=1e-12, do_assert=True,
         return_type=None, ref_mass=None, **kwargs):
    if engine == "fluka":
        if xc.fluka.engine.is_running():
            xc.fluka.engine.stop(clean=True)
        coll = xc.FlukaCollimator(length=0.4, material='MoGr')
        coll.jaw = 0.002
        xc.fluka.engine.particle_ref = particle_ref
        xc.fluka.engine.capacity = capacity
        xc.fluka.engine.relative_capacity = 20
        if return_type is not None:
            xc.fluka.engine.return_none = True
            setattr(xc.fluka.engine, f'return_{return_type}', True)
        else:
            xc.fluka.engine.return_all = True
        xc.fluka.engine.start(elements=coll, clean=True, verbose=True)
        particle_ref = xc.fluka.engine.particle_ref

    elif engine == "geant4":
        if xc.geant4.engine.is_running():
            xc.geant4.engine.stop(clean=True)
        coll = xc.Geant4Collimator(length=0.4, material='MoGr')
        coll.jaw = 0.002
        xc.geant4.engine.particle_ref = particle_ref
        if return_type is not None:
            xc.geant4.engine.return_none = True
            setattr(xc.geant4.engine, f'return_{return_type}', True)
        else:
            xc.geant4.engine.return_all = True
        xc.geant4.engine.start(elements=coll, clean=True, verbose=True)
        particle_ref = xc.geant4.engine.particle_ref

    if hit:
        part, part_init = _init_particles(num_part, particle_ref=particle_ref, capacity=capacity, ref_mass=ref_mass, **kwargs)
    else:
        part, part_init = _init_particles(num_part, particle_ref=particle_ref, x=0, px=0, capacity=capacity, ref_mass=ref_mass, **kwargs)

    print(f"Tracking {part._num_active_particles} {pdg.get_name_from_pdg_id(part.pdg_id[0])}"
         + "s...     ", flush=True)
    start = time.time()
    if engine == "fluka":  xc.fluka.engine.physics_settings()
    if engine == "geant4": xc.geant4.engine.physics_settings()
    coll.track(part)
    print(f"Done in {round(time.time()-start, 3)}s.", flush=True)

    if do_assert:
        E_ref = part_init.energy[0]
        if hit:
            coll_state = xcc.LOST_ON_FLUKA_COLL if engine == 'fluka' else xcc.LOST_ON_GEANT4_COLL
            _assert_hit(part, part_init, E_ref, coll, coll_state=coll_state, tol=tol)
        else:
            if engine == "fluka":
                _assert_missed(part, part_init, E_ref, coll, tol=tol, allow_spurious_lost=2)
            else:
                _assert_missed(part, part_init, E_ref, coll, tol=tol)

    if engine == 'fluka' and xc.fluka.engine.is_running():
        xc.fluka.engine.stop(clean=True)
    elif engine == 'geant4' and xc.geant4.engine.is_running():
        xc.geant4.engine.stop(clean=True)
    print()
    return part


def _init_particles(num_part, particle_ref, x=0.004, px=-1.e-5, capacity=None, ref_mass=None, **kwargs):
    x_init  = np.random.normal(loc=x, scale=0.2e-3, size=num_part)
    px_init = np.random.normal(loc=px, scale=5.e-6, size=num_part)
    y_init  = np.random.normal(loc=0., scale=1e-3, size=num_part)
    py_init = np.random.normal(loc=0., scale=5.e-6, size=num_part)
    p0c = particle_ref.p0c
    mass0 = particle_ref.mass0
    if ref_mass is not None:
        mratio = ref_mass / mass0
        E0 = np.sqrt(p0c**2 + mass0**2)
        E  = np.sqrt(p0c**2 + ref_mass**2)
        kwargs['ptau'] = (E/mratio - E0)/p0c
        kwargs['mass_ratio'] = mratio
    kwargs.setdefault('p0c', p0c)
    kwargs.setdefault('mass0', mass0)
    kwargs.setdefault('q0', particle_ref.q0)
    kwargs.setdefault('pdg_id', particle_ref.pdg_id)
    part_init = xt.Particles(x=x_init, px=px_init, y=y_init, py=py_init,
                             _capacity=capacity, **kwargs)
    # print(f"E0={part_init.energy0[0]} eV  P0c={part_init.p0c[0]} eV/c  E={part_init.energy[0]} eV  P={(1+part_init.delta[0])*part_init.p0c[0]*part_init.mass_ratio[0]} eV/c")
    print("Created particles with the following parameters:")
    E0 = particle_ref.energy[0]
    print(f"Reference energy: {E0} eV (in particles: {part_init.energy[0]} eV)")
    print(f"Mass {part_init.mass[0]} eV")
    print(f"Charge {part_init.charge_ratio[0]*part_init.q0}")
    print(f"PDG ID {part_init.pdg_id[0]}")
    return part_init.copy(), part_init


def _assert_missed(part, part_init, E0, coll, tol=1e-12, allow_spurious_lost=0):
    alive = part_init._num_active_particles
    assert part._num_active_particles == alive
    state, counts = np.unique(part.state, return_counts=True)
    states = dict(zip(state, counts))
    print(f"Particles alive after missing: {part._num_active_particles}. State counts: {states}")
    assert 1 in states
    assert LAST_INVALID_STATE in states
    invalid = set(state) - {1, LAST_INVALID_STATE}
    assert np.sum([states[s] for s in invalid]) <= allow_spurious_lost  # Allow max N spurious lost particles
    mask_alive = False*np.ones(part.x.shape, dtype=bool)
    mask_alive[:alive] = True
    _test_drift(part, part_init, E0, mask_alive, coll, tol=tol)


def _assert_hit(part, part_init, E0, coll, coll_state, tol=1e-12):
    alive_before = part_init._num_active_particles
    alive_after = part._num_active_particles
    part.sort(interleave_lost_particles=True)
    is_child = np.array([pid!=ppid for pid, ppid in zip(part.particle_id, part.parent_particle_id)])
    has_children = np.array([pid in part.parent_particle_id[is_child] for pid in part.particle_id])
    print(f"Particles alive before: {alive_before}, after: {alive_after}, children generated: {np.sum(is_child)}")

    mask_hit = part_init.state != LAST_INVALID_STATE

    # All that are supposed to hit and are still alive should have lost energy
    # Except those (which actually did not hit but barely scratched the collimator):
    mask_hit_but_no = False*np.ones(mask_hit.shape, dtype=bool)
    mask_hit_but_no[mask_hit & (part.state==1)] = [np.isclose(ee, E0, atol=100*tol, rtol=tol)
                                                   for ee in part.energy[mask_hit & (part.state==1)]]
    if np.any(mask_hit_but_no):
        _test_drift(part, part_init, E0, mask_hit_but_no, coll, tol=tol)
        assert not np.any(mask_hit_but_no & is_child)
        assert not np.any(mask_hit_but_no & has_children)
        mask_hit[mask_hit_but_no] = False

    # The following really should have hit
    assert np.all(part.energy[mask_hit & (part.state==1)] < E0)
    # All that are supposed to hit and died without leaving children should still have all their energy
    assert np.allclose(part.energy[mask_hit & (part.state==coll_state) & ~has_children], E0, atol=tol, rtol=tol)
    # All that are supposed to hit and died with children should have lost energy
    assert np.all(part.energy[mask_hit & has_children] < E0)
    # All children should have less energy than the initial
    assert np.all(part.energy[is_child] < E0)
    # All energies need to be positive
    assert np.all(part.energy[(part.state!=LAST_INVALID_STATE) & (part.state!=xcc.ACC_IONISATION_LOSS)] > -1.e-12)

    with flaky_assertions():
        Edead = part.energy[part.state==coll_state].sum()
        if np.any(part.state==coll_state):
            assert Edead > 0
        # Eacc = part.energy[part.state==xcc.ACC_IONISATION_LOSS].sum()
        Eacc = coll._acc_ionisation_loss
        assert Eacc >= 0
        Evirtual = part.energy[part.state==xcc.VIRTUAL_ENERGY].sum()
        if np.any(part.state==xcc.VIRTUAL_ENERGY):
            assert Evirtual > 0
        Emassless = part.energy[part.state==xcc.MASSLESS_OR_NEUTRAL].sum()
        if np.any(part.state==xcc.MASSLESS_OR_NEUTRAL):
            assert Emassless > 0
        Eout = part.energy[part.state==1].sum()
        if alive_after > 0:
            assert Eout > 0
        print(E0*alive_before, Edead + Eacc + Evirtual + Emassless + Eout)
        assert np.isclose(E0*alive_before, Edead + Eacc + Evirtual + Emassless + Eout, atol=100*tol, rtol=tol)


def _test_drift(part, part_init, E0, mask, coll, tol=1e-12):
    assert np.all(part.state[mask] == 1)
    assert len(np.unique(part.energy[mask] == 1))
    assert np.allclose(part.energy[mask], E0, atol=100*tol, rtol=tol)
    dri = xt.Drift(length=coll.length)
    dri.model = 'exact'
    this_part_init = part_init.copy()
    dri.track(this_part_init)
    assert np.allclose(part.x[mask], this_part_init.x[mask], atol=100*tol, rtol=tol)
    assert np.allclose(part.px[mask], this_part_init.px[mask], atol=100*tol, rtol=tol)
    assert np.allclose(part.y[mask], this_part_init.y[mask], atol=100*tol, rtol=tol)
    assert np.allclose(part.py[mask], this_part_init.py[mask], atol=100*tol, rtol=tol)
    assert np.allclose(part.zeta[mask], this_part_init.zeta[mask], atol=100*tol, rtol=tol)
    assert np.allclose(part.delta[mask], this_part_init.delta[mask], atol=100*tol, rtol=tol)
