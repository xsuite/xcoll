import time
import numpy as np
import pytest

import xtrack as xt
import xtrack.particles.pdg as pdg
from xtrack.particles import LAST_INVALID_STATE
from xtrack.particles import masses as xpm

import xcoll as xc
from  xcoll import particle_states as xcp
from xcoll import fluka_masses as xflm


FLUKA_PROTON_MASS_EV = xflm[2212][1]
FLUKA_ELECTRON_MASS_EV = xflm[11][1]
FLUKA_MUON_MASS_EV = xflm[13][1]
FLUKA_PION_MASS_EV = xflm[211][1]
FLUKA_LEAD_MASS_EV = xflm[1000822080][1]


@pytest.mark.parametrize('hit', [True, False], ids=['hit', 'miss'])
def test_protons(hit):
    print(f"Testing protons with hit={hit}")
    num_part = 500
    p0c = 6.8e12
    capacity = num_part*100
    _stop_engine()
    particle_ref = xt.Particles.reference_from_pdg_id(pdg_id='proton', p0c=p0c)
    coll = _init_fluka(particle_ref, capacity)
    if hit:
        part, part_init = _init_particles(num_part)
    else:
        part, part_init = _init_particles(num_part, x=0, px=0)
    _track_particles(coll, part)
    tol=1e-12
    if hit:
        _assert_hit(part, part_init, coll, tol=tol)
    else:
        _assert_missed(part, part_init, coll, tol=tol)
    _stop_engine()


@pytest.mark.parametrize('hit', [True, False], ids=['hit', 'miss'])
def test_lead(hit):
    print(f"Testing lead with hit={hit}.")
    num_part = 100
    p0c = 6.8e12*82
    capacity = 500_000
    _stop_engine()
    if hit:
        kwargs_ref = {}
    else:
        kwargs_ref = {'x': 0, 'px': 0}
    particle_ref = xt.Particles.reference_from_pdg_id(pdg_id='Pb208', p0c=p0c)
    coll = _init_fluka(particle_ref, capacity)
    part, part_init = _init_particles(num_part, **kwargs_ref)
    _track_particles(coll, part)
    tol=1e-12
    if hit:
        _assert_hit(part, part_init, coll, tol=tol)
    else:
        _assert_missed(part, part_init, coll, tol=tol)
    _stop_engine()


@pytest.mark.parametrize('hit', [True, False], ids=['hit', 'miss'])
@pytest.mark.parametrize('proton_ref', [True, False], ids=['proton_ref', 'antiproton_ref'])
def test_antiprotons(proton_ref, hit):
    print(f"Testing antiprotons with hit={hit} and proton_ref={proton_ref}")
    num_part = 100
    p0c = 6.8e12
    capacity = num_part*100
    _stop_engine()
    if hit:
        kwargs_ref = {}
    else:
        kwargs_ref = {'x': 0, 'px': 0}
    if proton_ref:
        particle_ref = xt.Particles.reference_from_pdg_id(pdg_id='proton', p0c=p0c)
        kwargs_ref.update({'mass_ratio': 1, 'charge_ratio': -1,
                           'pdg_id': -2212, 'ptau': 0})
        tol=3e-11
    else:
        particle_ref = xt.Particles.reference_from_pdg_id(pdg_id='antiproton', p0c=p0c)
        tol=1e-12
    coll = _init_fluka(particle_ref, capacity)
    part, part_init = _init_particles(num_part, **kwargs_ref)
    _track_particles(coll, part)
    if hit:
        _assert_hit(part, part_init, coll, tol=tol)
    else:
        _assert_missed(part, part_init, coll, tol=tol)
    _stop_engine()


@pytest.mark.parametrize('hit', [True, False], ids=['hit', 'miss'])
@pytest.mark.parametrize('proton_ref', [True, False], ids=['proton_ref', 'electron_ref'])
def test_electrons(proton_ref, hit):
    print(f"Testing electrons with hit={hit} and proton_ref={proton_ref}")
    num_part = 500
    p0c = 200e9
    capacity = num_part*100
    _stop_engine()
    if hit:
        kwargs_ref = {}
    else:
        kwargs_ref = {'x': 0, 'px': 0}
    if proton_ref:
        particle_ref = xt.Particles.reference_from_pdg_id(pdg_id='proton', p0c=p0c)
        mratio = FLUKA_ELECTRON_MASS_EV/FLUKA_PROTON_MASS_EV
        ptau = (1/mratio-1)*np.sqrt(p0c**2 + FLUKA_PROTON_MASS_EV**2)/p0c
        kwargs_ref.update({'mass_ratio': mratio, 'charge_ratio': -1,
                           'pdg_id': 11, 'ptau': ptau})
        tol=3e-11
    else:
        particle_ref = xt.Particles.reference_from_pdg_id(pdg_id='electron', p0c=p0c)
        tol=1e-12
    coll = _init_fluka(particle_ref, capacity)
    part, part_init = _init_particles(num_part, **kwargs_ref)
    _track_particles(coll, part)
    if hit:
        _assert_hit(part, part_init, coll, tol=tol)
    else:
        _assert_missed(part, part_init, coll, tol=tol)
    _stop_engine()


@pytest.mark.parametrize('hit', [True, False], ids=['hit', 'miss'])
@pytest.mark.parametrize('proton_ref', [True, False], ids=['proton_ref', 'positron_ref'])
def test_positrons(proton_ref, hit):
    print(f"Testing positrons with hit={hit} and proton_ref={proton_ref}")
    num_part = 500
    p0c = 200e9
    capacity = num_part*100
    _stop_engine()
    if hit:
        kwargs_ref = {}
    else:
        kwargs_ref = {'x': 0, 'px': 0}
    if proton_ref:
        particle_ref = xt.Particles.reference_from_pdg_id(pdg_id='proton', p0c=p0c)
        mratio = FLUKA_ELECTRON_MASS_EV/FLUKA_PROTON_MASS_EV
        ptau = (1/mratio-1)*np.sqrt(p0c**2 + FLUKA_PROTON_MASS_EV**2)/p0c
        kwargs_ref.update({'mass_ratio': mratio, 'charge_ratio': 1,
                           'pdg_id': -11, 'ptau': ptau})
        tol=3e-11
    else:
        particle_ref = xt.Particles.reference_from_pdg_id(pdg_id='positron', p0c=p0c)
        tol=1e-12
    coll = _init_fluka(particle_ref, capacity)
    part, part_init = _init_particles(num_part, **kwargs_ref)
    _track_particles(coll, part)
    if hit:
        _assert_hit(part, part_init, coll, tol=tol)
    else:
        _assert_missed(part, part_init, coll, tol=tol)
    _stop_engine()


@pytest.mark.parametrize('hit', [True, False], ids=['hit', 'miss'])
@pytest.mark.parametrize('proton_ref', [True, False], ids=['proton_ref', 'muon_ref'])
def test_muons(proton_ref, hit):
    print(f"Testing muons with hit={hit} and proton_ref={proton_ref}")
    num_part = 500
    p0c = 200e9
    capacity = num_part*100
    _stop_engine()
    if hit:
        kwargs_ref = {}
    else:
        kwargs_ref = {'x': 0, 'px': 0}
    if proton_ref:
        particle_ref = xt.Particles.reference_from_pdg_id(pdg_id='proton', p0c=p0c)
        mratio = FLUKA_MUON_MASS_EV/FLUKA_PROTON_MASS_EV
        ptau = (1/mratio-1)*np.sqrt(p0c**2 + FLUKA_PROTON_MASS_EV**2)/p0c
        kwargs_ref.update({'mass_ratio': mratio, 'charge_ratio': -1,
                           'pdg_id': 13, 'ptau': ptau})
        tol=3e-11
    else:
        particle_ref = xt.Particles.reference_from_pdg_id(pdg_id='muon', p0c=p0c)
        tol=1e-12
    coll = _init_fluka(particle_ref, capacity)
    part, part_init = _init_particles(num_part, **kwargs_ref)
    _track_particles(coll, part)
    if hit:
        _assert_hit(part, part_init, coll, tol=tol)
    else:
        _assert_missed(part, part_init, coll, tol=tol)
    _stop_engine()


@pytest.mark.parametrize('hit', [True, False], ids=['hit', 'miss'])
@pytest.mark.parametrize('proton_ref', [True, False], ids=['proton_ref', 'antimuon_ref'])
def test_antimuons(proton_ref, hit):
    print(f"Testing antimuons with hit={hit} and proton_ref={proton_ref}")
    num_part = 500
    p0c = 200e9
    capacity = num_part*100
    _stop_engine()
    if hit:
        kwargs_ref = {}
    else:
        kwargs_ref = {'x': 0, 'px': 0}
    if proton_ref:
        particle_ref = xt.Particles.reference_from_pdg_id(pdg_id='proton', p0c=p0c)
        mratio = FLUKA_MUON_MASS_EV/FLUKA_PROTON_MASS_EV
        ptau = (1/mratio-1)*np.sqrt(p0c**2 + FLUKA_PROTON_MASS_EV**2)/p0c
        kwargs_ref.update({'mass_ratio': mratio, 'charge_ratio': 1,
                           'pdg_id': -13, 'ptau': ptau})
        tol=3e-11
    else:
        particle_ref = xt.Particles.reference_from_pdg_id(pdg_id='antimuon', p0c=p0c)
        tol=1e-12
    coll = _init_fluka(particle_ref, capacity)
    part, part_init = _init_particles(num_part, **kwargs_ref)
    _track_particles(coll, part)
    if hit:
        _assert_hit(part, part_init, coll, tol=tol)
    else:
        _assert_missed(part, part_init, coll, tol=tol)
    _stop_engine()

@pytest.mark.parametrize('hit', [True, False], ids=['hit', 'miss'])
@pytest.mark.parametrize('proton_ref', [True, False], ids=['proton_ref', 'pion_ref'])
def test_positive_pions(proton_ref, hit):
    print(f"Testing positive pions with hit={hit} and proton_ref={proton_ref}")
    num_part = 500
    p0c = 450e9
    capacity = num_part*100
    _stop_engine()
    if hit:
        kwargs_ref = {}
    else:
        kwargs_ref = {'x': 0, 'px': 0}
    if proton_ref:
        particle_ref = xt.Particles.reference_from_pdg_id(pdg_id='proton', p0c=p0c)
        mratio = FLUKA_PION_MASS_EV/FLUKA_PROTON_MASS_EV
        ptau = (1/mratio-1)*np.sqrt(p0c**2 + FLUKA_PROTON_MASS_EV**2)/p0c
        kwargs_ref.update({'mass_ratio': mratio, 'charge_ratio': 1,
                           'pdg_id': 211, 'ptau': ptau})
        tol=3e-11
    else:
        particle_ref = xt.Particles.reference_from_pdg_id(pdg_id='pi+', p0c=p0c)
        tol=1e-12
    coll = _init_fluka(particle_ref, capacity)
    part, part_init = _init_particles(num_part, **kwargs_ref)
    _track_particles(coll, part)
    if hit:
        _assert_hit(part, part_init, coll, tol=tol)
    else:
        _assert_missed(part, part_init, coll, tol=tol)
    _stop_engine()


@pytest.mark.parametrize('hit', [True, False], ids=['hit', 'miss'])
@pytest.mark.parametrize('proton_ref', [True, False], ids=['proton_ref', 'pion_ref'])
def test_negative_pions(proton_ref, hit):
    print(f"Testing negative pions with hit={hit} and proton_ref={proton_ref}")
    num_part = 500
    p0c = 450e9
    capacity = num_part*100
    _stop_engine()
    if hit:
        kwargs_ref = {}
    else:
        kwargs_ref = {'x': 0, 'px': 0}
    if proton_ref:
        particle_ref = xt.Particles.reference_from_pdg_id(pdg_id='proton', p0c=p0c)
        mratio = FLUKA_PION_MASS_EV/FLUKA_PROTON_MASS_EV
        ptau = (1/mratio-1)*np.sqrt(p0c**2 + FLUKA_PROTON_MASS_EV**2)/p0c
        kwargs_ref.update({'mass_ratio': mratio, 'charge_ratio': -1,
                           'pdg_id': -211, 'ptau': ptau})
        tol=3e-11
    else:
        particle_ref = xt.Particles.reference_from_pdg_id(pdg_id='pi-', p0c=p0c)
        tol=1e-12
    coll = _init_fluka(particle_ref, capacity)
    part, part_init = _init_particles(num_part, **kwargs_ref)
    _track_particles(coll, part)
    if hit:
        _assert_hit(part, part_init, coll, tol=tol)
    else:
        _assert_missed(part, part_init, coll, tol=tol)
    _stop_engine()


def _stop_engine():
    if xc.FlukaEngine.is_running():
        xc.FlukaEngine.stop(clean=True)


def _init_fluka(particle_ref, capacity):
    coll = xc.FlukaCollimator(length=0.4, assembly='fcc_tcp')
    coll.jaw = 0.002
    xc.FlukaEngine.particle_ref = particle_ref
    xc.FlukaEngine.capacity = capacity
    xc.FlukaEngine.start(elements=coll, clean=True, verbose=True, return_all=True,
                        return_neutral=True, electron_lower_momentum_cut=1.e6)
    return coll


def _init_particles(num_part, x=0.0022, px=-1.e-5, **kwargs):
    x_init  = np.random.normal(loc=x, scale=0.2e-3, size=num_part)
    px_init = np.random.normal(loc=px, scale=5.e-6, size=num_part)
    y_init  = np.random.normal(loc=0., scale=1e-3, size=num_part)
    py_init = np.random.normal(loc=0., scale=5.e-6, size=num_part)
    pref    = xc.FlukaEngine.particle_ref
    kwargs.setdefault('p0c', pref.p0c)
    kwargs.setdefault('mass0', pref.mass0)
    kwargs.setdefault('q0', pref.q0)
    kwargs.setdefault('pdg_id', pref.pdg_id)
    part_init = xt.Particles(x=x_init, px=px_init, y=y_init, py=py_init,
                             _capacity=xc.FlukaEngine.capacity, **kwargs)
    print("Created particles with the following parameters:")
    E0 = xc.FlukaEngine.particle_ref.energy[0]
    print(f"Reference energy: {E0} eV (in particles: {part_init.energy[0]} eV)")
    print(f"Mass {part_init.mass[0]} eV")
    print(f"Charge {part_init.charge_ratio[0]*part_init.q0}")
    print(f"PDG ID {part_init.pdg_id[0]}")
    return part_init.copy(), part_init


def _track_particles(coll, part):
    print(f"Tracking {part._num_active_particles} {pdg.get_name_from_pdg_id(part.pdg_id[0])}"
         + "s...     ", end='', flush=True)
    start = time.time()
    coll.track(part)
    print(f"Done in {round(time.time()-start, 3)}s.", flush=True)
    print()


def _assert_missed(part, part_init, coll, tol=1e-12):
    E0 = xc.FlukaEngine.particle_ref.energy[0]
    alive = part_init._num_active_particles
    assert part._num_active_particles == alive
    assert set(np.unique(part.state)) == {1, LAST_INVALID_STATE}
    mask_alive = False*np.ones(part.x.shape, dtype=bool)
    mask_alive[:alive] = True
    _test_drift(part, part_init, mask_alive, coll, tol=tol)


def _test_drift(part, part_init, mask, coll, tol=1e-12):
    E0 = xc.FlukaEngine.particle_ref.energy[0]
    assert np.all(part.state[mask] == 1)
    assert len(np.unique(part.energy[mask] == 1))
    assert np.allclose(part.energy[mask], E0, atol=100*tol, rtol=tol)
    dri = xt.Drift(length=coll.length)
    this_part_init = part_init.copy()
    dri.track(this_part_init)
    assert np.allclose(part.x[mask], this_part_init.x[mask], atol=100*tol, rtol=tol)
    assert np.allclose(part.px[mask], this_part_init.px[mask], atol=100*tol, rtol=tol)
    assert np.allclose(part.y[mask], this_part_init.y[mask], atol=100*tol, rtol=tol)
    assert np.allclose(part.py[mask], this_part_init.py[mask], atol=100*tol, rtol=tol)
    assert np.allclose(part.zeta[mask], this_part_init.zeta[mask], atol=100*tol, rtol=tol)
    assert np.allclose(part.delta[mask], this_part_init.delta[mask], atol=100*tol, rtol=tol)


def _assert_hit(part, part_init, coll, tol=1e-12):
    E0 = xc.FlukaEngine.particle_ref.energy[0]
    alive_before = part_init._num_active_particles
    alive_after = part._num_active_particles
    part.sort(interleave_lost_particles=True)
    is_child = np.array([pid!=ppid for pid, ppid in zip(part.particle_id, part.parent_particle_id)])
    has_children = np.array([pid in part.parent_particle_id[is_child] for pid in part.particle_id])

    mask_hit = (part_init.x >= coll.jaw_L) | (part_init.x <= coll.jaw_R)
    mask_hit &= part_init.state != LAST_INVALID_STATE
    mask_not_hit = (part_init.x < coll.jaw_L) & (part_init.x > coll.jaw_R)
    mask_not_hit &= part_init.state != LAST_INVALID_STATE

    # All that have not hit should be alive, not have lost energy, and just have drifted
    _test_drift(part, part_init, mask_not_hit, coll, tol=tol)
    assert not np.any(mask_not_hit & is_child)
    assert not np.any(mask_not_hit & has_children)

    # All that are supposed to hit and are still alive should have lost energy
    # Except those (which actually did not hit but barely scratched the collimator):
    mask_hit_but_no = False*np.ones(mask_hit.shape, dtype=bool)
    mask_hit_but_no[mask_hit & (part.state==1)] = [np.isclose(ee, E0, atol=100*tol, rtol=tol)
                                                   for ee in part.energy[mask_hit & (part.state==1)]]
    if np.any(mask_hit_but_no):
        _test_drift(part, part_init, mask_hit_but_no, coll, tol=tol)
        assert not np.any(mask_hit_but_no & is_child)
        assert not np.any(mask_hit_but_no & has_children)
        mask_hit[mask_hit_but_no] = False
    # The following really should have hit
    assert np.all(part.energy[mask_hit & (part.state==1)] < E0)
    # All that are supposed to hit and died without leaving children should still have all their energy
    assert np.all(part.energy[mask_hit & (part.state==xcp.LOST_ON_FLUKA_COLL) & ~has_children] == E0)
    # All that are supposed to hit and died with children should have lost energy
    assert np.all(part.energy[mask_hit & has_children] < E0)
    # All children should have less energy than the initial
    assert np.all(part.energy[is_child] < E0)
    # All energies need to be positive
    assert np.all(part.energy[(part.state!=LAST_INVALID_STATE) & (part.state!=xcp.ACC_IONISATION_LOSS)] > 0)

    Edead = part.energy[part.state==xcp.LOST_ON_FLUKA_COLL].sum()
    if np.any(part.state==xcp.LOST_ON_FLUKA_COLL):
        assert Edead > 0
    # Eacc = part.energy[part.state==xcp.ACC_IONISATION_LOSS].sum()
    Eacc = coll._acc_ionisation_loss
    assert Eacc >= 0
    Evirtual = part.energy[part.state==xcp.VIRTUAL_ENERGY].sum()
    if np.any(part.state==xcp.VIRTUAL_ENERGY):
        assert Evirtual > 0
    Emassless = part.energy[part.state==xcp.MASSLESS_OR_NEUTRAL].sum()
    if np.any(part.state==xcp.MASSLESS_OR_NEUTRAL):
        assert Emassless > 0
    Eout = part.energy[part.state==1].sum()
    if alive_after > 0:
        assert Eout > 0
    assert np.isclose(E0*alive_before, Edead + Eacc + Evirtual + Emassless + Eout, atol=100*tol, rtol=tol)
