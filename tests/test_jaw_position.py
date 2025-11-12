# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

import pytest
import numpy as np
from pathlib import Path

import xpart as xp
import xtrack as xt
import xcoll as xc
try:
    import rpyc
except ImportError as e:
    rpyc = None

path = Path(__file__).parent / 'data'


# TODO: particle angles

jaws = [0.001, [0.0013, -0.002789], [-1.2e-6, -3.2e-3], [3.789e-3, 4.678e-7]]
jaw_ids = ['symmetric', 'asymmetric', 'negative', 'positive']
angles = [0, 90, 127.5]
tilts = [0, [2.2e-6, 1.3e-6], [1.9e-6, -2.7e-6]]
tilt_ids = ['no_tilt', 'positive_tilt', 'pos_neg_tilt']

particle_ref = xt.Particles('proton', p0c=6.8e12)


@pytest.mark.parametrize('tilt', tilts, ids=tilt_ids)
@pytest.mark.parametrize('angle', angles)
@pytest.mark.parametrize('jaw', jaws, ids=jaw_ids)
def test_everest(jaw, angle, tilt):
    num_part = 2_500_000
    length = 0.873
    coll = xc.EverestCollimator(length=length, jaw=jaw, angle=angle, tilt=tilt, material='MoGR')
    part_init, hit_ids, not_hit_ids = _generate_particles(coll, num_part=num_part, particle_ref=particle_ref,
                                                jaw_band=5e-9, angular_spread=1e-3, delta_spread=1e-3,
                                                zeta_spread=5e-2)
    part = part_init.copy()
    coll.track(part)
    # _plot_jaws(coll, part_init, part, hit_ids, not_hit_ids)
    _assert_valid_positions(part_init, part, hit_ids, not_hit_ids)


# TODO: why is BDSIM so imprecise? Most hits are very exact, but on the negative jaw, there are a few
# particles with positive kick and slightly inside the jaw (up to around 1e-9) that survive but should
# have hit, and similarily for the other jaw with negative kicks. Though most particles behave very well,
# also in these regions.
@pytest.mark.parametrize('tilt', tilts, ids=tilt_ids)
@pytest.mark.parametrize('angle', angles)
@pytest.mark.parametrize('jaw', jaws[:2], ids=jaw_ids[:2])  # BDSIM cannot handle fully positive or negative jaws
@pytest.mark.skipif(rpyc is None, reason="rpyc not installed")
@pytest.mark.skipif(not xc.geant4.environment.compiled, reason="BDSIM+Geant4 installation not found")
def test_geant4(jaw, angle, tilt):
    num_part = 50_000
    length = 0.873
    if xc.geant4.engine.is_running():
        xc.geant4.engine.stop(clean=True)
    coll = xc.Geant4Collimator(length=length, jaw=jaw, angle=angle, tilt=tilt, material='MoGR')
    xc.geant4.engine.particle_ref = particle_ref
    xc.geant4.engine.start(elements=coll, relative_energy_cut=0.1)
    part_init, hit_ids, not_hit_ids = _generate_particles(coll, num_part=num_part, particle_ref=particle_ref,
                                                jaw_band=1e-8, jaw_accuracy=5e-9, angular_spread=1e-3,
                                                delta_spread=1e-3, zeta_spread=5e-2, exact_drift=True,
                                                _capacity=2*num_part)
    part = part_init.copy()
    coll.track(part)
    # _plot_jaws(coll, part_init, part, hit_ids, not_hit_ids)
    _assert_valid_positions(part_init, part, hit_ids, not_hit_ids)
    xc.geant4.engine.stop(clean=True)


# TODO: lhc_tdi and fcc_tcdq still fail. Are they unrotatable? Side fixed or ill-defined?
# @pytest.mark.parametrize('tilt', tilts, ids=tilt_ids)
@pytest.mark.parametrize('assembly', ['sps_tcsm', 'lhc_tcp', 'lhc_tcsg', 'lhc_tcsp',
                                      'lhc_tcla', 'lhc_tct', 'lhc_tcl', 'lhc_tdi',
                                      'lhc_tclia', 'lhc_tclib', 'lhc_tcdqaa', 'lhc_tcdqab',
                                      'lhc_tcdqac', 'hilumi_tcppm', 'hilumi_tcspm', 'hilumi_tcspgrc',
                                      'hilumi_tcld', 'hilumi_tctx', 'hilumi_tcty', 'hilumi_tclx',
                                      'fcc_tcp', 'fcc_tcsg', 'fcc_tcdq'])
@pytest.mark.parametrize('angle', angles)
@pytest.mark.parametrize('jaw', jaws, ids=jaw_ids)
def test_fluka(jaw, angle, assembly):
    tilt = 0

    # If a previous test failed, stop the server manually
    if xc.fluka.engine.is_running():
        xc.fluka.engine.stop(clean=True)

    # Define collimator and start the FLUKA server
    coll = xc.FlukaCollimator(length=1, jaw=jaw, angle=angle, tilt=tilt, assembly=assembly)
    xc.fluka.engine.particle_ref = particle_ref
    xc.fluka.engine.start(elements=coll, capacity=10_000)

    part_init, hit_ids, not_hit_ids = _generate_particles(coll, num_part=5000, x_dim=0.015,
                            _capacity=xc.fluka.engine.capacity, particle_ref=xc.fluka.engine.particle_ref)
    part = part_init.copy()
    coll.track(part)
    _assert_valid_positions(part, hit_ids, not_hit_ids)

    # Stop the FLUKA server
    xc.fluka.engine.stop(clean=True)


def _generate_particles(coll, num_part, particle_ref, _capacity=None, jaw_band=1.e-6,
                        jaw_accuracy=1.e-12, x_dim=0.05,  y_dim=0.05, angular_spread=0,
                        delta_spread=0, zeta_spread=0, exact_drift=False):
    if _capacity is None:
        _capacity = num_part

    # Particle distribution (x and y are in the frame of the collimator)
    num_part_step = num_part//10
    num_part = 10*num_part_step

    x = np.random.uniform(-x_dim, x_dim, 2*num_part_step)

    if coll.side != 'both' and not coll.side is None:  # None is hack as FLUKA assemblies have no side yet TODO: update metadata
        num_part_step *= 2
    if coll.side == 'left' or coll.side == 'both' or coll.side is None:  # None is hack as FLUKA assemblies have no side yet TODO: update metadata
        x = np.concatenate([x, np.random.uniform(coll.jaw_LU - jaw_band, coll.jaw_LU - jaw_accuracy, num_part_step)])
        x = np.concatenate([x, np.random.uniform(coll.jaw_LU + jaw_accuracy, coll.jaw_LU + jaw_band, num_part_step)])
        x = np.concatenate([x, np.random.uniform(coll.jaw_LD - jaw_band, coll.jaw_LD - jaw_accuracy, num_part_step)])
        x = np.concatenate([x, np.random.uniform(coll.jaw_LD + jaw_accuracy, coll.jaw_LD + jaw_band, num_part_step)])
    if coll.side == 'right' or coll.side == 'both' or coll.side is None:  # None is hack as FLUKA assemblies have no side yet TODO: update metadata
        x = np.concatenate([x, np.random.uniform(coll.jaw_RU - jaw_band, coll.jaw_RU - jaw_accuracy, num_part_step)])
        x = np.concatenate([x, np.random.uniform(coll.jaw_RU + jaw_accuracy, coll.jaw_RU + jaw_band, num_part_step)])
        x = np.concatenate([x, np.random.uniform(coll.jaw_RD - jaw_band, coll.jaw_RD - jaw_accuracy, num_part_step)])
        x = np.concatenate([x, np.random.uniform(coll.jaw_RD + jaw_accuracy, coll.jaw_RD + jaw_band, num_part_step)])

    if y_dim > 0.:
        y = np.random.uniform(-y_dim, y_dim, num_part)
    else:
        y = np.zeros_like(x)

    x_new = np.cos(np.deg2rad(coll.angle))*x - np.sin(np.deg2rad(coll.angle))*y
    y_new = np.sin(np.deg2rad(coll.angle))*x + np.cos(np.deg2rad(coll.angle))*y
    if delta_spread > 0.:
        delta = np.random.uniform(-delta_spread, delta_spread, size=len(x_new))
    else:
        delta = np.zeros_like(x)
    if zeta_spread > 0.:
        zeta = np.random.uniform(-zeta_spread, zeta_spread, size=len(x_new))
    else:
        zeta = np.zeros_like(x)
    if angular_spread > 0.:
        px = np.random.uniform(-angular_spread, angular_spread, size=len(x_new))
        if y_dim > 0.:
            py = np.random.uniform(-angular_spread, angular_spread, size=len(y_new))
        else:
            py = np.zeros_like(px)
        if exact_drift:
            pz = np.sqrt((1. + delta)**2 - px**2 - py**2)
        else:
            pz = 1. + delta
        px_new = np.cos(np.deg2rad(coll.angle))*px - np.sin(np.deg2rad(coll.angle))*py
        py_new = np.sin(np.deg2rad(coll.angle))*px + np.cos(np.deg2rad(coll.angle))*py
    else:
        px = 0; py = 0; pz = 1; px_new = 0; py_new = 0
    part_init = xp.build_particles(x=x_new, y=y_new, px=px_new, py=py_new, delta=delta,
                                   zeta=zeta, particle_ref=particle_ref, _capacity=_capacity)
    mask  = (x + px/pz*coll.jaw_s_LU >= coll.jaw_LU) | (x + px/pz*coll.jaw_s_LD >= coll.jaw_LD)
    mask |= (x + px/pz*coll.jaw_s_RU <= coll.jaw_RU) | (x + px/pz*coll.jaw_s_RD <= coll.jaw_RD)
    mask = np.concatenate([mask, np.full(_capacity-num_part, False)])
    expected_hit_ids = part_init.particle_id[mask & (part_init.state > 0)]
    expected_not_hit_ids = part_init.particle_id[~mask & (part_init.state > 0)]

    return part_init, expected_hit_ids, expected_not_hit_ids


def _mask_hits(part_init, part, momentum_accuracy=1.e-12):
    assigned_particles = part.particle_id >= 0
    mask_hits  = (part.state[assigned_particles] < 1) # died
    parent_id = part.parent_particle_id[assigned_particles]
    mask_hits |= abs(part_init.px[parent_id] - part.px[assigned_particles]) > momentum_accuracy # kick in x
    mask_hits |= abs(part_init.py[parent_id] - part.py[assigned_particles]) > momentum_accuracy # kick in y
    mask_hits |= abs(part_init.delta[parent_id] - part.delta[assigned_particles]) > momentum_accuracy # energy loss
    mask_hits = np.concatenate([mask_hits, np.full(part._capacity-len(mask_hits), False)])
    return mask_hits

def _mask_not_hits(part_init, part, momentum_accuracy=1.e-12):
    mask_hits = _mask_hits(part_init, part, momentum_accuracy)
    return (~mask_hits) & (part.state > -999999)


def _assert_valid_positions(part_init, part, expected_hit_ids, expected_not_hit_ids, momentum_accuracy=1.e-12):
    expected_hit = np.isin(part.parent_particle_id, expected_hit_ids)
    expected_not_hit = np.isin(part.parent_particle_id, expected_not_hit_ids)

    hits = _mask_hits(part_init, part, momentum_accuracy)
    not_hits = _mask_not_hits(part_init, part, momentum_accuracy)

    print(f"{part.particle_id.max() - part_init.particle_id.max()} children generated.")
    print(f"{sum(expected_hit & not_hits)} particles were expected to hit but did not.")
    print(f"{sum(expected_not_hit & hits)} particles hit but were not expected.")

    # If one child particle has hit, no other child particle should be marked as not having hit
    assert set(part.parent_particle_id[hits]) & set(part.parent_particle_id[not_hits]) == set()
    assert set(part.parent_particle_id[hits]) | set(part.parent_particle_id[not_hits]) == set(range(part_init.particle_id.max() + 1))

    # Particles that are supposed to not have hit the collimator, but have a kick or are dead, are considered faulty
    assert sum(expected_not_hit & hits) == 0

    # Particles that are supposed to have hit the collimator, but are alive and have no kick, are considered faulty
    assert sum(expected_hit & not_hits) <= 1 # We allow for a small margin of error


def _plot_jaws(coll, part_init, part, expected_hit_ids, expected_not_hit_ids, momentum_accuracy=1.e-12):
    import matplotlib.pyplot as plt

    hit_ids = part.parent_particle_id[_mask_hits(part_init, part, momentum_accuracy)]
    not_hit_ids = part.parent_particle_id[_mask_not_hits(part_init, part, momentum_accuracy)]

    expected_but_not_hit = set(expected_hit_ids) - set(hit_ids)
    hit_but_not_expected = set(hit_ids) - set(expected_hit_ids)
    not_expected_but_hit = set(expected_not_hit_ids) - set(not_hit_ids)
    not_hit_but_expected = set(not_hit_ids) - set(expected_not_hit_ids)
    assert expected_but_not_hit == not_hit_but_expected # sanity check
    assert hit_but_not_expected == not_expected_but_hit # sanity check
    expected_but_not_hit = list(expected_but_not_hit)
    hit_but_not_expected = list(hit_but_not_expected)

    expected_and_hit = list(set(expected_hit_ids) & set(hit_ids))
    not_expected_and_not_hit = list(set(expected_not_hit_ids) & set(not_hit_ids))

    x = np.cos(np.deg2rad(coll.angle))*part_init.x + np.sin(np.deg2rad(coll.angle))*part_init.y
    y = - np.sin(np.deg2rad(coll.angle))*part_init.x + np.cos(np.deg2rad(coll.angle))*part_init.y
    px = np.cos(np.deg2rad(coll.angle))*part_init.px + np.sin(np.deg2rad(coll.angle))*part_init.py
    py = - np.sin(np.deg2rad(coll.angle))*part_init.px + np.cos(np.deg2rad(coll.angle))*part_init.py

    fig, ax = plt.subplots(2, 2, figsize=(14, 8))
    for i, arr in enumerate([[x, y], [x, px], [y, py], [part_init.zeta, part_init.delta]]):
        idx = i % 2
        idy = i // 2
        ax[idx][idy].scatter(arr[0][expected_and_hit],         arr[1][expected_and_hit], c='tab:blue', s=4, label='Hit (as expected)')
        ax[idx][idy].scatter(arr[0][not_expected_and_not_hit], arr[1][not_expected_and_not_hit], c='tab:green', s=4, label='Not hit (as expected)')
        ax[idx][idy].scatter(arr[0][expected_but_not_hit],     arr[1][expected_but_not_hit], c='tab:red', s=4, label='Not hit (but expected to hit)')
        ax[idx][idy].scatter(arr[0][hit_but_not_expected],     arr[1][hit_but_not_expected], c='tab:brown', s=4, label='Hit (but not expected)')
        if idy == 0:
            ax[idx][idy].axvline(coll.jaw_LU, c='k', linestyle='--')
            ax[idx][idy].axvline(coll.jaw_LD, c='k', linestyle='--')
            ax[idx][idy].axvline(coll.jaw_RU, c='k', linestyle='--')
            ax[idx][idy].axvline(coll.jaw_RD, c='k', linestyle='--')
        ax[idx][idy].set_xlabel(['x [mm]', 'x [mm]', 'y [mm]', 'zeta [mm]'][i])
        ax[idx][idy].set_ylabel(['y [mm]', 'px [mrad]', 'py [mrad]', 'delta [1e-3]'][i])
        ax[idx][idy].legend(loc='upper right')
    fig.tight_layout()
    plt.show()
