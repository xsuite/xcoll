# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2026.                 #
# ######################################### #

import time
import pytest
import numpy as np
from pathlib import Path
from warnings import warn
from scipy.stats import landau

import xpart as xp
import xtrack as xt
import xcoll as xc
import xcoll.constants as xcc

from _common_api import old_bdsim, engine_params


@pytest.mark.fluka
def test_ionisation_loss():
    num_part = 20_000
    mat = xc.materials.db['MG6403Fc']
    coll = xc.FlukaCollimator(length=0.6, angle=0, jaw=0.001, material=mat)
    particle_ref = xt.Particles('proton', p0c=6.8e12)
    xc.fluka.engine.particle_ref = particle_ref
    xc.fluka.engine.capacity = 2*num_part
    xc.fluka.engine.include_elastic = False
    xc.fluka.engine.include_inelastic = False
    xc.fluka.engine.include_showers = False
    xc.fluka.engine.include_single_coulomb = False
    xc.fluka.engine.include_multiple_coulomb = False
    xc.fluka.engine.seed = 3456543
    part_init = xp.build_particles(
        x=np.random.uniform(0.002-1e-6, 0.002+1e-6, num_part),
        px=np.random.uniform(-1e-6, 1e-6, num_part),
        y=np.random.uniform(-1e-6, 1e-6, num_part),
        py=np.random.uniform(-1e-6, 1e-6, num_part),
        particle_ref=xc.fluka.engine.particle_ref,
        _capacity=xc.fluka.engine.capacity
    )

    # Ionisation parameters: analytic expressions to test against
    rho = mat.density
    Z_A = mat._ZA_mean
    I = mat.excitation_energy
    beta0 = particle_ref.beta0[0]
    gamma0 = particle_ref.gamma0[0]
    K = 0.307075e6         # eV cm^2 / g
    me = 510.99895e3       # electron mass
    mp = 938.27208816e6    # proton mass
    plasma = np.sqrt(rho * Z_A) * 28.816
    Tmax = 2 * me * beta0**2 * gamma0**2 / (1 + 2*gamma0*me/mp + (me/mp)**2)
    delta = (2 * np.log(plasma/I) + 2 * np.log(beta0 * gamma0) - 1)
    dEdx_mass = K*Z_A/beta0**2 * (0.5 * np.log(2*me*beta0**2*gamma0**2*Tmax/I**2) - beta0**2 - 0.5 * delta)
    dEdx = dEdx_mass * rho
    dE = dEdx * (100 * coll.length)
    xi = K/2 * Z_A * rho * 100*coll.length
    Empv = xi * np.log(2 * me * xi / plasma**2) + 0.2*xi   # Vavilov high-energy approximation
    # Scipy parameters for Landau distribution
    scipy_Empv = Empv + xi * (1 - np.euler_gamma - 0.20005183774398613) + xi * np.log(np.pi / 2)
    scipy_xi = xi * np.pi / 2

    # Check the mean stopping power with all physics deactivated
    xc.fluka.engine.include_ionisation_fluctuations = False
    xc.fluka.engine.include_pair_production = False
    xc.fluka.engine.include_bremsstrahlung = False
    xc.fluka.engine.start(elements=coll, clean=True, verbose=False, fortran_debug_level=1)
    xc.fluka.engine.physics_settings()
    part = part_init.copy()
    coll.track(part)
    xc.fluka.engine.stop(clean=True)
    E_diff = part.energy0[part.state > 0] - part.energy[part.state > 0]
    mean_stopping_power = np.unique(E_diff)
    assert len(mean_stopping_power) == 1
    print(f"Mean stopping power: {mean_stopping_power[0]/1e6} MeV")
    assert np.isclose(mean_stopping_power[0], dE, rtol=1e-2, atol=0)
    # We expect 410.12 MeV, which corresponds to 1/rho dE/dx of 2.68 MeV cm2/g
    assert np.isclose(mean_stopping_power[0], 410.12e6, rtol=1e-7, atol=0)

    # Check the regular ionisation losses
    xc.fluka.engine.include_ionisation_fluctuations = True
    xc.fluka.engine.include_pair_production = False
    xc.fluka.engine.include_bremsstrahlung = False
    # Higher precision ionisation loss
    extra_card = "IONFLUCT         1.0       0.0       4.0  BLCKHOLE  @LASTMAT"
    xc.fluka.engine.start(elements=coll, clean=True, verbose=False, fortran_debug_level=1, extra_cards=[extra_card])
    xc.fluka.engine.physics_settings()
    part = part_init.copy()
    coll.track(part)
    xc.fluka.engine.stop(clean=True)
    E_diff = part.energy0[part.state > 0] - part.energy[part.state > 0]
    assert np.isclose(np.median(E_diff), landau.median(scipy_Empv, scipy_xi), rtol=2e-3, atol=0)
    # Landau/Vavilov upper tail: P(dE > xi) = xi/T  =>  T_threshold = xi/P
    #       1 particle  above this threshold  =>  P = 1/N
    #      10 particles above this threshold  =>  P = 10/N
    #     100 particles above this threshold  =>  P = 100/N
    dE_threshold_1         = landau.ppf(1 - 1/num_part,   scipy_Empv, scipy_xi)
    dE_threshold_10        = landau.ppf(1 - 10/num_part,  scipy_Empv, scipy_xi)
    dE_threshold_100       = landau.ppf(1 - 100/num_part, scipy_Empv, scipy_xi)
    dE_threshold_lower_200 = landau.ppf(200/num_part, scipy_Empv, scipy_xi)
    num_outliers_1         = np.sum(E_diff > dE_threshold_1)
    num_outliers_10        = np.sum(E_diff > dE_threshold_10)
    num_outliers_100       = np.sum(E_diff > dE_threshold_100)
    num_outliers_lower_200 = np.sum(E_diff < dE_threshold_lower_200)
    # Bounds are Poissonian; only 1e-4 tests should fail
    assert num_outliers_1 < 6      # Expect 1 outlier
    assert num_outliers_10 < 24    # Expect 10 outliers
    assert num_outliers_100 < 139  # Expect 100 outliers
    assert num_outliers_lower_200 < 256 # Expect 200 outliers

    # Check the full energy losses (including pair production and bremsstrahlung)
    xc.fluka.engine.include_ionisation_fluctuations = True
    xc.fluka.engine.include_pair_production = True
    xc.fluka.engine.include_bremsstrahlung = True
    # Higher precision ionisation loss
    extra_card = "IONFLUCT         1.0       0.0       4.0  BLCKHOLE  @LASTMAT"
    xc.fluka.engine.start(elements=coll, clean=True, verbose=False, fortran_debug_level=1, extra_cards=[extra_card])
    xc.fluka.engine.physics_settings()
    part = part_init.copy()
    coll.track(part)
    xc.fluka.engine.stop(clean=True)
    E_diff = part.energy0[part.state > 0] - part.energy[part.state > 0]
    # Thresholds from file
    data = xc.json.json_load('data/fluka_calibration_data.json')
    dE_threshold_1         = data['full']['thresholds']['1']
    dE_threshold_10        = data['full']['thresholds']['10']
    dE_threshold_100       = data['full']['thresholds']['100']
    dE_threshold_lower_200 = data['full']['1st percentile']
    num_outliers_1         = np.sum(E_diff > dE_threshold_1)
    num_outliers_10        = np.sum(E_diff > dE_threshold_10)
    num_outliers_100       = np.sum(E_diff > dE_threshold_100)
    num_outliers_lower_200 = np.sum(E_diff < dE_threshold_lower_200)
    # Bounds are Poissonian; only 1e-4 tests should fail
    assert num_outliers_1 < 6      # Expect 1 outlier
    assert num_outliers_10 < 24    # Expect 10 outliers
    assert num_outliers_100 < 139  # Expect 100 outliers
    assert num_outliers_lower_200 < 256 # Expect 200 outliers


@pytest.mark.parametrize("engine", engine_params)
@pytest.mark.parametrize("log_impacts, mark_scattered_particles", [
                            [False, False],
                            [False, True],
                            [True,  False],
                            [True,  True]
                         ], ids=["default", "mark", "impacts", "impacts_mark"])
def test_deep_physics_check(engine, log_impacts, mark_scattered_particles, running_with_xdist):
    num_part = 24_000      # When this is changed, need to re-generate input distribution (rm data/{engine}_part_init.json)
    capacity = 4*num_part  # When this is changed, need to re-generate input distribution
    particle_ref = xt.Particles('proton', p0c=6.8e12)

    # Prepare collimators, engine, and initial particles
    if engine == 'fluka':
        xc_engine = xc.fluka.engine
        coll_has_flanges = True
        coll1 = xc.FlukaCollimator(length=0.6, angle=0,   jaw=0.001,  assembly='hilumi_tcppm')
        coll2 = xc.FlukaCollimator(length=0.6, angle=123, jaw=0.0005, assembly='hilumi_tcppm')
    elif engine == 'geant4':
        xc_engine = xc.geant4.engine
        coll_has_flanges = False
        coll1 = xc.Geant4Collimator(length=0.6, angle=0,   jaw=0.001,  material='MG6403Fc')
        coll2 = xc.Geant4Collimator(length=0.6, angle=123, jaw=0.0005, material='MG6403Fc')
    if mark_scattered_particles:
        coll1.mark_scattered_particles = True
        coll2.mark_scattered_particles = True
    else:
        coll1.mark_scattered_particles = False
        coll2.mark_scattered_particles = False

    xc_engine.particle_ref = particle_ref
    if engine == 'fluka':
        xc_engine.capacity = capacity
    xc_engine.seed = 453532
    xc_engine.reset_physics_settings()
    xc_engine.return_baryons = True   # To get some massless particles as well
    xc_engine.start(elements=[coll1, coll2], clean=True, verbose=False)
    xc_engine.physics_settings()

    # Create black absorbers after starting engine so length_front and length_back are known
    black1 = xc.BlackAbsorber(length=0.6 + coll1.length_front + coll1.length_back, angle=0,   jaw=0.001)
    black2 = xc.BlackAbsorber(length=0.6 + coll1.length_front + coll1.length_back, angle=123, jaw=0.0005)

    if log_impacts:
        impacts = xc.InteractionRecord(elements=[coll1, coll2], record_impacts=True)

    part_init, mask_miss, mask_hitbox_but_miss, mask_hit, mask_sec = _create_masked_particles(
        num_part,
        capacity,
        engine,
        coll_has_flanges
    )
    num_part = len(part_init.x[part_init.state > 0])
    part = part_init.copy()
    part_black = part_init.copy()

    # Track
    coll1.track(part)
    part.at_element[part.state > 0] += 1  # Need to do this manually for the impact table
    part_mid = part.copy()
    part_mid.sort(interleave_lost_particles=True)

    coll2.track(part)
    part.at_element[part.state > 0] += 1
    part.sort(interleave_lost_particles=True)

    coll1._drift(part_black, -coll1.length_front)
    black1.track(part_black)
    coll1._drift(part_black, -coll1.length_back)
    part_black.at_element[part_black.state > 0] += 1
    part_black_mid = part_black.copy()
    part_black_mid.sort(interleave_lost_particles=True)

    coll2._drift(part_black, -coll2.length_front)
    black2.track(part_black)
    coll2._drift(part_black, -coll2.length_back)
    part_black.at_element[part_black.state > 0] += 1
    part_black.sort(interleave_lost_particles=True)

    xc_engine.stop(clean=True)

    if log_impacts:
        df = impacts.to_pandas(frame='lattice')

    print(f"Particle types generated: ")
    for pdg_id, count in zip(*np.unique(part.pdg_id[(part.state > -99999) & (part.particle_id >= num_part)], return_counts=True)):
        print(f"    PDG ID {pdg_id}: {count} particles")
    print()

    # =======================================
    # === CHECKS AFTER FIRST PASS (coll1) ===
    # =======================================

    # Preliminary info
    print(f"Checks after first pass (coll1):")
    print(f"  Primary hit states: {np.unique(part_mid.state[~mask_sec])}")
    print(f"  Secondary hit states: {np.unique(part_mid.state[mask_sec])}")
    print(f"  Children generated: {((part_mid.state > -9999999) & (part_mid.particle_id >= num_part)).sum()}")
    if xcc.LOST_ON_MATERIAL not in part_mid.state:
        raise ValueError("No particles lost on material. Choose a different seed.")
    if xcc.LOST_ON_MATERIAL_SEC not in part_mid.state:
        raise ValueError("No secondary particles lost on material. Choose a different seed.")
    if xcc.VIRTUAL_ENERGY not in part_mid.state:
        raise ValueError("No virtual energy particles created. Choose a different seed.")
    if xcc.VIRTUAL_ENERGY_SEC not in part_mid.state:
        raise ValueError("No secondary virtual energy particles created. Choose a different seed.")
    if xcc.MASSLESS_OR_NEUTRAL not in part_mid.state:
        if old_bdsim:
            warn("No massless or neutral particles created.")
        else:
            raise ValueError("No massless or neutral particles created. Choose a different seed.")

    # Compare to previous result; should be independent of logging impacts or marking scattered particles
    _compare_particles(part_mid, f'temp/{engine}_part_mid.json', running_with_xdist)

    # Verify all state flags.
    # Do not use USE_IN_LOSSMAP_PRIM/SEC to check the states directly.
    primary_lost_states   = [xcc.LOST_ON_MATERIAL, xcc.VIRTUAL_ENERGY]
    kill_states           = [xcc.MASSLESS_OR_NEUTRAL, xcc.EXCITED_ION_STATE]
    secondary_lost_states = [xcc.LOST_ON_MATERIAL_SEC, xcc.VIRTUAL_ENERGY_SEC] + kill_states

    # Verify that there are no leftover hit states
    assert not np.any(part_mid.state == xcc.HIT_ON_FLUKA)
    assert not np.any(part_mid.state == xcc.HIT_ON_FLUKA_SEC)
    assert not np.any(part_mid.state == xcc.HIT_ON_GEANT4)
    assert not np.any(part_mid.state == xcc.HIT_ON_GEANT4_SEC)

    # Verify that all missed particles survived
    assert np.all(np.unique(part_mid.state[mask_miss & ~mask_sec]) == [1])
    assert np.all(np.unique(part_mid.state[mask_miss & mask_sec]) == [xcc.SECONDARY_PARTICLE])
    if coll_has_flanges:
        assert np.all(np.unique(part_mid.state[mask_hitbox_but_miss & ~mask_sec]) == [1])
        assert np.all(np.unique(part_mid.state[mask_hitbox_but_miss & mask_sec]) == [xcc.SECONDARY_PARTICLE])

    # Verify that primary particles die as primary, and survive as secondary
    # if marked (otherwise survive as primary)
    if mark_scattered_particles:
        assert np.all([ss in primary_lost_states + [xcc.SECONDARY_PARTICLE]
                        for ss in part_mid.state[mask_hit & ~mask_sec]])
    else:
        assert np.all([ss in primary_lost_states + [1]
                        for ss in part_mid.state[mask_hit & ~mask_sec]])

    # Verify that secondary particles die as secondary, and survive as secondary
    assert np.all([ss in secondary_lost_states + [xcc.SECONDARY_PARTICLE]
                    for ss in part_mid.state[mask_hit & mask_sec]])

    # Verify that parents are not flagged as massless/neutral.
    # They can be absorbed, virtual, or survived.
    mask_parents  = (part_mid.particle_id < num_part) & (part_mid.state > -99999)
    if mark_scattered_particles:
        assert np.all([ss in primary_lost_states + [1, xcc.SECONDARY_PARTICLE]
                    for ss in part_mid.state[mask_parents & ~mask_sec]])
    else:
        assert np.all([ss in primary_lost_states + [1]
                    for ss in part_mid.state[mask_parents & ~mask_sec]])
    assert np.all([ss in [xcc.VIRTUAL_ENERGY_SEC, xcc.LOST_ON_MATERIAL_SEC, xcc.SECONDARY_PARTICLE]
                for ss in part_mid.state[mask_parents & mask_sec]])

    # Verify that children always survive, unless massless/neutral
    mask_children = part_mid.particle_id >= num_part
    if mark_scattered_particles:
        assert np.all([ss in kill_states + [xcc.SECONDARY_PARTICLE]
                    for ss in part_mid.state[mask_children]])
    else:
        secondary_parents = part_init.particle_id[mask_sec]
        secondary_parents_mask = np.isin(part_mid.parent_particle_id, secondary_parents)
        assert np.all([ss in kill_states + [1]
                    for ss in part_mid.state[mask_children & ~secondary_parents_mask]])
        assert np.all([ss in kill_states + [xcc.SECONDARY_PARTICLE]
                    for ss in part_mid.state[mask_children & secondary_parents_mask]])

    # Verify the final positions
    assert np.allclose(part_mid.s[mask_miss], coll1.length)
    if coll_has_flanges:
        assert np.allclose(part_mid.s[mask_hitbox_but_miss], coll1.length)
    assert np.allclose(part_mid.s[np.isin(part_mid.state, kill_states)], coll1.length + coll1.length_back)

    # The energy of missed particles should not have changed
    energy0 = xc_engine.particle_ref.energy0[0]
    assert np.allclose(part_mid.energy[mask_miss], energy0)
    if coll_has_flanges:
        assert np.allclose(part_mid.energy[mask_hitbox_but_miss], energy0)

    # Check total summed energy
    Etot  = part_mid.energy[part_mid.state > -99999].sum()  # Should be close to initial total energy
    Etot += coll1._acc_ionisation_loss
    Etot += coll1._acc_ionisation_loss_sec
    assert np.isclose(Etot, part_init.energy[part_init.state > -99999].sum())
    assert coll1._acc_ionisation_loss > 0
    assert coll1._acc_ionisation_loss_sec > 0
    assert not np.isclose(coll1._acc_ionisation_loss, coll1._acc_ionisation_loss_sec)

    # Check the energy sum per parent->children chain, and the number of
    # outliers in the Landau/Vavilov tail.
    _check_energy_sum(part_mid, mask_hit)

    # Check the impacts
    if log_impacts:
        df_mid = df[df.collimator == coll1.name]
        assert not np.any([pid in df_mid.particle_id_before.values for pid in part_mid.particle_id[mask_miss]])
        if coll_has_flanges:
            assert np.all([pid in df_mid.particle_id_before.values for pid in part_mid.particle_id[mask_hitbox_but_miss]])
        assert np.all([pid in df_mid.particle_id_before.values for pid in part_mid.particle_id[mask_hit]])
        df_mid = df_mid.sort_values("particle_id_before")
        if coll_has_flanges:
            assert np.all(part_black_mid.particle_id[mask_hit | mask_hitbox_but_miss] == df_mid.particle_id_before.values)
            assert np.allclose(part_black_mid.s[mask_hit | mask_hitbox_but_miss],     df_mid.s_before.values)
            assert np.allclose(part_black_mid.x[mask_hit | mask_hitbox_but_miss],     df_mid.x_before.values)
            assert np.allclose(part_black_mid.px[mask_hit | mask_hitbox_but_miss],    df_mid.px_before.values)
            assert np.allclose(part_black_mid.y[mask_hit | mask_hitbox_but_miss],     df_mid.y_before.values)
            assert np.allclose(part_black_mid.py[mask_hit | mask_hitbox_but_miss],    df_mid.py_before.values)
            assert np.allclose(part_black_mid.zeta[mask_hit | mask_hitbox_but_miss],  df_mid.zeta_before.values)
            assert np.allclose(part_black_mid.delta[mask_hit | mask_hitbox_but_miss], df_mid.delta_before.values)
        else:
            assert np.all(part_black_mid.particle_id[mask_hit] == df_mid.particle_id_before.values)
            assert np.allclose(part_black_mid.s[mask_hit],     df_mid.s_before.values)
            assert np.allclose(part_black_mid.x[mask_hit],     df_mid.x_before.values)
            assert np.allclose(part_black_mid.px[mask_hit],    df_mid.px_before.values)
            assert np.allclose(part_black_mid.y[mask_hit],     df_mid.y_before.values)
            assert np.allclose(part_black_mid.py[mask_hit],    df_mid.py_before.values)
            assert np.allclose(part_black_mid.zeta[mask_hit],  df_mid.zeta_before.values)
            assert np.allclose(part_black_mid.delta[mask_hit], df_mid.delta_before.values)
    print()

    # ========================================
    # === CHECKS AFTER SECOND PASS (coll2) ===
    # ========================================

    # Preliminary info
    print(f"Checks after second pass (coll2):")
    print(f"  Primary hit states: {np.unique(part.state[~mask_sec])}")
    print(f"  Secondary hit states: {np.unique(part.state[mask_sec])}")
    print(f"  Children generated: {((part.state > -9999999) & (part.particle_id >= num_part)).sum()}")

    # Compare to previous result; should be independent of logging impacts or marking scattered particles
    _compare_particles(part, f'temp/{engine}_part.json', running_with_xdist)

    # Verify that there are no leftover hit states
    assert not np.any(part_mid.state == xcc.HIT_ON_FLUKA)
    assert not np.any(part_mid.state == xcc.HIT_ON_FLUKA_SEC)
    assert not np.any(part_mid.state == xcc.HIT_ON_GEANT4)
    assert not np.any(part_mid.state == xcc.HIT_ON_GEANT4_SEC)

    # Verify state flags on the different collimators
    mask_coll1 = part.at_element == 0
    mask_coll2 = part.at_element == 1
    mask_surv  = part.at_element == 2
    mask_original = part.particle_id < num_part
    _ang = np.radians(coll2.angle)
    mask_hit_coll2  = part_mid.y >= np.tan(_ang - np.pi/2) * (part_mid.x - coll2.jaw_L*np.cos(_ang)) + coll2.jaw_L*np.sin(_ang)
    mask_hit_coll2 |= part_mid.y <= np.tan(_ang - np.pi/2) * (part_mid.x - coll2.jaw_R*np.cos(_ang)) + coll2.jaw_R*np.sin(_ang)
    mask_hit_coll2 |= part_mid.y + 0.6*part_mid.py >= np.tan(_ang - np.pi/2) * (part_mid.x + 0.6*part_mid.px - coll2.jaw_L*np.cos(_ang)) + coll2.jaw_L*np.sin(_ang)
    mask_hit_coll2 |= part_mid.y + 0.6*part_mid.py <= np.tan(_ang - np.pi/2) * (part_mid.x + 0.6*part_mid.px - coll2.jaw_R*np.cos(_ang)) + coll2.jaw_R*np.sin(_ang)

    # Particles flagged as primary that died on the first collimator should have primary states or massless/neutral
    assert np.all([ss in primary_lost_states for ss in part.state[mask_original & mask_coll1 & ~mask_sec]])
    # Particles flagged as secondary that died on the first collimator should have secondary states
    assert np.all([ss in secondary_lost_states for ss in part.state[mask_original & mask_coll1 & mask_sec]])

    # Particles flagged as primary that did not hit the first collimator but died on the second should have primary states or massless/neutral
    assert np.all([ss in primary_lost_states for ss in part.state[mask_original & mask_coll2 & ~mask_hit & ~mask_sec]])
    # Particles flagged as secondary that did not hit the first collimator but died on the second should have secondary states
    assert np.all([ss in secondary_lost_states for ss in part.state[mask_original & mask_coll2 & ~mask_hit & mask_sec]])
    # Particles flagged as primary that scattered on the first collimator should be flagged as secondary only if requested
    if mark_scattered_particles:
        assert np.all([ss in secondary_lost_states for ss in part.state[mask_coll2 & mask_hit & ~mask_sec]])
    else:
        assert np.all([ss in primary_lost_states + kill_states for ss in part.state[mask_coll2 & mask_hit & ~mask_sec]])
    # Particles flagged as secondary that scattered on the first collimator should be flagged as secondary
    assert np.all([ss in secondary_lost_states for ss in part.state[mask_coll2 & mask_hit & mask_sec]])

    # Particles flagged as primary that survived without hitting any collimator should be flagged as primary
    assert np.unique(part.state[mask_surv & mask_miss & ~mask_hit_coll2 & ~mask_sec]) == [1]
    # Particles flagged as secondary that survived without hitting any collimator should be flagged as secondary
    assert np.unique(part.state[mask_surv & mask_miss & ~mask_hit_coll2 & mask_sec]) == [xcc.SECONDARY_PARTICLE]
    # Particles flagged as primary that survived after hitting any collimator should be flagged as secondary only if requested
    if mark_scattered_particles:
        assert np.unique(part.state[mask_surv & (mask_hit | mask_hit_coll2) & ~mask_sec]) == [xcc.SECONDARY_PARTICLE]
    else:
        # Children of particles flagged as secondary will ALWAYS be flagged as secondary,
        # even if marking is not requested.
        secondary_parents = part_init.particle_id[mask_sec]
        secondary_parents_mask = np.isin(part.parent_particle_id, secondary_parents)
        secondary_parents_mask |= np.isin(part.parent_particle_id, part.particle_id[secondary_parents_mask])  # Also mark grandchildren
        assert np.unique(part.state[mask_surv & (mask_hit | mask_hit_coll2) & secondary_parents_mask]) == [xcc.SECONDARY_PARTICLE]
        assert np.unique(part.state[mask_surv & (mask_hit | mask_hit_coll2) & ~mask_sec & ~secondary_parents_mask]) == [1]
    # Particles flagged as secondary that survived after hitting any collimator should be flagged as secondary
    assert np.unique(part.state[mask_surv & (mask_hit | mask_hit_coll2) & mask_sec]) == [xcc.SECONDARY_PARTICLE]

    # Verify the final positions
    assert np.allclose(part.s[mask_surv], coll1.length + coll2.length)
    assert np.allclose(part.s[np.isin(part.state, kill_states) & (part.at_element == 0)], coll1.length + coll1.length_back)
    assert np.allclose(part.s[np.isin(part.state, kill_states) & (part.at_element == 1)], coll1.length + coll2.length + coll2.length_back)

    # The energy of missed particles should not have changed
    assert np.allclose(part.energy[mask_miss & ~mask_hit_coll2], energy0)

    # Check total summed energy
    Etot  = part.energy[part.state > -99999].sum()  # Should be close to initial total energy
    Etot += coll1._acc_ionisation_loss
    Etot += coll1._acc_ionisation_loss_sec
    Etot += coll2._acc_ionisation_loss
    Etot += coll2._acc_ionisation_loss_sec
    assert np.isclose(Etot, part_init.energy[part_init.state > -99999].sum(), atol=1e-12)
    assert coll1._acc_ionisation_loss > 0
    assert coll1._acc_ionisation_loss_sec > 0
    assert not np.isclose(coll1._acc_ionisation_loss, coll1._acc_ionisation_loss_sec)
    assert coll2._acc_ionisation_loss > 0
    assert coll2._acc_ionisation_loss_sec > 0
    assert not np.isclose(coll2._acc_ionisation_loss, coll2._acc_ionisation_loss_sec)

    # Check the energy sum per parent->children chain, and the number of
    # outliers in the Landau/Vavilov tail.
    _check_energy_sum(part, mask_hit)

    # Check the impacts
    if log_impacts:
        df_end = df[df.collimator == coll2.name]
        df_end = df_end.sort_values("particle_id_before")
        # Only compare particles that missed the first collimator and hit the second one
        mask_end = (part_black.at_element == 1) & (part_black.state < 0)
        mask_df = np.isin(df_end.particle_id_before.values, part_black.particle_id[mask_end])
        assert np.all(part_black.particle_id[mask_end] == df_end.particle_id_before.values[mask_df])
        assert np.allclose(part_black.s[mask_end],     df_end.s_before.values[mask_df] + coll1.length)
        assert np.allclose(part_black.x[mask_end],     df_end.x_before.values[mask_df])
        assert np.allclose(part_black.px[mask_end],    df_end.px_before.values[mask_df])
        assert np.allclose(part_black.y[mask_end],     df_end.y_before.values[mask_df])
        assert np.allclose(part_black.py[mask_end],    df_end.py_before.values[mask_df])
        assert np.allclose(part_black.zeta[mask_end],  df_end.zeta_before.values[mask_df])
        assert np.allclose(part_black.delta[mask_end], df_end.delta_before.values[mask_df])


def _compare_particles(part, file, running_with_xdist):
    if running_with_xdist:
        warn("Not comparing to previous result (running with xdist).")
        return

    file = Path(file)
    file.parent.mkdir(parents=True, exist_ok=True)
    if file.exists():
        dct = xc.json.json_load(file)
        if time.time() - dct['time'] < 600:
            # Only use recent results for comparison (to increase reproducibility)
            print("Comparing to previous result...")
            part2 = xt.Particles.from_dict(dct['data'])
            part2.sort(interleave_lost_particles=True)
            assert np.allclose(part.s,  part2.s)
            assert np.allclose(part.x,  part2.x)
            assert np.allclose(part.px, part2.px)
            assert np.allclose(part.y,  part2.y)
            assert np.allclose(part.py, part2.py)
            assert np.allclose(part.zeta,  part2.zeta)
            assert np.allclose(part.delta, part2.delta)
            assert np.array_equal(np.isin(part.state,  [1, xcc.SECONDARY_PARTICLE]),
                                  np.isin(part2.state, [1, xcc.SECONDARY_PARTICLE]))
            assert np.array_equal(np.isin(part.state,  [xcc.LOST_ON_MATERIAL, xcc.LOST_ON_MATERIAL_SEC]),
                                  np.isin(part2.state, [xcc.LOST_ON_MATERIAL, xcc.LOST_ON_MATERIAL_SEC]))
            assert np.array_equal(np.isin(part.state,  [xcc.VIRTUAL_ENERGY, xcc.VIRTUAL_ENERGY_SEC]),
                                  np.isin(part2.state, [xcc.VIRTUAL_ENERGY, xcc.VIRTUAL_ENERGY_SEC]))
            assert np.array_equal(np.isin(part.state,  [xcc.MASSLESS_OR_NEUTRAL]),
                                  np.isin(part2.state, [xcc.MASSLESS_OR_NEUTRAL]))
            assert np.array_equal(np.isin(part.state,  [xcc.EXCITED_ION_STATE]),
                                  np.isin(part2.state, [xcc.EXCITED_ION_STATE]))
            return

    print("Storing result for future comparison")
    dct = {'time': time.time(), 'data': part.to_dict()}
    xc.json.json_dump(dct, file)


def _create_masked_particles(num_part, capacity, engine, coll_has_flanges):
    # When this is changed, need to re-generate input distribution
    if engine == 'fluka':
        xc_engine = xc.fluka.engine
    elif engine == 'geant4':
        xc_engine = xc.geant4.engine
    num_steps_miss = 2   # *2 for primary/secondary
    if coll_has_flanges:
        num_steps_hitbox_but_miss = 2   # *2 for primary/secondary
    else:
        num_steps_hitbox_but_miss = 0
    num_steps_hit = 4   # *2 for primary/secondary
    num_steps = 2*(num_steps_miss + num_steps_hitbox_but_miss + num_steps_hit)
    step_size = num_part//num_steps
    print(f"{num_part=}  {capacity=}")
    print(f"{num_steps} steps of {step_size} particles each, for a total of {num_steps*step_size} particles.")

    # Masks
    mask_miss = np.concatenate([
        np.full(2*num_steps_miss*step_size, True),
        np.full(2*num_steps_hitbox_but_miss*step_size, False),
        np.full(2*num_steps_hit*step_size, False),
        np.full(capacity - num_steps*step_size, False)
    ])
    mask_hitbox_but_miss = np.concatenate([
        np.full(2*num_steps_miss*step_size, False),
        np.full(2*num_steps_hitbox_but_miss*step_size, True),
        np.full(2*num_steps_hit*step_size, False),
        np.full(capacity - num_steps*step_size, False)
    ])
    mask_hit = np.concatenate([
        np.full(2*num_steps_miss*step_size, False),
        np.full(2*num_steps_hitbox_but_miss*step_size, False),
        np.full(2*num_steps_hit*step_size, True),
        np.full(capacity - num_steps*step_size, False)
    ])
    mask_sec = np.concatenate([
        np.full(num_steps_miss*step_size, False),
        np.full(num_steps_miss*step_size, True),
        np.full(num_steps_hitbox_but_miss*step_size, False),
        np.full(num_steps_hitbox_but_miss*step_size, True),
        np.full(num_steps_hit*step_size, False),
        np.full(num_steps_hit*step_size, True),
        np.full(capacity - num_steps*step_size, False)
    ])

    # Coordinates
    init_file = Path(f'data/{engine}_part_init.json')
    if init_file.exists():
        part = xt.Particles.from_dict(xc.json.json_load(init_file))

    else:
        x_miss = np.linspace(-0.999e-3, 0.999e-3, num_steps_miss*step_size)
        px_miss = np.zeros(num_steps_miss*step_size)
        if coll_has_flanges:
            x_hitbox_but_miss = np.concatenate([
                -0.99e-3*np.ones(step_size),
                0.99e-3*np.ones(step_size)
            ])
            px_hitbox_but_miss = np.concatenate([
                np.linspace(3e-5, 4e-5, step_size),
                np.linspace(-4e-5, -3e-5, step_size)
            ])
        else:
            x_hitbox_but_miss = np.array([])
            px_hitbox_but_miss = np.array([])
        x_hit  = np.concatenate([
            np.linspace(1.001e-3, 2e-3, step_size),
            np.linspace(-2e-3, -1.001e-3, step_size),
            np.linspace(0.9e-3, 0.99e-3, step_size),
            np.linspace(-0.99e-3, -0.9e-3, step_size)
        ])
        px_hit = np.concatenate([
            np.zeros(2*step_size),
            np.linspace(2e-4, 1e-3, step_size),
            np.linspace(-1e-3, -2e-4, step_size)
        ])

        # Sanity checks
        assert len(x_miss) == num_steps_miss*step_size
        assert len(x_hitbox_but_miss) == num_steps_hitbox_but_miss*step_size
        assert len(x_hit) == num_steps_hit*step_size

        part = xp.build_particles(
            x=np.concatenate([
                x_miss, x_miss,
                x_hitbox_but_miss, x_hitbox_but_miss,
                x_hit, x_hit
            ]),
            px=np.concatenate([
                px_miss, px_miss,
                px_hitbox_but_miss, px_hitbox_but_miss,
                px_hit, px_hit
            ]),
            y=np.linspace(-1e-6, 1e-6, step_size*num_steps),
            py=np.linspace(-1e-7, 1e-7, step_size*num_steps),
            particle_ref=xc_engine.particle_ref,
            _capacity=capacity)
        part.state[mask_sec] = xcc.SECONDARY_PARTICLE  # Mark secondary particles in initial distribution
        xc.json.json_dump(part.to_dict(), init_file)

    return part, mask_miss, mask_hitbox_but_miss, mask_hit, mask_sec


def _get_num_coll_traversed(part, pids):
    # Check how many collimators are traversed by all particles
    # (parent and children) without dying. This gives us an estimate
    # of the total traversed length, and hence for the amount of
    # "missing" energy allowed.
    num_coll_traversed = 0
    for this_pid in pids:
        at_element = part.at_element[part.particle_id == this_pid][0]
        this_ppid = part.parent_particle_id[part.particle_id == this_pid][0]
        if this_ppid == this_pid:
            # Primary
            num_coll_traversed += at_element
        else:
            # Child
            parent_at_element = part.at_element[part.particle_id == this_ppid][0]
            if parent_at_element < at_element:
                num_coll_traversed += at_element - parent_at_element - 1
    return num_coll_traversed


def _check_energy_sum(part, mask):
    # Check the sum of the children energy and leftover energy of the parent.
    # This sum will not match exactly the initial energy, as ionisation losses
    # are not accounted for on a particle-by-particle basis (only accumulated
    # in the collimator). When a particle dies in a collimator, all ionisation
    # losses in that collimator are automatically added to the dead particle's
    # energy. So ionisation losses are not individually accounted for ONLY when
    # a particle survives a collimator.
    print("Checking energy sum for each parent->children chain...")
    tot_part = len(part.particle_id[mask])
    energy0 = part.energy0[0]
    tree = xc.ParticlesTree(part)
    # For the first loop, we only check how many primary particles actually
    # caused ionisation losses, to be able to get our statistics right.
    for pid, e_parent in zip(part.particle_id[mask], part.energy[mask]):
        des = tree.descendants_ids(pid)
        num_coll_traversed = _get_num_coll_traversed(part, [pid, *des])
        if num_coll_traversed == 0:
            # This parent and its children did not traverse+survive any
            # collimator, so remove it from the statistics.
            tot_part -= 1
    print(f"Total number of primary particles that traversed at least one collimator: {tot_part}")
    # For the second loop, we check for each parent->children chain that the
    # missing energy is statistically compatible with the Landau/Vavilov
    # distribution of ionisation losses, plus extras. This is benchmarked
    # earlier for this material (see first test above). In particular, we check
    # that it is lower than the Landau/Vavilov tail, estimated probabilistically
    # in such a way that we expect a given number of outliers out of all particles.
    # Note that we cannot check the lower tail, as not all particles traverse
    # the full collimator length; there might hence be zero ionisation losses.
    # Tail: P(dE > xi) = xi/T  =>  T_threshold = xi/P
    #       1 particle  above this threshold  =>  P = 1/N
    #      10 particles above this threshold  =>  P = 10/N
    #     100 particles above this threshold  =>  P = 100/N
    # Thresholds from file
    data = xc.json.json_load('data/fluka_calibration_data.json')
    dE_threshold_1   = data['full']['thresholds']['1']
    dE_threshold_10  = data['full']['thresholds']['10']
    dE_threshold_100 = data['full']['thresholds']['100']
    num_outliers_1   = 0
    num_outliers_10  = 0
    num_outliers_100 = 0
    num_zero_ionisation_losses = 0
    loss_per_coll = {1: [], 2: []}
    for pid, e_parent in zip(part.particle_id[mask], part.energy[mask]):
        des = tree.descendants_ids(pid)
        energy_children = part.energy[np.isin(part.particle_id, des)]
        this_diff = (energy0 - energy_children.sum() - e_parent)
        num_coll_traversed = _get_num_coll_traversed(part, [pid, *des])
        if num_coll_traversed == 0:
            # In this case the energy sum needs to match the initial energy exactly
            # if len(des) > 0:
            assert np.isclose(this_diff/energy0, 0, atol=1.e-12)
        else:
            # In this case, we have missing energy due to ionisation losses,
            # which should be compatible with the Landau/Vavilov distribution
            assert np.all(this_diff/energy0 >= -1.e-12)  # No energy should be created
            if np.isclose(this_diff/energy0, 0, atol=1.e-12):
                num_zero_ionisation_losses += 1
            if num_coll_traversed not in loss_per_coll:
                loss_per_coll[num_coll_traversed] = []
            loss_per_coll[num_coll_traversed].append(this_diff)
            if this_diff > dE_threshold_1*num_coll_traversed:
                num_outliers_1 += 1
            if this_diff > dE_threshold_10*num_coll_traversed:
                num_outliers_10 += 1
            if this_diff > dE_threshold_100*num_coll_traversed:
                num_outliers_100 += 1
    print(f"Traversed particles without ionisation losses: {num_zero_ionisation_losses}")
    print(f"Ionisation losses statistics:")
    print(f"    Number of outliers above 1 particle threshold: {num_outliers_1}")
    print(f"    Number of outliers above 10 particle threshold: {num_outliers_10}")
    print(f"    Number of outliers above 100 particle threshold: {num_outliers_100}")
    print("Ionisation losses quantiles:")
    for kkk, lll in loss_per_coll.items():
        if len(lll) > 0:
            print(f"  num_coll_traversed = {kkk}:")
            for n in [1, 10, 100]:
                q = 1 - n / tot_part
                print(f"    Quantile for {n} particles: {np.quantile(np.asarray(lll), q) / 1e6} MeV")
    # _plot_ionisation_losses(loss_per_coll)
    assert num_outliers_1 < 3      # Expect maximally 1 outlier
    assert num_outliers_10 < 15    # Expect maximally 10 outliers
    assert num_outliers_100 < 120  # Expect maximally 100 outliers
    assert num_zero_ionisation_losses < 100 * len(part.x[mask])/10_000  # These are corner hits


def _plot_ionisation_losses(loss_per_coll):
    import matplotlib.pyplot as plt
    nbins = 250
    pos_loss_per_coll = {
        ncoll: np.asarray(losses)[np.asarray(losses) > 0]
        for ncoll, losses in loss_per_coll.items()
    }
    E_min  = min([min(ll) for ll in pos_loss_per_coll.values() if len(ll) > 0])
    E_high = max([max(ll) for ll in pos_loss_per_coll.values() if len(ll) > 0])
    bins = np.logspace(np.log10(E_min), np.log10(E_high), nbins + 1)
    bin_centres = np.sqrt(bins[:-1] * bins[1:])
    dlog = np.diff(np.log10(bins))
    _, ax = plt.subplots(2, 1, figsize=(6, 8))
    for ii, (ncoll, losses) in enumerate(loss_per_coll.items()):
        losses = np.asarray(losses)
        n_zero = np.count_nonzero(losses <= 0)
        pos_losses = losses[losses > 0]
        if len(pos_losses) > 0:
            counts, _ = np.histogram(pos_losses, bins=bins)
            dNdlogE = counts / (len(losses) * dlog)
            ax[ii].step(bin_centres, dNdlogE, where='mid')
            ax[ii].set_xscale('log')
            ax[ii].set_yscale('log')
            ax[ii].set_xlabel('Energy [eV]')
            ax[ii].set_ylabel(r'Normalised frequency $\frac{dN}{d\log E}$')
            ax[ii].set_title(
                f"{ncoll} collimator(s), {n_zero} zero-loss chains"
            )
            ax[ii].grid(True, which='both', linestyle='--', linewidth=0.5)
    plt.tight_layout()
    plt.show()
