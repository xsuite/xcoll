# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2026.                 #
# ######################################### #

import time
import pytest
import numpy as np
from pathlib import Path
from warnings import warn

import xpart as xp
import xtrack as xt
import xcoll as xc
import xcoll.constants as xcc


@pytest.mark.geant4
@pytest.mark.parametrize("log_impacts, mark_scattered_particles", [
                            [False, False],
                            [False, True],
                            [True,  False],
                            [True,  True]
                         ], ids=["default", "mark", "impacts", "impacts_mark"])
def test_geant4_deep_check(log_impacts, mark_scattered_particles, running_with_xdist):
    num_part = 50000       # When this is changed, need to re-generate input distribution
    capacity = 4*num_part  # When this is changed, need to re-generate input distribution
    particle_ref = xt.Particles('proton', p0c=6.8e12)

    # Prepare collimators, engine, and initial particles
    coll1 = xc.Geant4Collimator(length=0.6, angle=0,   jaw=0.001,  material=xc.materials.MolybdenumGraphite)
    coll2 = xc.Geant4Collimator(length=0.6, angle=123, jaw=0.0005, material=xc.materials.MolybdenumGraphite)
    if mark_scattered_particles:
        coll1.mark_scattered_particles = True
        coll2.mark_scattered_particles = True
    else:
        coll1.mark_scattered_particles = False
        coll2.mark_scattered_particles = False

    xc.geant4.engine.particle_ref = particle_ref
    xc.geant4.engine.seed = 453532
    xc.geant4.engine.return_baryons = True   # To get some massless particles as well
    xc.geant4.engine.start(elements=[coll1, coll2], clean=True, verbose=False)

    # Create black absorbers after starting engine so length_front and length_back are known
    black1 = xc.BlackAbsorber(length=0.6 + coll1.length_front + coll1.length_back, angle=0,   jaw=0.001)
    black2 = xc.BlackAbsorber(length=0.6 + coll1.length_front + coll1.length_back, angle=123, jaw=0.0005)

    if log_impacts:
        impacts = xc.InteractionRecord.start(elements=[coll1, coll2], record_impacts=True)

    part_init, mask_miss, mask_hit, mask_sec = _create_masked_particles(num_part, capacity)
    part = part_init.copy()
    part_black = part_init.copy()

    # Track
    coll1.track(part)
    part_mid = part.copy()
    part_mid.sort(interleave_lost_particles=True)

    part.at_element[part.state > 0] += 1  # Need to do this manually for the impact table
    coll2.track(part)
    part.sort(interleave_lost_particles=True)
    part.at_element[part.state > 0] += 1

    coll1._drift(part_black, -coll1.length_front)
    black1.track(part_black)
    coll1._drift(part_black, -coll1.length_back)
    part_black_mid = part_black.copy()
    part_black_mid.sort(interleave_lost_particles=True)

    coll2._drift(part_black, -coll2.length_front)
    part_black.at_element[part_black.state > 0] += 1
    black2.track(part_black)
    coll2._drift(part_black, -coll2.length_back)
    part_black.sort(interleave_lost_particles=True)
    part_black.at_element[part_black.state > 0] += 1

    xc.geant4.engine.stop(clean=True)

    if log_impacts:
        df = impacts.to_pandas(frame='lattice')

    print(f"Particle types generated: ")
    for pdg_id, count in zip(*np.unique(part.pdg_id[(part.state > -99999) & (part.particle_id >= num_part)], return_counts=True)):
        print(f"    PDG ID {pdg_id}: {count} particles")

    # =======================================
    # === CHECKS AFTER FIRST PASS (coll1) ===
    # =======================================

    # Preliminary info
    print(f"Primary hit states: {np.unique(part_mid.state[~mask_sec])}")
    print(f"Secondary hit states: {np.unique(part_mid.state[mask_sec])}")
    print(f"Children generated: {((part_mid.state > -9999999) & (part_mid.particle_id >= num_part)).sum()}")
    if xcc.LOST_ON_MATERIAL not in part_mid.state:
        raise ValueError("No particles lost on material. Choose a different seed.")
    if xcc.LOST_ON_MATERIAL_SEC not in part_mid.state:
        raise ValueError("No secondary particles lost on material. Choose a different seed.")
    if xcc.VIRTUAL_ENERGY not in part_mid.state:
        raise ValueError("No virtual energy particles created. Choose a different seed.")
    if xcc.VIRTUAL_ENERGY_SEC not in part_mid.state:
        raise ValueError("No secondary virtual energy particles created. Choose a different seed.")
    if xcc.MASSLESS_OR_NEUTRAL not in part_mid.state:
        warn("No massless or neutral particles created.")

    # Compare to previous result; should be independent of logging impacts or marking scattered particles
    _compare_particles(part_mid, 'temp_geant4_part_mid.json', running_with_xdist)

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
    assert np.allclose(part_mid.s[np.isin(part_mid.state, kill_states)], coll1.length + coll1.length_back)

    # The energy of missed particles should not have changed
    energy0 = xc.geant4.engine.particle_ref.energy0[0]
    assert np.allclose(part_mid.energy[mask_miss], energy0)

    # Check total summed energy
    Etot  = part_mid.energy[part_mid.state > -99999].sum()  # Should be close to initial total energy
    Etot += coll1._acc_ionisation_loss
    Etot += coll1._acc_ionisation_loss_sec
    assert np.isclose(Etot, part_init.energy[part_init.state > -99999].sum())
    assert coll1._acc_ionisation_loss > 0
    assert coll1._acc_ionisation_loss_sec > 0
    assert not np.isclose(coll1._acc_ionisation_loss, coll1._acc_ionisation_loss_sec)

    # Check that the sum of the children energy and leftover energy of the parent is close to the initial energy
    tree = xc.ParticlesTree(part_mid)
    for pid, e_parent in zip(part_mid.particle_id[mask_hit], part_mid.energy[mask_hit]):
        des = tree.descendants_ids(pid)
        if len(des) > 0:
            energy_children = part_mid.energy[np.isin(part_mid.particle_id, des)]
            assert np.isclose(energy_children.sum() + e_parent - energy0, 0., atol=1e-2)

    # Check the impacts
    if log_impacts:
        df_mid = df[df.collimator == coll1.name]
        assert not np.any([pid in df_mid.id_before.values for pid in part_mid.particle_id[mask_miss]])
        assert np.all([pid in df_mid.id_before.values for pid in part_mid.particle_id[mask_hit]])
        df_mid = df_mid.sort_values("id_before")
        assert np.all(part_black_mid.particle_id[mask_hit] == df_mid.id_before.values)
        assert np.allclose(part_black_mid.s[mask_hit],     df_mid.s_before.values)
        assert np.allclose(part_black_mid.x[mask_hit],     df_mid.x_before.values)
        assert np.allclose(part_black_mid.px[mask_hit],    df_mid.px_before.values)
        assert np.allclose(part_black_mid.y[mask_hit],     df_mid.y_before.values)
        assert np.allclose(part_black_mid.py[mask_hit],    df_mid.py_before.values)
        assert np.allclose(part_black_mid.zeta[mask_hit],  df_mid.zeta_before.values)
        assert np.allclose(part_black_mid.delta[mask_hit], df_mid.delta_before.values)

    # ========================================
    # === CHECKS AFTER SECOND PASS (coll2) ===
    # ========================================

    # Preliminary info
    print(f"Primary hit states: {np.unique(part.state[~mask_sec])}")
    print(f"Secondary hit states: {np.unique(part.state[mask_sec])}")
    print(f"Children generated: {((part.state > -9999999) & (part.particle_id >= num_part)).sum()}")

    # Compare to previous result; should be independent of logging impacts or marking scattered particles
    _compare_particles(part, 'temp_geant4_part.json', running_with_xdist)

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
    assert np.isclose(Etot, part_init.energy[part_init.state > -99999].sum())
    assert coll1._acc_ionisation_loss > 0
    assert coll1._acc_ionisation_loss_sec > 0
    assert not np.isclose(coll1._acc_ionisation_loss, coll1._acc_ionisation_loss_sec)
    assert coll2._acc_ionisation_loss > 0
    assert coll2._acc_ionisation_loss_sec > 0
    assert not np.isclose(coll2._acc_ionisation_loss, coll2._acc_ionisation_loss_sec)

    # Check that the sum of the children energy and leftover energy of the parent is close to the initial energy
    # Difference of 1e-2 is allowed (empirical) to allow for ionisation losses
    tree = xc.ParticlesTree(part)
    for pid, e_parent in zip(part.particle_id[mask_hit], part.energy[mask_hit]):
        des = tree.descendants_ids(pid)
        if len(des) > 0:
            energy_children = part.energy[np.isin(part.particle_id, des)]
            assert np.isclose((energy0 - energy_children.sum() - e_parent)/energy0, 0., atol=1e-2)

    # Check the impacts
    if log_impacts:
        df_end = df[df.collimator == coll2.name]
        df_end = df_end.sort_values("id_before")
        # Only compare particles that missed the first collimator and hit the second one
        mask_end = (part_black.at_element == 1) & (part_black.state < 0)
        mask_df = np.isin(df_end.id_before.values, part_black.particle_id[mask_end])
        assert np.all(part_black.particle_id[mask_end] == df_end.id_before.values[mask_df])
        assert np.allclose(part_black.s[mask_end],     df_end.s_before.values[mask_df] + coll1.length)
        assert np.allclose(part_black.x[mask_end],     df_end.x_before.values[mask_df])
        assert np.allclose(part_black.px[mask_end],    df_end.px_before.values[mask_df])
        assert np.allclose(part_black.y[mask_end],     df_end.y_before.values[mask_df])
        assert np.allclose(part_black.py[mask_end],    df_end.py_before.values[mask_df])
        assert np.allclose(part_black.zeta[mask_end],  df_end.zeta_before.values[mask_df])
        assert np.allclose(part_black.delta[mask_end], df_end.delta_before.values[mask_df])


def _compare_particles(part, file, running_with_xdist):
    if running_with_xdist:
        warn("Not comparing to previous result since running with xdist.")
        return

    file = Path(file)
    if file.exists():
        dct = xc.json.json_load(file)
        if time.time() - dct['time'] < 7200:
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

    # Otherwise, save the current result for future comparison
    dct = {'time': time.time(),
           'data': part.to_dict()}
    xc.json.json_dump(dct, file)


def _create_masked_particles(num_part, capacity):
    # When this is changed, need to re-generate input distribution
    step_size = num_part//10
    mask_miss = np.concat([np.full(step_size, True),    np.full(step_size, True),
                           np.full(4*step_size, False), np.full(4*step_size, False),
                           np.full(capacity - 10*step_size, False)])
    mask_hit = np.concat([np.full(step_size, False),   np.full(step_size, False),
                          np.full(4*step_size, True),  np.full(4*step_size, True),
                          np.full(capacity - 10*step_size, False)])
    mask_sec = np.concat([np.full(step_size, False),   np.full(step_size, True),
                          np.full(4*step_size, False), np.full(4*step_size, True),
                          np.full(capacity - 10*step_size, False)])

    init_file = Path('data/geant4_part_init.json')
    if init_file.exists():
        part_init = xt.Particles.from_dict(xc.json.json_load(init_file))
    else:
        x_miss  = np.linspace(-0.999e-3, 0.999e-3, step_size)
        px_miss = np.zeros(step_size)
        x_hit  = np.concat([np.linspace(1.001e-3, 2e-3, step_size),
                            np.linspace(-2e-3, -1.001e-3, step_size),
                            np.linspace(0.9e-3, 0.99e-3, step_size),
                            np.linspace(-0.99e-3, -0.9e-3, step_size)])
        px_hit = np.concat([np.zeros(2*step_size),
                            np.linspace(2e-4, 1e-3, step_size),
                            np.linspace(-1e-3, -2e-4, step_size)])
        x_miss_sec  = x_miss
        px_miss_sec = px_miss
        x_hit_sec  = x_hit
        px_hit_sec = px_hit
        part_init = xp.build_particles(
            x=np.concat([x_miss, x_miss_sec, x_hit, x_hit_sec]),
            px=np.concat([px_miss, px_miss_sec, px_hit, px_hit_sec]),
            y=np.linspace(-1e-6, 1e-6, step_size*10),
            py=np.linspace(-1e-7, 1e-7, step_size*10),
            particle_ref=xc.geant4.engine.particle_ref,
            _capacity=capacity)
        part_init.state[mask_sec] = xcc.SECONDARY_PARTICLE  # Mark secondary particles in initial distribution
        xc.json.json_dump(part_init.to_dict(), init_file)

    return part_init, mask_miss, mask_hit, mask_sec
