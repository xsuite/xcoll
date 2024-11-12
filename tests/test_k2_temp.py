# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

import numpy as np

import xpart as xp
import xcoll as xc
from xcoll.scattering_routines.k2 import K2Engine


_ACCURACY = 1e-8  # Anything in this region around the jaw might or might not hit; we can't be sure
num_part_step = 50000

def test_k2_position():
    coll = xc.beam_elements._K2Collimator(length=1, jaw=0.01, material='C', emittance=3.5e-6)
    K2Engine().start(elements=coll, particle_ref='proton', p0c=400e9, _capacity=num_part_step*5)
    x = np.random.uniform(-0.1, 0.1, num_part_step)
    jaw_step = 1.e-7
    x = np.concatenate([x, np.random.uniform(coll.jaw_L - jaw_step, coll.jaw_L -_ACCURACY, num_part_step)])
    x = np.concatenate([x, np.random.uniform(coll.jaw_L +_ACCURACY, coll.jaw_L + jaw_step, num_part_step)])
    x = np.concatenate([x, np.random.uniform(coll.jaw_R - jaw_step, coll.jaw_R -_ACCURACY, num_part_step)])
    x = np.concatenate([x, np.random.uniform(coll.jaw_R +_ACCURACY, coll.jaw_R + jaw_step, num_part_step)])
    y = np.random.uniform(-0.1, 0.1, 5*num_part_step)
    part = xp.build_particles(x=x, y=y, particle_ref=K2Engine().particle_ref)
    part_init = part.copy()
    mask = (part_init.x > coll.jaw_R) & (part_init.x < coll.jaw_L)
    not_hit_ids = part_init.particle_id[mask]
    hit_ids = part_init.particle_id[~mask]
    coll.track(part)
    mask_hit = np.isin(part.particle_id, hit_ids)
    mask_not_hit = np.isin(part.particle_id, not_hit_ids)
    print(part.px[mask_not_hit][abs(part.px[mask_not_hit]) > 1.e-12])
    print(len(part.px[mask_not_hit][abs(part.px[mask_not_hit]) > 1.e-12]))
    print(len(part.px[mask_not_hit]))
    assert np.allclose(part.px[mask_not_hit], 0., atol=1e-12)
    assert not np.any(part.state[mask_not_hit] < 1)
    
    px_hit = part.px[[pp in hit_ids for pp in part.particle_id]]
    print(len(px_hit[np.abs(px_hit) < 1e-12]))
    K2Engine().stop()



import numpy as np
import xpart as xp
import xcoll as xc
from xcoll.scattering_routines.k2 import K2Engine
import matplotlib.pyplot as plt
_ACCURACY = 0
num_part_step = 50000
jaws = []
spread_outside_L = []
spread_outside_R = []
spread_inside_L = []
spread_inside_R = []
for jaw in [10**i for i in np.linspace(-1, -6, 50)]:
    for ran in np.random.uniform(-5.e-10, 5e-10, 2):
        jaws.append(jaw+ran)
        coll = xc.beam_elements._K2Collimator(length=1.234578623, jaw=jaws[-1], material='C', emittance=3.5e-6)
        K2Engine().start(elements=coll, particle_ref='proton', p0c=400e9, _capacity=num_part_step*5)
        x = np.random.uniform(-0.2, 0.2, num_part_step)
        jaw_step = 1.e-7
        x = np.concatenate([x, np.random.uniform(coll.jaw_L - jaw_step, coll.jaw_L -_ACCURACY, num_part_step)])
        x = np.concatenate([x, np.random.uniform(coll.jaw_L +_ACCURACY, coll.jaw_L + jaw_step, num_part_step)])
        x = np.concatenate([x, np.random.uniform(coll.jaw_R - jaw_step, coll.jaw_R -_ACCURACY, num_part_step)])
        x = np.concatenate([x, np.random.uniform(coll.jaw_R +_ACCURACY, coll.jaw_R + jaw_step, num_part_step)])
        y = np.random.uniform(-0.2, 0.2, 5*num_part_step)
        part = xp.build_particles(x=x, y=y, particle_ref=K2Engine().particle_ref)
        part_init = part.copy()
        mask = (part_init.x >= coll.jaw_L) | (part_init.x <= coll.jaw_R)
        hit_ids = part_init.particle_id[mask]
        not_hit_ids = part_init.particle_id[~mask]
        coll.track(part)
        mask_hit = np.isin(part.particle_id, hit_ids)
        mask_not_hit = np.isin(part.particle_id, not_hit_ids)
        # Particles that are supposed to not have hit the collimator, but have a kick, are considered faulty
        faulty = abs(part.px[mask_not_hit]) > 1.e-12
        fids = part.particle_id[mask_not_hit][faulty]
        faulty_outside = part_init.x[np.isin(part_init.particle_id, fids)]
        faulty_outside_y = part_init.y[np.isin(part_init.particle_id, fids)]
        if len(faulty_outside[faulty_outside > 0]) > 0:
            spread_outside_L.append(coll.jaw_L - faulty_outside[faulty_outside > 0].min())
        else:
            spread_outside_L.append(0)
        if len(faulty_outside[faulty_outside < 0]) > 0:
            spread_outside_R.append(faulty_outside[faulty_outside < 0].max() - coll.jaw_R)
        else:
            spread_outside_R.append(0)
        # Particles that are supposed to have hit the collimator, but have no kick, are considered faulty
        faulty = (abs(part.px[mask_hit]) < 1.e-12) & (part.state[mask_hit] > 0)
        fids = part.particle_id[mask_hit][faulty]
        faulty_inside = part_init.x[np.isin(part_init.particle_id, fids)]
        faulty_inside_y = part_init.y[np.isin(part_init.particle_id, fids)]
        clean_mask = (faulty_inside > 0) & (faulty_inside < 8.e-2)
        if len(faulty_inside[clean_mask]) > 0:
            spread_inside_L.append(faulty_inside[clean_mask].max() - coll.jaw_L)
        else:
            spread_inside_L.append(0)
        clean_mask = (faulty_inside < 0) & (faulty_inside > -8.e-2)
        if len(faulty_inside[clean_mask]) > 0:
            spread_inside_R.append(coll.jaw_R - faulty_inside[clean_mask].min())
        else:
            spread_inside_R.append(0)
        # plt.scatter(part_init.x, part_init.y, c='b', s=0.5)
        # plt.scatter(faulty_outside, faulty_outside_y, c='r', s=0.5)
        # plt.scatter(faulty_inside, faulty_inside_y, c='magenta', s=0.5)
        # plt.axvline(x=coll.jaw_L, c='g')
        # plt.axvline(x=coll.jaw_R, c='g')
        # plt.show()
        K2Engine().stop()
for j, iL, iR, oL, oR in zip(jaws, spread_inside_L, spread_inside_R, spread_outside_L, spread_outside_R):
    iL = f"{iL:.4e}" if iL > 0 else '          '
    iR = f"{iR:.4e}" if iR > 0 else '          '
    oL = f"{oL:.4e}" if oL > 0 else '          '
    oR = f"{oR:.4e}" if oR > 0 else '          '
    print(f"Jaw: {j}   outside: [{oL}, {oR}]   inside: [{iL}, {iR}]")


# TODO TODO TODO: issue with inside (and only when inside): outliers at x = +- 0.1
#                 seems to be related to the gap size (fixed large beta gives small gaps, and online then do we get inside vetoes)
plt.plot(jaws, spread_inside_L)
plt.plot(jaws, spread_inside_R)
plt.plot(jaws, spread_outside_L)
plt.plot(jaws, spread_outside_R)
plt.xscale('log')
plt.yscale('log')
plt.show()

plt.scatter(part_init.x, part_init.y, c='b', s=0.1)
plt.scatter(faulty_outside, faulty_outside_y, c='r', s=1.5)
plt.scatter(faulty_inside, faulty_inside_y, c='magenta', s=1.5)
plt.axvline(x=coll.jaw_L, c='g')
plt.axvline(x=coll.jaw_R, c='g')
plt.show()