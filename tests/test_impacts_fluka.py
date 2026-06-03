# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

import time
import pytest
import numpy as np

import xpart as xp
from xpart.test_helpers import flaky_assertions, retry
import xtrack as xt
import xcoll as xc
import pandas as pd  
import os


jaws = [0.001, [0.0013, -0.002789], [-1.2e-6, -3.2e-3], [3.789e-3, 4.678e-7]]
jaw_ids = ['symmetric', 'asymmetric', 'negative', 'positive']
angles = [0, 90, 127.5]
tilts = [0, [2.2e-6, 1.3e-6], [1.9e-6, -2.7e-6]]
tilt_ids = ['no_tilt', 'positive_tilt', 'pos_neg_tilt']


@pytest.mark.parametrize('tilt', tilts, ids=tilt_ids)
@pytest.mark.parametrize('angle', angles)
@pytest.mark.parametrize('jaw', jaws, ids=jaw_ids)
@retry()
# XXX Maybe test all collimator assemblies?
def test_impacts(jaw, angle, tilt):
    length = 0.6
    material = xc.materials.MolybdenumGraphite
    particle_ref = xt.Particles('proton', p0c=6.8e12)
    assembly = "lhc_tcp"

    num_part = 1000
    capacity = 2*num_part
    jaw_band = 5e-9
    jaw_accuracy = 1.e-9
    x_dim = 0.015
    y_dim = 0.015
    exact_drift = True
    if xc.fluka.engine.is_running():
        xc.fluka.engine.stop(clean=True)
    # coll = xc.FlukaCollimator(length=length, jaw=jaw, angle=angle, tilt=tilt,  material=material)
    coll = xc.FlukaCollimator(length=length, jaw=jaw, angle=angle, tilt=tilt,  assembly=assembly)
    impacts = xc.InteractionRecord.start(names=[coll.name], elements=[coll])

    xc.fluka.engine.particle_ref = particle_ref
    xc.fluka.engine.start(elements=coll, capacity=capacity, verbose=True, touches=True)
    particle_ref = xc.fluka.engine.particle_ref

    part_init = _generate_particles(coll, num_part=num_part, particle_ref=particle_ref,
                                                jaw_band=jaw_band, jaw_accuracy=jaw_accuracy, angular_spread=1e-3,
                                                delta_spread=1e-3, zeta_spread=5e-2, exact_drift=exact_drift,
                                                x_dim=x_dim, y_dim=y_dim, _capacity=capacity)

    part = part_init.copy()
    t1 = time.time()
    coll.track(part)

    fluka_path = xc.fluka.engine.cwd
    file_path = os.path.join(fluka_path, "fluka_input001_toucMap.dat")

    xc.fluka.engine.stop(clean=False)

    impacts.stop(names=[coll.name], elements=[coll])

    df = pd.read_csv(file_path, sep=r"\s+", comment='*')

    # FLUKA ids start in 1. Changed to 0
    fluka_ids = df.iloc[:, 7].values -1

    xcoll_ids = np.unique(impacts.id_before)

    # 4. Compare: find which xcoll_ids are present in IDP (and vice versa)
    in_fluka_not_in_xcoll = np.setdiff1d(fluka_ids, xcoll_ids)
    in_xcoll_not_in_fluka = np.setdiff1d(xcoll_ids, fluka_ids)
    common = np.intersect1d(fluka_ids, xcoll_ids)

    assert len(in_fluka_not_in_xcoll) == 0

    # print(f"FLUKA values:           {fluka_ids}")
    # print(f"Xcoll: {xcoll_ids}")
    # print(f"Common:               {common}")
    print(f"In FLUKA but not Xcoll: {in_fluka_not_in_xcoll}")
    print(f"In Xcoll but not FLUKA: {in_xcoll_not_in_fluka}")

    import subprocess

    # subprocess.run(["rm", "-rf", str(fluka_path)], check=True)
    path_str = str(fluka_path)
    subprocess.run(["find", path_str, "-type", "f", "-delete"], check=False)
    subprocess.run(["find", path_str, "-type", "d", "-empty", "-delete"], check=False)


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

    return part_init
