# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

import pytest
import numpy as np
import pandas as pd
from shutil import rmtree

import xtrack as xt
import xcoll as xc

# from test_jaw_position import  _generate_particles


# jaws = [0.001, [0.0013, -0.002789], [-1.2e-6, -3.2e-3], [3.789e-3, 4.678e-7]]
# jaw_ids = ['symmetric', 'asymmetric', 'negative', 'positive']
# angles = [0, 90, 127.5]
# tilts = [0, [2.2e-6, 1.3e-6], [1.9e-6, -2.7e-6]]
# tilt_ids = ['no_tilt', 'positive_tilt', 'pos_neg_tilt']


# @pytest.mark.parametrize('tilt', tilts, ids=tilt_ids)
# @pytest.mark.parametrize('angle', angles)
# @pytest.mark.parametrize('jaw', jaws, ids=jaw_ids)
# def test_impacts(jaw, angle, tilt):
#     length = 0.6
#     particle_ref = xt.Particles('proton', p0c=6.8e12)
#     assembly = "lhc_tcp"

#     num_part = 1000
#     capacity = 2*num_part
#     jaw_band = 5e-9
#     jaw_accuracy = 1.e-9
#     x_dim = 0.015
#     y_dim = 0.015
#     exact_drift = True
#     if xc.fluka.engine.is_running():
#         xc.fluka.engine.stop(clean=True)
#     coll = xc.FlukaCollimator(length=length, jaw=jaw, angle=angle, tilt=tilt,  assembly=assembly)
#     impacts = xc.InteractionRecord(names=[coll.name], elements=[coll])

#     xc.fluka.engine.particle_ref = particle_ref
#     xc.fluka.engine.start(elements=coll, capacity=capacity, verbose=True, touches=True)
#     particle_ref = xc.fluka.engine.particle_ref

#     part_init, _, _ = _generate_particles(coll, num_part=num_part, particle_ref=particle_ref,
#                                                 jaw_band=jaw_band, jaw_accuracy=jaw_accuracy, angular_spread=1e-3,
#                                                 delta_spread=1e-3, zeta_spread=5e-2, exact_drift=exact_drift,
#                                                 x_dim=x_dim, y_dim=y_dim, _capacity=capacity)
#     part = part_init.copy()
#     coll.track(part)

#     fluka_path = xc.fluka.engine.cwd
#     file_path = fluka_path / "fluka_input001_toucMap.dat"

#     xc.fluka.engine.stop(clean=False)
#     impacts.stop(elements=[coll])

#     df = pd.read_csv(file_path, sep=r"\s+", comment='*')
#     file_path.unlink()
#     rmtree(fluka_path)

#     # FLUKA ids start in 1. Changed to 0
#     fluka_ids = df.iloc[:, 7].values -1
#     xcoll_ids = np.unique(impacts.id_before)

#     # Compare: find which xcoll_ids are present in IDP (and vice versa)
#     in_fluka_not_in_xcoll = np.setdiff1d(fluka_ids, xcoll_ids)
#     in_xcoll_not_in_fluka = np.setdiff1d(xcoll_ids, fluka_ids)

#     assert len(in_fluka_not_in_xcoll) == 0
#     print(f"In FLUKA but not Xcoll: {in_fluka_not_in_xcoll}")
#     print(f"In Xcoll but not FLUKA: {in_xcoll_not_in_fluka}")
