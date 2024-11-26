# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

import numpy as np
import pytest

import xobjects as xo
import xpart as xp
import xtrack as xt
import xcoll as xc
import time

import matplotlib.pyplot as plt
from scipy.stats import ks_2samp

def fluka_collimator():
    length = 0.6 # XXX TO BE REMOVED
    coll = xc.FlukaCollimator(length=length)
    coll_name = 'tcp.c6l7.b1'
    coll.jaw = 0.001
    coll.gap = 1 # XXX STILL NEEDED?
    return coll

def particle_distribution(num_part, xi, xf):
    # Distributed equally in distance between xi and xf
    x_init   = np.linspace(xi, xf, num_part)
    # The rest is 0
    px_init  = np.zeros(num_part)
    y_init   = np.zeros(num_part)
    py_init  = np.zeros(num_part)

    return x_init, px_init, y_init, py_init

def test_fluka_connection():
    num_part = 200
    _capacity = 10000
    coll = fluka_collimator()
    # # Create a FLUKA collimator

    # Connect to FLUKA
    coll_name = 'tcp.c6l7.b1'
    xc.FlukaEngine.start(elements=coll, names=coll_name, debug_level=1, _capacity=_capacity)

    part_1 = particle_distribution(num_part=200, xi=-0.002, xf=0.002)
    part_2 = particle_distribution(num_part=2400, xi=-0.0002-0.001, xf=-0.001+0.0002)
    part_3 = particle_distribution(num_part=2400, xi=-0.0002-0.001, xf=-0.001+0.0002)

    x_init = np.append(np.append(part_1[0], part_2[0]), part_3[0])
    px_init = np.append(np.append(part_1[1], part_2[1]), part_3[1])
    y_init = np.append(np.append(part_1[2], part_2[2]), part_3[2])
    py_init = np.append(np.append(part_1[3], part_2[3]), part_3[3])

    particle_ref = xp.Particles.reference_from_pdg_id(pdg_id='proton', p0c=6.8e12)
    xc.FlukaEngine.set_particle_ref(particle_ref=particle_ref)
    part_init = xp.build_particles(x=x_init, px=px_init, y=y_init, py=py_init, particle_ref=particle_ref, _capacity=_capacity)

    part = part_init.copy()

    coll.track(part)

    mask = part.state > 0

    # survival_particles_id = part.particle_id[mask]
    survival_particles_parents_id = part.parent_particle_id[mask]

    for pt in range(len(part_init.x[0:_capacity])):
        if part_init.x[pt]<-0.001: # Hitting left jaw
            if part_init.particle_id[pt] in survival_particles_parents_id:
                index = np.where(survival_particles_parents_id == part_init.particle_id[pt])
                # print("left jaw:")
                if len(part.px[index])>1: # Case of multiple particles produced
                    for i in range(len(part.px[index])):
                        assert abs(part.px[index][i])>1e-12
                        continue
                # print(part.px[index])
                assert abs(part.px[index][0])>1e-12
        elif part_init.x[pt]>-0.001 and part_init.x[pt]<0.001: # Hitting centre
            if part_init.particle_id[pt] in survival_particles_parents_id:
                index = np.where(survival_particles_parents_id == part_init.particle_id[pt])
                # print("centre:")
                if len(part.px[index])>1:
                    for i in range(len(part.px[index])):
                        assert abs(part.px[index][i])<1e-12
                        continue
                # print(part.px[index])
                assert abs(part.px[index][0])<1e-12
        elif part_init.x[pt]>0.001: # Hitting right jaw
            if part_init.particle_id[pt] in survival_particles_parents_id:
                index = np.where(survival_particles_parents_id == part_init.particle_id[pt])
                # print("right jaw:")
                if len(part.px[index])>1:
                    for i in range(len(part.px[index])):
                        assert abs(part.px[index][i])>1e-12
                        continue
                # print(part.px[index])
                assert abs(part.px[index][0])>1e-12
