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
from collections.abc import Iterable

global particles_absorbed_in_coll_1, particles_absorbed_in_coll_2, particles_absorbed_in_coll_3, particles_absorbed_in_coll_4
@pytest.fixture
def fluka_collimator():
    length = 0.6 # XXX TO BE REMOVED
    coll = xc.FlukaCollimator(length=length)
    coll_name = 'tcp.c6l7.b1'
    coll.jaw = 0.001
    coll.gap = 1 # XXX STILL NEEDED?
    return coll

@pytest.fixture
def particle_distribution(num_part):
    x_init   = np.random.normal(loc=0.003, scale=1e-3, size=num_part)
    px_init  = np.random.normal(loc=0., scale=5.e-6, size=num_part)
    y_init   = np.random.normal(loc=0., scale=1e-3, size=num_part)
    py_init  = np.random.normal(loc=0., scale=5.e-6, size=num_part)
    particle_ref = xp.Particles.reference_from_pdg_id(pdg_id='proton', p0c=6.8e12)
    # particle_ref = xp.Particles.reference_from_pdg_id(pdg_id='Pb208', p0c=522.34e12)
    xc.FlukaEngine.set_particle_ref(particle_ref=particle_ref)
    return xp.build_particles(x=x_init, px=px_init, y=y_init, py=py_init, particle_ref=particle_ref, _capacity=num_part*2)

@pytest.mark.parametrize("num_part", [10])
def test_track_id(fluka_collimator, particle_distribution, num_part):
    class ParticleTree:
        def __init__(self):
            self.children = {}
            for n in range(num_part):
                self.children[n] = []

        def add_particles(self, particle_id_list, parent_id_list):
            for pid, parent in zip(particle_id_list, parent_id_list):
                if parent not in self.children:
                    self.children[parent] = []
                if pid not in self.children[parent]:
                    self.children[parent].append(pid)
        
        def print_childs(self, file, childs, level):
            prefix = "  " * level 
            print(childs)

            if not isinstance(childs, Iterable):
                childs = [childs]  # Wrap in a list if it's a scalar
 
            for child in childs:
                if child in particles_absorbed_in_coll_1:
                    file.write(f"{prefix}|-- {child} absorbed in first collimator\n")
                elif child in particles_absorbed_in_coll_2:
                    file.write(f"{prefix}|-- {child} absorbed in second collimator\n")
                elif child in particles_absorbed_in_coll_3:
                    file.write(f"{prefix}|-- {child} absorbed in third collimator\n")
                elif child in particles_absorbed_in_coll_4:
                    file.write(f"{prefix}|-- {child} absorbed in fourth collimator\n")
                else:
                    file.write(f"{prefix}|-- {child}\n")
                if child in self.children.keys():
                    self.print_childs(file, self.children[child], level + 1)


        def dump_tree_to_file(self, file):
            # Write children recursively using the correct child lookup
            for p_id in range(num_part):
                level = 0
                prefix = "  " * level  # Indentation based on depth
                #if self.children[p_id] == []:
                if p_id in particles_absorbed_in_coll_1:
                    file.write(f"{prefix}|-- {p_id} absorbed in first collimator\n")
                elif p_id in particles_absorbed_in_coll_2:
                    file.write(f"{prefix}|-- {p_id} absorbed in second collimator\n")
                elif p_id in particles_absorbed_in_coll_3:
                    file.write(f"{prefix}|-- {p_id} absorbed in third collimator\n")    
                elif p_id in particles_absorbed_in_coll_4:
                    file.write(f"{prefix}|-- {p_id} absorbed in fourth collimator\n")
                else:
                    file.write(f"{prefix}|-- {p_id}\n")

                # import pdb; pdb.set_trace()
                if self.children[p_id]:
                    for childs in self.children[p_id]:
                        if childs == p_id: 
                            if childs in particles_absorbed_in_coll_1 or childs in particles_absorbed_in_coll_2 or childs in particles_absorbed_in_coll_3 or childs in particles_absorbed_in_coll_4:
                                continue
                            file.write(f"{prefix}|-- {childs} is alive!\n")
                            continue
                        self.print_childs(file, childs, level + 1)

    _capacity = num_part*2
    coll = fluka_collimator
    # Create a FLUKA collimator

    # Connect to FLUKA
    coll_name = 'tcp.c6l7.b1'
    xc.FlukaEngine.start(elements=coll, names=coll_name, debug_level=1, _capacity=_capacity)

    part_init = particle_distribution

    part = part_init.copy()

    particle_tree = ParticleTree()

    # Do the tracking in first TCP
    start = time.time()
    coll.track(part)
    fluka_time = round(time.time()-start, 3)

    mask = part.state > 0
    maks_absorbed = part.state == -334
    print(maks_absorbed)
    print(part.particle_id[maks_absorbed])

    particles_alive_id_coll_1 = part.particle_id[mask]
    particles_absorbed_in_coll_1 = part.particle_id[part.state == -334]
    parents_1 = part.parent_particle_id[mask]
    print(particles_absorbed_in_coll_1)

    particle_tree.add_particles(particles_alive_id_coll_1, parents_1)
    # Dump the tree structure to a file
    print("First TCP, particles alive and their parents:")
    print(particles_alive_id_coll_1)
    print(parents_1)

    # Do the tracking in second TCP

    part2 = part.copy()

    start = time.time()
    coll.track(part2)
    fluka_time = round(time.time()-start, 3)

    mask = part2.state > 0

    particles_alive_id_coll_2 = part2.particle_id[mask]
    particles_absorbed_in_coll_2 = part2.particle_id[part2.state == -334]
    parents_2 = part2.parent_particle_id[mask]
    particle_tree.add_particles(particles_alive_id_coll_2, parents_2)

    print("Second TCP, particles alive and their parents:")
    print(particles_alive_id_coll_2)
    print(parents_2)

    # Do the tracking in third TCP

    part3 = part2.copy()

    start = time.time()
    coll.track(part3)
    fluka_time = round(time.time()-start, 3)

    mask = part3.state > 0

    particles_alive_id_coll_3 = part3.particle_id[mask]
    particles_absorbed_in_coll_3 = part3.particle_id[part3.state == -334]
    parents_3 = part3.parent_particle_id[mask]
    particle_tree.add_particles(particles_alive_id_coll_3, parents_3)

    print("Third TCP, particles alive and their parents:")
    print(particles_alive_id_coll_3)
    print(parents_3)

    print(particle_tree.children)

    # Do the tracking in fourth TCP

    part4 = part3.copy()

    start = time.time()
    coll.track(part4)
    fluka_time = round(time.time()-start, 3)

    mask = part4.state > 0

    particles_alive_id_coll_4 = part4.particle_id[mask]
    particles_absorbed_in_coll_4 = part4.particle_id[part4.state == -334]
    parents_4 = part4.parent_particle_id[mask]
    particle_tree.add_particles(particles_alive_id_coll_4, parents_4)

    

    with open("particle_tree.txt", "w") as file:
        particle_tree.dump_tree_to_file(file)
    # Stop the FLUKA server
    xc.FlukaEngine.stop()

    # Compare distribution using Kolmogorov-Smirnov test


#     print(f"Survived in FLUKA: {len(part.state[part