# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #


import numpy as np
from pathlib import Path
import pytest
from scipy.stats import skew, kurtosis, moment

import xobjects as xo
import xpart as xp
import xtrack as xt
import xcoll as xc
from xcoll.scattering_routines.k2 import K2Engine
from xcoll.beam_elements.k2 import _K2Collimator

def _get_moment(hist0, histpos, histneg, histrand):
    mean  = np.array([np.mean(hist0), np.mean(histpos), np.mean(histneg), np.mean(histrand)])
    var   = np.array([np.var(hist0), np.var(histpos), np.var(histneg), np.var(histrand)])
    skewn = np.array([skew(hist0), skew(histpos), skew(histneg), skew(histrand)])
    kurt  = np.array([kurtosis(hist0), kurtosis(histpos), kurtosis(histneg), kurtosis(histrand)])
    return np.array([mean, var, skewn, kurt])

def _are_numbers_equal_within_tolerance(num1, num2, tolerance_percent): # should this be relative 
    if not hasattr(num1,'__iter__') and not hasattr(num2, '__iter__'):
        difference = abs(num1 - num2)
        relative_error = difference / num1 * 100 if num1 != 0 else 0
        print(f" difference {relative_error}, {num1}, {num2}")
        return relative_error <= tolerance_percent
    else:
        for i in range(len(num1)):
            difference = abs(num1[i] - num2[i])
            relative_error = (difference / abs(num1[i])) * 100 if num1[i] != 0 else 0
            if relative_error > tolerance_percent:
                print(f" relative error {relative_error}, {i}")
                print(f" num1 {num1[i]}, num2 {num2[i]}")
                return False
        return True

def _create_4_particles(line, pos, num_particles, plane):
    angle0     = np.zeros(num_particles)
    angle_pos  = np.random.uniform(0,1e-6, num_particles)
    angle_neg  = np.random.uniform(-1e-6,0, num_particles)
    angle_rand = np.random.uniform(-1e-6,1e-6, num_particles)
    pos_distr  = np.random.normal(pos, 1e-9, int(num_particles))
    if plane == 'H':
        part_0    = line.build_particles(x = pos_distr, px = angle0)
        part_pos  = line.build_particles(x = pos_distr, px = angle_pos)
        part_neg  = line.build_particles(x = pos_distr, px = angle_neg)
        part_rand = line.build_particles(x = pos_distr, px = angle_rand)
    else:
        part_0    = line.build_particles(y = pos_distr, py = angle0)
        part_pos  = line.build_particles(y = pos_distr, py = angle_pos)
        part_neg  = line.build_particles(y = pos_distr, py = angle_neg)
        part_rand = line.build_particles(y = pos_distr, py = angle_rand)
    return part_0, part_pos, part_neg, part_rand

def _track_with_angles(beam, plane, pos, everest=False, material=False, K2coll=None, Ecoll=None):
    path = Path(__file__).parent / 'data'
    num_particles = int(1e6)

    line = xt.Line.from_json(path / f'sequence_lhc_run3_b{beam}.json')

    colldb = xc.CollimatorDatabase.from_yaml(path / f'colldb_lhc_run3_ir7.yaml',
                                    beam=beam, ignore_crystals=False)
    if everest and not material:
        colldb.install_everest_collimators(line=line, verbose=True)
        coll  = line[f"tcp.{'c' if plane=='H' else 'd'}6{'l' if beam==1 else 'r'}7.b{beam}"]
        line.build_tracker()
        xc.assign_optics_to_collimators(line=line)
    elif not everest and not material:
        colldb._install_k2_collimators(line=line, verbose=True)
        coll  = line[f"tcp.{'c' if plane=='H' else 'd'}6{'l' if beam==1 else 'r'}7.b{beam}"]
        line.build_tracker()
        xc.assign_optics_to_collimators(line=line)
        K2Engine.start(line=line, cwd='run_1',_capacity=num_particles)
    elif material and everest:
        xc.install_elements(line, [Ecoll[1]], [Ecoll[0]], need_apertures=False)
        line.build_tracker()
        xc.assign_optics_to_collimators(line=line)
        coll = Ecoll[0]
    else:
        xc.install_elements(line, [K2coll[1]], [K2coll[0]], need_apertures=False)
        line.build_tracker()
        xc.assign_optics_to_collimators(line=line)
        K2Engine.start(line=line, cwd='run_1',_capacity=num_particles)
        coll = K2coll[0]

    part_zero_init, part_pos_init, part_neg_init, part_rand_init = _create_4_particles(line=line, pos=pos, num_particles=num_particles, plane=plane)
    part_zero = part_zero_init.copy()
    part_pos = part_pos_init.copy()
    part_neg = part_neg_init.copy()
    part_rand = part_rand_init.copy()

    xc.enable_scattering(line=line)
    coll.track(part_zero)
    coll.track(part_pos)
    coll.track(part_neg)
    coll.track(part_rand)
    xc.disable_scattering(line=line)
    line.discard_tracker()

    part_zero.sort(interleave_lost_particles=True)
    part_pos.sort(interleave_lost_particles=True)
    part_neg.sort(interleave_lost_particles=True)    
    part_rand.sort(interleave_lost_particles=True)

    return part_zero, part_pos, part_neg, part_rand, part_zero_init, part_pos_init, part_neg_init, part_rand_init, line, coll

############################################### TESTS #######################################################################################

@pytest.mark.parametrize("beam, plane",[[1,'V'],[2,'H'], [1,'H'], [2,'V']])
def test_everest_and_K2_angles(beam, plane):
    if plane == 'H':
        pos = [0.00098, 0.00093, 0.0009198, 0.00091, 0.00088] 
    else:
        pos = [0.0015, 0.00131, 0.0013092, 0.0013082, 0.00012]

    for idx,i in enumerate(pos):
        part_0, part_pos, part_neg, part_rand, part_0_init, part_pos_init, part_neg_init, part_rand_init,_,_ = _track_with_angles(beam, plane, pos=i)
        part_0_E, part_pos_E, part_neg_E, part_rand_E, part_0_init_E, \
        part_pos_init_E, part_neg_init_E, part_rand_init_E, _, coll = _track_with_angles(beam, plane, pos=i, everest=True)
        print(f"pos {i}")
        # checks that same number of particles are alive within a tolerance of 1 %
        assert _are_numbers_equal_within_tolerance(np.sum(part_0.state < 1),np.sum(part_0_E.state < 1),1)
        assert _are_numbers_equal_within_tolerance(np.sum(part_pos.state < 1),np.sum(part_pos_E.state < 1),1)
        assert _are_numbers_equal_within_tolerance(np.sum(part_neg.state < 1),np.sum(part_neg_E.state < 1),1)
        assert _are_numbers_equal_within_tolerance(np.sum(part_rand.state < 1),np.sum(part_rand_E.state < 1),1)

        # check states before/after jaw
        if plane == 'H':
            print(np.unique(part_0.state))
            print(np.unique(part_0_E.state))
            print(np.unique(part_0.x[part_0.state>0]))
            print(np.unique(part_0_E.x[part_0_E.state>0]))
            print(np.unique(part_0.s))
            print(np.unique(part_0_E.s))
            assert _are_numbers_equal_within_tolerance(np.sum( part_0.state[part_0_init.x >= (coll.jaw_L)] < 1), np.sum(part_0_E.state[part_0_init_E.x >= (coll.jaw_L)])< 1,1) 
            assert _are_numbers_equal_within_tolerance(np.sum( part_0.state[part_0_init.x < (coll.jaw_L)]<1), np.sum(part_0_E.state[part_0_init_E.x < (coll.jaw_L)]<1),1)
            assert _are_numbers_equal_within_tolerance(np.sum( part_pos.state[part_pos_init.x >= (coll.jaw_L)]<1), np.sum(part_pos_E.state[part_pos_init_E.x >= (coll.jaw_L)]<1),1)
            assert _are_numbers_equal_within_tolerance(np.sum( part_pos.state[part_pos_init.x < (coll.jaw_L)]<1), np.sum(part_pos_E.state[part_pos_init_E.x < (coll.jaw_L)]<1),1)
            assert _are_numbers_equal_within_tolerance(np.sum( part_neg.state[part_neg_init.x >= (coll.jaw_L)]<1), np.sum(part_neg_E.state[part_neg_init_E.x >= (coll.jaw_L)]<1),1)
            assert _are_numbers_equal_within_tolerance(np.sum( part_neg.state[part_neg_init.x < (coll.jaw_L)]<1), np.sum(part_neg_E.state[part_neg_init_E.x < (coll.jaw_L)]<1),1)
            assert _are_numbers_equal_within_tolerance(np.sum( part_rand.state[part_rand_init.x >= (coll.jaw_L)]<1), np.sum(part_rand_E.state[part_rand_init_E.x >= (coll.jaw_L)]<1),1)
            assert _are_numbers_equal_within_tolerance(np.sum( part_rand.state[part_rand_init.x < (coll.jaw_L)]<1), np.sum(part_rand_E.state[part_rand_init_E.x < (coll.jaw_L)]<1),1)
            if idx != 0 or idx != 4: # tacky 
                assert np.sum(part_pos.state < 1) <= np.sum(part_neg.state < 1)
        else:
            print(np.unique(part_0.state))
            print(np.unique(part_0_E.state))
            print(np.unique(part_0.y[part_0.state>0]))
            print(np.unique(part_0_E.y[part_0_E.state>0]))
            print(np.unique(part_0.s))
            print(np.unique(part_0_E.s))
            assert _are_numbers_equal_within_tolerance(np.sum(part_0.state[part_0_init.y >= (coll.jaw_R)]<1), np.sum(part_0_E.state[part_0_init_E.y >= (coll.jaw_R)]<1),1)
            assert _are_numbers_equal_within_tolerance(np.sum(part_0.state[part_0_init.y < (coll.jaw_R)]<1), np.sum(part_0_E.state[part_0_init_E.y < (coll.jaw_R)]<1),1)
            assert _are_numbers_equal_within_tolerance(np.sum(part_pos.state[part_pos_init.y >= (coll.jaw_R)]<1), np.sum(part_pos_E.state[part_pos_init_E.y >= (coll.jaw_R)]<1),1)
            assert _are_numbers_equal_within_tolerance(np.sum(part_pos.state[part_pos_init.y < (coll.jaw_R)]<1), np.sum(part_pos_E.state[part_pos_init_E.y < (coll.jaw_R)]<1),1)
            assert _are_numbers_equal_within_tolerance(np.sum(part_neg.state[part_neg_init.y >= (coll.jaw_R)]<1), np.sum(part_neg_E.state[part_neg_init_E.y >= (coll.jaw_R)]<1),1)
            assert _are_numbers_equal_within_tolerance(np.sum(part_neg.state[part_neg_init.y < (coll.jaw_R)]<1), np.sum(part_neg_E.state[part_neg_init_E.y < (coll.jaw_R)]<1),1)
            assert _are_numbers_equal_within_tolerance(np.sum(part_rand.state[part_rand_init.y >= (coll.jaw_R)]<1), np.sum(part_rand_E.state[part_rand_init_E.y >= (coll.jaw_R)]<1),1)
            assert _are_numbers_equal_within_tolerance(np.sum(part_rand.state[part_rand_init.y < (coll.jaw_R)]<1), np.sum(part_rand_E.state[part_rand_init_E.y < (coll.jaw_R)]<1),1)
            if idx != 0 or idx != 4:
                assert np.sum(part_pos.state < 1) >= np.sum(part_neg.state < 1)

        # garbage collection for memory
        del part_0, part_pos, part_neg, part_rand
        del part_0_init, part_pos_init, part_neg_init, part_rand_init
        del part_0_E, part_pos_E, part_neg_E, part_rand_E
        del part_0_init_E, part_pos_init_E, part_neg_init_E, part_rand_init_E

@pytest.mark.parametrize("beam, plane",[[1,'V']])  
def test_everest_and_K2_materials(beam, plane):
    pos    = [0.0015, 0.00131, 0.0013, 0.00129, 0.0011]

    light  = _K2Collimator(length=0.6, jaw=0.0013, material='C', angle=90, emittance=3.5e-6) # 1.67
    middle = _K2Collimator(length=0.6, jaw=0.0013, material='MoGR', angle=90, emittance=3.5e-6) # 10.22 
    heavy  = _K2Collimator(length=0.6, jaw=0.0013, material='Iner', angle=90, emittance=3.5e-6) # 18 
    light_E  = xc.EverestCollimator(length=0.6, jaw=0.0013, material=xc.materials.Carbon, angle=90, emittance=3.5e-6) # 1.67 
    middle_E = xc.EverestCollimator(length=0.6, jaw=0.0013, material=xc.materials.MolybdenumGraphite, angle=90, emittance=3.5e-6) # 10.22 
    heavy_E  = xc.EverestCollimator(length=0.6, jaw=0.0013, material=xc.materials.Inermet, angle=90, emittance=3.5e-6) # 18 

    K2      = np.array([[light, 'tcp.b6l7.b1'], [middle,'tcp.c6l7.b1'],[heavy,'tcp.d6l7.b1']])
    Everest = np.array([[light_E, 'tcp.b6l7.b1'], [middle_E,'tcp.c6l7.b1'],[heavy_E,'tcp.d6l7.b1']])

    abs_light_inside  = np.array([0,0]) # K2, EVEREST
    abs_light_corner  = np.array([0,0])
    abs_light_drift   = np.array([0,0])
    abs_middle_inside = np.array([0,0])  
    abs_middle_corner = np.array([0,0])
    abs_middle_drift  = np.array([0,0])
    abs_heavy_inside  = np.array([0,0])
    abs_heavy_corner  = np.array([0,0])
    abs_heavy_drift   = np.array([0,0])

    for j in range(len(K2)): 
        for idx, i in enumerate(pos):
            part_0, part_pos, part_neg, part_rand,_, _, _, _,_,_ = _track_with_angles(beam=beam, plane=plane, pos=i, material=True, K2coll=K2[j])
            part_0_E, part_pos_E, part_neg_E, part_rand_E,_,_,_,_,_,_ = _track_with_angles(beam, plane, pos=i, material=True, everest=True, Ecoll=Everest[j])

            assert _are_numbers_equal_within_tolerance(np.sum(part_0.state < 1), np.sum(part_0_E.state < 1),1.5), f"{np.sum(part_0.state < 1)}, {np.sum(part_0_E.state<1)}"
            assert _are_numbers_equal_within_tolerance(np.sum(part_pos.state < 1), np.sum(part_pos_E.state < 1),1.5), f"{np.sum(part_pos.state<1)}, {np.sum(part_pos_E.state<1)}"
            assert _are_numbers_equal_within_tolerance(np.sum(part_neg.state < 1), np.sum(part_neg_E.state < 1),1.5), f"{np.sum(part_neg.state<1)}, {np.sum(part_neg_E.state<1)}"
            assert _are_numbers_equal_within_tolerance(np.sum(part_rand.state < 1), np.sum(part_rand_E.state < 1),1.5),  f"{np.sum(part_rand.state<1)}, {np.sum(part_rand_E.state<1)}"
            if j == 0:
                if idx == 0:
                    abs_light_inside[0] += np.sum(part_neg.state < 1) + np.sum(part_pos.state < 1) + np.sum(part_rand.state < 1) + np.sum(part_0.state < 1)
                    abs_light_inside[1] += np.sum(part_neg_E.state < 1) + np.sum(part_pos_E.state < 1) + np.sum(part_rand_E.state < 1) + np.sum(part_0_E.state < 1)
                    print(f" abs light inside 0: {(np.sum(part_0.state < 1)/int(1e6))*100} %, {(np.sum(part_0_E.state<1)/int(1e6))*100} %")
                    print(f" abs light inside pos {(np.sum(part_pos.state < 1)/int(1e6)*100)} %, {(np.sum(part_pos_E.state<1)/int(1e6))*100} %")
                    print(f" abs light inside neg {np.sum(part_neg.state < 1)/int(1e6)} %, {np.sum(part_neg_E.state<1)/int(1e6)} %")
                    print(f" abs light inside rand {np.sum(part_rand.state < 1)/int(1e6)} %, {np.sum(part_rand_E.state<1)/int(1e6)} %")
                    print(f"abs_light_inside {(np.sum(part_neg.state < 1) + np.sum(part_pos.state < 1) + np.sum(part_rand.state < 1) + (np.sum(part_0.state < 1))/(4*int(1e6)))*100} %")
                    print(f"abs light inside E {(np.sum(part_neg_E.state < 1) + np.sum(part_pos_E.state < 1) + np.sum(part_rand_E.state < 1) + (np.sum(part_0_E.state < 1))/(4*int(1e6)))*100} %")
                elif idx == 2:
                    abs_light_corner[0] += np.sum(part_neg.state < 1) + np.sum(part_pos.state < 1) + np.sum(part_rand.state < 1) + np.sum(part_0.state < 1)
                    abs_light_corner[1] += np.sum(part_neg_E.state < 1) + np.sum(part_pos_E.state < 1) + np.sum(part_rand_E.state < 1) + np.sum(part_0_E.state < 1)
                    print(f"abs light corner 0 {(np.sum(part_0.state < 1)/int(1e6))*100} %,{(np.sum(part_0_E.state<1)/int(1e6))*100}%")
                    print(f"abs light corner pos {(np.sum(part_pos.state < 1)/int(1e6)*100)} %, {(np.sum(part_pos_E.state<1)/int(1e6))*100} %")
                    print(f"abs light corner neg {np.sum(part_neg.state < 1)/int(1e6)} %, {np.sum(part_neg_E.state<1)/int(1e6)} %")
                    print(f"abs light corner rand {np.sum(part_rand.state < 1)/int(1e6)} %, {np.sum(part_rand_E.state<1)/int(1e6)} %")
                    print(f"abs_light_corner {(np.sum(part_neg.state < 1) + np.sum(part_pos.state < 1) + np.sum(part_rand.state < 1) + (np.sum(part_0.state < 1))/(4*int(1e6)))*100} %")
                    print(f"abs light corner E {(np.sum(part_neg_E.state < 1) + np.sum(part_pos_E.state < 1) + np.sum(part_rand_E.state < 1) + (np.sum(part_0_E.state < 1))/(4*int(1e6)))*100} %")
                elif idx == 4:
                    abs_light_drift[0] += np.sum(part_neg.state < 1) + np.sum(part_pos.state < 1) + np.sum(part_rand.state < 1) + np.sum(part_0.state < 1)
                    abs_light_drift[1] += np.sum(part_neg_E.state < 1) + np.sum(part_pos_E.state < 1) + np.sum(part_rand_E.state < 1) + np.sum(part_0_E.state < 1)
                    print(f"abs light drift 0 {(np.sum(part_0.state < 1)/int(1e6))*100} %,{(np.sum(part_0_E.state<1)/int(1e6))*100}%")
                    print(f"abs light drift pos {(np.sum(part_pos.state < 1)/int(1e6)*100)} %, {(np.sum(part_pos_E.state<1)/int(1e6))*100} %")
                    print(f"abs light drift neg {np.sum(part_neg.state < 1)/int(1e6)} %, {np.sum(part_neg_E.state<1)/int(1e6)} %")
                    print(f"abs light drift rand {np.sum(part_rand.state < 1)/int(1e6)} %, {np.sum(part_rand_E.state<1)/int(1e6)} %")
                    print(f"abs_light_drift {(np.sum(part_neg.state < 1) + np.sum(part_pos.state < 1) + np.sum(part_rand.state < 1) + (np.sum(part_0.state < 1))/(4*int(1e6)))*100} %")
                    print(f"abs light drift E {(np.sum(part_neg_E.state < 1) + np.sum(part_pos_E.state < 1) + np.sum(part_rand_E.state < 1) + (np.sum(part_0_E.state < 1))/(4*int(1e6)))*100} %")
            if j == 1:
                if idx == 0:
                    abs_middle_inside[0] += np.sum(part_neg.state < 1) + np.sum(part_pos.state < 1) + np.sum(part_rand.state < 1) + np.sum(part_0.state < 1)
                    abs_middle_inside[1] += np.sum(part_neg_E.state < 1) + np.sum(part_pos_E.state < 1) + np.sum(part_rand_E.state < 1) + np.sum(part_0_E.state < 1)
                    print(f" abs middle inside 0 {(np.sum(part_0.state < 1)/int(1e6))*100}%, {(np.sum(part_0_E.state<1)/int(1e6))*100} %")
                    print(f" abs middle inside pos {(np.sum(part_pos.state < 1)/int(1e6)*100)} %, {(np.sum(part_pos_E.state<1)/int(1e6))*100} %")
                    print(f" abs middle inside neg {np.sum(part_neg.state < 1)/int(1e6)} %, {np.sum(part_neg_E.state<1)/int(1e6)} %")
                    print(f" abs middle inside rand {np.sum(part_rand.state < 1)/int(1e6)} %, {np.sum(part_rand_E.state<1)/int(1e6)} %")
                    print(f"abs_middle_inside {(np.sum(part_neg.state < 1) + np.sum(part_pos.state < 1) + np.sum(part_rand.state < 1) + (np.sum(part_0.state < 1))/(4*int(1e6)))*100} %")
                    print(f"abs middle inside E {(np.sum(part_neg_E.state < 1) + np.sum(part_pos_E.state < 1) + np.sum(part_rand_E.state < 1) + (np.sum(part_0_E.state < 1))/(4*int(1e6)))*100} %")
                elif idx == 2:
                    abs_middle_corner[0] += np.sum(part_neg.state < 1) + np.sum(part_pos.state < 1) + np.sum(part_rand.state < 1) + np.sum(part_0.state < 1)
                    abs_middle_corner[1] += np.sum(part_neg_E.state < 1) + np.sum(part_pos_E.state < 1) + np.sum(part_rand_E.state < 1) + np.sum(part_0_E.state < 1)
                    print(f"abs_middle_corner 0 {(np.sum(part_0.state < 1)/int(1e6))*100} %,{(np.sum(part_0_E.state<1)/int(1e6))*100} %")
                    print(f"abs_middle_corner pos {(np.sum(part_pos.state < 1)/int(1e6)*100)} %, {(np.sum(part_pos_E.state<1)/int(1e6))*100} %")
                    print(f"abs_middle_corner neg {np.sum(part_neg.state < 1)/int(1e6)} %, {np.sum(part_neg_E.state<1)/int(1e6)} %")
                    print(f"abs_middle_corner rand {np.sum(part_rand.state < 1)/int(1e6)} %, {np.sum(part_rand_E.state<1)/int(1e6)} %")
                    print(f"abs_middle_corner {(np.sum(part_neg.state < 1) + np.sum(part_pos.state < 1) + np.sum(part_rand.state < 1) + (np.sum(part_0.state < 1))/(4*int(1e6)))*100} %")
                    print(f"abs middle corner E {(np.sum(part_neg_E.state < 1) + np.sum(part_pos_E.state < 1) + np.sum(part_rand_E.state < 1) + (np.sum(part_0_E.state < 1))/(4*int(1e6)))*100} %")
                elif idx == 4:
                    abs_middle_drift[0] += np.sum(part_neg.state < 1) + np.sum(part_pos.state < 1) + np.sum(part_rand.state < 1) + np.sum(part_0.state < 1)
                    abs_middle_drift[1] += np.sum(part_neg_E.state < 1) + np.sum(part_pos_E.state < 1) + np.sum(part_rand_E.state < 1) + np.sum(part_0_E.state < 1)
                    print(f"abs_middle_drift 0 {(np.sum(part_0.state < 1)/int(1e6))*100}%, {(np.sum(part_0_E.state<1)/int(1e6))*100} %")
                    print(f"abs_middle_drift pos {(np.sum(part_pos.state < 1)/int(1e6)*100)} %, {(np.sum(part_pos_E.state<1)/int(1e6))*100} %")
                    print(f"abs_middle_drift neg {np.sum(part_neg.state < 1)/int(1e6)} %, {np.sum(part_neg_E.state<1)/int(1e6)} %") 
                    print(f"abs_middle_drift rand {np.sum(part_rand.state < 1)/int(1e6)} %, {np.sum(part_rand_E.state<1)/int(1e6)} %")
                    print(f"abs_middle_drift {(np.sum(part_neg.state < 1) + np.sum(part_pos.state < 1) + np.sum(part_rand.state < 1) + (np.sum(part_0.state < 1))/(4*int(1e6)))*100} %")
                    print(f"abs middle drift E {(np.sum(part_neg_E.state < 1) + np.sum(part_pos_E.state < 1) + np.sum(part_rand_E.state < 1) + (np.sum(part_0_E.state < 1))/(4*int(1e6)))*100} %")
            if j == 2:
                if idx == 0:
                    abs_heavy_inside[0] += np.sum(part_neg.state < 1) + np.sum(part_pos.state < 1) + np.sum(part_rand.state < 1) + np.sum(part_0.state < 1)
                    abs_heavy_inside[1] += np.sum(part_neg_E.state < 1) + np.sum(part_pos_E.state < 1) + np.sum(part_rand_E.state < 1) + np.sum(part_0_E.state < 1)
                    print(f"abs_heavy_inside 0 {(np.sum(part_0.state < 1)/int(1e6)*100)} %, {(np.sum(part_0_E.state<1)/int(1e6))*100} %")
                    print(f"abs_heavy_inside pos {(np.sum(part_pos.state < 1)/int(1e6)*100)} %, {(np.sum(part_pos_E.state<1)/int(1e6))*100} %")
                    print(f"abs_heavy_inside neg {np.sum(part_neg.state < 1)/int(1e6)} %, {np.sum(part_neg_E.state<1)/int(1e6)} %")
                    print(f"abs_heavy_inside rand {np.sum(part_rand.state < 1)/int(1e6)} %, {np.sum(part_rand_E.state<1)/int(1e6)} %")
                    print(f"abs_heavy_inside {(np.sum(part_neg.state < 1) + np.sum(part_pos.state < 1) + np.sum(part_rand.state < 1) + (np.sum(part_0.state < 1))/(4*int(1e6)))*100} %")
                    print(f"abs heavy inside E {(np.sum(part_neg_E.state < 1) + np.sum(part_pos_E.state < 1) + np.sum(part_rand_E.state < 1) + (np.sum(part_0_E.state < 1))/(4*int(1e6)))*100} %")
                elif idx == 2:
                    abs_heavy_corner[0] += np.sum(part_neg.state < 1) + np.sum(part_pos.state < 1) + np.sum(part_rand.state < 1) + np.sum(part_0.state < 1)
                    abs_heavy_corner[1] += np.sum(part_neg_E.state < 1) + np.sum(part_pos_E.state < 1) + np.sum(part_rand_E.state < 1) + np.sum(part_0_E.state < 1)
                    print(f"abs_heavy_corner 0 {(np.sum(part_0.state < 1)/int(1e6))*100} %,  {(np.sum(part_0_E.state<1)/int(1e6))*100} %")
                    print(f"abs_heavy_corner pos {(np.sum(part_pos.state < 1)/int(1e6)*100)} %, {(np.sum(part_pos_E.state<1)/int(1e6))*100} %")
                    print(f"abs_heavy_corner neg {np.sum(part_neg.state < 1)/int(1e6)} %, {np.sum(part_neg_E.state<1)/int(1e6)} %")
                    print(f"abs_heavy_corner rand {np.sum(part_rand.state < 1)/int(1e6)} %, {np.sum(part_rand_E.state<1)/int(1e6)} %")
                    print(f"abs_heavy_corner {(np.sum(part_neg.state < 1) + np.sum(part_pos.state < 1) + np.sum(part_rand.state < 1) + (np.sum(part_0.state < 1))/(4*int(1e6)))*100} %")
                    print(f"abs heavy corner E {(np.sum(part_neg_E.state < 1) + np.sum(part_pos_E.state < 1) + np.sum(part_rand_E.state < 1) + (np.sum(part_0_E.state < 1))/(4*int(1e6)))*100} %")
                elif idx == 4:
                    abs_heavy_drift[0] += np.sum(part_neg.state < 1) + np.sum(part_pos.state < 1) + np.sum(part_rand.state < 1) + np.sum(part_0.state < 1)
                    abs_heavy_drift[1] += np.sum(part_neg_E.state < 1) + np.sum(part_pos_E.state < 1) + np.sum(part_rand_E.state < 1) + np.sum(part_0_E.state < 1)
                    print(f"abs_heavy_drift 0 {(np.sum(part_0.state < 1)/int(1e6))*100} %,  {(np.sum(part_0_E.state<1)/int(1e6))*100} %")
                    print(f"abs_heavy_drift pos {(np.sum(part_pos.state < 1)/int(1e6)*100)} %, {(np.sum(part_pos_E.state<1)/int(1e6))*100} %")
                    print(f"abs_heavy_drift neg {(np.sum(part_neg.state < 1)/int(1e6)*100)} %, {(np.sum(part_neg_E.state<1)/int(1e6)*100)} %")
                    print(f"abs_heavy_drift rand {np.sum(part_rand.state < 1)/int(1e6)} %, {np.sum(part_rand_E.state<1)/int(1e6)} %")
                    print(f"abs_heavy_drift {(np.sum(part_neg.state < 1) + np.sum(part_pos.state < 1) + np.sum(part_rand.state < 1) + (np.sum(part_0.state < 1))/(4*int(1e6)))*100} %")
                    print(f"abs heavy drift E {(np.sum(part_neg_E.state < 1) + np.sum(part_pos_E.state < 1) + np.sum(part_rand_E.state < 1) + (np.sum(part_0_E.state < 1))/(4*int(1e6)))*100} %")

            del part_0, part_pos, part_neg, part_rand
            del part_0_E, part_pos_E, part_neg_E, part_rand_E
    
    assert np.all(np.greater(abs_light_inside, abs_light_drift))
    assert np.all(np.greater(abs_light_inside , abs_light_corner))
    assert np.all(np.greater(abs_middle_inside, abs_middle_drift))
    assert np.all(np.greater(abs_middle_inside, abs_middle_corner))
    assert np.all(np.greater(abs_heavy_inside , abs_heavy_drift))
    assert np.all(np.greater(abs_heavy_inside , abs_heavy_corner))

    assert np.all(np.greater(abs_heavy_inside , abs_middle_inside))
    assert np.all(np.greater(abs_middle_inside, abs_light_inside))
    assert np.all(np.greater(abs_heavy_corner , abs_light_corner))


@pytest.mark.parametrize("beam, plane",[[1,'V']]) 
def test_everest_and_K2_histogram(beam, plane):
    pos    = [0.0015, 0.00131, 0.001305, 0.00129, 0.0011] 

    K2Coll = _K2Collimator(length=0.6, jaw=0.0013, material='MoGR', angle=90.0, emittance=3.5e-6)
    EverestColl = xc.EverestCollimator(length=0.6, jaw=0.0013, material=xc.materials.MolybdenumGraphite,angle=90.0, emittance=3.5e-6) 
    for idx, i in enumerate(pos):
        part_0, part_pos, part_neg, part_rand,zero_init,pos_init,neg_init, rand_init,_,_ = _track_with_angles(beam=beam, plane=plane, pos=i, material=True, K2coll=[K2Coll,'tcp.b6l7.b1'])
        part_0_E, part_pos_E, part_neg_E, part_rand_E,zero_initE,pos_initE,neg_initE, rand_initE,_,_ = _track_with_angles(beam, plane, pos=i, material=True, everest=True, Ecoll=[EverestColl,'tcp.b6l7.b1'])

        mask_0     = part_0.state > 0
        ids0       = part_0.particle_id[mask_0]
        mask_pos   = part_pos.state > 0
        ids_pos    = part_pos.particle_id[mask_pos]
        mask_neg   = part_neg.state > 0
        ids_neg    = part_neg.particle_id[mask_neg]
        mask_rand  = part_rand.state > 0
        ids_rand   = part_rand.particle_id[mask_rand]
        mask_0E    = part_0_E.state > 0
        ids0E      = part_0_E.particle_id[mask_0E]
        mask_posE  = part_pos_E.state > 0
        ids_posE   = part_pos_E.particle_id[mask_posE]
        mask_negE  = part_neg_E.state > 0
        ids_negE   = part_neg_E.particle_id[mask_negE]
        mask_randE = part_rand_E.state > 0
        ids_randE  = part_rand_E.particle_id[mask_randE]

        # histograms angular distribution
        hist_part_0,_      = np.histogram(part_0.kin_yprime[mask_0] - zero_init.kin_yprime[ids0], bins=500)
        hist_part_pos,_    = np.histogram(part_pos.kin_yprime[mask_pos] - pos_init.kin_yprime[ids_pos], bins=500)
        hist_part_neg,_    = np.histogram(part_neg.kin_yprime[mask_neg] - neg_init.kin_yprime[ids_neg], bins=500)
        hist_part_rand,_   = np.histogram(part_rand.kin_yprime[mask_rand] - rand_init.kin_yprime[ids_rand], bins=500)
        hist_part_0_E,_    = np.histogram(part_0_E.kin_yprime [mask_0E]- zero_initE.kin_yprime[ids0E], bins=500)
        hist_part_pos_E,_  = np.histogram(part_pos_E.kin_yprime[mask_posE]- pos_initE.kin_yprime[ids_posE], bins=500)
        hist_part_neg_E,_  = np.histogram(part_neg_E.kin_yprime[mask_negE]- neg_initE.kin_yprime[ids_negE], bins=500) 
        hist_part_rand_E,_ = np.histogram(part_rand_E.kin_yprime[mask_randE] - rand_initE.kin_yprime[ids_randE], bins=500)

        # histograms energy distribution
        hist_part_0_energy,_     = np.histogram(part_0.energy[mask_0], bins=500)
        hist_part_pos_energy,_   = np.histogram(part_pos.energy[mask_pos], bins=500)
        hist_part_neg_energy,_   = np.histogram(part_neg.energy[mask_neg], bins=500)
        hist_part_rand_energy,_  = np.histogram(part_rand.energy[mask_rand], bins=500)
        hist_part_0_energy_E,_   = np.histogram(part_0_E.energy[mask_0E], bins=500)
        hist_part_pos_energy_E,_ = np.histogram(part_pos_E.energy[mask_posE], bins=500)
        hist_part_neg_energy_E,_ = np.histogram(part_neg_E.energy[mask_negE], bins=500)
        hist_part_rand_energy_E,_= np.histogram(part_rand_E.energy[mask_randE], bins=500)

        # moments 
        moments          = _get_moment(hist_part_0, hist_part_pos, hist_part_neg, hist_part_rand)
        moments_E        = _get_moment(hist_part_0_E, hist_part_pos_E, hist_part_neg_E, hist_part_rand_E)
        moments_energy   = _get_moment(hist_part_0_energy, hist_part_pos_energy, hist_part_neg_energy, hist_part_rand_energy)
        moments_energy_E = _get_moment(hist_part_0_energy_E, hist_part_pos_energy_E, hist_part_neg_energy_E, hist_part_rand_energy_E)
    
        # check if they are the same
        print(f" K2 moments {moments}, Everest moments {moments_E}")
        print(f" K2 energy moments {moments_energy}, Everest energy moments {moments_energy_E}")
       
        assert _are_numbers_equal_within_tolerance(moments[0], moments_E[0], 1), f"{np.mean(hist_part_0), np.mean(hist_part_0_E), np.var(hist_part_0), np.var(hist_part_0_E)}"
        assert _are_numbers_equal_within_tolerance(moments[1], moments_E[1], 1), f"angle pos, turn {i}, {moments[1], moments_E[1]}"
        assert _are_numbers_equal_within_tolerance(moments[2], moments_E[2], 1), f"angle neg, turn {i}, {moments[2], moments_E[2]}"
        assert _are_numbers_equal_within_tolerance(moments[3], moments_E[3], 1), f"angle rand, turn {i}, {moments[3], moments_E[3]}"

        assert _are_numbers_equal_within_tolerance(moments_energy[0], moments_energy_E[0], 1), f"energy, turn {i}, {moments_energy[0], moments_energy_E[0]}"
        assert _are_numbers_equal_within_tolerance(moments_energy[1], moments_energy_E[1], 1), f"energy, turn {i}, {moments_energy[1], moments_energy_E[1]}"
        assert _are_numbers_equal_within_tolerance(moments_energy[2], moments_energy_E[2], 1), f"energy, turn {i}, {moments_energy[2], moments_energy_E[2]}"
        assert _are_numbers_equal_within_tolerance(moments_energy[3], moments_energy_E[3], 1), f"energy, turn {i}, {moments_energy[3], moments_energy_E[3]}"

        assert _are_numbers_equal_within_tolerance(part_0.energy, part_0_E.energy,3), f"energy, turn {i}, {part_0.energy, part_0_E.energy}"
        assert _are_numbers_equal_within_tolerance(part_pos.energy, part_pos_E.energy,3), f"energy, turn {i}, {part_pos.energy, part_pos_E.energy}"
        assert _are_numbers_equal_within_tolerance(part_neg.energy, part_neg_E.energy,3), f"energy, turn {i}, {part_neg.energy, part_neg_E.energy}"
        assert _are_numbers_equal_within_tolerance(part_rand.energy, part_rand_E.energy, 3), f"energy, turn {i}, {part_rand.energy, part_rand_E.energy}"