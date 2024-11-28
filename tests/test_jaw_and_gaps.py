# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

from pathlib import Path
import pytest
import numpy as np

import xobjects as xo
import xtrack as xt
import xcoll as xc
from xpart.test_helpers import flaky_assertions, retry
from xobjects.test_helpers import for_all_test_contexts

path = Path(__file__).parent / 'data'


@for_all_test_contexts(
    excluding=('ContextCupy', 'ContextPyopencl')  # Rutherford RNG not on GPU
)
@pytest.mark.parametrize("beam", [1, 2], ids=["B1", "B2"])
def test_gaps(beam, test_context):
    line = xt.Line.from_json(path / f'sequence_lhc_run3_b{beam}.json')
    coll = xc.BlackAbsorber(length=1.738, angle=127.5)
    name = 'tcp.b6l7.b1' if beam == 1 else 'tcp.b6r7.b2'
    line.collimators.install(name, coll, need_apertures=True)
    line.build_tracker()
    tw = line.twiss()
    beta_gamma_rel = line.particle_ref._xobject.gamma0[0]*line.particle_ref._xobject.beta0[0]
    coll.assign_optics(name=name, nemitt_x=3.5e-6, nemitt_y=2.5e-6, twiss=tw, beta_gamma_rel=beta_gamma_rel)
    coll.jaw_L = None
    coll.jaw_R = None
    coll.gap_L = None
    coll.gap_R = None
    coll.angle = 127.5
    coll.tilt = 0.0

    # gap L
    newval_gapL = np.array([[2.3, -0.0013],[2.3, 0.0016],[2.3, -2.4],[2.3, 98],[2.3, np.pi/5]]) 
    second_gapL = np.array(['jaw_R', 'jaw_L', 'gap_R', 'angle', 'tilt']) 
    for r, i in enumerate(second_gapL):
        check_3_times(coll, newval_gapL, np.array(['gap_L', i]), r)

    # gap R
    newval_gapR = np.array([[-1.3, -0.00024],[-1.3, 0.0033],[-1.3, 3.3],[-1.3, 122],[-1.3, np.pi/6]]) 
    second_gapR = np.array(['jaw_R', 'jaw_L', 'gap_L', 'angle', 'tilt'])
    for r, i in enumerate(second_gapR):
        check_3_times(coll, newval_gapR, np.array(['gap_R', i]), r)

    # jaw L
    newval_jawL = np.array([[0.0036, -0.0028],[0.0036, -4.4],[0.0036, 2.1],[0.0036, 130],[0.0036, np.pi/9]]) 
    second_jawL = np.array(['jaw_R', 'gap_R', 'gap_L', 'angle', 'tilt'])
    for r, i in enumerate(second_jawL):
        check_3_times(coll, newval_jawL, np.array(['jaw_L', i]), r)

    # jaw R
    newval_jawR = np.array([[-0.0011, -4.6],[-0.0011, 0.0015],[-0.0011, 4.5],[-0.0011, 89],[-0.0011, np.pi/5]]) 
    second_jawR = np.array(['gap_R', 'jaw_L', 'gap_L', 'angle', 'tilt'])
    for r, i in enumerate(second_jawR):
        check_3_times(coll, newval_jawR, np.array(['jaw_R', i]), r)

    # Tilt
    newval_tilt = np.array([[np.pi/10, -0.0038],[np.pi/10, 0.0030],[np.pi/10, 3.1],[np.pi/10, 111],[np.pi/10, -3.0]]) 
    second_tilt = np.array(['jaw_R', 'jaw_L', 'gap_L', 'angle', 'gap_R'])
    for r, i in enumerate(second_tilt):
        check_3_times(coll, newval_tilt, np.array(['tilt', i]), r)

    # Angle 
    newval_angle = np.array([[114, -0.0018],[114, 0.0018],[114, 3.7],[114, -3.6],[114, np.pi/7]]) 
    second_angle = np.array(['jaw_R', 'jaw_L', 'gap_L', 'gap_R', 'tilt'])
    for r, i in enumerate(second_angle):
        check_3_times(coll, newval_angle, np.array(['angle', i]), r)

def change_and_check(coll, parameter, val=None, is_third=False): # this is ugly
    check = False
    if is_third == True:
        if parameter == 'gap_L':
            coll.gap_L = 4.8
            check = np.isclose(coll.gap_L, 4.8, atol=1e-9)
        if parameter == 'gap_R': 
            coll.gap_R = -4.9
            check = np.isclose(coll.gap_R, -4.9, atol=1e-9)
        if parameter == 'jaw_L':
            coll.jaw_L = 0.0032
            check = np.isclose(coll.jaw_L, 0.0032, atol=1e-9)
        if parameter == 'jaw_R':
            coll.jaw_R = -0.0014
            check = np.isclose(coll.jaw_R, -0.0014, atol=1e-9)
        if parameter == 'angle':
            coll.angle = 86
            check = np.isclose(coll.angle, 86, atol=1e-9)
        if parameter == 'tilt':
            coll.tilt = np.pi/5
            check = np.isclose(coll.tilt, np.pi/5, atol=1e-9)
    else:
        if parameter == 'gap_L':
            coll.gap_L = val
            check = np.isclose(coll.gap_L, val, atol=1e-9)
        if parameter == 'gap_R': 
            coll.gap_R = val
            check = np.isclose(coll.gap_R, val, atol=1e-9)
        if parameter == 'jaw_L':
            coll.jaw_L = val
            check = np.isclose(coll.jaw_L, val, atol=1e-9)
        if parameter == 'jaw_R':
            coll.jaw_R = val
            check = np.isclose(coll.jaw_R, val, atol=1e-9)
        if parameter == 'angle':
            coll.angle = val
            check = np.isclose(coll.angle, val, atol=1e-9)
        if parameter == 'tilt':
            coll.tilt = val
            check = np.isclose(coll.tilt, val, atol=1e-9)
    return check

def check_3_times(coll, newval, succession, round):
    """
    newval = np.array([[origin, second],[orgin,second2],...]) 
    round = starts from 0, and describes where in we are in the sequence (each 2nd change has 5 possibilities/rounds)
    succession = describes the order (1st and 2nd) of the sequence (1st, 2nd, 3rd)
    third = describes what to look at for the 3rd change
    """
    third = np.array(['jaw_L', 'jaw_R', 'gap_L', 'gap_R', 'angle', 'tilt'])
    oldv = np.array([coll.jaw_L, coll.jaw_R, coll.gap_L, coll.gap_R, coll.angle, coll.tilt])                # TODO: it got stuck, so this is a quick fix
    only_R = False
    only_L = False

    # first change =================================================================================
    check = change_and_check(coll, succession[0], newval[round][0])
    assert check == True

    if succession[0] == 'gap_L' or succession[0] == 'jaw_L':
        assert np.isclose(coll.gap_L, (coll.jaw_L - coll.co[0][0]) / coll.sigma[0][0], atol=1e-9)
        assert np.isclose(coll.jaw_L, (coll.gap_L * coll.sigma[0][0]) + coll.co[0][0], atol=1e-9)

    elif succession[0] == 'gap_R' or succession[0] == 'jaw_R':
        assert np.isclose(coll.gap_R,(coll.jaw_R - coll.co[0][1]) / coll.sigma[0][1], atol=1e-9)
        assert np.isclose(coll.jaw_R, (coll.gap_R * coll.sigma[0][1]) + coll.co[0][1], atol=1e-9)

    if not succession[0] == 'tilt':
        assert np.isclose(coll.tilt, oldv[5], atol=1e-9)
    if not succession[0] == 'angle':
        assert np.isclose(coll.angle, oldv[4], atol=1e-9)

    third = np.delete(third, np.where(third == succession[0]))
    oldv[:] = [coll.jaw_L, coll.jaw_R, coll.gap_L, coll.gap_R, coll.angle, coll.tilt]

    # 2nd change =================================================================================
    check = change_and_check(coll, succession[1], newval[round][1])
    assert check == True

    if not succession[1] == 'tilt':                                             # check if tilt/angle not changed
        assert np.isclose(coll.tilt, oldv[5], atol=1e-9)
    if not succession[1] == 'angle':
        assert np.isclose(coll.angle, oldv[4], atol=1e-9)

    third = np.delete(third, np.where(third == succession[1]))
    if all(x not in ['angle', 'tilt'] for x in succession[:2]):                  # check if we ever had tilt/angle
        if (succession[0].endswith('L') and succession[1].endswith('L')):
            only_L = True
            assert coll.gap_R == None
            assert coll.jaw_R == None
            assert np.isclose(coll.jaw_L, (coll.gap_L * coll.sigma[0][0]) + coll.co[0][0], atol=1e-6)

        elif (succession[0].endswith('R') and succession[1].endswith('R')):
            only_R = True
            assert coll.gap_L == None
            assert coll.jaw_L == None
            assert np.isclose(coll.jaw_R, (coll.gap_R * coll.sigma[0][1]) + coll.co[0][1], atol=1e-6)
        else:
            assert coll.jaw_L > coll.jaw_R
            assert np.isclose(coll.jaw_R, (coll.gap_R * coll.sigma[0][1]) + coll.co[0][1], atol=1e-6)
            assert np.isclose(coll.jaw_L, (coll.gap_L * coll.sigma[0][0]) + coll.co[0][0], atol=1e-6)
    else:
        if succession[0].endswith('L') or succession[1].endswith('L'):
            only_L = True
        elif succession[0].endswith('R') or succession[1].endswith('R'):
            only_R = True
    oldv = [coll.jaw_L, coll.jaw_R, coll.gap_L, coll.gap_R, coll.angle, coll.tilt]

    # 3rd change =================================================================================
    check = change_and_check(coll, third[0], is_third=True)
    assert check == True
    if only_R == False and only_L == False:
        if not any(x in ['jaw_L', 'jaw_R', 'gap_R', 'gap_L'] for x in succession[:2]):
            if (third[0].endswith('L')):
                assert coll.gap_R == None
                assert coll.jaw_R == None
                assert np.isclose(coll.jaw_L, (coll.gap_L * coll.sigma[0][0]) + coll.co[0][0], atol=1e-6)
            elif (third[0].endswith('R')):
                assert coll.gap_L == None
                assert coll.jaw_L == None
                assert np.isclose(coll.jaw_R, (coll.gap_R * coll.sigma[0][1]) + coll.co[0][1], atol=1e-6)
            else:
                assert coll.jaw_L == None
                assert coll.jaw_R == None
                assert coll.gap_L == None
                assert coll.gap_R == None
        else:
            assert coll.jaw_L > coll.jaw_R
            assert np.isclose(coll.jaw_R, (coll.gap_R * coll.sigma[0][1]) + coll.co[0][1], atol=1e-6)
            assert np.isclose(coll.jaw_L, (coll.gap_L * coll.sigma[0][0]) + coll.co[0][0], atol=1e-6)
    if only_L == True:
        if third[0] == 'angle' or third[0] == 'tilt':
            assert coll.gap_R == None
            assert coll.jaw_R == None
            assert np.isclose(coll.jaw_L, (coll.gap_L * coll.sigma[0][0]) + coll.co[0][0], atol=1e-6)
        else:
            if third[0].endswith('L'):
                assert coll.gap_R == None
                assert coll.jaw_R == None
                assert np.isclose(coll.jaw_L, (coll.gap_L * coll.sigma[0][0]) + coll.co[0][0], atol=1e-6)
            else:
                assert coll.gap_R != None
                assert coll.jaw_R != None
                assert np.isclose(coll.jaw_R, (coll.gap_R * coll.sigma[0][1]) + coll.co[0][1], atol=1e-6)
    if only_R == True:
        if third[0] == 'angle' or third[0] == 'tilt':
            assert coll.gap_L == None
            assert coll.jaw_L == None
            assert np.isclose(coll.jaw_R, (coll.gap_R * coll.sigma[0][1]) + coll.co[0][1], atol=1e-6)
        else:
            if third[0].endswith('R'):
                assert coll.gap_L == None
                assert coll.jaw_L == None
                assert np.isclose(coll.jaw_R, (coll.gap_R * coll.sigma[0][1]) + coll.co[0][1], atol=1e-6)
            else:
                assert coll.gap_L != None
                assert coll.jaw_L != None
                assert np.isclose(coll.jaw_L, (coll.gap_L * coll.sigma[0][0]) + coll.co[0][0], atol=1e-6)
    if third[0] != 'tilt':
        assert coll.tilt == oldv[5]
    if third[0] != 'angle':
        assert coll.angle == oldv[4]
    coll.jaw_L = oldv[0]
    coll.jaw_R = oldv[1]
    coll.gap_L = oldv[2]
    coll.gap_R = oldv[3]
    coll.angle = oldv[4]
    coll.tilt = oldv[5]

    check = change_and_check(coll, third[1], is_third=True) # 2 
    assert check == True
    if only_R == False and only_L == False:
        if not any(x in ['jaw_L', 'jaw_R', 'gap_R', 'gap_L'] for x in succession[:2]):
            if (third[1].endswith('L')):
                assert coll.gap_R == None
                assert coll.jaw_R == None
                assert np.isclose(coll.jaw_L, (coll.gap_L * coll.sigma[0][0]) + coll.co[0][0], atol=1e-6)
            elif (third[1].endswith('R')):
                assert coll.gap_L == None
                assert coll.jaw_L == None
                assert np.isclose(coll.jaw_R, (coll.gap_R * coll.sigma[0][1]) + coll.co[0][1], atol=1e-6)
            else:
                assert coll.jaw_L == None
                assert coll.jaw_R == None
                assert coll.gap_L == None
                assert coll.gap_R == None
        else:
            assert coll.jaw_L > coll.jaw_R
            assert np.isclose(coll.jaw_R, (coll.gap_R * coll.sigma[0][1]) + coll.co[0][1], atol=1e-6)
            assert np.isclose(coll.jaw_L, (coll.gap_L * coll.sigma[0][0]) + coll.co[0][0], atol=1e-6)
    if only_L == True:
        if third[1] == 'angle' or third[1] == 'tilt':
            assert coll.gap_R == None
            assert coll.jaw_R == None
            assert np.isclose(coll.jaw_L, (coll.gap_L * coll.sigma[0][0]) + coll.co[0][0], atol=1e-6)
        else:
            if third[1].endswith('L'):
                assert coll.gap_R == None
                assert coll.jaw_R == None
                assert np.isclose(coll.jaw_L, (coll.gap_L * coll.sigma[0][0]) + coll.co[0][0], atol=1e-6)
            else:
                assert coll.gap_R != None
                assert coll.jaw_R != None
                assert np.isclose(coll.jaw_R, (coll.gap_R * coll.sigma[0][1]) + coll.co[0][1], atol=1e-6)
    if only_R == True:
        if third[1] == 'angle' or third[1] == 'tilt':
            assert coll.gap_L == None
            assert coll.jaw_L == None
            assert np.isclose(coll.jaw_R, (coll.gap_R * coll.sigma[0][1]) + coll.co[0][1], atol=1e-6)
        else:
            if third[1].endswith('R'):
                assert coll.gap_L == None
                assert coll.jaw_L == None
                assert np.isclose(coll.jaw_R, (coll.gap_R * coll.sigma[0][1]) + coll.co[0][1], atol=1e-6)
            else:
                assert coll.gap_L != None
                assert coll.jaw_L != None
                assert np.isclose(coll.jaw_L, (coll.gap_L * coll.sigma[0][0]) + coll.co[0][0], atol=1e-6)
    if third[1] != 'tilt':
        assert coll.tilt == oldv[5]
    if third[1] != 'angle':
        assert coll.angle == oldv[4]
    coll.jaw_L = oldv[0]
    coll.jaw_R = oldv[1]
    coll.gap_L = oldv[2]
    coll.gap_R = oldv[3]
    coll.angle = oldv[4]
    coll.tilt = oldv[5]

    check = change_and_check(coll, third[2], is_third=True) # 3 
    assert check == True
    if only_R == False and only_L == False:
        if not any(x in ['jaw_L', 'jaw_R', 'gap_R', 'gap_L'] for x in succession[:2]):
            if (third[2].endswith('L')):
                assert coll.gap_R == None
                assert coll.jaw_R == None
                assert np.isclose(coll.jaw_L, (coll.gap_L * coll.sigma[0][0]) + coll.co[0][0], atol=1e-6)
            elif (third[2].endswith('R')):
                assert coll.gap_L == None
                assert coll.jaw_L == None
                assert np.isclose(coll.jaw_R, (coll.gap_R * coll.sigma[0][1]) + coll.co[0][1], atol=1e-6)
            else:
                assert coll.jaw_L == None
                assert coll.jaw_R == None
                assert coll.gap_L == None
                assert coll.gap_R == None
        else:
            assert coll.jaw_L > coll.jaw_R
            assert np.isclose(coll.jaw_R, (coll.gap_R * coll.sigma[0][1]) + coll.co[0][1], atol=1e-6)
            assert np.isclose(coll.jaw_L, (coll.gap_L * coll.sigma[0][0]) + coll.co[0][0], atol=1e-6)
    if only_L == True:
        if third[2] == 'angle' or third[2] == 'tilt':
            assert coll.gap_R == None
            assert coll.jaw_R == None
            assert np.isclose(coll.jaw_L, (coll.gap_L * coll.sigma[0][0]) + coll.co[0][0], atol=1e-6)
        else:
            if third[2].endswith('L'):
                assert coll.gap_R == None
                assert coll.jaw_R == None
                assert np.isclose(coll.jaw_L, (coll.gap_L * coll.sigma[0][0]) + coll.co[0][0], atol=1e-6)
            else:
                assert coll.gap_R != None
                assert coll.jaw_R != None
                assert np.isclose(coll.jaw_R, (coll.gap_R * coll.sigma[0][1]) + coll.co[0][1], atol=1e-6)
    if only_R == True:
        if third[2] == 'angle' or third[2] == 'tilt':
            assert coll.gap_L == None
            assert coll.jaw_L == None
            assert np.isclose(coll.jaw_R, (coll.gap_R * coll.sigma[0][1]) + coll.co[0][1], atol=1e-6)
        else:
            if third[2].endswith('R'):
                assert coll.gap_L == None
                assert coll.jaw_L == None
                assert np.isclose(coll.jaw_R, (coll.gap_R * coll.sigma[0][1]) + coll.co[0][1], atol=1e-6)
            else:
                assert coll.gap_L != None
                assert coll.jaw_L != None
                assert np.isclose(coll.jaw_L, (coll.gap_L * coll.sigma[0][0]) + coll.co[0][0], atol=1e-6)
    if third[2] != 'tilt':
        assert coll.tilt == oldv[5]
    if third[2] != 'angle':
        assert coll.angle == oldv[4]
    coll.jaw_L = oldv[0]
    coll.jaw_R = oldv[1]
    coll.gap_L = oldv[2]
    coll.gap_R = oldv[3]
    coll.angle = oldv[4]
    coll.tilt = oldv[5]

    check = change_and_check(coll, third[3], is_third=True) # 3 
    assert check == True
    if only_R == False and only_L == False:
        if not any(x in ['jaw_L', 'jaw_R', 'gap_R', 'gap_L'] for x in succession[:2]):
            if (third[3].endswith('L')):
                assert coll.gap_R == None
                assert coll.jaw_R == None
                assert np.isclose(coll.jaw_L, (coll.gap_L * coll.sigma[0][0]) + coll.co[0][0], atol=1e-6)
            elif (third[3].endswith('R')):
                assert coll.gap_L == None
                assert coll.jaw_L == None
                assert np.isclose(coll.jaw_R, (coll.gap_R * coll.sigma[0][1]) + coll.co[0][1], atol=1e-6)
            else:
                assert coll.jaw_L == None
                assert coll.jaw_R == None
                assert coll.gap_L == None
                assert coll.gap_R == None
        else:
            assert coll.jaw_L > coll.jaw_R
            assert np.isclose(coll.jaw_R, (coll.gap_R * coll.sigma[0][1]) + coll.co[0][1], atol=1e-6)
            assert np.isclose(coll.jaw_L, (coll.gap_L * coll.sigma[0][0]) + coll.co[0][0], atol=1e-6)
    if only_L == True:
        if third[3] == 'angle' or third[3] == 'tilt':
            assert coll.gap_R == None
            assert coll.jaw_R == None
            assert np.isclose(coll.jaw_L, (coll.gap_L * coll.sigma[0][0]) + coll.co[0][0], atol=1e-6)
        else:
            if third[3].endswith('L'):
                assert coll.gap_R == None
                assert coll.jaw_R == None
                assert np.isclose(coll.jaw_L, (coll.gap_L * coll.sigma[0][0]) + coll.co[0][0], atol=1e-6)
            else:
                assert coll.gap_R != None
                assert coll.jaw_R != None
                assert np.isclose(coll.jaw_R, (coll.gap_R * coll.sigma[0][1]) + coll.co[0][1], atol=1e-6)
    if only_R == True:
        if third[3] == 'angle' or third[3] == 'tilt':
            assert coll.gap_L == None
            assert coll.jaw_L == None
            assert np.isclose(coll.jaw_R, (coll.gap_R * coll.sigma[0][1]) + coll.co[0][1], atol=1e-6)
        else:
            if third[3].endswith('R'):
                assert coll.gap_L == None
                assert coll.jaw_L == None
                assert np.isclose(coll.jaw_R, (coll.gap_R * coll.sigma[0][1]) + coll.co[0][1], atol=1e-6)
            else:
                assert coll.gap_L != None
                assert coll.jaw_L != None
                assert np.isclose(coll.jaw_L, (coll.gap_L * coll.sigma[0][0]) + coll.co[0][0], atol=1e-6)
    if third[3] != 'tilt':
        assert coll.tilt == oldv[5]
    if third[3] != 'angle':
        assert coll.angle == oldv[4]
    coll.jaw_L = None
    coll.jaw_R = None
    coll.gap_L = None
    coll.gap_R = None
    coll.angle = 127.5
    coll.tilt = 0.0