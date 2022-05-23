import numpy as np
from pathlib import Path
from xcoll import _pkg_root

try:
    import xcoll.beam_elements.pyk2 as pyk2
except ImportError:
    raise Exception("Error: Failed importing pyK2 (did you compile?). Cannot track.")
    

__pseudo_rands__ = np.load(Path(_pkg_root,'beam_elements','k2','randoms.npy'))
__index__ = -1
__n_ruths__ = 0



def get_random():
    global __index__
    __index__ += 1
    # a = pyk2.pyk2_rand()
    thisran = __pseudo_rands__[__index__]
    # if abs(a-thisran) > 0:
    #     print(__index__, a, thisran)
    return thisran


def get_random_ruth(cgen):
    
    global __n_ruths__
    __n_ruths__ += 1

    ran = get_random()
    gap=0.010101010091602802 #1./99.
    gapinv=99.
    tleft= 2./99
    bright=97./99.
    gaps=tleft/49.
    gapins=1./gaps
    j = int(ran*gapinv) + 1
    if j < 3:
        j1 = int(ran *gapins)
        j = j1 + 101
        j = max(j,102)
        j = min(j,148)
        p = (ran -gaps*(j1-1)) * gapins
        a = (p+1.0) * cgen[j+1] - (p-2.0)*cgen[j-2]
        b = (p-1.0) * cgen[j-1] - p * cgen[j]
       
    elif j > 97:
        j1 = int((ran-bright)*gapins)
        j = j1 + 151
        j = max(j,152)
        j = min(j,198)
        p = ((ran-bright) -gaps*(j1-1)) * gapins
        a = (p+1.0) * cgen[j+1] - (p-2.0)*cgen[j-2]
        b = (p-1.0) * cgen[j-1] - p * cgen[j]
       
    else:
        p = (ran -gap*(j-1)) * gapinv
        a = (p+1.) * cgen[j+1] - (p-2.)*cgen[j-2]
        b = (p-1.) * cgen[j-1] - p * cgen[j]
    ran = ((a*p)*(p-1.))*0.16666667 + ((b*(p+1.))*(p-2.))*0.5
       
    return ran




