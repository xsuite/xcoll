import numpy as np


def calculate_scattering(p0,anuc,rho,zatom,emr,csref0,csref1,csref5,bnref):
    cprob = np.array([0,0,0,0,0,0], dtype=np.float64)
    xintl = 0
    bn = 0
    return cprob, xintl, bn
