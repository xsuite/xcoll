# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

import numpy as np
import scipy.constants as sc


# Be careful! Changing those default values will change all materials that do not set them explicitly!
# Hence, when loading jsons saved before the change, the materials will change!

_default_excitation_energies = {
    1: 19.2, 2: 41.8, 3: 40, 4: 63.7, 5: 76, 6: 81, 7: 82, 8: 95, 9: 115, 10: 137, 11: 149, 12: 156, 13: 166,
    14: 173, 15: 173, 16: 180, 17: 174, 18: 188, 19: 190, 20: 191, 21: 216, 22: 233, 23: 245, 24: 257, 25: 272,
    26: 286, 27: 297, 28: 311, 29: 322, 30: 330, 31: 334, 32: 350, 33: 347, 34: 348, 35: 343, 36: 352, 37: 363,
    38: 366, 39: 379, 40: 393, 41: 417, 42: 424, 43: 428, 44: 441, 45: 449, 46: 470, 47: 470, 48: 469, 49: 488,
    50: 488, 51: 487, 52: 485, 53: 491, 54: 482, 55: 488, 56: 491, 57: 501, 58: 523, 59: 535, 60: 546, 61: 560,
    62: 574, 63: 580, 64: 591, 65: 614, 66: 628, 67: 650, 68: 658, 69: 674, 70: 684, 71: 694, 72: 705, 73: 718,
    74: 727, 75: 736, 76: 746, 77: 757, 78: 790, 79: 790, 80: 800, 81: 810, 82: 823, 83: 823, 84: 830, 85: 825,
    86: 794, 87: 827, 88: 826, 89: 841, 90: 847, 91: 878, 92: 890, 93: 902, 94: 921, 95: 934, 96: 939, 97: 952,
    98: 966, 100: 994, 116: 1213
}

def _approximate_radiation_length(Z, A, density):
    # Tsai's more accurate expression (PDG'24 Sec 34.4.2)
    match int(round(Z)):
        case 1:  # Hydrogen
            L_rad = 5.31
            L_rad_prime = 6.144
        case 2:  # Helium
            L_rad = 4.79
            L_rad_prime = 5.621
        case 3:  # Lithium
            L_rad = 4.74
            L_rad_prime = 5.805
        case 4:  # Beryllium
            L_rad = 4.71
            L_rad_prime = 5.924
        case _:
            L_rad = np.log(184.15 * Z ** (-1.0 / 3.0))
            L_rad_prime = np.log(1194.0 * Z ** (-2.0 / 3.0))
    denom = Z**2 * (L_rad - _f_coulomb(Z)) + Z * L_rad_prime
    X0_mass = 716.405 * A / denom  # in g/cm^2
    return X0_mass / density / 100  # in meters

def _f_coulomb(Z):
    # Coulomb correction f(Z) from PDG.
    a = sc.alpha * Z
    a2 = a*a
    return a2 * (
        1/(1+a2)
        + 0.20206
        - 0.0369*a2
        + 0.0083*a2*a2
        - 0.0020*a2*a2*a2
    )


# --- Mixture rules ------------------------------------------------------------

def _average_Z_over_A(components, mass_fractions):
    return sum([fr * el.Z / el.A for el, fr in zip(components,
                                                   mass_fractions)])

def _effective_Z2(components, molar_fractions):
    return np.sum([mf * el.Z**2 for el, mf in zip(components,
                                                  molar_fractions)])

def _combine_radiation_lengths(components, mass_fractions, density):
    # radiation_length in m
    X0_mass = [el.radiation_length * 100 * el.density for el in components]  # g/cm^2
    X0_mass_mix = _inverse_weighted_mean(X0_mass, mass_fractions)            # g/cm^2
    return X0_mass_mix / density / 100                                       # m

def _combine_excitation_energies(components, mass_fractions):
    fractions = np.array([fr * el.Z / el.A for el, fr in zip(components,
                                                             mass_fractions)])
    return _logarithmic_mean([el.excitation_energy for el in components],
                              fractions)



def _inverse_weighted_mean(vals, weights):
    return 1 / (np.array(weights) / np.array(vals)).sum()

def _logarithmic_mean(vals, weights):
    res = np.sum(np.array(weights) * np.log(np.array(vals))) / np.sum(weights)
    return np.exp(res)
