# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2026.                 #
# ######################################### #

import math
import numpy as np
from .isotopes import ISOTOPES
from .glauber_gribov import (
    load_all_splines, make_nucleon_cs,
    glauber_element_single, GG_KEYS,
)

# Map element symbol -> Z
SYMBOL_TO_Z = {
    "H":1,"He":2,"Li":3,"Be":4,"B":5,"C":6,"N":7,"O":8,"F":9,"Ne":10,
    "Na":11,"Mg":12,"Al":13,"Si":14,"P":15,"S":16,"Cl":17,"Ar":18,
    "K":19,"Ca":20,"Sc":21,"Ti":22,"V":23,"Cr":24,"Mn":25,"Fe":26,
    "Co":27,"Ni":28,"Cu":29,"Zn":30,"Ga":31,"Ge":32,"As":33,"Se":34,
    "Br":35,"Kr":36,"Rb":37,"Sr":38,"Y":39,"Zr":40,"Nb":41,"Mo":42,
    "Tc":43,"Ru":44,"Rh":45,"Pd":46,"Ag":47,"Cd":48,"In":49,"Sn":50,
    "Sb":51,"Te":52,"I":53,"Xe":54,"Cs":55,"Ba":56,"La":57,"Ce":58,
    "Pr":59,"Nd":60,"Pm":61,"Sm":62,"Eu":63,"Gd":64,"Tb":65,"Dy":66,
    "Ho":67,"Er":68,"Tm":69,"Yb":70,"Lu":71,"Hf":72,"Ta":73,"W":74,
    "Re":75,"Os":76,"Ir":77,"Pt":78,"Au":79,"Hg":80,"Tl":81,"Pb":82,
    "Bi":83,"Po":84,"At":85,"Rn":86,"Fr":87,"Ra":88,"Ac":89,"Th":90,
    "Pa":91,"U":92,"Np":93,"Pu":94,
}


def compute_gg_arrays(A, Z, sqrt_s_arr, cs_hN):
    rows = [glauber_element_single(A, Z, s, cs_hN) for s in sqrt_s_arr]
    return {key: np.array([r[key] for r in rows]) for key in GG_KEYS}


def isotope_weighted_gg(symbol, Z, sqrt_s_arr, cs_hN):
    """
    Compute the isotope-abundance-weighted GG cross sections for an element.
    Only stable isotopes (abundance is not None) are included.
    Returns dict {key: np.array}, or None if no stable isotopes found.
    """
    isotopes = [iso for iso in ISOTOPES[symbol]["isotopes"]
                if iso["abundance"] is not None]
    if not isotopes:
        return None

    # Normalise abundances (they should sum to 1 but do so defensively)
    total = sum(iso["abundance"] for iso in isotopes)
    combined = {k: np.zeros(len(sqrt_s_arr)) for k in GG_KEYS}
    for iso in isotopes:
        frac = iso["abundance"] / total
        gg   = compute_gg_arrays(iso["atomic_mass"], Z, sqrt_s_arr, cs_hN)
        for k in GG_KEYS:
            combined[k] += frac * gg[k]

    return combined


def relative_difference(arr_a, arr_b):
    denom = (np.abs(arr_a) + np.abs(arr_b)) / 2.0
    with np.errstate(invalid='ignore', divide='ignore'):
        rel = np.where(denom > 0, np.abs(arr_a - arr_b) / denom, 0.0)
    return rel

def check_isotope_approximation(threshold=0.10, n_points=10000, cs_key="cs_inel_hA"):
    """
    For each element in the ISOTOPES table, compare the GG cross section
    computed from the standard atomic weight A against the isotope-abundance-
    weighted sum of per-isotope GG cross sections.
 
    Parameters
    ----------
    threshold : float
        Flag elements whose max relative difference exceeds this value.
        Default 0.10 (10%).
    n_points : int
        Number of energy grid points. Default 1000.
    cs_key : str
        Which GG cross section to compare. One of GG_KEYS.
        Default 'cs_inel_hA'.
 
    Returns
    -------
    results : list of dict
        One dict per element with keys:
            symbol, Z, A_std, max_rel_diff, mean_rel_diff, flagged
    flagged : list of dict
        Subset of results where flagged is True.
    """
    splines, grid_min, grid_max = load_all_splines()
    sqrt_s_arr = np.logspace(math.log10(grid_min), math.log10(grid_max), n_points)
    _, _, _, _, cs_hN = make_nucleon_cs(splines)
 
    print(f"\nComparing GG cross sections: standard A vs isotope-weighted sum")
    print(f"Cross section: {cs_key}")
    print(f"Flag threshold: {threshold*100:.1f}%")
 
    print(f"{'Symbol':<6}  {'Z':>4}  {'A_std':>10}  "
          f"{'Max rel diff':>14}  {'Mean rel diff':>14}  {'Flag'}")
    print("-" * 65)
 
    results = []
    flagged = []
 
    for symbol, data in ISOTOPES.items():
        Z = SYMBOL_TO_Z.get(symbol)
        if Z is None:
            continue
 
        A_std = data["standard_atomic_weight"]
 
        # Method 1: standard atomic weight
        gg_std = compute_gg_arrays(A_std, Z, sqrt_s_arr, cs_hN)
 
        # Method 2: isotope-weighted sum
        gg_iso = isotope_weighted_gg(symbol, Z, sqrt_s_arr, cs_hN)
        if gg_iso is None:
            print(f"{symbol:<6}  {Z:>4}  {A_std:>10.4f}  "
                  f"{'no stable isotopes':>30}")
            continue
 
        rel      = relative_difference(gg_std[cs_key], gg_iso[cs_key])
        max_rel  = float(rel.max())
        mean_rel = float(rel.mean())
        is_flagged = max_rel > threshold
 
        print(f"{symbol:<6}  {Z:>4}  {A_std:>10.4f}  "
              f"{max_rel*100:>13.4f}%  {mean_rel*100:>13.4f}%  "
              f"{'*** FLAG ***' if is_flagged else ''}")
 
        row = {
            "symbol":        symbol,
            "Z":             Z,
            "A_std":         A_std,
            "max_rel_diff":  max_rel,
            "mean_rel_diff": mean_rel,
            "flagged":       is_flagged,
        }
        results.append(row)
        if is_flagged:
            flagged.append(row)
 
    print()
    if flagged:
        print(f"Flagged elements (max rel diff > {threshold*100:.1f}%):")
        for row in flagged:
            print(f"  {row['symbol']} (Z={row['Z']}, A={row['A_std']:.4f}): "
                  f"{row['max_rel_diff']*100:.4f}%")
    else:
        print(f"No elements flagged above {threshold*100:.1f}% threshold.")
 
    return results, flagged