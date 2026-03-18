# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2026.                 #
# ######################################### #

"""
materials_cross_section.py
--------------------------
Stage 2: combine elemental cross section arrays (from stage1_elements.py)
into material cross sections using the Material class's own molar_fractions
and components properties.

This module is a pure library — no file I/O, no main.  It is imported by
the Material class and by stage1_elements.py.

The Material class already handles all composition logic internally:
  - Elements    : mat.is_elemental  →  use mat.Z, mat.A directly
  - Compounds   : mat.is_compound   →  mat.molar_fractions, mat.components
  - Mixtures    : mat.is_mixture    →  mat.molar_fractions, mat.components
  - Nesting     : the Material class fully flattens nested compounds in
                  _resolve_n_atoms / _resolve_mass_fractions at construction
                  time, so mat.components always contains pure elements.

We therefore read mat.molar_fractions and mat.components directly —
no custom flattening code needed here.

Typical usage (from the Material class)
----------------------------------------
    from .materials_cross_section import load_element_cache, get_material_cs

    element_cs, sqrt_s_arr = load_element_cache("elements_cs.npz")
    cs = get_material_cs(self, element_cs, sqrt_s_arr)

    # cs["cs_inel_hA"]  ->  np.array [mb], ready to pass to C
"""

import numpy as np

from .glauber_gribov import GG_KEYS, NUCLEON_KEYS


def load_element_cache(filename="elements_cs.npz"):
    """
    Load the NPZ produced by stage1_elements.py.

    Returns
    -------
    element_cs : dict
        element_cs[0]  = {nucleon_key: np.array, ...}   (pp/pn, shared)
        element_cs[Z]  = {"A": float, gg_key: np.array, ...}
    sqrt_s_arr : np.array
        The energy grid [GeV].
    """
    data       = np.load(filename)
    sqrt_s_arr = data["sqrt_s"]
    Zs         = data["element_Z"].astype(int)
    As         = data["element_A"]

    element_cs = {}

    # Nucleon-level arrays — element-independent, stored under Z=0
    element_cs[0] = {key: data[f"nucleon_{key}"] for key in NUCLEON_KEYS}

    # Per-element GG arrays
    for Z, A in zip(Zs, As):
        element_cs[Z] = {"A": float(A)}
        for key in GG_KEYS:
            element_cs[Z][key] = data[f"Z{Z}_{key}"]

    return element_cs, sqrt_s_arr


def get_material_cs(mat, element_cs, sqrt_s_arr):
    """
    Combine pre-computed elemental GG arrays for one Material instance.

    Uses mat.molar_fractions and mat.components directly — the Material
    class already flattens compounds and mixtures to pure elements at
    construction time, so no custom recursion is needed here.

    For a pure element, mat.Z and mat.A are used directly.

    Parameters
    ----------
    mat : Material
        An xcoll Material instance (element, compound, or mixture).
    element_cs : dict
        Cache returned by load_element_cache().
    sqrt_s_arr : np.array
        Energy grid [GeV], also from load_element_cache().

    Returns
    -------
    dict with keys:
        sqrt_s        np.array [GeV]   energy grid
        cs_tot_hA     np.array [mb]    total nucleus cross section
        cs_inel_hA    np.array [mb]    inelastic
        cs_el_hA      np.array [mb]    elastic
        cs_prod_hA    np.array [mb]    production
        cs_sd_hA      np.array [mb]    single diffractive
        cs_qel_hA     np.array [mb]    quasi-elastic
        cs_tot_pp     np.array [mb]    pp total   (nucleon level, same for all)
        cs_el_pp      np.array [mb]    pp elastic
        cs_inel_pp    np.array [mb]    pp inelastic
        cs_tot_pn     np.array [mb]    pn total

    Raises
    ------
    KeyError    if a required element Z is missing from element_cs
                (re-run stage1_elements.py to rebuild the cache)
    """
    if mat.is_elemental:
        # Single element: one component, weight = 1
        pairs = [(int(mat.Z), 1.0)]
    else:
        # Compound or mixture: Material class gives us molar_fractions
        # and components already flattened to pure elements
        pairs = [
            (int(el.Z), float(frac))
            for el, frac in zip(mat.components, mat.molar_fractions)
        ]

    # Weighted sum of elemental GG arrays
    combined = {k: np.zeros(len(sqrt_s_arr)) for k in GG_KEYS}
    for Z, frac in pairs:
        if Z not in element_cs:
            raise KeyError(
                f"Z={Z} not found in element cache. "
                "Re-run stage1_elements.py to rebuild the cache."
            )
        for k in GG_KEYS:
            combined[k] += frac * element_cs[Z][k]

    return {
        "sqrt_s": sqrt_s_arr,
        **combined,
        **{k: element_cs[0][k] for k in NUCLEON_KEYS},
    }