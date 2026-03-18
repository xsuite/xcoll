# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2026.                 #
# ######################################### #

"""
xcoll_cs_common.py
------------------
Shared physics used by both stage1_elements.py and stage2_materials.py:

  - PDG data loading and spline fitting
  - Nucleon-level cross sections (pp, pn)
  - Glauber-Gribov nucleus cross sections for a single element
  - Material composition flattening (element / compound / mixture)
"""

import math
from collections import defaultdict

import numpy as np
import xcoll as xc
from scipy.interpolate import UnivariateSpline


# ===========================================================================
# Physical constants
# ===========================================================================
MP = 0.938270   # proton mass [GeV]

GG_KEYS      = ["cs_tot_hA", "cs_inel_hA", "cs_el_hA",
                "cs_prod_hA", "cs_sd_hA",  "cs_qel_hA"]
NUCLEON_KEYS = ["cs_tot_pp", "cs_el_pp", "cs_inel_pp", "cs_tot_pn"]


# ===========================================================================
# PDG data loading and spline fitting
# ===========================================================================

def plab_to_sqrts(plab):
    """Convert lab-frame momentum [GeV/c] to sqrt(s) [GeV]."""
    plab = np.asarray(plab, dtype=float)
    E    = np.sqrt(plab**2 + MP**2)
    return np.sqrt(2*MP**2 + 2*MP*E)


def load_pdg_table(filename):
    """
    Read a PDG-format reaction data file.
    Column 1 (0-indexed) = plab [GeV/c],  column 4 = sigma [mb].
    Lines with fewer than 9 columns or non-numeric values are skipped.
    """
    plab, sigma = [], []
    with open(filename) as f:
        for line in f:
            parts = line.split()
            if len(parts) < 9:
                continue
            try:
                nums = [float(parts[i]) for i in range(9)]
            except ValueError:
                continue
            plab.append(nums[1])
            sigma.append(nums[4])
    if not plab:
        raise ValueError(f"No numeric data found in {filename!r}")
    return np.array(plab), np.array(sigma)


def build_spline(filename, smoothing):
    """
    Load a PDG file, convert plab -> sqrt(s), sort, fit a smoothing spline
    in log(sqrt_s) space.  Returns (spline, sqrts_array).
    """
    plab, sigma = load_pdg_table(filename)
    sqrts       = plab_to_sqrts(plab)
    idx         = np.argsort(sqrts)
    sqrts, sigma = sqrts[idx], sigma[idx]
    spline      = UnivariateSpline(np.log(sqrts), sigma, s=smoothing)
    return spline, sqrts


def load_all_splines():
    """
    Load the four PDG splines, print their data ranges, and return
    (splines_dict, grid_min, grid_max) where the grid bounds are the
    intersection of all data ranges.
    """
    datasets = [
        ("pp tot", "ppcrosstot.dat", 2000),
        ("pp el",  "ppcross.dat",    150),
        ("pn tot", "pncrosstot.dat", 140),
        ("pd tot", "pdcrosstot.dat", 200),
    ]

    print("Loading PDG spline tables ...")
    print(f"  {'Dataset':<10}  {'Points':>6}  {'sqrt(s) min [GeV]':>18}  "
          f"{'sqrt(s) max [GeV]':>18}  File")
    print(f"  {'-'*10}  {'-'*6}  {'-'*18}  {'-'*18}  {'-'*20}")

    splines      = {}
    sqrts_ranges = []
    for label, fname, s in datasets:
        spline, sqrts = build_spline(fname, s)
        splines[label] = spline
        sqrts_ranges.append((sqrts.min(), sqrts.max()))
        print(f"  {label:<10}  {len(sqrts):>6}  {sqrts.min():>18.6g}  "
              f"{sqrts.max():>18.6g}  {fname}")

    grid_min = max(lo for lo, hi in sqrts_ranges)
    grid_max = min(hi for lo, hi in sqrts_ranges)
    print(f"\n  Auto grid: sqrt(s) = {grid_min:.6g} to {grid_max:.6g} GeV")
    print("  (intersection of all spline data ranges)\n")

    return splines, grid_min, grid_max


# ===========================================================================
# Nucleon-level cross sections
# ===========================================================================

def make_nucleon_cs_fns(splines):
    """
    Return a bundle of functions that evaluate nucleon-level cross sections
    at a given sqrt(s), bound to the provided spline dict.
    """
    def cs_tot_pp(sqrt_s):
        return float(splines["pp tot"](np.log(sqrt_s)))

    def cs_el_pp(sqrt_s):
        return float(splines["pp el"](np.log(sqrt_s)))

    def cs_inel_pp(sqrt_s):
        return max(0.0, cs_tot_pp(sqrt_s) - cs_el_pp(sqrt_s))

    def cs_tot_pn(sqrt_s):
        """pn spline below 30 GeV; fall back to pp above (no world data)."""
        if sqrt_s < 30.0:
            return float(splines["pn tot"](np.log(sqrt_s)))
        return cs_tot_pp(sqrt_s)

    def cs_inel_pn(sqrt_s):
        return cs_inel_pp(sqrt_s)   # no pn inelastic data

    def cs_el_pn(sqrt_s):
        return cs_el_pp(sqrt_s)     # no pn elastic data

    def cs_hN(A, Z, sqrt_s):
        """
        Hadron-nucleon cross sections combining pp and pn with nuclear weights.
        Returns (cs_tot_hN, cs_inel_hN, cs_el_hN).
        """
        tot  = Z * cs_tot_pp(sqrt_s)  + (A - Z) * cs_tot_pn(sqrt_s)
        inel = Z * cs_inel_pp(sqrt_s) + (A - Z) * cs_inel_pn(sqrt_s)
        el   = Z * cs_el_pp(sqrt_s)   + (A - Z) * cs_el_pn(sqrt_s)
        return tot, inel, el

    return cs_tot_pp, cs_el_pp, cs_inel_pp, cs_tot_pn, cs_hN


# ===========================================================================
# Glauber-Gribov nucleus cross sections for one element at one energy
# ===========================================================================

def get_R(A):
    """Nuclear radius [m]."""
    if A > 21:
        return 1.1 * A**(1.0/3.0) * 1e-15 * 0.9
    return 1.1 * A**(1.0/3.0) * 1e-15 * 1.05

def pi_R2_mb(A):
    """pi*R^2 in mb  (1 fm^2 = 10 mb)."""
    R_fm = get_R(A) * 1e15
    return math.pi * R_fm**2 * 10.0

def glauber_element_single(A, Z, sqrt_s, cs_hN_fn):
    """
    GG nucleus cross sections [mb] for one element at one energy.

    Returns dict with keys: cs_tot_hA, cs_inel_hA, cs_el_hA,
                             cs_prod_hA, cs_sd_hA, cs_qel_hA.

    For A < 4, GG shadowing is not applied (H, He, Li have no nuclear
    shadowing — the nucleus is too small for multiple scattering).
    """
    piR2                      = pi_R2_mb(A)
    sig_tot, sig_inel, sig_el = cs_hN_fn(A, Z, sqrt_s)

    if A < 4:
        return {
            "cs_tot_hA":  sig_tot,
            "cs_inel_hA": sig_inel,
            "cs_el_hA":   sig_el,
            "cs_prod_hA": sig_inel,
            "cs_sd_hA":   0.0,
            "cs_qel_hA":  0.0,
        }

    cs_tot_hA  = 2 * piR2 * math.log(1.0 + sig_tot  / (2 * piR2))
    cs_inel_hA =     piR2 * math.log(1.0 + cs_tot_hA /      piR2)
    cs_el_hA   = max(0.0, cs_tot_hA - cs_inel_hA)
    cs_prod_hA =     piR2 * math.log(1.0 + sig_inel  /      piR2)
    alpha      = cs_tot_hA / (2 * piR2 + cs_tot_hA)
    cs_sd_hA   =     piR2 * (alpha - math.log(1.0 + alpha))
    cs_qel_hA  = max(0.0, cs_inel_hA - cs_prod_hA)

    return {
        "cs_tot_hA":  cs_tot_hA,
        "cs_inel_hA": cs_inel_hA,
        "cs_el_hA":   cs_el_hA,
        "cs_prod_hA": cs_prod_hA,
        "cs_sd_hA":   cs_sd_hA,
        "cs_qel_hA":  cs_qel_hA,
    }


# ===========================================================================
# Material composition flattening
# ===========================================================================

def _is_element(mat) -> bool:
    return (
        hasattr(mat, 'Z') and mat.Z is not None
        and not (hasattr(mat, 'components') and mat.components)
    )

def _effective_A(mat) -> float:
    """Mean atomic mass of a material, molar-fraction weighted."""
    return sum(A * f for (Z, A), f in flatten_material(mat).items())

def flatten_material(mat) -> dict:
    """
    Recursively reduce any xcoll Material to {(Z, A): molar_fraction}.
    All fractions sum to 1.

    Compound    (n_atoms=[1,2])         f_j = n_j / sum(n)
    Molar mix   (molar_fractions=[...]) used directly, normalised
    Mass mix    (mass_fractions=[...])  converted to molar via w_j/A_j_eff
    """
    if _is_element(mat):
        return {(int(mat.Z), float(mat.A)): 1.0}

    components = mat.components

    if hasattr(mat, 'n_atoms') and mat.n_atoms is not None:
        n = np.array(mat.n_atoms, dtype=float)
        mol_fracs = n / n.sum()

    elif hasattr(mat, 'molar_fractions') and mat.molar_fractions is not None:
        mol_fracs = np.array(mat.molar_fractions, dtype=float)
        mol_fracs = mol_fracs / mol_fracs.sum()

    elif hasattr(mat, 'mass_fractions') and mat.mass_fractions is not None:
        w      = np.array(mat.mass_fractions, dtype=float)
        A_effs = np.array([_effective_A(c) for c in components])
        num    = w / A_effs
        mol_fracs = num / num.sum()

    else:
        raise ValueError(
            f"Cannot determine composition of {mat!r}: "
            "no n_atoms, molar_fractions, or mass_fractions found."
        )

    combined: dict = defaultdict(float)
    for comp, frac in zip(components, mol_fracs):
        for (Z, A), child_frac in flatten_material(comp).items():
            combined[(Z, A)] += frac * child_frac

    total = sum(combined.values())
    return {k: v / total for k, v in combined.items()}


def material_composition(mat) -> list:
    """
    Return [(Z, A, molar_fraction), ...] sorted by Z.
    This is the only thing Stage 2 needs from the composition step.
    """
    return sorted(
        [(Z, A, f) for (Z, A), f in flatten_material(mat).items()],
        key=lambda x: x[0]
    )