# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2026.                 #
# ######################################### #

"""
stage1_elements.py
------------------
Stage 1: compute Glauber-Gribov cross sections for every element (Z=1..118)
on a log-spaced sqrt(s) grid and save to disk.

This is the physics-heavy step and only needs to be run once.  The output
NPZ file is then read by stage2_materials.py to build compound/mixture
cross sections without recomputing any GG physics.

Output
------
  elements_cs.npz
    Arrays stored under keys:
      "sqrt_s"              shape (n_points,)   [GeV]
      "Z{z}_cs_tot_hA"      shape (n_points,)   [mb]
      "Z{z}_cs_inel_hA"     shape (n_points,)   [mb]
      "Z{z}_cs_el_hA"       shape (n_points,)   [mb]
      "Z{z}_cs_prod_hA"     shape (n_points,)   [mb]
      "Z{z}_cs_sd_hA"       shape (n_points,)   [mb]
      "Z{z}_cs_qel_hA"      shape (n_points,)   [mb]
      "nucleon_cs_tot_pp"   shape (n_points,)   [mb]
      "nucleon_cs_el_pp"    shape (n_points,)   [mb]
      "nucleon_cs_inel_pp"  shape (n_points,)   [mb]
      "nucleon_cs_tot_pn"   shape (n_points,)   [mb]
      "element_Z"           shape (n_elements,) Z values present
      "element_A"           shape (n_elements,) corresponding A values

Usage
-----
  python stage1_elements.py
  python stage1_elements.py --n-points 500 --out elements_cs.npz
"""

import numpy as np
import xcoll as xc

from .glauber_gribov import (
    GG_KEYS, NUCLEON_KEYS,
    load_all_splines,
    make_nucleon_cs,
    glauber_element_single,
)

def find_all_elements():
    """
    Scan xc.materials for pure elements (have Z and A, no components).
    Returns a dict {Z: material_object} sorted by Z.
    """
    elements = {}
    for name in dir(xc.materials):
        mat = getattr(xc.materials, name)
        try:
            if (hasattr(mat, 'Z') and mat.Z is not None
                    and hasattr(mat, 'A') and mat.A is not None
                    and not (hasattr(mat, 'components') and mat.components)):
                Z = int(mat.Z)
                if Z not in elements:
                    elements[Z] = mat
        except Exception:
            continue
    return dict(sorted(elements.items()))
 
 
def resolve_materials(only=None):
    """
    Return {Z: material_object} for the requested materials.
 
    Parameters
    ----------
    only : list of str or None
        Material names e.g. ['Carbon', 'Tungsten'].
        If None, all elements in xc.materials are used.
    """
    if only is None:
        return find_all_elements()
 
    elements = {}
    for name in only:
        mat = getattr(xc.materials, name, None)
        if mat is None:
            raise ValueError(f"Material {name!r} not found in xc.materials.")
        if not (hasattr(mat, 'Z') and mat.Z is not None
                and hasattr(mat, 'A') and mat.A is not None
                and not (hasattr(mat, 'components') and mat.components)):
            raise ValueError(f"{name!r} is not a pure element.")
        Z = int(mat.Z)
        elements[Z] = mat
 
    return dict(sorted(elements.items()))
 
 
# ---------------------------------------------------------------------------
# Core computation
# ---------------------------------------------------------------------------
 
def _gg_isotope_weighted(Z, mat, sqrt_s_arr, cs_hN):
    """
    Compute isotope-abundance-weighted GG cross sections for one element.
 
    If the material has isotope data with known abundances, the result is:
        sigma = sum_i  f_i * GG(A_i, Z, sqrt_s)
    where f_i = abundance_i / sum(abundances).
 
    Falls back to mean A if no isotope data is available.
 
    Returns dict {key: np.array} for all GG_KEYS.
    """
    # Get stable isotopes (abundance is not None)
    stable = []
    if hasattr(mat, 'isotopes') and mat.isotopes:
        stable = [iso for iso in mat.isotopes if iso.abundance is not None]
 
    if stable:
        total    = sum(iso.abundance for iso in stable)
        combined = {k: np.zeros(len(sqrt_s_arr)) for k in GG_KEYS}
        for iso in stable:
            frac = iso.abundance / total
            rows = [glauber_element_single(iso.atomic_mass, Z, s, cs_hN)
                    for s in sqrt_s_arr]
            for k in GG_KEYS:
                combined[k] += frac * np.array([r[k] for r in rows])
        return combined
    else:
        # Fallback: use mean atomic weight
        A    = float(mat.A)
        rows = [glauber_element_single(A, Z, s, cs_hN) for s in sqrt_s_arr]
        return {k: np.array([r[k] for r in rows]) for k in GG_KEYS}
 
 
def compute_elements(sqrt_s_arr, splines, only=None, testing=False):
    """
    Compute isotope-weighted GG cross sections for elements over the energy grid.
 
    Parameters
    ----------
    sqrt_s_arr : np.array
        Energy grid [GeV].
    splines : dict
        Spline dict from load_all_splines().
    only : list of str or None
        Material names to compute, e.g. ['Carbon'].
        If None, all elements in xc.materials are computed.
 
    Returns
    -------
    element_cs : dict
        element_cs[0]  = nucleon-level pp/pn arrays (shared)
        element_cs[Z]  = {"A": float, "name": str, cs_key: np.array, ...}
    """
    cs_tot_pp, cs_el_pp, cs_inel_pp, cs_tot_pn, cs_hN = make_nucleon_cs(splines)
 
    elements = resolve_materials(only)
 
    print(f"Stage 1: computing GG cross sections for {len(elements)} element(s) "
          f"over {len(sqrt_s_arr)} SQRT(S) points ...")
 
    element_cs = {}
    element_cs[0] = {
        "cs_tot_pp":  np.array([cs_tot_pp(s)  for s in sqrt_s_arr]),
        "cs_el_pp":   np.array([cs_el_pp(s)   for s in sqrt_s_arr]),
        "cs_inel_pp": np.array([cs_inel_pp(s) for s in sqrt_s_arr]),
        "cs_tot_pn":  np.array([cs_tot_pn(s)  for s in sqrt_s_arr]),
    }
    if testing:
        # Print nucleon-level spot check
        print(f"\n  Nucleon-level cross sections at selected energies [mb]:")
        print(f"  {'sqrt_s [GeV]':>14}  {'cs_tot_pp':>12}  {'cs_el_pp':>12}  "
            f"{'cs_inel_pp':>12}  {'cs_tot_pn':>12}")
        for s in [10, 20, 50, 100]:
            if s >= sqrt_s_arr.min() and s <= sqrt_s_arr.max():
                print(f"  {s:>14}  {cs_tot_pp(s):>12.4f}  {cs_el_pp(s):>12.4f}  "
                    f"{cs_inel_pp(s):>12.4f}  {cs_tot_pn(s):>12.4f}")
 
    for Z, mat in elements.items():
        name = mat.name or f"Z{Z}"
        A    = float(mat.A)
 
        gg = _gg_isotope_weighted(Z, mat, sqrt_s_arr, cs_hN)
 
        element_cs[Z] = {"A": A, "name": name, **gg}
    if testing:
        # Spot check
        print(f"\n  GG cross sections for {name} (Z={Z}, A={A}) [mb]:")
        print(f"  {'sqrt_s [GeV]':>14}  {'cs_tot_hA':>12}  {'cs_inel_hA':>12}  "
              f"{'cs_el_hA':>12}  {'cs_prod_hA':>12}  {'cs_sd_hA':>12}  {'cs_qel_hA':>12}")
        for s in [10, 20, 50, 100]:
            if s >= sqrt_s_arr.min() and s <= sqrt_s_arr.max():
                r = glauber_element_single(A, Z, s, cs_hN)
                print(f"  {s:>14}  {r['cs_tot_hA']:>12.4f}  {r['cs_inel_hA']:>12.4f}  "
                      f"{r['cs_el_hA']:>12.4f}  {r['cs_prod_hA']:>12.4f}  "
                      f"{r['cs_sd_hA']:>12.4f}  {r['cs_qel_hA']:>12.4f}")
 
    print("\n  Done.\n")
    return element_cs
 
 
# # ---------------------------------------------------------------------------
# # Save to NPZ
# # ---------------------------------------------------------------------------
 
# def save_elements_npz(element_cs, sqrt_s_arr, filename="elements_cs.npz"):
#     """Save all elemental cross section arrays to a single NPZ file."""
#     arrays = {"sqrt_s": sqrt_s_arr}
 
#     for key in NUCLEON_KEYS:
#         arrays[f"nucleon_{key}"] = element_cs[0][key]
 
#     Zs, As = [], []
#     for Z in sorted(k for k in element_cs if k > 0):
#         for key in GG_KEYS:
#             arrays[f"Z{Z}_{key}"] = element_cs[Z][key]
#         Zs.append(Z)
#         As.append(element_cs[Z]["A"])
 
#     arrays["element_Z"] = np.array(Zs)
#     arrays["element_A"] = np.array(As)
 
#     np.savez(filename, **arrays)
#     print(f"Saved {len(Zs)} element(s) -> {filename}")
#     print(f"  Keys: sqrt_s, nucleon_cs_*, Z{{z}}_cs_*  (Z = {Zs[0]} .. {Zs[-1]})")
 
 

# def main():
#     parser = argparse.ArgumentParser(
#         description="Stage 1: compute GG cross sections for all elements."
#     )
#     parser.add_argument("--sqrt-s-min", type=float, default=None,
#                         help="Min sqrt(s) [GeV]  (default: auto from spline data)")
#     parser.add_argument("--sqrt-s-max", type=float, default=None,
#                         help="Max sqrt(s) [GeV]  (default: auto from spline data)")
#     parser.add_argument("--n-points",   type=int,   default=1000,
#                         help="Number of energy grid points")
#     parser.add_argument("--out",        default="elements_cs.npz",
#                         help="Output NPZ file")
#     args = parser.parse_args()

#     # Load splines and determine grid
#     splines, grid_min, grid_max = load_all_splines()

#     sqrt_s_min = args.sqrt_s_min if args.sqrt_s_min is not None else grid_min
#     sqrt_s_max = args.sqrt_s_max if args.sqrt_s_max is not None else grid_max

#     sqrt_s_arr = (
#         np.logspace(math.log10(sqrt_s_min), math.log10(sqrt_s_max), args.n_points)
#     )

#     print(f"Grid: {args.n_points} points, "
#           f"sqrt(s) = {sqrt_s_min:.4g} to {sqrt_s_max:.4g} GeV\n")

#     # Run Stage 1
#     element_cs = compute_all_elements(sqrt_s_arr, splines)
#     save_elements_npz(element_cs, sqrt_s_arr, args.out)


# if __name__ == "__main__":
#     main()