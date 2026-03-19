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
    Returns a dict  {Z: (name, A)}  sorted by Z.
    """
    elements = {}
    for name in dir(xc.materials):
        mat = getattr(xc.materials, name)
        try:
            if (hasattr(mat, 'Z') and mat.Z is not None
                    and hasattr(mat, 'A') and mat.A is not None
                    and not (hasattr(mat, 'components') and mat.components)):
                Z = int(mat.Z)
                A = float(mat.A)
                if Z not in elements:
                    elements[Z] = (name, A)
        except Exception:
            continue
    return dict(sorted(elements.items()))


def compute_elements(sqrt_s_arr, splines, only=None):
    """
    Compute GG cross sections for all elements over the energy grid.

    Parameters
    ----------
    only : list of (Z, A) tuples or None
        If provided, only compute for these specific elements.
        e.g. only=[(6, 12.011), (74, 183.84)]
        If None, compute for all elements found in xc.materials.
    """
    cs_tot_pp, cs_el_pp, cs_inel_pp, cs_tot_pn, cs_hN = \
        make_nucleon_cs(splines)

    if only is not None:
        elements = {Z: (f"Z{Z}", A) for Z, A in only}
    else:
        elements = find_all_elements()

    print(f"Stage 1: computing GG cross sections for {len(elements)} elements "
          f"over {len(sqrt_s_arr)} energy points ...")

    element_cs = {}
    element_cs[0] = {
        "cs_tot_pp":  np.array([cs_tot_pp(s)  for s in sqrt_s_arr]),
        "cs_el_pp":   np.array([cs_el_pp(s)   for s in sqrt_s_arr]),
        "cs_inel_pp": np.array([cs_inel_pp(s) for s in sqrt_s_arr]),
        "cs_tot_pn":  np.array([cs_tot_pn(s)  for s in sqrt_s_arr]),
    }

    # Print nucleon-level spot check
    print("\nNucleon-level cross sections at selected energies [mb]:")
    print(f"  {'sqrt_s [GeV]':>14}  {'cs_tot_pp':>12}  {'cs_el_pp':>12}  {'cs_inel_pp':>12}  {'cs_tot_pn':>12}")
    for s in [10, 20, 50, 100, 115]:
        print(f"  {s:>14}  {cs_tot_pp(s):>12.4f}  {cs_el_pp(s):>12.4f}  "
            f"{cs_inel_pp(s):>12.4f}  {cs_tot_pn(s):>12.4f}")

    for Z, (name, A) in elements.items():
        rows = [glauber_element_single(A, Z, s, cs_hN) for s in sqrt_s_arr]
        element_cs[Z] = {"A": A, "name": name}
        for key in GG_KEYS:
            element_cs[Z][key] = np.array([r[key] for r in rows])

        # Print GG spot check for this element
        print(f"\nGG cross sections for {name} (Z={Z}, A={A}) at selected energies [mb]:")
        print(f"  {'sqrt_s [GeV]':>14}  {'cs_tot_hA':>12}  {'cs_inel_hA':>12}  "
            f"{'cs_el_hA':>12}  {'cs_prod_hA':>12}  {'cs_sd_hA':>12}  {'cs_qel_hA':>12}")
        for s in [10, 20, 50, 100, 115]:
            r = glauber_element_single(A, Z, s, cs_hN)
            print(f"  {s:>14}  {r['cs_tot_hA']:>12.10f}  {r['cs_inel_hA']:>12.10f}  "
                f"{r['cs_el_hA']:>12.10f}  {r['cs_prod_hA']:>12.10f}  "
                f"{r['cs_sd_hA']:>12.10f}  {r['cs_qel_hA']:>12.10f}")
    print("  Done.\n")
    return element_cs


def save_elements_npz(element_cs, sqrt_s_arr, filename="elements_cs.npz"):
    """
    Save all elemental cross section arrays to a single NPZ file.

    Keys follow the pattern  "Z{z}_cs_inel_hA"  for GG arrays and
    "nucleon_cs_tot_pp" for nucleon-level arrays.  Metadata arrays
    "element_Z" and "element_A" list which elements are present.
    """
    arrays = {"sqrt_s": sqrt_s_arr}

    # Nucleon-level arrays
    for key in NUCLEON_KEYS:
        arrays[f"nucleon_{key}"] = element_cs[0][key]

    # Per-element GG arrays + metadata
    Zs, As = [], []
    for Z in sorted(k for k in element_cs if k > 0):
        for key in GG_KEYS:
            arrays[f"Z{Z}_{key}"] = element_cs[Z][key]
        Zs.append(Z)
        As.append(element_cs[Z]["A"])

    arrays["element_Z"] = np.array(Zs)
    arrays["element_A"] = np.array(As)

    np.savez(filename, **arrays)
    print(f"Saved {len(Zs)} elements -> {filename}")
    print(f"  Keys: sqrt_s, nucleon_cs_*, Z{{z}}_cs_*  "
          f"(Z = {Zs[0]} .. {Zs[-1]})")


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