# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2026.                 #
# ######################################### #
import numpy as np
from pathlib import Path
from scipy.interpolate import CubicSpline
from .isotopes import ISOTOPES
from .glauber_gribov import make_nucleon_cs, load_all_splines, SPECIES
N_CS_POINTS = 1000
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
def get_R(A):
    """Nuclear radius [m]."""
    A = np.asarray(A)  # Ensure A is an array
    R = np.where(A > 21,
                 1.1 * A**(1.0/3.0) * 1e-15 * 0.9,
                 1.1 * A**(1.0/3.0) * 1e-15 * 1.05)
    return R

def pi_R2_mb(A):
    """pi*R^2 in mb  (1 fm^2 = 10 mb)."""
    R_fm = get_R(A=A) * 1e15
    return np.pi * R_fm**2 * 10.0

def _compute_glauber_cs(A_sig_tot, A_sig_inel, piR2):
        cs_tot   = 2 * piR2 * np.log(1.0 + A_sig_tot / (2 * piR2))
        cs_inel  = piR2 * np.log(1.0 + A_sig_tot / piR2)
        cs_el    = np.maximum(0.0, cs_tot - cs_inel)

        cs_prod  = piR2 * np.log(1.0 +A_sig_inel / piR2)

        alpha    = A_sig_tot / (2 * piR2 + A_sig_tot)
        cs_sd    = piR2 * (alpha - np.log(1.0 + alpha))

        cs_el_nucleon = np.maximum(0.0, cs_inel - cs_prod)

        return cs_tot, cs_prod, cs_el, cs_el_nucleon, cs_sd

def _glauber_isotope(A, Z, sqrt_s, cs_hN, species):

    sqrt_s = np.asarray(sqrt_s)
    piR2 = pi_R2_mb(A=A)

    # storage only for ONE species now
    A_sig = {
        "tot":  np.zeros_like(sqrt_s),
        "inel": np.zeros_like(sqrt_s),
        "el":   np.zeros_like(sqrt_s),
    }

    # -----------------------------
    # fill nucleon-level inputs
    # -----------------------------
    for i, s in enumerate(sqrt_s):
        values = cs_hN(A, Z, s)

        # unpack ONLY the requested species
        if species == "pp":
            A_sig["tot"][i], A_sig["inel"][i], A_sig["el"][i] = values[0:3]

        elif species == "kmin":
            A_sig["tot"][i], A_sig["inel"][i], A_sig["el"][i] = values[3:6]

        elif species == "kplus":
            A_sig["tot"][i], A_sig["inel"][i], A_sig["el"][i] = values[6:9]

        elif species == "pimin":
            A_sig["tot"][i], A_sig["inel"][i], A_sig["el"][i] = values[9:12]

        elif species == "piplus":
            A_sig["tot"][i], A_sig["inel"][i], A_sig["el"][i] = values[12:15]

        else:
            raise ValueError(f"Unknown species: {species}")

    # -----------------------------
    # small-A limit (no shadowing)
    # -----------------------------
    if A < 4:
        return {
            f"cs_tot_{species}_GG": A_sig["tot"],
            f"cs_el_{species}_GG": A_sig["el"],
            f"cs_prod_{species}_GG": np.zeros_like(sqrt_s),
            f"cs_sd_{species}_GG": np.zeros_like(sqrt_s),
            f"cs_el_nucleon_{species}_GG": np.zeros_like(sqrt_s),
        }

    # -----------------------------
    # Glauber computation
    # -----------------------------
    cs_tot, cs_prod, cs_el, cs_el_nucl, cs_sd = _compute_glauber_cs(
        A_sig["tot"],
        A_sig["inel"],
        piR2
    )

    return {
        f"cs_tot_{species}_GG": cs_tot,
        f"cs_el_{species}_GG": cs_el,
        f"cs_prod_{species}_GG": cs_prod,
        f"cs_sd_{species}_GG": cs_sd,
        f"cs_el_nucleon_{species}_GG": cs_el_nucl,
    }

# def _glauber_isotope(A, Z, sqrt_s, cs_hN):
#     """
#     GG nucleus cross sections [mb] for isotopes.

#     Returns dict with keys: cs_tot_hA, cs_inel_hA, cs_el_hA,
#                             cs_prod_hA, cs_sd_hA, cs_qel_hA.

#     For A < 4, GG shadowing is not applied.
#     """
#     sqrt_s = np.asarray(sqrt_s)  # Ensure sqrt_s is an array
#     piR2 = pi_R2_mb(A=A)
#     A_sig_tot = np.zeros_like(sqrt_s)
#     A_sig_inel = np.zeros_like(sqrt_s)
#     A_sig_el = np.zeros_like(sqrt_s)
#     for i, s in enumerate(sqrt_s):
#         Atot, Ainel, Ael = cs_hN(A, Z, s)
#         A_sig_tot[i] = Atot
#         A_sig_inel[i] = Ainel
#         A_sig_el[i] = Ael

#     if A < 4:
#         return {
#             "cs_tot_hA":  A_sig_tot,
#             # "cs_inel_hA": A_sig_inel,
#             "cs_el_hA":   A_sig_el,
#             "cs_prod_hA": np.zeros_like(A_sig_tot),
#             "cs_sd_hA":   np.zeros_like(A_sig_tot),
#             "cs_el_nucleon":  np.zeros_like(A_sig_tot),
#         }

#     cs_tot_hA  = 2 * piR2 * np.log(1.0 + A_sig_tot / (2 * piR2))
#     cs_inel_hA = piR2 * np.log(1.0 + A_sig_tot / piR2)
#     cs_el_hA   = np.maximum(0.0, cs_tot_hA - cs_inel_hA)
#     cs_prod_hA = piR2 * np.log(1.0 + A_sig_inel / piR2)
#     alpha      = A_sig_tot / (2 * piR2 + A_sig_tot)
#     cs_sd_hA   = piR2 * (alpha - np.log(1.0 + alpha))
#     cs_el_nucleon  = np.maximum(0.0, cs_inel_hA - cs_prod_hA)

#     return {
#         "cs_tot_hA":  cs_tot_hA,
#         # "cs_inel_hA": cs_inel_hA,
#         "cs_el_hA":   cs_el_hA,
#         "cs_prod_hA": cs_prod_hA,
#         "cs_sd_hA":   cs_sd_hA,
#         "cs_el_nucleon":  cs_el_nucleon,
#     }
def build_element_file(
    only_element=None,
    n_points=None,
    sqrt_s_min=None,
    sqrt_s_max=None
):
    splines, grid_min, grid_max = load_all_splines()

    # NOTE: cs_hN is assumed to return ALL species in fixed order
    _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, cs_hN = make_nucleon_cs(splines)

    n_points = n_points or N_CS_POINTS
    # -----------------------------
    # choose elements
    # -----------------------------
    if only_element is not None:
        elements = [(only_element, ISOTOPES[only_element])]
    else:
        elements = ISOTOPES.items()

    # -----------------------------
    # main loop over elements
    # -----------------------------
    for element, iso_data in elements:

        Z = SYMBOL_TO_Z.get(element)
        if Z is None:
            print(f"[WARN] missing Z for {element}")
            continue

        isotopes = [
            iso for iso in iso_data["isotopes"]
            if iso["abundance"] is not None
        ]

        if not isotopes:
            continue

        # output container (species × GG_KEYS)
        out = {sp: {f"cs_{ch}_{sp}_GG" if ch != "el_nucleon" else f"cs_el_nucleon_{sp}_GG": [] 
               for ch in ["tot", "el", "prod", "sd", "el_nucleon"]}
               for sp in SPECIES}
        print(out.keys())
        A_list = []

        for iso in isotopes:
            A = iso["mass_number"]
            A_list.append(A)
            for sp in SPECIES:
                print(f"Processing {element} ({sp})...")
                # energy grid per species
                smin = sqrt_s_min or grid_min[sp]
                smax = sqrt_s_max or grid_max[sp]

                sqrt_s_arr = np.logspace(
                    np.log10(smin),
                    np.log10(smax),
                    n_points
                )

                log_sqrt_s = np.log(sqrt_s_arr)
                gg = _glauber_isotope(A, Z, sqrt_s_arr, cs_hN, species=sp)
                gg = _glauber_isotope(A, Z, sqrt_s_arr, cs_hN, species=sp)
                for key in gg.keys():
                    if sp in key:  # only keep keys belonging to this species
                        print(key)
                        out[sp][key].append({
                            "y_values": gg[key],
                            "knots": log_sqrt_s
                        })
            
        # convert lists → object arrays
        for sp in SPECIES:
            print(out[sp].keys())
            print(out.keys())
            for key in out[sp].keys():
                print(f"Converting {element} ({sp}) {key} to array...")
                out[sp][key] = np.array(out[sp][key], dtype=object)

        out["A"] = np.array(A_list)

        filename = (
            Path(__file__).parent
            / "isotopes"
            / f"{element}_isotopes.npz"
        )

        np.savez(filename, **out)

        print(f"[OK] wrote {filename}")
# def build_element_file():
#     splines, grid_min, grid_max = load_all_splines()
#     n_points   = n_points   or N_CS_POINTS
#     sqrt_s_min = sqrt_s_min or grid_min
#     sqrt_s_max = sqrt_s_max or grid_max
#     sqrt_s_arr = np.logspace(np.log10(sqrt_s_min), np.log10(sqrt_s_max), n_points)
#     log_sqrt_s = np.log(sqrt_s_arr)
#     _,_,_,_,cs_hN = make_nucleon_cs(splines)

#     for element, iso_data in ISOTOPES.items():

#         Z = SYMBOL_TO_Z.get(element)
#         if Z is None:
#             print(f"[WARN] missing Z for {element}")
#             continue

#         isotopes = [
#             iso for iso in iso_data["isotopes"]
#             if iso["abundance"] is not None
#         ]

#         if not isotopes:
#             continue

#         out = {key: [] for key in GG_KEYS}
#         A_list = []
#         abundance_list = []

#         for iso in isotopes:

#             A = iso["mass_number"]

#             gg = _glauber_isotope(A, Z, sqrt_s_arr, cs_hN)

#             A_list.append(A)

#             for key in GG_KEYS:
#                 spline = CubicSpline(log_sqrt_s, gg[key])
#                 out[key].append({"coeffs": spline.c, "knots": log_sqrt_s})

#         out = {k: np.array(v) for k, v in out.items()}
#         out["A"] = np.array(A_list)
#         filename = Path(__file__).parent / "isotopes" / f"{element}_isotopes.npz"
#         np.savez(filename, **out)
#         print(f"[OK] wrote {filename}")

# def isotope_weighted_gg(symbol, Z, sqrt_s_arr, cs_hN):
#     """
#     Compute the isotope-abundance-weighted GG cross sections for an element.
#     Only stable isotopes (abundance is not None) are included.
#     Returns dict {key: np.array}, or None if no stable isotopes found.
#     """
#     isotopes = [iso for iso in ISOTOPES[symbol]["isotopes"]
#                 if iso["abundance"] is not None]
#     if not isotopes:
#         return None

#     # Normalise abundances (they should sum to 1 but do so defensively)
#     total = sum(iso["abundance"] for iso in isotopes)
#     combined = {k: np.zeros(len(sqrt_s_arr)) for k in GG_KEYS}
#     for iso in isotopes:
#         frac = iso["abundance"] / total
#         gg   = _glauber_isotope(iso["atomic_mass"], Z, sqrt_s_arr, cs_hN)
#         for k in GG_KEYS:
#             combined[k] += frac * gg[k]
#     return combined


# def relative_difference(arr_a, arr_b):
#     denom = (np.abs(arr_a) + np.abs(arr_b)) / 2.0
#     with np.errstate(invalid='ignore', divide='ignore'):
#         rel = np.where(denom > 0, np.abs(arr_a - arr_b) / denom, 0.0)
#     return rel

# def check_isotope_approximation(threshold=0.10, n_points=10000, cs_key="cs_inel_hA"):
#     """
#     For each element in the ISOTOPES table, compare the GG cross section
#     computed from the standard atomic weight A against the isotope-abundance-
#     weighted sum of per-isotope GG cross sections.
 
#     Parameters
#     ----------
#     threshold : float
#         Flag elements whose max relative difference exceeds this value.
#         Default 0.10 (10%).
#     n_points : int
#         Number of energy grid points. Default 1000.
#     cs_key : str
#         Which GG cross section to compare. One of GG_KEYS.
#         Default 'cs_inel_hA'.
 
#     Returns
#     -------
#     results : list of dict
#         One dict per element with keys:
#             symbol, Z, A_std, max_rel_diff, mean_rel_diff, flagged
#     flagged : list of dict
#         Subset of results where flagged is True.
#     """
#     splines, grid_min, grid_max = load_all_splines()
#     sqrt_s_arr = np.logspace(math.log10(grid_min), math.log10(grid_max), n_points)
#     _, _, _, _, cs_hN = make_nucleon_cs(splines)
 
#     print(f"\nComparing GG cross sections: standard A vs isotope-weighted sum")
#     print(f"Cross section: {cs_key}")
#     print(f"Flag threshold: {threshold*100:.1f}%")
 
#     print(f"{'Symbol':<6}  {'Z':>4}  {'A_std':>10}  "
#           f"{'Max rel diff':>14}  {'Mean rel diff':>14}  {'Flag'}")
#     print("-" * 65)
 
#     results = []
#     flagged = []
 
#     for symbol, data in ISOTOPES.items():
#         Z = SYMBOL_TO_Z.get(symbol)
#         if Z is None:
#             continue
 
#         A_std = data["standard_atomic_weight"]
 
#         # Method 1: standard atomic weight
#         gg_std = compute_gg_arrays(A_std, Z, sqrt_s_arr, cs_hN)
 
#         # Method 2: isotope-weighted sum
#         gg_iso = isotope_weighted_gg(symbol, Z, sqrt_s_arr, cs_hN)
#         if gg_iso is None:
#             print(f"{symbol:<6}  {Z:>4}  {A_std:>10.4f}  "
#                   f"{'no stable isotopes':>30}")
#             continue
 
#         rel      = relative_difference(gg_std[cs_key], gg_iso[cs_key])
#         max_rel  = float(rel.max())
#         mean_rel = float(rel.mean())
#         is_flagged = max_rel > threshold
 
#         print(f"{symbol:<6}  {Z:>4}  {A_std:>10.4f}  "
#               f"{max_rel*100:>13.4f}%  {mean_rel*100:>13.4f}%  "
#               f"{'*** FLAG ***' if is_flagged else ''}")
 
#         row = {
#             "symbol":        symbol,
#             "Z":             Z,
#             "A_std":         A_std,
#             "max_rel_diff":  max_rel,
#             "mean_rel_diff": mean_rel,
#             "flagged":       is_flagged,
#         }
#         results.append(row)
#         if is_flagged:
#             flagged.append(row)
 
#     print()
#     if flagged:
#         print(f"Flagged elements (max rel diff > {threshold*100:.1f}%):")
#         for row in flagged:
#             print(f"  {row['symbol']} (Z={row['Z']}, A={row['A_std']:.4f}): "
#                   f"{row['max_rel_diff']*100:.4f}%")
#     else:
#         print(f"No elements flagged above {threshold*100:.1f}% threshold.")
 
#     return results, flagged