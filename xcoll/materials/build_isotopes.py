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

        # output container (species x GG_KEYS)
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
            
        # convert lists to object arrays
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