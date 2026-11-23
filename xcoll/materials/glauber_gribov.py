# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2026.                 #
# ######################################### #
import math
import numpy as np
from scipy.interpolate import UnivariateSpline
from pathlib import Path

# ===========================================================================
# Physical constants
# ===========================================================================
MP = 0.938270   # proton mass [GeV]

# GG_KEYS      = ["cs_tot_pp_GG", "cs_inel_pp_GG", "cs_el_pp_GG",
#                 "cs_prod_pp_GG", "cs_sd_pp_GG",  "cs_el_nucleon_pp_GG",
#                 "cs_tot_kmin_GG", "cs_inel_kmin_GG", "cs_el_kmin_GG", 
#                 "cs_prod_kmin_GG", "cs_sd_kmin_GG",  "cs_el_nucleon_kmin_GG",
#                 "cs_tot_kplus_GG", "cs_inel_kplus_GG", "cs_el_kplus_GG",
#                 "cs_prod_kplus_GG", "cs_sd_kplus_GG",  "cs_el_nucleon_kplus",
#                 "cs_tot_pimin_GG", "cs_inel_pimin_GG", "cs_el_pimin_GG",
#                 "cs_prod_pimin_GG", "cs_sd_pimin_GG",  "cs_el_nucleon_pimin_GG",
#                 "cs_tot_piplus_GG", "cs_inel_piplus_GG", "cs_el_piplus_GG",
#                 "cs_prod_piplus_GG", "cs_sd_piplus_GG",  "cs_el_nucleon_piplus_GG"]
SPECIES = ["pp", "kmin", "kplus", "pimin", "piplus"]
# NUCLEON_KEYS = ["cs_tot_pp", "cs_el_pp", "cs_inel_pp", "cs_tot_pn"] # Maybe call spline keys


# ===========================================================================
# PDG data saving and spline fitting
# ===========================================================================

def plab_to_sqrts(plab):
    """Convert lab-frame momentum [GeV/c] to sqrt(s) [GeV]."""
    plab = np.asarray(plab, dtype=float)
    E    = np.sqrt(plab**2 + MP**2)
    return np.sqrt(2*MP**2 + 2*MP*E)


def load_pdg_table(filepath):
    """
    Read a PDG-format reaction data file.
    Column 1 (0-indexed) = plab [GeV/c],  column 4 = sigma [mb].
    Lines with fewer than 9 columns or non-numeric values are skipped.
    """
    plab, sigma = [], []
    with open(filepath) as f:
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
        raise ValueError(f"No numeric data found in {filepath!r}")
    return np.array(plab), np.array(sigma)


def build_spline(filename, smoothing):
    """
    Load a PDG file, convert plab -> sqrt(s), sort, fit a smoothing spline
    in log(sqrt_s) space.  Returns (spline, sqrts_array).
    """
    filepath     = Path(__file__).parent / "data" / filename
    plab, sigma  = load_pdg_table(filepath)
    sqrts        = plab_to_sqrts(plab)
    idx          = np.argsort(sqrts)
    sqrts, sigma = sqrts[idx], sigma[idx]
    spline       = UnivariateSpline(np.log(sqrts), sigma, s=smoothing)
    return spline, sqrts, sigma   # also return sigma for saving


def save_splines_npz(out_path,
                     fname_pp_tot   ="ppcrosstot.dat",
                     fname_pp_el    ="ppcross.dat",
                     fname_pn_tot   ="pncrosstot.dat",
                     fname_kmin_tot ="k_min_p_tot.dat",
                     fname_kmin_el  ="k_min_p_el.dat",
                     fname_kplu_tot ="k_plu_p_tot.dat",
                     fname_kplu_el  ="k_plu_p_el.dat",
                     fname_pimin_tot="pi_min_p_tot.dat",
                     fname_pimin_el ="pi_min_p_el.dat",
                     fname_piplu_tot="pi_plu_p_tot.dat",
                     fname_piplu_el ="pi_plu_p_el.dat"):
    """
    Build splines from raw PDG dat files and save the raw data points
    (in log(sqrt_s) space) to an NPZ file.  Run this once to generate
    splines_data.npz.

    Parameters
    ----------
    out_path : str or Path
        Where to write the NPZ, e.g. "xcoll/materials/data/splines_data.npz"
    fname_tot, fname_el, fname_pn : str
        Paths to the three PDG dat files.
    """
    _, sqrts_tot, sigma_tot             = build_spline(fname_pp_tot, smoothing=2000)
    _, sqrts_el,  sigma_el              = build_spline(fname_pp_el,  smoothing=150)
    _, sqrts_pn,  sigma_pn              = build_spline(fname_pn_tot,  smoothing=140)
    _, sqrts_kmin_tot, sigma_kmin_tot   = build_spline(fname_kmin_tot, smoothing=1200)
    _, sqrts_kmin_el,  sigma_kmin_el    = build_spline(fname_kmin_el,  smoothing=900)
    _, sqrts_kplu_tot, sigma_kplu_tot   = build_spline(fname_kplu_tot, smoothing=500)
    _, sqrts_kplu_el,  sigma_kplu_el    = build_spline(fname_kplu_el,  smoothing=315)
    _, sqrts_pimin_tot, sigma_pimin_tot = build_spline(fname_pimin_tot, smoothing=1300)
    _, sqrts_pimin_el,  sigma_pimin_el  = build_spline(fname_pimin_el,  smoothing=1300)
    _, sqrts_piplu_tot, sigma_piplu_tot = build_spline(fname_piplu_tot, smoothing=1300)
    _, sqrts_piplu_el,  sigma_piplu_el  = build_spline(fname_piplu_el,  smoothing=500)

    np.savez(out_path,
             x_tot      = np.log(sqrts_tot),       y_tot=sigma_tot,
             x_el       = np.log(sqrts_el),        y_el=sigma_el,
             x_pn       = np.log(sqrts_pn),        y_pn=sigma_pn,
             x_kmin_tot = np.log(sqrts_kmin_tot),  y_kmin_tot=sigma_kmin_tot,
             x_kmin_el  = np.log(sqrts_kmin_el),   y_kmin_el=sigma_kmin_el,
             x_kplu_tot = np.log(sqrts_kplu_tot),  y_kplu_tot=sigma_kplu_tot,
             x_kplu_el  = np.log(sqrts_kplu_el),   y_kplu_el=sigma_kplu_el,
             x_pimin_tot= np.log(sqrts_pimin_tot), y_pimin_tot=sigma_pimin_tot,
             x_pimin_el = np.log(sqrts_pimin_el),  y_pimin_el=sigma_pimin_el,
             x_piplu_tot= np.log(sqrts_piplu_tot), y_piplu_tot=sigma_piplu_tot,
             x_piplu_el = np.log(sqrts_piplu_el),  y_piplu_el=sigma_piplu_el,
    )
    print(f"Saved splines data -> {out_path}")

# ===========================================================================
# Loading splines
# ==========================================================================
def load_all_splines():
    """
    Load the PDG splines from splines_data.npz bundled with the package.
    Returns splines_dict with keys 'pp tot', 'pp el', 'pn tot'.
    """
    data = np.load(Path(__file__).parent / "data" / "splines_data.npz")

    splines = {
        "pp tot": UnivariateSpline(data["x_tot"], data["y_tot"], s=900),
        "pp el" : UnivariateSpline(data["x_el"],  data["y_el"],  s=500),
        "pn tot": UnivariateSpline(data["x_pn"],  data["y_pn"],  s=200),
        "kmin tot": UnivariateSpline(data["x_k_min_tot"], data["y_k_min_tot"], s=1200),
        "kmin el":  UnivariateSpline(data["x_k_min_el"],  data["y_k_min_el"],  s=900),
        "kplus tot": UnivariateSpline(data["x_k_plu_tot"], data["y_k_plu_tot"], s=500),
        "kplus el":  UnivariateSpline(data["x_k_plu_el"],  data["y_k_plu_el"],  s=315),
        "pimin tot": UnivariateSpline(data["x_pi_min_tot"], data["y_pi_min_tot"], s=1300),
        "pimin el":  UnivariateSpline(data["x_pi_min_el"],  data["y_pi_min_el"],  s=1300),
        "piplus tot": UnivariateSpline(data["x_pi_plu_tot"], data["y_pi_plu_tot"], s=1300),
        "piplus el":  UnivariateSpline(data["x_pi_plu_el"],  data["y_pi_plu_el"],  s=500),
    }
    grid_min_pp = max(np.exp(data[x].min()) for x in ["x_tot", "x_el"])
    grid_max_pp = min(np.exp(data[x].max()) for x in ["x_tot", "x_el"])
    grid_min_kmin = max(np.exp(data[x].min()) for x in ["x_k_min_tot", "x_k_min_el"])
    grid_max_kmin = min(np.exp(data[x].max()) for x in ["x_k_min_tot", "x_k_min_el"])
    grid_min_kplu = max(np.exp(data[x].min()) for x in ["x_k_plu_tot", "x_k_plu_el"])
    grid_max_kplu = min(np.exp(data[x].max()) for x in ["x_k_plu_tot", "x_k_plu_el"])
    grid_min_pimin = max(np.exp(data[x].min()) for x in ["x_pi_min_tot", "x_pi_min_el"])
    grid_max_pimin = min(np.exp(data[x].max()) for x in ["x_pi_min_tot", "x_pi_min_el"])
    grid_min_piplu = max(np.exp(data[x].min()) for x in ["x_pi_plu_tot", "x_pi_plu_el"])
    grid_max_piplu = min(np.exp(data[x].max()) for x in ["x_pi_plu_tot", "x_pi_plu_el"])
    grid_min = {
        "pp": grid_min_pp,
        "kmin": grid_min_kmin,
        "kplus": grid_min_kplu,
        "pimin": grid_min_pimin,
        "piplus": grid_min_piplu,
    }
    grid_max = {
        "pp": grid_max_pp,
        "kmin": grid_max_kmin,
        "kplus": grid_max_kplu,
        "pimin": grid_max_pimin,
        "piplus": grid_max_piplu,
    }

    return splines, grid_min, grid_max

# ===========================================================================
# Nucleon-level cross sections
# ===========================================================================

def make_nucleon_cs(splines):
    """
    Return a bundle of functions that evaluate nucleon-level cross sections
    at a given sqrt(s), bound to the provided spline dict.
    Returns (cs_tot_pp, cs_el_pp, cs_inel_pp, cs_tot_pn, cs_hN).
    """
    def cs_tot_pp(sqrt_s):
        return (splines["pp tot"](np.log(sqrt_s)))

    def cs_el_pp(sqrt_s):
        return (splines["pp el"](np.log(sqrt_s)))

    def cs_inel_pp(sqrt_s):
        return max(0.0, cs_tot_pp(sqrt_s) - cs_el_pp(sqrt_s))

    def cs_tot_pn(sqrt_s):
        pn = splines["pn tot"](np.log(sqrt_s))
        pp = splines["pp tot"](np.log(sqrt_s))
        lo, hi = 25.0, 35.0
        if sqrt_s <= lo:
            return pn
        if sqrt_s >= hi:
            return pp
        t = (sqrt_s - lo) / (hi - lo)
        return (1 - t) * pn + t * pp

    def cs_inel_pn(sqrt_s):
        return cs_inel_pp(sqrt_s)   # no pn inelastic data

    def cs_el_pn(sqrt_s):
        return cs_el_pp(sqrt_s)     # no pn elastic data

    def cs_tot_kmin_p(sqrt_s):
        return splines["kmin tot"](np.log(sqrt_s))
    
    def cs_el_kmin_p(sqrt_s):
        return splines["kmin el"](np.log(sqrt_s))
    
    def cs_inel_kmin_p(sqrt_s):
        return max(0.0, cs_tot_kmin_p(sqrt_s) - cs_el_kmin_p(sqrt_s))
    
    def cs_tot_kplu_p(sqrt_s):
        return splines["kplus tot"](np.log(sqrt_s))
    
    def cs_el_kplu_p(sqrt_s):
        return splines["kplus el"](np.log(sqrt_s))
    
    def cs_inel_kplu_p(sqrt_s):
        return max(0.0, cs_tot_kplu_p(sqrt_s) - cs_el_kplu_p(sqrt_s))

    def cs_tot_pimin_p(sqrt_s):
        return splines["pimin tot"](np.log(sqrt_s))
    
    def cs_el_pimin_p(sqrt_s):
        return splines["pimin el"](np.log(sqrt_s))
    
    def cs_inel_pimin_p(sqrt_s):
        return max(0.0, cs_tot_pimin_p(sqrt_s) - cs_el_pimin_p(sqrt_s))

    def cs_tot_piplu_p(sqrt_s):
        return splines["piplus tot"](np.log(sqrt_s))
    
    def cs_el_piplu_p(sqrt_s):
        return splines["piplus el"](np.log(sqrt_s))

    def cs_inel_piplu_p(sqrt_s):
        return max(0.0, cs_tot_piplu_p(sqrt_s) - cs_el_piplu_p(sqrt_s))

    def cs_hN(A, Z, sqrt_s):
        """
        Hadron-nucleon cross sections: A*sigma = Z*sigma_pp + (A-Z)*sigma_pn.
        Returns (cs_tot_hN, cs_inel_hN, cs_el_hN).
        """
        tot_pp  = Z * cs_tot_pp(sqrt_s)  + (A - Z) * cs_tot_pn(sqrt_s)
        inel_pp = Z * cs_inel_pp(sqrt_s) + (A - Z) * cs_inel_pn(sqrt_s)
        el_pp   = Z * cs_el_pp(sqrt_s)   + (A - Z) * cs_el_pn(sqrt_s)

        tot_kmin_p  = Z * cs_tot_kmin_p(sqrt_s)  + (A - Z) * cs_tot_kmin_p(sqrt_s)
        inel_kmin_p = Z * cs_inel_kmin_p(sqrt_s) + (A - Z) * cs_inel_kmin_p(sqrt_s)
        el_kmin_p   = Z * cs_el_kmin_p(sqrt_s)   + (A - Z) * cs_el_kmin_p(sqrt_s)
        
        tot_kplu_p = Z * cs_tot_kplu_p(sqrt_s)  + (A - Z) * cs_tot_kplu_p(sqrt_s)
        inel_kplu_p = Z * cs_inel_kplu_p(sqrt_s) + (A - Z) * cs_inel_kplu_p(sqrt_s)
        el_kplu_p   = Z * cs_el_kplu_p(sqrt_s)   + (A - Z) * cs_el_kplu_p(sqrt_s)

        tot_pimin_p = Z * cs_tot_pimin_p(sqrt_s)  + (A - Z) * cs_tot_pimin_p(sqrt_s)
        inel_pimin_p = Z * cs_inel_pimin_p(sqrt_s) + (A - Z) * cs_inel_pimin_p(sqrt_s)
        el_pimin_p   = Z * cs_el_pimin_p(sqrt_s)   + (A - Z) * cs_el_pimin_p(sqrt_s)

        tot_piplu_p = Z * cs_tot_piplu_p(sqrt_s)  + (A - Z) * cs_tot_piplu_p(sqrt_s)
        inel_piplu_p = Z * cs_inel_piplu_p(sqrt_s) + (A - Z) * cs_inel_piplu_p(sqrt_s)
        el_piplu_p   = Z * cs_el_piplu_p(sqrt_s)   + (A - Z) * cs_el_piplu_p(sqrt_s)

        return tot_pp, inel_pp, el_pp, tot_kmin_p, inel_kmin_p, el_kmin_p, \
               tot_kplu_p, inel_kplu_p, el_kplu_p, tot_pimin_p, inel_pimin_p, \
               el_pimin_p, tot_piplu_p, inel_piplu_p, el_piplu_p

    return  cs_tot_pp, cs_el_pp, cs_inel_pp, cs_tot_pn, \
            cs_tot_kmin_p, cs_el_kmin_p, cs_inel_kmin_p, \
            cs_tot_kplu_p, cs_el_kplu_p, cs_inel_kplu_p, \
            cs_tot_pimin_p, cs_el_pimin_p, cs_inel_pimin_p, \
            cs_tot_piplu_p, cs_el_piplu_p, cs_inel_piplu_p, \
            cs_hN
 