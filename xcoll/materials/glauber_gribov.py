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

GG_KEYS      = ["cs_tot_hA", "cs_inel_hA", "cs_el_hA",
                "cs_prod_hA", "cs_sd_hA",  "cs_el_nucleon"]
NUCLEON_KEYS = ["cs_tot_pp", "cs_el_pp", "cs_inel_pp", "cs_tot_pn"] # Maybe call spline keys


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
                     fname_tot="ppcrosstot.dat",
                     fname_el="ppcross.dat",
                     fname_pn="pncrosstot.dat"):
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
    _, sqrts_tot, sigma_tot = build_spline(fname_tot, smoothing=2000)
    _, sqrts_el,  sigma_el  = build_spline(fname_el,  smoothing=150)
    _, sqrts_pn,  sigma_pn  = build_spline(fname_pn,  smoothing=140)

    np.savez(out_path,
             x_tot=np.log(sqrts_tot), y_tot=sigma_tot,
             x_el=np.log(sqrts_el),   y_el=sigma_el,
             x_pn=np.log(sqrts_pn),   y_pn=sigma_pn)
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
    }
    grid_min = max(np.exp(data[x].min()) for x in ["x_tot", "x_el"])
    grid_max = min(np.exp(data[x].max()) for x in ["x_tot", "x_el"])

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

    def cs_hN(A, Z, sqrt_s):
        """
        Hadron-nucleon cross sections: A*sigma = Z*sigma_pp + (A-Z)*sigma_pn.
        Returns (cs_tot_hN, cs_inel_hN, cs_el_hN).
        """
        tot  = Z * cs_tot_pp(sqrt_s)  + (A - Z) * cs_tot_pn(sqrt_s)
        inel = Z * cs_inel_pp(sqrt_s) + (A - Z) * cs_inel_pn(sqrt_s)
        el   = Z * cs_el_pp(sqrt_s)   + (A - Z) * cs_el_pn(sqrt_s)
        return tot, inel, el

    return cs_tot_pp, cs_el_pp, cs_inel_pp, cs_tot_pn, cs_hN
 