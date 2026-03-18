# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2026.                 #
# ######################################### #

"""
Cross section sampler for hadron-nucleus interactions.
=======================================================

WHAT THIS SCRIPT DOES — STEP BY STEP
--------------------------------------

Step 1: Load PDG data and build splines
    Four PDG data files are read (pp total, pp elastic, pn total, pd total).
    Each file contains lab-frame momenta (plab) and cross sections (sigma).
    plab is converted to sqrt(s) via kinematics, then a smoothing spline is
    fitted in log(sqrt_s) space.  The data range of each spline is printed
    so you can see exactly what energy coverage you have.

Step 2: Determine the sampling grid
    The grid runs from sqrt_s_min to sqrt_s_max on a log scale with n_points
    steps.  By default, sqrt_s_min and sqrt_s_max are set automatically to
    the overlap of all spline data ranges — i.e. the highest lower bound and
    the lowest upper bound across all four splines.  This ensures we never
    sample outside any spline's data.  You can override these with CLI flags.

Step 3: Nucleon-level cross sections
    At each sqrt_s grid point, the pp and pn cross sections are evaluated
    from the splines:
      - cs_tot_pp, cs_el_pp  : directly from splines
      - cs_inel_pp            : tot_pp - el_pp
      - cs_tot_pn             : pn spline below 30 GeV, falls back to pp above
      - cs_inel_pn, cs_el_pn  : assumed equal to pp (no pn data, same as C code)

Step 4: Glauber-Gribov nucleus cross sections per element
    For each element (Z, A), the nucleon cross sections are first combined as:
      sigma_hN = Z * sigma_pp + (A-Z) * sigma_pn
    Then Glauber-Gribov shadowing gives the nucleus cross sections:
      cs_tot_hA  = 2*pi*R^2 * ln(1 + sigma_tot_hN  / (2*pi*R^2))
      cs_inel_hA =   pi*R^2 * ln(1 + cs_tot_hA     /    pi*R^2 )
      cs_el_hA   = cs_tot_hA - cs_inel_hA
      cs_prod_hA =   pi*R^2 * ln(1 + sigma_inel_hN /    pi*R^2 )
      cs_sd_hA   =   pi*R^2 * (alpha - ln(1 + alpha)),
                     alpha = cs_tot_hA / (2*pi*R^2 + cs_tot_hA)
      cs_qel_hA  = cs_inel_hA - cs_prod_hA
    where R is the nuclear radius from the standard parametrisation.
    For A < 4 (very light nuclei) GG shadowing is not applied.

Step 5: Material averaging
    Pure elements are used directly.  Compounds and mixtures are first
    flattened recursively to a dict of {(Z, A): molar_fraction}, then the
    GG cross sections for each constituent element are averaged with those
    molar fractions:
      sigma_material = sum_j  f_j * sigma_GG_j(sqrt_s)

    The three composition types in xcoll are handled as follows:
      Compound     (n_atoms=[...])         f_j = n_j / sum(n)
      Molar mix    (molar_fractions=[...]) f_j = molar_fraction_j  (normalised)
      Mass mix     (mass_fractions=[...])  f_j = (w_j/A_j) / sum(w_i/A_i)
                                           i.e. mass fractions converted to
                                           molar fractions via effective A.
    Nesting is supported: a mixture whose components are themselves compounds
    (e.g. BoronNitride5000 containing BoronNitride) is fully flattened.

Step 6: Coulomb cross section
    The Coulomb cross section is computed once per material (it does not
    depend on sqrt_s, only on the beam momentum and geometry):
      cs_coulomb = -C * (R^2 * B * E1(R^2*B*t_cut) - exp(-R^2*B*t_cut)/t_cut)
    where C = 4*pi*Z^2*alpha^2*hbar_c^2, B is the hadron-nucleus slope,
    R = 2*hbar_c*sqrt(B), t_cut = (pc * 2.325 * theta_init)^2,
    and E1 is the exponential integral.

Step 7: Save outputs
    - cross_sections.csv / .npz : one row per (material, sqrt_s) with all
      cross sections in mb.
    - cross_sections_splines.h  : C header with one flat double[] per
      (material, cs_type).  Since the grid is fixed log-spaced, the C code
      computes the array index in O(1):
        i = (int)((log(sqrt_s) - LOG_SQRT_S_MIN) / LOG_STEP)
      Grid parameters are written as #defines.

Usage
-----
  python xcoll_cross_sections.py
  python xcoll_cross_sections.py --materials Carbon CopperDiamond Air
  python xcoll_cross_sections.py --sqrt-s-min 5 --sqrt-s-max 1000 --n-points 500
"""

import argparse
import csv
import math
import textwrap
from collections import defaultdict

import numpy as np
import xcoll as xc
from scipy.interpolate import UnivariateSpline
from scipy.special import exp1


# ===========================================================================
# Step 0: Physical constants
# ===========================================================================
XC_AVOGADRO = 6.02214076e23     # mol^-1
HBAR_C      = math.sqrt(0.389)  # sqrt(mb)*GeV ~ 0.6239 GeV*fm
ALPHA_EM    = 1.0 / 137.0
MP          = 0.938270           # proton mass [GeV]


# ===========================================================================
# Step 1: Load PDG data files and build smoothing splines
# ===========================================================================

def plab_to_sqrts(plab):
    """Convert lab-frame momentum [GeV/c] to centre-of-mass energy sqrt(s) [GeV]."""
    plab = np.asarray(plab, dtype=float)
    E    = np.sqrt(plab**2 + MP**2)        # total energy of beam proton
    s    = 2*MP**2 + 2*MP*E                # Mandelstam s for p+p_target
    return np.sqrt(s)


def load_pdg_table(filename):
    """
    Read a PDG-format reaction data file.
    Expected columns (space-separated):  idx  plab  ...  sigma  ...
    Column 1 (0-indexed) = plab [GeV/c]
    Column 4 (0-indexed) = sigma [mb]
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


def build_spline(filename, smoothing, label=""):
    """
    Load a PDG file, convert plab -> sqrt(s), sort, and fit a smoothing
    spline in log(sqrt_s) space.  Returns (spline, sqrts_array, sigma_array)
    where sqrts_array and sigma_array are the raw sorted data points.
    """
    plab, sigma = load_pdg_table(filename)
    sqrts       = plab_to_sqrts(plab)
    idx         = np.argsort(sqrts)
    sqrts, sigma = sqrts[idx], sigma[idx]
    spline      = UnivariateSpline(np.log(sqrts), sigma, s=smoothing)
    return spline, sqrts, sigma


# Load all four splines and print their data ranges
print("Loading PDG spline tables ...")
print(f"  {'Dataset':<12}  {'Points':>6}  {'sqrt(s) min':>12}  {'sqrt(s) max':>12}  File")
print(f"  {'-'*12}  {'-'*6}  {'-'*12}  {'-'*12}  {'-'*20}")

_spline_tot_pp,  _sqrts_tot_pp,  _ = build_spline("ppcrosstot.dat", s=2000, label="pp tot")
_spline_el_pp,   _sqrts_el_pp,   _ = build_spline("ppcross.dat",    s=150,  label="pp el")
_spline_tot_pn,  _sqrts_tot_pn,  _ = build_spline("pncrosstot.dat", s=140,  label="pn tot")
_spline_tot_pd,  _sqrts_tot_pd,  _ = build_spline("pdcrosstot.dat", s=200,  label="pd tot")

_all_spline_data = [
    ("pp tot",  _sqrts_tot_pp,  "ppcrosstot.dat"),
    ("pp el",   _sqrts_el_pp,   "ppcross.dat"),
    ("pn tot",  _sqrts_tot_pn,  "pncrosstot.dat"),
    ("pd tot",  _sqrts_tot_pd,  "pdcrosstot.dat"),
]
for label, sqrts, fname in _all_spline_data:
    print(f"  {label:<12}  {len(sqrts):>6}  {sqrts.min():>12.4g}  {sqrts.max():>12.4g}  {fname}")

# Auto grid range: intersection of all spline data ranges
_GRID_SQRT_S_MIN = max(s.min() for _, s, _ in _all_spline_data)
_GRID_SQRT_S_MAX = min(s.max() for _, s, _ in _all_spline_data)
print(f"\nAuto grid range: sqrt(s) = {_GRID_SQRT_S_MIN:.4g} to {_GRID_SQRT_S_MAX:.4g} GeV")
print("(intersection of all spline data ranges — override with --sqrt-s-min / --sqrt-s-max)\n")


# ===========================================================================
# Step 2 (helpers): Evaluate splines at a given sqrt(s)
# ===========================================================================

def cs_tot_pp(sqrt_s):
    return float(_spline_tot_pp(np.log(sqrt_s)))

def cs_el_pp(sqrt_s):
    return float(_spline_el_pp(np.log(sqrt_s)))

def cs_inel_pp(sqrt_s):
    return max(0.0, cs_tot_pp(sqrt_s) - cs_el_pp(sqrt_s))

def cs_tot_pn(sqrt_s):
    """pn spline below 30 GeV; fall back to pp above (no world data)."""
    if sqrt_s < 30.0:
        return float(_spline_tot_pn(np.log(sqrt_s)))
    return cs_tot_pp(sqrt_s)

def cs_inel_pn(sqrt_s):
    return cs_inel_pp(sqrt_s)   # no pn inelastic data — same as C code

def cs_el_pn(sqrt_s):
    return cs_el_pp(sqrt_s)     # no pn elastic data  — same as C code


# ===========================================================================
# Step 3 (helpers): Hadron-nucleon combinations  Z*sigma_pp + (A-Z)*sigma_pn
# ===========================================================================

def cs_tot_hN(A, Z, sqrt_s):
    return Z * cs_tot_pp(sqrt_s) + (A - Z) * cs_tot_pn(sqrt_s)

def cs_inel_hN(A, Z, sqrt_s):
    return Z * cs_inel_pp(sqrt_s) + (A - Z) * cs_inel_pn(sqrt_s)

def cs_el_hN(A, Z, sqrt_s):
    return Z * cs_el_pp(sqrt_s) + (A - Z) * cs_el_pn(sqrt_s)


# ===========================================================================
# Step 4: Glauber-Gribov nucleus cross sections for a single element
# ===========================================================================

def get_R(A):
    """Nuclear radius [m], standard parametrisation."""
    if A > 21:
        return 1.1 * A**(1.0/3.0) * 1e-15 * 0.9
    return 1.1 * A**(1.0/3.0) * 1e-15 * 1.05

def _pi_R2_mb(A):
    """pi*R^2 in mb  (unit conversion: 1 fm^2 = 10 mb)."""
    R_fm = get_R(A) * 1e15
    return math.pi * R_fm**2 * 10.0

def glauber_element(A, Z, sqrt_s):
    """
    Compute all GG nucleus cross sections [mb] for one element at one energy.
    Returns a dict with keys: cs_tot_hA, cs_inel_hA, cs_el_hA,
                               cs_prod_hA, cs_sd_hA, cs_qel_hA.
    For A < 4 the GG approximation is not applied (no shadowing).
    """
    piR2     = _pi_R2_mb(A)
    sig_tot  = cs_tot_hN(A, Z, sqrt_s)
    sig_inel = cs_inel_hN(A, Z, sqrt_s)

    if A < 4:
        return {
            "cs_tot_hA":  sig_tot,
            "cs_inel_hA": sig_inel,
            "cs_el_hA":   cs_el_hN(A, Z, sqrt_s),
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
# Step 6 (helper): Coulomb cross section  [mb]
# ===========================================================================

def _get_slope_hadron_nucleus(A):
    """Forward elastic B-slope [GeV^-2] for hadron-nucleus, from GEANT4."""
    return 14.5 * A**(2.0/3.0) if A <= 62 else 60.0 * A**(1.0/3.0)

def coulomb_cross_section(Z_eff, A_eff, pc, theta_init):
    """
    Coulomb cross section [mb] for a material with effective Z and A.
    Does not depend on sqrt(s), only on beam momentum and geometry.
    """
    b_c   = _get_slope_hadron_nucleus(A_eff)
    R_c   = 2.0 * HBAR_C * math.sqrt(b_c)
    t_cut = (pc * 2.325 * theta_init)**2
    arg   = R_c**2 * b_c * t_cut
    E1    = float(exp1(arg))
    const = 4 * math.pi * Z_eff**2 * ALPHA_EM**2 * HBAR_C**2
    return -const * (R_c**2 * b_c * E1 - math.exp(-arg) / t_cut)


# ===========================================================================
# Step 5: Material flattening — recursively resolve to {(Z, A): molar_fraction}
# ===========================================================================

def _is_element(mat) -> bool:
    """True if mat is a pure element (has Z and A, no components list)."""
    return (
        hasattr(mat, 'Z') and mat.Z is not None
        and not (hasattr(mat, 'components') and mat.components)
    )

def _effective_A(mat) -> float:
    """Molar-fraction weighted mean atomic mass of a (possibly composite) material."""
    return sum(A * f for (Z, A), f in flatten_to_elements(mat).items())

def flatten_to_elements(mat) -> dict:
    """
    Recursively reduce any xcoll Material to {(Z, A): molar_fraction}.
    All molar fractions sum to 1.

    Compound    (n_atoms)         f_j = n_j / sum(n)
    Molar mix   (molar_fractions) f_j = molar_fraction_j  (normalised)
    Mass mix    (mass_fractions)  f_j = (w_j/A_j_eff) / sum(w_i/A_i_eff)
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
        # Convert mass fractions to molar fractions via effective A of each component.
        # For a pure element A_eff is just its atomic mass.
        # For a compound/mixture component (e.g. BoronNitride) it recurses.
        w      = np.array(mat.mass_fractions, dtype=float)
        A_effs = np.array([_effective_A(c) for c in components])
        num    = w / A_effs
        mol_fracs = num / num.sum()

    else:
        raise ValueError(
            f"Cannot determine composition of {mat!r}: "
            "no n_atoms, molar_fractions, or mass_fractions found."
        )

    # Recurse into each component and accumulate weighted elemental fractions
    combined: dict = defaultdict(float)
    for comp, frac in zip(components, mol_fracs):
        for (Z, A), child_frac in flatten_to_elements(comp).items():
            combined[(Z, A)] += frac * child_frac

    total = sum(combined.values())
    return {k: v / total for k, v in combined.items()}


def resolve_material(name: str):
    """Return an xcoll material object from its string name."""
    mat = getattr(xc.materials, name, None)
    if mat is None:
        raise ValueError(
            f"Material {name!r} not found in xc.materials. "
            f"Browse dir(xc.materials) for available names."
        )
    return mat


# ===========================================================================
# Steps 3-6 combined: sample one material over the full energy grid
# ===========================================================================

GG_KEYS = ["cs_tot_hA", "cs_inel_hA", "cs_el_hA",
           "cs_prod_hA", "cs_sd_hA",  "cs_qel_hA"]

def sample_material(
    material_name: str,
    pc_GeV:     float,
    theta_init: float,
    sqrt_s_min: float,
    sqrt_s_max: float,
    n_points:   int,
    log_scale:  bool = True,
) -> list[dict]:
    """
    Sample all cross sections for one material over the sqrt_s grid.
    Returns a list of row dicts (one per energy point).
    """
    mat        = resolve_material(material_name)
    elem_fracs = flatten_to_elements(mat)   # {(Z, A): molar_fraction}

    # Effective Z and A (molar-fraction weighted) used for Coulomb
    Z_eff = sum(Z * f for (Z, A), f in elem_fracs.items())
    A_eff = sum(A * f for (Z, A), f in elem_fracs.items())

    sqrt_s_arr = (
        np.logspace(math.log10(sqrt_s_min), math.log10(sqrt_s_max), n_points)
        if log_scale else
        np.linspace(sqrt_s_min, sqrt_s_max, n_points)
    )

    # Coulomb does not depend on sqrt_s — compute once
    cs_coul = coulomb_cross_section(Z_eff, A_eff, pc_GeV, theta_init)

    rows = []
    for sqrt_s in sqrt_s_arr:

        # Nucleon spline values at this energy
        pp_tot  = cs_tot_pp(sqrt_s)
        pp_el   = cs_el_pp(sqrt_s)
        pp_inel = cs_inel_pp(sqrt_s)
        pn_tot  = cs_tot_pn(sqrt_s)

        # Molar-fraction weighted GG cross sections
        gg = {k: 0.0 for k in GG_KEYS}
        for (Z, A), frac in elem_fracs.items():
            el = glauber_element(A, Z, sqrt_s)
            for k in GG_KEYS:
                gg[k] += frac * el[k]

        rows.append({
            # identity
            "material":         material_name,
            "A_eff":            A_eff,
            "Z_eff":            Z_eff,
            "sqrt_s_GeV":       sqrt_s,
            # nucleon spline cross sections [mb]
            "cs_tot_pp_mb":     pp_tot,
            "cs_el_pp_mb":      pp_el,
            "cs_inel_pp_mb":    pp_inel,
            "cs_tot_pn_mb":     pn_tot,
            # Glauber-Gribov nucleus cross sections [mb]
            "cs_tot_hA_mb":     gg["cs_tot_hA"],
            "cs_inel_hA_mb":    gg["cs_inel_hA"],
            "cs_el_hA_mb":      gg["cs_el_hA"],
            "cs_prod_hA_mb":    gg["cs_prod_hA"],
            "cs_sd_hA_mb":      gg["cs_sd_hA"],
            "cs_qel_hA_mb":     gg["cs_qel_hA"],
            # Coulomb [mb]
            "cs_coulomb_mb":    cs_coul,
        })

    return rows


# ===========================================================================
# Step 7a: Save CSV and NPZ
# ===========================================================================

def save_csv(rows: list[dict], filename: str = "cross_sections.csv") -> None:
    if not rows:
        return
    with open(filename, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    print(f"Saved {len(rows)} rows -> {filename}")


def save_npz(rows: list[dict], filename: str = "cross_sections.npz") -> None:
    if not rows:
        return
    arrays = {}
    for k in rows[0]:
        vals = [r[k] for r in rows]
        try:
            arrays[k] = np.array(vals, dtype=float)
        except (ValueError, TypeError):
            arrays[k] = np.array(vals, dtype=object)
    np.savez(filename, **arrays)
    print(f"Saved numpy arrays -> {filename}")


# ===========================================================================
# Step 7b: Save C header
# ===========================================================================

# Maps row dict key -> C array name suffix
_CS_FIELDS = [
    ("cs_tot_hA_mb",  "tot_hA"),
    ("cs_inel_hA_mb", "inel_hA"),
    ("cs_el_hA_mb",   "el_hA"),
    ("cs_prod_hA_mb", "prod_hA"),
    ("cs_sd_hA_mb",   "sd_hA"),
    ("cs_qel_hA_mb",  "qel_hA"),
    ("cs_coulomb_mb", "coulomb"),
    ("cs_tot_pp_mb",  "tot_pp"),
    ("cs_inel_pp_mb", "inel_pp"),
    ("cs_el_pp_mb",   "el_pp"),
]

def _c_safe(name: str) -> str:
    """Make a material name safe for use as a C identifier."""
    return name.replace(" ", "_").replace("-", "_")


def save_c_header(
    all_rows:   list[dict],
    sqrt_s_min: float,
    sqrt_s_max: float,
    filename:   str = "cross_sections_splines.h",
) -> None:
    """
    Write a C header with flat double[] sigma arrays for every material and
    cross section type.

    Because the grid is fixed log-spaced, the C code finds the right array
    index in O(1):
        i = (int)((log(sqrt_s) - XCOLL_CS_LOG_SQRT_S_MIN) / XCOLL_CS_LOG_STEP)
        if (i < 0)                  i = 0;
        if (i >= XCOLL_CS_N_POINTS) i = XCOLL_CS_N_POINTS - 1;
        double sigma = XCOLL_CS_Carbon_inel_hA[i];

    Grid parameters are written as #defines so the C code never needs to
    hard-code them.
    """
    # Group rows by material, preserving insertion order
    by_mat: dict[str, list[dict]] = {}
    for row in all_rows:
        by_mat.setdefault(row["material"], []).append(row)

    n_points = len(next(iter(by_mat.values())))
    log_min  = math.log(sqrt_s_min)
    log_max  = math.log(sqrt_s_max)
    log_step = (log_max - log_min) / (n_points - 1)

    lines = []

    # --- file header + grid #defines ---
    lines.append(textwrap.dedent(f"""\
        // copyright ############################### //
        // This file is part of the Xcoll package.   //
        // Copyright (c) CERN, 2026.                 //
        // ######################################### //
        //
        // AUTO-GENERATED by xcoll_cross_sections.py -- do not edit by hand.
        //
        // All sigma values are in mb.
        // Grid is log-spaced in sqrt_s [GeV].
        //
        // O(1) lookup:
        //   int i = (int)((log(sqrt_s) - XCOLL_CS_LOG_SQRT_S_MIN) / XCOLL_CS_LOG_STEP);
        //   if (i < 0)                  i = 0;
        //   if (i >= XCOLL_CS_N_POINTS) i = XCOLL_CS_N_POINTS - 1;
        //   double sigma = XCOLL_CS_Carbon_inel_hA[i];

        #ifndef XCOLL_CROSS_SECTIONS_SPLINES_H
        #define XCOLL_CROSS_SECTIONS_SPLINES_H

        #define XCOLL_CS_N_POINTS       {n_points}
        #define XCOLL_CS_SQRT_S_MIN     {sqrt_s_min}       /* [GeV] */
        #define XCOLL_CS_SQRT_S_MAX     {sqrt_s_max}    /* [GeV] */
        #define XCOLL_CS_LOG_SQRT_S_MIN {log_min:.15e}
        #define XCOLL_CS_LOG_SQRT_S_MAX {log_max:.15e}
        #define XCOLL_CS_LOG_STEP       {log_step:.15e}

    """))

    # --- one flat double[] per (material, cs_type) ---
    for mat_name, rows in by_mat.items():
        cname = _c_safe(mat_name)
        lines.append(f"/* {'='*62} */")
        lines.append(f"/* {mat_name}  "
                     f"(A_eff={rows[0]['A_eff']:.4g}, "
                     f"Z_eff={rows[0]['Z_eff']:.4g}) */")
        for row_key, c_suffix in _CS_FIELDS:
            arr_name = f"XCOLL_CS_{cname}_{c_suffix}"
            vals     = [f"{row[row_key]:.6e}" for row in rows]
            # 8 values per line for readability
            chunks   = [", ".join(vals[i:i+8]) for i in range(0, len(vals), 8)]
            body     = ",\n  ".join(chunks)
            lines.append(
                f"static const double {arr_name}[{n_points}] = {{\n"
                f"  {body}\n"
                f"}};\n"
            )

    # --- lookup table struct ---
    field_decls = "\n".join(
        f"    const double* {suf};"
        for _, suf in _CS_FIELDS
    )
    lines.append(textwrap.dedent(f"""\
        /* Lookup table entry -- one per material */
        typedef struct {{
            const char* name;
            double      A_eff;
            double      Z_eff;
        {field_decls}
        }} XcollCSEntry;

    """))

    # --- table initialiser ---
    entries = []
    for mat_name, rows in by_mat.items():
        cname = _c_safe(mat_name)
        field_inits = "\n    ".join(
            f"XCOLL_CS_{cname}_{suf},"
            for _, suf in _CS_FIELDS
        )
        entries.append(
            f'  {{ "{mat_name}", {rows[0]["A_eff"]:.6g}, {rows[0]["Z_eff"]:.6g},\n'
            f'    {field_inits}\n  }}'
        )

    n_mats = len(by_mat)
    lines.append(
        f"static const XcollCSEntry XCOLL_CS_TABLE[{n_mats}] = {{\n"
        + ",\n".join(entries)
        + "\n};\n"
    )
    lines.append(f"#define XCOLL_CS_N_MATERIALS {n_mats}\n")
    lines.append("#endif /* XCOLL_CROSS_SECTIONS_SPLINES_H */\n")

    with open(filename, "w") as f:
        f.write("\n".join(lines))
    print(f"Saved C header -> {filename}")


# ===========================================================================
# Default material list
# ===========================================================================
DEFAULT_MATERIALS = [
    # Pure elements
    "Beryllium", "Carbon", "Aluminium", "Copper", "Tungsten", "Lead",
    # Compounds  (n_atoms)
    "CarbonDioxide", "Water",
    # Molar mixtures
    "CopperDiamond", "MolybdenumGraphite",
    # Mass mixtures
    "Air", "Inermet180", "StainLessSteel316L",
]


# ===========================================================================
# Main
# ===========================================================================
def main():
    parser = argparse.ArgumentParser(
        description="Sample Xcoll cross sections vs sqrt(s) for a list of materials."
    )
    parser.add_argument(
        "--materials", nargs="+", default=DEFAULT_MATERIALS, metavar="NAME",
        help="xcoll material names -- elements, compounds, or mixtures"
    )
    parser.add_argument("--pc",         type=float, default=7000.0,
                        help="Beam momentum [GeV/c]")
    parser.add_argument("--theta",      type=float, default=1e-6,
                        help="Initial angle [rad]")
    parser.add_argument("--sqrt-s-min", type=float, default=None,
                        help="Min sqrt(s) [GeV]  (default: auto from spline data)")
    parser.add_argument("--sqrt-s-max", type=float, default=None,
                        help="Max sqrt(s) [GeV]  (default: auto from spline data)")
    parser.add_argument("--n-points",   type=int,   default=1000,
                        help="Number of energy grid points per material")
    parser.add_argument("--linear",     action="store_true",
                        help="Use linear spacing instead of log spacing")
    parser.add_argument("--out-csv",    default="cross_sections.csv")
    parser.add_argument("--out-npz",    default="cross_sections.npz")
    parser.add_argument("--out-h",      default="cross_sections_splines.h")
    args = parser.parse_args()

    # Use auto range unless overridden
    sqrt_s_min = args.sqrt_s_min if args.sqrt_s_min is not None else _GRID_SQRT_S_MIN
    sqrt_s_max = args.sqrt_s_max if args.sqrt_s_max is not None else _GRID_SQRT_S_MAX

    print(f"Sampling grid: {args.n_points} points, "
          f"sqrt(s) = {sqrt_s_min:.4g} to {sqrt_s_max:.4g} GeV\n")

    all_rows = []
    for name in args.materials:
        print(f"Sampling {name} ...")
        try:
            rows = sample_material(
                material_name=name,
                pc_GeV=args.pc,
                theta_init=args.theta,
                sqrt_s_min=sqrt_s_min,
                sqrt_s_max=sqrt_s_max,
                n_points=args.n_points,
                log_scale=not args.linear,
            )
            all_rows.extend(rows)
            print(f"  -> {len(rows)} points")
        except Exception as e:
            print(f"  WARNING: skipping {name!r}: {e}")

    if not all_rows:
        print("No data collected -- check material names and PDG data files.")
        return

    save_csv(all_rows, args.out_csv)
    save_npz(all_rows, args.out_npz)
    save_c_header(all_rows, sqrt_s_min, sqrt_s_max, args.out_h)
    print(f"\nDone.  {len(all_rows)} total rows across {len(args.materials)} materials.")


if __name__ == "__main__":
    main()