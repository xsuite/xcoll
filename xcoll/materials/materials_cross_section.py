# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2026.                 #
# ######################################### #
import numpy as np

from .glauber_gribov import GG_KEYS, NUCLEON_KEYS


# def load_element_cache(filename="elements_cs.npz"):
#     """
#     Load the NPZ produced by stage1_elements.py.
#     """
#     data       = np.load(filename)
#     sqrt_s_arr = data["sqrt_s"]
#     Zs         = data["element_Z"].astype(int)
#     As         = data["element_A"]

#     element_cs = {}

#     # Nucleon-level arrays — element-independent, stored under Z=0
#     element_cs[0] = {key: data[f"nucleon_{key}"] for key in NUCLEON_KEYS}

#     # Per-element GG arrays
#     for Z, A in zip(Zs, As):
#         element_cs[Z] = {"A": float(A)}
#         for key in GG_KEYS:
#             element_cs[Z][key] = data[f"Z{Z}_{key}"]

#     return element_cs, sqrt_s_arr


# def get_material_cs(mat, element_cs, sqrt_s_arr):
#     """
#     Combine pre-computed elemental GG arrays for one Material instance.

#     Uses mat.molar_fractions and mat.components directly — the Material
#     class already flattens compounds and mixtures to pure elements at
#     construction time, so no custom recursion is needed here.
#     """
#     if mat.is_elemental:
#         # Single element: one component, weight = 1
#         pairs = [(int(mat.Z), 1.0)]
#     else:
#         # Compound or mixture: Material class gives us molar_fractions
#         # and components already flattened to pure elements
#         pairs = [
#             (int(el.Z), float(frac))
#             for el, frac in zip(mat.components, mat.molar_fractions)
#         ]

#     # Weighted sum of elemental GG arrays
#     combined = {k: np.zeros(len(sqrt_s_arr)) for k in GG_KEYS}
#     for Z, frac in pairs:
#         if Z not in element_cs:
#             raise KeyError(
#                 f"Z={Z} not found in element cache. "
#                 "Re-run stage1_elements.py to rebuild the cache."
#             )
#         for k in GG_KEYS:
#             combined[k] += frac * element_cs[Z][k]

#     return {
#         "sqrt_s": sqrt_s_arr,
#         **combined,
#         **{k: element_cs[0][k] for k in NUCLEON_KEYS},
#     }