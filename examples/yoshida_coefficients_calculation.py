#!/usr/bin/env python3
"""
Yoshida symplectic integrator coefficients (triple-jump recursion).

We build a 2k-th order symmetric composition from a (2k-2)-th order symmetric method:
    S_{2k}(h) = S_{2k-2}(chi1*h)  S_{2k-2}(chi0*h)  S_{2k-2}(chi1*h)

Finally, we flatten everything into weights w_i such that:
    S_{2k}(h) = Π_i S_2(w_i h)

This script fixes the common mistake in the recursion: the chi-coefficients must be
computed using the order of the method being composed (order-2), not the target order.
(Your original code computed chi from the target order; see yoshida_weights() in
yoshida_coefficients_calculation.py :contentReference[oaicite:0]{index=0})
"""

import math
import csv
import os
from functools import lru_cache


# ------------------------------------------------------------
# Yoshida chi coefficients
# ------------------------------------------------------------
def yoshida_chi(base_order: int):
    """
    Returns (chi1, chi0) for composing a symmetric method of even order `base_order`
    to obtain a method of order base_order + 2 via the triple-jump.

    The formula uses 2^(1/(base_order+1)).
    """
    if base_order % 2 != 0 or base_order < 2:
        raise ValueError("base_order must be an even integer >= 2.")

    two_p = 2.0 ** (1.0 / (base_order + 1.0))

    chi1 = 1.0 / (2.0 - two_p)
    chi0 = -two_p / (2.0 - two_p)

    return chi1, chi0


# ------------------------------------------------------------
# Recursive to flat Yoshida weights
# ------------------------------------------------------------
@lru_cache(None)
def yoshida_weights(order: int):
    """
    Returns tuple w_i such that:
        S_{order}(h) = Π_i S_2(w_i h)
    for even order in {2,4,6,...} constructed via repeated triple-jump recursion.
    """
    if order % 2 != 0 or order < 2:
        raise ValueError("order must be an even integer >= 2.")

    if order == 2:
        return (1.0,)

    # Compose the previous (order-2) method; chi must be computed from that base order.
    base_order = order - 2
    chi1, chi0 = yoshida_chi(base_order)

    w_prev = yoshida_weights(base_order)

    w = (
        tuple(chi1 * wi for wi in w_prev)
        + tuple(chi0 * wi for wi in w_prev)
        + tuple(chi1 * wi for wi in w_prev)
    )
    return w


# ------------------------------------------------------------
# Printing
# ------------------------------------------------------------
def print_weights(order: int, precision: int = 17):
    w = yoshida_weights(order)
    fmt = f"{{:.{precision}g}}"

    expected = 3 ** ((order // 2) - 1)

    print(f"\nYoshida order {order}")
    print(f"Number of weights: {len(w)} (expected {expected})")
    print(f"Sum of weights: {sum(w):.15e}\n")

    for i, wi in enumerate(w):
        print(f"{fmt.format(wi)},", end="")
        if (i + 1) % 3 == 0:
            print()
    print()


# ------------------------------------------------------------
# Save to CSV
# ------------------------------------------------------------
def save_weights_csv(order: int, outdir: str = "yoshida_coefficients"):
    os.makedirs(outdir, exist_ok=True)

    w = yoshida_weights(order)
    filename = os.path.join(outdir, f"yoshida_weights_order_{order}.csv")

    with open(filename, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["index", "weight"])
        for i, wi in enumerate(w):
            writer.writerow([i, f"{wi:.18e}"])

    print(f"Saved: {filename}")


if __name__ == "__main__":

    orders = [2, 4, 6, 8, 10, 12]

    for order in orders:
        print_weights(order)
        save_weights_csv(order)

