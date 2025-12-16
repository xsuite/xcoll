#!/usr/bin/env python3
"""
Yoshida symplectic integrator coefficients.



    S_{2k}(h) = Π_i S_2(w_i h)

"""

import math
import csv
import os
from functools import lru_cache


# ------------------------------------------------------------
# Yoshida chi coefficients
# ------------------------------------------------------------
def yoshida_chi(order: int):
    p = order - 1
    two_p = 2.0 ** (1.0 / p)

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
        S_{2k}(h) = Π S_2(w_i h)
    """
    if order == 2:
        return (1.0,)

    chi1, chi0 = yoshida_chi(order)
    w_prev = yoshida_weights(order - 2)

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

    print(f"\nYoshida order {order}")
    print(f"Number of weights: {len(w)} (expected {3**((order//2)-1)})")
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
