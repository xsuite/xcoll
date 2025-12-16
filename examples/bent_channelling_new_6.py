#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import time
import xtrack as xt
import xcoll as xc
import os

# Only for qualitive and not quantitive analysis,
# because dont use the in-class slicing, but we slice more in order to visualise what's happeing.
# -----------------------------------------------------------
# Global settings
# -----------------------------------------------------------
L_tot = 4.0e-4  # total crystal length [m]

# Case 05
cases = {
    0: {'x0': 8.5e-11, 'theta0': 4.4e-7},
}

# -----------------------------------------------------------
# Load Mathematica reference trajectory for Case 05
# -----------------------------------------------------------
def load_reference_case05():
    filepath = "bent_channelling_reference_functions/solution_case05.csv"

    if not os.path.exists(filepath):
        raise FileNotFoundError(f"Reference file not found: {filepath}")

    try:
        arr = np.loadtxt(filepath, delimiter=",", skiprows=1)
    except ValueError:
        arr = np.loadtxt(filepath, delimiter=",")

    return arr[:, 0], arr[:, 1], arr[:, 2]  # s, x, px


# -----------------------------------------------------------
# Run one simulation and return trajectory
# -----------------------------------------------------------
def run_case(case_id, model, order, variant, n_samples):

    x0 = cases[case_id]['x0']
    px0 = cases[case_id]['theta0']

    ds = L_tot / n_samples

    ch = xc.BentChannellingDev(
        length=ds,
        method=model,
        order=order,
        variant=variant,
    )
    
    
   # Very important to set p0c and delta, because we worked
   # with bpc = 150e9 eV for Mathematica Reference functions.

    part = xt.Particles(p0c=1e9, x=x0, px=px0, delta=0.0)

    xs = np.zeros(n_samples + 1)
    pxs = np.zeros(n_samples + 1)

    xs[0] = x0
    pxs[0] = px0

    ch.track(part.copy())  # warm-up compile

    s_axis = np.linspace(0, L_tot, n_samples + 1)

    t0 = time.perf_counter()
    for i in range(n_samples):
        ch.track(part)
        xs[i+1]  = part.x[0]
        pxs[i+1] = part.px[0]
    elapsed = time.perf_counter() - t0

    return s_axis, xs, pxs, elapsed

# -----------------------------------------------------------
# Main execution
# -----------------------------------------------------------
if __name__ == "__main__":

    n_samples = 20

    # Run models
    s_M2, x_M2, px_M2, T_M2 = run_case(0, 2, 4, 1, n_samples)
    s_M3, x_M3, px_M3, T_M3 = run_case(0, 3, 4, 1, n_samples)
    s_M4, x_M4, px_M4, T_M4 = run_case(0, 4, 4, 1, n_samples)

    print(f"M2 time: {T_M2*1e3:.3f} ms")
    print(f"M3 time: {T_M3*1e3:.3f} ms")
    print(f"M4 time: {T_M4*1e3:.3f} ms")

    # Load Mathematica reference
    s_ref, x_ref, px_ref = load_reference_case05()

    # ======================================================================
    # (A) THREE SEPARATE SCATTER PLOTS (one for each model)
    # ======================================================================

    # ---- MODEL 2 ----
    fig2, ax2 = plt.subplots(1, 2, figsize=(12, 4))
    ax2[0].scatter(s_M2, x_M2, s=15, alpha=0.8, label="M2")
    ax2[0].plot(s_ref, x_ref, 'k-', lw=2, label="Mathematica")
    ax2[0].set_title("M2 vs Mathematica — x(s)")
    ax2[0].set_xlabel("s [m]")
    ax2[0].set_ylabel("x [m]")
    ax2[0].grid(True, ls=":")
    ax2[0].legend()

    ax2[1].scatter(s_M2, px_M2, s=15, alpha=0.8, label="M2")
    ax2[1].plot(s_ref, px_ref, 'k-', lw=2, label="Mathematica")
    ax2[1].set_title("M2 vs Mathematica — px(s)")
    ax2[1].set_xlabel("s [m]")
    ax2[1].set_ylabel("px")
    ax2[1].grid(True, ls=":")
    ax2[1].legend()

    plt.tight_layout()
    plt.show()

    # ---- MODEL 3 ----
    fig3, ax3 = plt.subplots(1, 2, figsize=(12, 4))
    ax3[0].scatter(s_M3, x_M3, s=15, alpha=0.8, label="M3")
    ax3[0].plot(s_ref, x_ref, 'k-', lw=2, label="Mathematica")
    ax3[0].set_title("M3 vs Mathematica — x(s)")
    ax3[0].set_xlabel("s [m]")
    ax3[0].set_ylabel("x [m]")
    ax3[0].grid(True, ls=":")
    ax3[0].legend()

    ax3[1].scatter(s_M3, px_M3, s=15, alpha=0.8, label="M3")
    ax3[1].plot(s_ref, px_ref, 'k-', lw=2, label="Mathematica")
    ax3[1].set_title("M3 vs Mathematica — px(s)")
    ax3[1].set_xlabel("s [m]")
    ax3[1].set_ylabel("px")
    ax3[1].grid(True, ls=":")
    ax3[1].legend()

    plt.tight_layout()
    plt.show()

    # ---- MODEL 4 ----
    fig4, ax4 = plt.subplots(1, 2, figsize=(12, 4))
    ax4[0].scatter(s_M4, x_M4, s=15, alpha=0.8, label="M4")
    ax4[0].plot(s_ref, x_ref, 'k-', lw=2, label="Mathematica")
    ax4[0].set_title("M4 vs Mathematica — x(s)")
    ax4[0].set_xlabel("s [m]")
    ax4[0].set_ylabel("x [m]")
    ax4[0].grid(True, ls=":")
    ax4[0].legend()

    ax4[1].scatter(s_M4, px_M4, s=15, alpha=0.8, label="M4")
    ax4[1].plot(s_ref, px_ref, 'k-', lw=2, label="Mathematica")
    ax4[1].set_title("M4 vs Mathematica — px(s)")
    ax4[1].set_xlabel("s [m]")
    ax4[1].set_ylabel("px")
    ax4[1].grid(True, ls=":")
    ax4[1].legend()

    plt.tight_layout()
    plt.show()


    # ======================================================================
    # (B) ONE COMBINED COMPARISON PLOT (all models together)
    # ======================================================================

    figC, axC = plt.subplots(1, 2, figsize=(14, 5))

    # x(s)
    axC[0].scatter(s_M2, x_M2, s=12, label="M2")
    axC[0].scatter(s_M3, x_M3, s=10, label="M3")
    axC[0].scatter(s_M4, x_M4, s=8, label="M4")
    axC[0].plot(s_ref, x_ref, 'k-', lw=2, label="Mathematica")

    axC[0].set_title("Comparison: M2, M3, M4 vs Mathematica — x(s)")
    axC[0].set_xlabel("s [m]")
    axC[0].set_ylabel("x [m]")
    axC[0].grid(True, ls=":")
    axC[0].legend()

    # px(s)
    axC[1].scatter(s_M2, px_M2, s=12, label="M2")
    axC[1].scatter(s_M3, px_M3, s=12, label="M3")
    axC[1].scatter(s_M4, px_M4, s=12, label="M4")
    axC[1].plot(s_ref, px_ref, 'k-', lw=2, label="Mathematica")

    axC[1].set_title("Comparison: M2, M3, M4 vs Mathematica — px(s)")
    axC[1].set_xlabel("s [m]")
    axC[1].set_ylabel("px")
    axC[1].grid(True, ls=":")
    axC[1].legend()

    plt.tight_layout()
    plt.show()

