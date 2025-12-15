#!/usr/bin/env python3

import numpy as np
import time
import xtrack as xt
import xcoll as xc
import os
import itertools
import pandas as pd

# ============================================================
# Global settings
# ============================================================

L_tot = 4.0e-4  # crystal length [m]

# ============================================================
# This script loads all Mathematica reference functions and takes:
#   - initial conditions from the first row
#   - reference final point from last row
# ============================================================

def load_all_references():
    refs = {}
    base = "bent_channelling_reference_functions"

    for case_id in range(13):  # case00 → case12
        fname = f"solution_case{case_id:02d}.csv"
        path = os.path.join(base, fname)

        if not os.path.exists(path):
            raise FileNotFoundError(path)

        arr = np.loadtxt(path, delimiter=",", skiprows=1)

        refs[case_id] = dict(
            # initial conditions 
            x0=arr[0, 1],
            px0=arr[0, 2],

            # references 
            s_ref=arr[:, 0],
            x_ref=arr[:, 1],
            px_ref=arr[:, 2],

            # final reference point
            x_last=arr[-1, 1],
            px_last=arr[-1, 2],
        )

    return refs


references = load_all_references()

# ============================================================
# Scan parameters
# ============================================================

methods   = [2, 3, 4]
orders    = [2, 4, 6, 8, 10, 12]
variants  = [1, 2]
# I found it easier to caclculate the steps and then find the steps per period, as long as 
# we take the harmonic period as a reference (constant - independent of x, xp).
# All of our cases have ~10-13 periods for L = 4e-4 m.
n_steps_list = list(range(10, 601, 10))

# ============================================================
# Run ONE configuration
# ============================================================

def run_case(case_id, method, order, variant, n_steps):

    ref = references[case_id]

    crystal = xc.BentChannellingDev(
        length=L_tot,
        method=method,
        order=order,
        variant=variant,
        n_steps=n_steps,  # explicit override
    )
#Particles(p0c=150e9, x=x0, px=px0, delta=0.0)
    part = xt.Particles(p0c=150e9,
        x=ref['x0'],
        px=ref['px0'], delta=0.0
    )

    # warm-up compilation
    crystal.track(part.copy())

    # timing
    t0 = time.perf_counter()
    crystal.track(part)
    runtime = time.perf_counter() - t0

    x, px = part.x[0], part.px[0]

    abs_err_x  = abs(x  - ref['x_last'])
    abs_err_px = abs(px - ref['px_last'])
    # avoiding devision by zero
    rel_err_x  = abs_err_x  / max(abs(ref['x_last']), 1e-30)
    rel_err_px = abs_err_px / max(abs(ref['px_last']), 1e-30)

    return dict(
        case=case_id,
        method=method,
        order=order,
        variant=variant,
        n_steps=n_steps,
        runtime_s=runtime,
        abs_err_x=abs_err_x,
        abs_err_px=abs_err_px,
        rel_err_x=rel_err_x,
        rel_err_px=rel_err_px,
    )

# ============================================================
# Main loop
# ============================================================

rows = []

for combo in itertools.product(
        range(13), methods, orders, variants, n_steps_list):

    case_id, method, order, variant, n_steps = combo

    try:
        row = run_case(case_id, method, order, variant, n_steps)
        rows.append(row)

        print(
            f"OK | case{case_id:02d} "
            f"M{method} O{order} v{variant} n={n_steps}"
        )

    except Exception as err:
        print(
            f"FAIL | case{case_id:02d} "
            f"M{method} O{order} v{variant} n={n_steps}"
        )
        print(err)

# ============================================================
# Save results
# ============================================================

df = pd.DataFrame(rows)
df.to_csv("bent_channelling_all_cases_scan.csv", index=False)

print("\nSaved: bent_channelling_all_cases_scan.csv")


