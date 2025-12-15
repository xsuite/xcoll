#!/usr/bin/env python3

import pandas as pd
import numpy as np

# ============================================================
# Load results
# ============================================================

df = pd.read_csv("bent_channelling_all_cases_scan.csv")

# Convenience: combined error metric
df["rel_err_max"] = df[["rel_err_x", "rel_err_px"]].max(axis=1)

# ============================================================
# Containers
# ============================================================

best_rows = []
failed_rows = []

# ============================================================
# Loop over (case, method)
# ============================================================

grouped = df.groupby(["case", "method"])

for (case_id, method), g in grouped:

    # --------------------------------------------------------
    # (a) MIN STEPS with acceptable accuracy (<= 1 percent)
    # --------------------------------------------------------
    acceptable = g[
        (g["rel_err_x"] <= 1e-2) &
        (g["rel_err_px"] <= 1e-2)
    ]

    if not acceptable.empty:
        best_fast = acceptable.loc[
            acceptable["n_steps"].idxmin()
        ]
        best_fast = best_fast.copy()
        best_fast["selection"] = "min_steps_under_1pct"
        best_rows.append(best_fast)

    # --------------------------------------------------------
    # (b) MINIMUM error overall (no matter the steps and the runtime)
    # --------------------------------------------------------
    finite = g[np.isfinite(g["rel_err_max"])]

    if not finite.empty:
        best_acc = finite.loc[
            finite["rel_err_max"].idxmin()
        ]
        best_acc = best_acc.copy()
        best_acc["selection"] = "min_error"
        best_rows.append(best_acc)

    # --------------------------------------------------------
    # FAILED entries: NaNs with MAX n_steps 
    # Did that to see what would be a safe number for the default n steps.
    # Turns out that the higher the order tthe more steps we need for M2/M3 to be avoid infinity.
    # --------------------------------------------------------
    failed = g[
        (~np.isfinite(g["rel_err_x"])) |
        (~np.isfinite(g["rel_err_px"]))
    ]

    if not failed.empty:
        max_steps = failed["n_steps"].max()
        worst = failed[failed["n_steps"] == max_steps]

        for _, row in worst.iterrows():
            row = row.copy()
            row["selection"] = "nan_max_steps"
            failed_rows.append(row)

# ============================================================
# Save outputs
# ============================================================

df_best = pd.DataFrame(best_rows)
df_failed = pd.DataFrame(failed_rows)

df_best.to_csv(
    "bent_channelling_all_cases_best_entries.csv",
    index=False
)

df_failed.to_csv(
    "bent_channelling_all_cases_failed_entries.csv",
    index=False
)

print("Saved:")

