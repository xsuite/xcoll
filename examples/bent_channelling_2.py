#!/bin/env python

import matplotlib.pyplot as plt
import numpy as np
import xtrack as xt
import xcoll as xc
import time


L_tot   = 0.0004
n_steps = 10_000
ds      = L_tot / n_steps
x0      = -0.95e-10
theta0  = 12e-6
s_axis  = np.linspace(0.0, L_tot, n_steps+1)

cases = {0: {'x0': -0.95e-10, 'theta0': 12e-6}}

def run_case(case_id : int , model : int , order : int , version : int , n : int):
    x0 = cases[case_id]['x0']
    theta0 = cases[case_id]['theta0']
    ds = L_tot / n
    ch = xc.BentChannellingDev(length=ds, method=model, variant=version, order=order)
    part = xt.Particles(x=x0, px=theta0)
    xs = np.zeros(n+1, dtype=np.float64)
    xs[0] = x0
    xps = np.zeros(n+1, dtype=np.float64)
    xps[0] = theta0
    s_axis = np.linspace(0.0, L_tot, n+1)
    part_pre = part.copy()
    ch.track(part_pre)
    t0 = time.perf_counter()
    for i in range(n):
        ch.track(part)
        xs[i+1] = part.x[0]
        xps[i+1] = part.px[0]
    elapsed = time.perf_counter() - t0
    return s_axis, xs, xps, elapsed


s_M2, x_M2, px_M2, elapsed_M2 = run_case(0, 2, 4, 1, 1000)
s_M3, x_M3, px_M3, elapsed_M3 = run_case(0, 3, 4, 1, 1000)
s_M4, x_M4, px_M4, elapsed_M4 = run_case(0, 4, 4, 1, 1000)

print(f"M2 took {1000*elapsed_M2} ms")
print(f"M3 took {1000*elapsed_M3} ms")
print(f"M4 took {1000*elapsed_M4} ms")

fig, ax =plt.subplots(3, 2)
ax[0][0].plot(s_M2, x_M2)
ax[1][0].plot(s_M3, x_M3)
ax[2][0].plot(s_M4, x_M4)
ax[0][1].plot(s_M2, px_M2)
ax[1][1].plot(s_M3, px_M3)
ax[2][1].plot(s_M4, px_M4)
plt.show()
