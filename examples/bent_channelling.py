#!/bin/env python

import matplotlib.pyplot as plt
import numpy as np
import xtrack as xt
import xcoll as xc

L_tot   = 0.0004
n_steps = 10_000
ds      = L_tot / n_steps
x0      = -0.95e-10
theta0  = 12e-6
s_axis  = np.linspace(0.0, L_tot, n_steps+1)

def run_with(ch, method, variant=2, order=12):
    ch.method  = method
    ch.variant = variant
    ch.order   = order

    part = xt.Particles(x=x0, px=theta0)
    xs  = [x0]
    xps = [theta0]
    for _ in range(n_steps):
        ch.track(part)           # ένα βήμα μήκους ds
        xs.append(part.x[0])
        xps.append(part.px[0])
        part = xt.Particles(x=xs[-1], px=xps[-1])  # όπως έχεις στο παράδειγμα
    return np.array(xs), np.array(xps)

# ίδιο element, runtime επιλογή
ch = xc.BentChannellingDev(length=ds, method=4, variant=2, order=12)

x_M2, px_M2 = run_with(ch, 2)
x_M3, px_M3 = run_with(ch, 3)
x_M4, px_M4 = run_with(ch, 4)

fig, ax = plt.subplots(2, 1, figsize=(12, 8), sharex=True)
ax[0].plot(s_axis, x_M2, label='M2')
ax[0].plot(s_axis, x_M3, label='M3')
ax[0].plot(s_axis, x_M4, label='M4')
ax[0].set_ylabel('x [m]')
ax[0].grid(True); ax[0].legend()

ax[1].plot(s_axis, px_M2, label='M2')
ax[1].plot(s_axis, px_M3, label='M3')
ax[1].plot(s_axis, px_M4, label='M4')
ax[1].set_xlabel('s [m]')
ax[1].set_ylabel('theta [rad]')
ax[1].grid(True); ax[1].legend()

plt.tight_layout(); plt.show()

















'''import time
import matplotlib.pyplot as plt
import xtrack as xt
import numpy as np
import xcoll as xc


n_steps = 10000
ch = xc.BentChannellingDev(length=0.0004/n_steps)
#x0 = -0.95 # in Angstrom
x0 = -0.95e-10 # [m]
#x0 = 0
#theta0 = 12 # in urad
theta0 = 12e-6 # [rad]

# TODO: if x0 = 0 and theta0 = 0 we get NaNs

# Particles
part = xt.Particles(x=x0, px=theta0)
x = [x0]
xp = [theta0]

# Track
for i in range(n_steps):
  ch.track(part)
  x.append(part.x[0])
  xp.append(part.px[0])
  part = xt.Particles(x=x[-1], px=xp[-1])


# Test performance
num_part = 100_000
ch = xc.BentChannellingDev(length=0.0004)
particles = xt.Particles(x=np.random.normal(0, 1e-12, size=num_part), px=np.random.normal(0, 1e-12, size=num_part))
t0_cpu = time.process_time()
t0_wall = time.perf_counter()
ch.track(particles)
t1_cpu = time.process_time()
t1_wall = time.perf_counter()

print(f"CPU time passed: {t1_cpu - t0_cpu:.6f} s for {num_part} particles")
print(f"Real time passed: {t1_wall - t0_wall:.6f} s for {num_part} particles")
s = np.linspace(0, 0.0004, n_steps+1)
fig, ax = plt.subplots(2, 1, figsize=(12, 8))
ax[0].plot(s, x)
ax[1].plot(s, xp)
ax[1].set_xlabel('s [m]')
ax[0].set_ylabel('x [m]')
ax[1].set_ylabel('theta [rad]')
ax[0].grid(True)
ax[1].grid(True)
plt.show()'''