#!/bin/env python
import time
import matplotlib.pyplot as plt
import xtrack as xt
import numpy as np
import xcoll as xc


n_steps = 1000
ch = xc.ChannellingDev(length=0.0004/n_steps)
x0 = -0.95 # in Angstrom
# x0 = 0
theta0 = 12 # in urad
# theta0 = 0

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
ch = xc.ChannellingDev(length=0.0004)
particles = xt.Particles(x=np.random.normal(0, 0.2, size=num_part), px=np.random.normal(0, 17, size=num_part))
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
ax[0].set_ylabel('x [Angstrom]')
ax[1].set_ylabel('theta [urad]')
ax[0].grid(True)
ax[1].grid(True)
plt.show()

