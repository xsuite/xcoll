#!/bin/env python
import matplotlib.pyplot as plt
import xtrack as xt
import numpy as np
import xcoll as xc

n_steps = 1000
ch = xc.ChannellingDev(length=0.0004/n_steps)
x0 = 0.95 # in Angstrom
theta0 = -1e-10 # in urad
theta0 = 0

part = xt.Particles(x=x0, px=theta0)
x = [x0]
xp = [theta0]
for i in range(n_steps):
  ch.track(part)
  x.append(part.x[0])
  xp.append(part.px[0])
  part = xt.Particles(x=x[-1], px=xp[-1])

s = np.linspace(0, 0.0004, n_steps+1)
fig, ax = plt.subplots(2, 1, figsize=(12, 8))
ax[0].plot(s, x, label ='x')
ax[1].plot(s, xp, label ='theta')
ax[1].set_xlabel('s [m]')
ax[0].set_ylabel('x')
ax[0].grid(True)
ax[1].grid(True)
ax[0].legend(loc='upper right')
ax[1].legend(loc='upper right')
plt.show()

