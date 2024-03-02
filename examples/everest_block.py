import numpy as np
import matplotlib.pyplot as plt

import xpart as xp
import xcoll as xc


block = xc.EverestBlock(length=1., material=xc.materials.Tungsten)

part = xp.Particles(x=np.zeros(1000000), energy0=450.e9)

block.track(part)

plt.hist(part.x, bins=200, density=True)
plt.hist(part.px*part.rpp, bins=200, density=True)
plt.show()
