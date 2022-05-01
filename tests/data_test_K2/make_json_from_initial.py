from pathlib import Path
import json
import numpy as np

import sixtrack_files as sf
import xpart as xp
import xobjects as xo

p0c = np.sqrt(6.8e12**2 - xp.PROTON_MASS_EV**2)
file = Path(Path.cwd(),'initial.dat')
part = sf.sixtrack_initial_to_particles(file, p0c=p0c)
print(part.x[0])
print(part.y[0])
print(part.px[0]/(1+part.delta[0])) # xp
print(part.py[0]/(1+part.delta[0])) # yp
print(part.zeta[0]/part.rvv[0])     # sigma
print(part.delta[0])
print((part.delta[0]+1)*p0c)       # p
print(part.energy[0])

file = Path(Path.cwd(),'initial.json')
with open(file, 'w') as fid:
    json.dump(part.to_dict(), fid, cls=xo.JEncoder)