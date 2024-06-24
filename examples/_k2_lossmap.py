import json
import numpy as np
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt
import xobjects as xo
import xpart as xp
import xtrack as xt
import xcoll as xc
from xcoll.scattering_routines.k2 import K2Engine

beam          =  1
plane         = 'V'
num_turns     = 20
num_particles = 20000

line = xt.Line.from_json(xc._pkg_root.parent / 'examples' / 'machines' / 'lhc_run3_b1.json')

colldb = xc.CollimatorDatabase.from_yaml(xc._pkg_root.parent / 'examples' / 'colldb' / 'lhc_run3_crystals.yaml',
                                         beam=beam, ignore_crystals=False)

colldb._install_k2_collimators(line=line, verbose=True)

line.build_tracker()

xc.assign_optics_to_collimators(line=line)

K2Engine.start(line=line, cwd='run_1')

part_init = xc.generate_pencil_on_collimator(line, 'tcp.d6l7.b1', num_particles, side='+')
part = part_init.copy()

xc.enable_scattering(line=line)

print(np.unique(part.s))
line['tcp.d6l7.b1'].track(part)
print(np.unique(part.s))
print(np.unique(part.state))


dri = xt.Drift(length=line['tcp.d6l7.b1'].length)
part2 = part_init.copy()
dri.track(part2)

plt.plot(part_init.kin_yprime, part.kin_yprime - part2.kin_yprime, 'o', markersize=1)
plt.show()
