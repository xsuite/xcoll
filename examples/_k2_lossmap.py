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
num_particles = 50000

line = xt.Line.from_json(xc._pkg_root.parent / 'examples' / 'machines' / 'lhc_run3_b1.json')

colldb = xc.CollimatorDatabase.from_yaml(xc._pkg_root.parent / 'examples' / 'colldb' / 'lhc_run3_crystals.yaml',
                                         beam=beam, ignore_crystals=False)

colldb._install_k2_collimators(line=line, verbose=True)

line.build_tracker()

xc.assign_optics_to_collimators(line=line)

K2Engine.start(line=line, cwd='run_1')

part_init = xc.generate_pencil_on_collimator(line, 'tcp.d6l7.b1', num_particles, side='+', impact_parameter=1e-3)
part_init_corner = xc.generate_pencil_on_collimator(line, 'tcp.d6l7.b1', num_particles, side='+', impact_parameter=-0.4e-6)
part_init_drift = xc.generate_pencil_on_collimator(line, 'tcp.d6l7.b1', num_particles, side='+', impact_parameter=-9e-4)

part = part_init.copy()
part_corner = part_init_corner.copy()
part_drift = part_init_drift.copy()

xc.enable_scattering(line=line)
line['tcp.d6l7.b1'].track(part)
line['tcp.d6l7.b1'].track(part_corner)
line['tcp.d6l7.b1'].track(part_drift)

mask_alive = part.state>0
mask_alive_corner = part_corner.state>0
mask_alive_drift = part_drift.state>0

dri = xt.Drift(length=line['tcp.d6l7.b1'].length)
part_copy = part_init.copy()
part_copy_corner= part_init_corner.copy()
part_copy_drift = part_init_drift.copy()
dri.track(part_copy)
dri.track(part_copy_corner)
dri.track(part_copy_drift)

# checks before plotting
print(np.unique(part.state))
print(np.unique(part_corner.state))
print(np.unique(part_drift.state))
print(np.sum(part.state < 1))
print(np.sum(part_corner.state < 1))
print(np.sum(part_drift.state < 1))

ids = part.particle_id[mask_alive]
ids_corner = part_corner.particle_id[mask_alive_corner]
ids_drift = part_drift.particle_id[mask_alive_drift]

plt.figure(1)
plt.scatter(part_init.kin_yprime[ids],part.kin_yprime[mask_alive]-part_copy.kin_yprime[mask_alive], s=1)
plt.xlabel("Initial particles y'")
plt.ylabel("y' scattered - y' drifted" )
plt.title('Scatter vs drift inside collimator')
plt.savefig('/home/ssolstra/Documents/figures/K2_drift_vs_scatter_in_coll.png')
plt.show()

plt.figure(2)
plt.scatter(part_init_corner.kin_yprime[ids_corner],part_corner.kin_yprime[mask_alive_corner]-part_copy_corner.kin_yprime[mask_alive_corner], s=1)
plt.xlabel("Initial particles y'")
plt.ylabel("y' scattered - y' drifted" )
plt.title('Corner vs scatter on corner')
plt.savefig('/home/ssolstra/Documents/figures/K2_drift_vs_scatter_on_corner.png')
plt.show()

plt.figure(3)
plt.scatter(part_init_drift.kin_yprime[ids_drift],part_drift.kin_yprime[mask_alive_drift]-part_copy_drift.kin_yprime[mask_alive_drift], s=1)
plt.xlabel("Initial particles y'")
plt.ylabel("y' tracked - y' drifted" )
plt.title('Drift vs drift in drift')
plt.ylim(-0.0001,0.0001)
plt.savefig('/home/ssolstra/Documents/figures/K2_drift_vs_drift_in_drift.png')
plt.show()

# plt.figure(4)
# plt.plot(part2.y, part2.py, 'o', markersize=1)
# plt.xlabel("Initial particles y'")
# plt.ylabel("y' scattered - y' drifted" )
# plt.title('Scatter')
# plt.savefig('/home/ssolstra/Documents/figures/K2_initial.png')
# plt.show()

# plt.figure(5)
# plt.plot(part.y[mask_alive], part.py[mask_alive], 'o', markersize=1, color='g')
# plt.plot(part.y[~mask_alive], part.py[~mask_alive], 'o', markersize=1, color='r')
# plt.xlabel("y")
# plt.ylabel("py" )
# plt.axvline(x=line['tcp.d6l7.b1'].jaw_L, color='c')
# plt.title('Dead particles inside')
# plt.savefig('/home/ssolstra/Documents/figures/K2_dead.png')
# plt.show()

# plt.figure(6)
# plt.plot(part_corner.y[mask_alive_corner], part_corner.py[mask_alive_corner], 'o', markersize=1, color = 'r')
# plt.plot(part_corner.y[~mask_alive_corner], part_corner.py[~mask_alive_corner], 'o', markersize=1, color = 'g')
# plt.xlabel("Initial particles y'")
# plt.axvline(x=line['tcp.d6l7.b1'].jaw_L, color='c')
# plt.ylabel("y' scattered - y' drifted" )
# plt.title('Corner')
# plt.savefig('/home/ssolstra/Documents/figures/K2_dead_corner.png')
# plt.show()

# plt.figure(7)
# plt.plot(part_drift.y[mask_alive_drift], part_drift.py[mask_alive_drift], 'o', markersize=1, color = 'r')
# plt.plot(part_drift.y[~mask_alive_drift], part_drift.py[~mask_alive_drift], 'o', markersize=1, color = 'g')
# plt.xlabel("y")
# plt.axvline(x=line['tcp.d6l7.b1'].jaw_L, color='c')
# plt.ylabel("py" )
# plt.title('Drift part')
# plt.savefig('/home/ssolstra/Documents/figures/K2_dead_drift.png')
# plt.show()