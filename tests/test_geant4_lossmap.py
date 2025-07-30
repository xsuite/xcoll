import json
import matplotlib.colors
import numpy as np
from pathlib import Path
import xpart as xp
import xtrack as xt
import xcoll as xc
import collimasim as cs
import time


path = xc._pkg_root.parent / 'tests' / 'data'
beam = 1
plane = 'H'
npart = int(5e5)
interpolation = 0.1

line = xt.Line.from_json(path / f'sequence_lhc_run3_b{beam}.json')
colldb = xc.CollimatorDatabase.from_yaml(path / f'colldb_lhc_run3.yaml', beam=beam, ignore_crystals=True)
colldb.install_geant4_collimators(verbose=False, line=line, random_seed=1336, names=colldb.collimator_names,
                                  bdsim_config_file=str(path / f'geant4_protons.gmad'))

df_with_coll = line.check_aperture()
assert not np.any(df_with_coll.has_aperture_problem)

line['tdisa.a4l2.b1'].active = 0
line['tdisb.a4l2.b1'].active = 0
line['tdisc.a4l2.b1'].active = 0
line['tclia.4r2'].active = 0
line['tclib.6r2.b1'].active = 0
line['tcld.a11r2.b1'].active = 0
line['tcpcv.a6l7.b1'].active = 0
line['tcpch.a4l7.b1'].active = 0

line.build_tracker()
xc.assign_optics_to_collimators(line=line)
tcp = f"tcp.{'c' if plane=='H' else 'd'}6{'l' if beam==1 else 'r'}7.b{beam}"

part = xc.generate_pencil_on_collimator(line, tcp, num_particles=npart, impact_parameter=4e-6, _capacity=2*npart)

xc.Geant4Engine.start(line=line, seed=1336, particle_ref='proton', p0c=7e12, relative_energy_cut=0.1,
                      bdsim_config_file=str(path / 'geant4_protons.gmad'))

t00 = time.time()
xc.enable_scattering(line)
line.track(part, num_turns=1)
xc.disable_scattering(line)
t11 = time.time()

print(f'total time: {t11-t00}')

assert np.sum(part.state == -333) > 80000
print(np.sum(part.state == -333))
assert np.sum(part.state == 0) > 10000
print(np.sum(part.state == 0))
assert len(part.s[(part.state==0) & (part.s > 20290) & (part.s < 20538)]) > 10
