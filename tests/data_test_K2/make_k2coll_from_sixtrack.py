from pathlib import Path
import json
import numpy as np

import sixtrack_files as sf
import sixtracktools as st
import xobjects as xo
import xtrack as xt
import xpart as xp
import xcoll as xc


print("Loading beam 1...")
path = Path(Path.cwd(),'SixTrack_B1')
line = st.SixInput(path).generate_xtrack_line()
line.particle_ref = xp.Particles(mass0=xp.PROTON_MASS_EV, p0c=np.sqrt(6.8e12**2 - xp.PROTON_MASS_EV**2) )
coll_manager = xc.CollimatorManager(
                    line=line,
                    colldb=xc.load_SixTrack_colldb(path / 'CollDB-Testing.dat', emit=3.5e-6)
                )

print("Installing collimators...")
coll_manager.install_k2_collimators(verbose=False, seed=6574654)
coll_manager.align_collimators_to('front')
tracker = coll_manager.build_tracker()

print("Setting openings...")
coll_manager.set_openings()

print("Dumping json's...")
for name in coll_manager.collimator_names:
    file = Path(Path.cwd(),'Collimators',name + '.json')
    with open(file, 'w') as fid:
        json.dump(line[name].to_dict(), fid, cls=xo.JEncoder)



print("Loading beam 2...")
path = Path(Path.cwd(),'SixTrack_B2')
line = st.SixInput(path).generate_xtrack_line()
line.particle_ref = xp.Particles(mass0 = xp.PROTON_MASS_EV, p0c=np.sqrt(6.8e12**2 - xp.PROTON_MASS_EV**2)
coll_manager = xc.CollimatorManager(
                    line=line,
                    colldb=xc.load_SixTrack_colldb(path / 'CollDB-Testing.dat', emit=3.5e-6)
                )

print("Installing collimators...")
coll_manager.install_k2_collimators(verbose=False, seed=6574654)
coll_manager.align_collimators_to('front')
tracker = coll_manager.build_tracker()

print("Setting openings...")
coll_manager.set_openings()

print("Dumping json's...")
for name in coll_manager.collimator_names:
    file = Path(Path.cwd(),'Collimators',name + '.json')
    with open(file, 'w') as fid:
        json.dump(line[name].to_dict(), fid, cls=xo.JEncoder)

print("Done.")

                                 
                                 
                                 