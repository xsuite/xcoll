import json
import numpy as np

import xtrack as xt
import xpart as xp
import xcoll as xc


print('Load file...')
with open('./HL_LHC_v1p5_clean_feb2022/HL_LHC_v1p5_line.json') as fid:
    dct = json.load(fid)
print('Build line...')
line = xt.Line.from_dict(dct)

# Attach reference particle (a proton a 7 TeV)
line.particle_ref = xp.Particles(mass0 = xp.PROTON_MASS_EV, p0c=7e12)

# Switch on RF (needed to twiss)
line['acsca.a5l4.b1'].voltage = 16e6
line['acsca.a5l4.b1'].frequency = 1e6

coll_manager = xc.CollimatorManager(
    coll_db_txt_file='HL_LHC_v1p5_clean_feb2022/CollDB_HL_relaxed_b1.data',
    nemitt_ref_x=2.5e-6, nemitt_ref_y=2.5e-6)

coll_manager.install_collimators_in_line(line=line) # They are left open

# Build tracker
tracker = line.build_tracker()

# Compute half-gaps and close collimators (black absorbers)
coll_manager.set_collimator_openings(collimator_names='all')

# Characterize machine aperture
n_sigmas = 30
n_part = 10000
x_norm = np.random.uniform(-n_sigmas, n_sigmas, n_part)
y_norm = np.random.uniform(-n_sigmas, n_sigmas, n_part)
part = xp.build_particles(tracker=tracker, x_norm=x_norm, y_norm=y_norm,
                          scale_with_transverse_norm_emitt=(2.5e-6, 2.5e-6),
                         )
tracker.track(part, num_turns=5)



# Serializeable

# coll_manager.initialize_geant4_engine(.....)
# coll_manager_initialize_k2_engine(....)
# 
# coll_manager.set_engine(engine='k2', names='all')
# coll_manager.set_engine(engine='k2', names=('tcp.d6l7.b1', 'tcp.c6l7.b1', 'tcp.b6l7.b1',))
# coll_manager.set_engine(engine='geant4',
#                         names=('tcsg.a6l7.b1',
#                             'tcsg.b5l7.b1',
#                             'tcsg.a5l7.b1',
#                             'tcsg.d4l7.b1',
#                             'tcspm.b4l7.b1',
#                             'tcsg.a4l7.b1',
#                             'tcsg.a4r7.b1',
#                             'tcsg.b5r7.b1',
#                             'tcsg.d5r7.b1',
#                             'tcspm.e5r7.b1',
#                             'tcspm.6r7.b1',)
#                         )



state_sorted = part.state.copy()
state_sorted[part.particle_id] = part.state
import matplotlib.pyplot as plt
plt.close('all')
plt.figure(1)
plt.plot(x_norm, y_norm, '.', color='red')
plt.plot(x_norm[state_sorted>0], y_norm[state_sorted>0], '.', color='green')
plt.axis('equal')
plt.show()



