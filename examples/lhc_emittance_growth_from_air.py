import json
import numpy as np
from pathlib import Path
import time
start_time = time.time()
import sys, os, contextlib

import xobjects as xo
import xtrack as xt
import xcoll as xc


# On a modern CPU, this will take around 1h30 (5000 part*turn/s)
# ==============================================================
num_turns = 500
num_part  = 50000


# Load the line
# =============
path_in  = xc._pkg_root.parent / 'examples'
line = xt.Line.from_json(path_in / 'machines' / f'lhc_run3_b1_no_aper.json')


# Create a 2m block of air and insert it at IP5
# =============================================
P = 1.01325 # Standard air pressure at sea level in Bar
P_Torr = P*750.062 # Standard air pressure at sea level in Torr
L_rad0 = 301 # For air. Table with radiation lengths: https://cds.cern.ch/record/941314/files/p245.pdf
L_rad = L_rad0/(P_Torr/760)

air = xc.Material(radiation_length=L_rad, name="Air (1 atm 20C)")
airblock = xc.EverestBlock(length=2, material=air)

s_IP5 = 13329.289
line.insert_element(element=airblock, name="Air IP5", at_s=s_IP5-1)
line.build_tracker()


# Generate an initial distribution of particles ad twiss
# ======================================================
nemitt_x = 3.5e-6
nemitt_y = 3.5e-6
sigma_z  = 7.61e-2
bunch_intensity = 1e11

# Scattering need to be disabled to be able to twiss
line["Air IP5"]._tracking = False
part = xp.generate_matched_gaussian_bunch(
         num_particles=num_part, total_intensity_particles=bunch_intensity,
         nemitt_x=nemitt_x, nemitt_y=nemitt_y, sigma_z=sigma_z, line=line)
part._init_random_number_generator()
tw = line.twiss()
line["Air IP5"]._tracking = True


# Track and store emittance at every turn
# =======================================
eps_x = [np.std(np.pi*line.particle_ref.gamma0*(
        tw['gamx'][0]*part.x*part.x + 2*tw['alfx'][0]*part.x*part.px*part.rpp
        + tw['betx'][0]*part.px*part.px*part.rpp*part.rpp))]
for _ in range(200):
    line.track(part)    
    eps_x = [*eps_x, np.std(np.pi*line.particle_ref.gamma0*(
            tw['gamx'][0]*part.x*part.x + 2*tw['alfx'][0]*part.x*part.px*part.rpp
            + tw['betx'][0]*part.px*part.px*part.rpp*part.rpp))]

# Plot the result
# ===============
_ = plt.plot(eps_x)
