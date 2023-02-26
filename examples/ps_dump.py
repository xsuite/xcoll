import json
import numpy    as np
import matplotlib.pyplot as plt

import xobjects as xo
import xtrack   as xt
import xpart    as xp
import xcoll    as xc


## Define beam parameters
emitt_x=2.5e-6
emitt_y=2.5e-6
sigma_z=9.0e-2
BRHO   =2.794987e9
betx0=21.03871984
bety0=11.41817476



# Load from json
with open('machines/ps_b1.json', 'r') as fid:
    loaded_dct = json.load(fid)
line = xt.Line.from_dict(loaded_dct)



# Initialise collimator
coll_manager = xc.CollimatorManager(
    line=line,
    colldb=xc.load_SixTrack_colldb('colldb/ps_b1.dat', emit=emitt_x)
    )

coll_manager.colldb.parking = 1
coll_manager.install_everest_collimators(verbose=True, seed=6554)
coll_manager.align_collimators_to('front')

coll_manager.build_tracker()
# coll_manager.set_openings()



# We set the openings etc manually
coll = 'pr.tdi47'
line[coll].jaw_LU = 0.001
line[coll].jaw_LD = 0.001
line[coll].ref_xU = 0
line[coll].ref_yU = 0
line[coll].dpx = 0
line[coll].dpy = 0
line[coll].material = 'BE'
line[coll].is_active = True
line[coll].onesided = True

coll = 'pr.tdi48'
line[coll].jaw_LU = 0.001
line[coll].jaw_LD = 0.001
line[coll].ref_xU = 0
line[coll].ref_yU = 0
line[coll].dpx = 0
line[coll].dpy = 0
line[coll].material = 'BE'
line[coll].is_active = True
line[coll].onesided = True



# Make initial particles
n_part = 10000
particles = xp.Particles(p0c=BRHO,
                         x=np.random.normal(loc=0.0, scale=np.sqrt(emitt_x*betx0), size=n_part),
                         y=np.random.normal(loc=0.0, scale=np.sqrt(emitt_y*bety0), size=n_part)
                        )

part_init = particles.copy()



# Track one turn, starting at the first dump block
start = 'pr.tdi47'
particles.at_element = line.element_names.index(start)
particles.s = line.get_s_position(start)
particles.start_tracking_at_element = line.element_names.index(start)
coll_manager.track(particles, num_turns=1)

# The survival flags are sorted as surviving particles first,
# hence we need to 'unsort' them using their IDs
surv = particles.state.copy()
surv[:] = surv[np.argsort(particles.particle_id)]

n = 0.025
# Plot the surviving particles as green, scattered particles as orange (which are lost a bit further in the aperture)
plt.figure(1,figsize=(20,20))
plt.plot(part_init.x, part_init.y, '.', color='red')
plt.plot(part_init.x[surv>0], part_init.y[surv>0], '.', color='green')
plt.plot(part_init.x[surv==0], part_init.y[surv==0], '.', color='orange')
plt.axis('equal')
plt.axis([n, -n, -n, n])
plt.show()



## Plot histogram of losses along the acceleratores
## ------------------------------------------------------------------

wdth=0.5;
particles_l=particles.filter(particles.state<=0)  # Lost particles
S, count = np.unique(np.floor(particles_l.s/wdth)*wdth, return_counts = True);

fig, ax=plt.subplots(1,1,figsize=(30,7));
plt.bar(S+wdth*.5, count, width = wdth);
ax.set_xlim(-10, 650);
ax.set_xlabel('S [m]');
ax.set_ylabel('No. of lost particles');
plt.show()




# ----------------------------------------------------------------------------------------------------
# Track, and after every turn move the dump block more inwards
coll = 'pr.tdi47'
line[coll].jaw_LU = 0.06
line[coll].jaw_LD = 0.06
coll = 'pr.tdi48'
line[coll].jaw_LU = 0.06
line[coll].jaw_LD = 0.06
part = {}

step = 0.005
steps = int(0.06*2/step +1)

for turn in range(1,steps+1):
    coll = 'pr.tdi47'
    line[coll].jaw_LU -= step
    line[coll].jaw_LD -= step
    coll = 'pr.tdi48'
    line[coll].jaw_LU -= step
    line[coll].jaw_LD -= step
    coll_manager.track(particles, num_turns=1)
    part[turn] = particles.copy()
    
    
    
