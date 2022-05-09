import json
import numpy    as np
import matplotlib.pyplot as plt

import xobjects as xo
import xtrack   as xt
import xpart    as xp
import xcoll    as xc



# Make a context and get a buffer
context = xo.ContextCpu()         # For CPU
# context = xo.ContextCupy()      # For CUDA GPUs
# context = xo.ContextPyopencl()  # For OpenCL GPUs
buffer = context.new_buffer()



# Load from json
with open('machines/lhc_run3_b1.json', 'r') as fid:
    loaded_dct = json.load(fid)
line = xt.Line.from_dict(loaded_dct)



# Initialise collmanager,on the specified buffer
coll_manager = xc.CollimatorManager(
    line=line,
    colldb=xc.load_SixTrack_colldb('colldb/lhc_run3_b1.dat', emit=3.5e-6),
    _buffer=buffer
    )

# Install collimators in line as black absorbers
coll_manager.install_k2_collimators(verbose=True)

# Build the tracker
tracker = coll_manager.build_tracker()

# Align the collimators
coll_manager.align_collimators_to('front')

# Set the collimator openings based on the colldb,
# or manually override with the option gaps={collname: gap}
coll_manager.set_openings()





# Horizontal loss map
num_particles = 50000
coll = 'tcp.c6l7.b1'

# Collimator plane: generate pencil distribution in normalized coordinates
x_norm, py_norm, _, _ = xp.generate_2D_pencil(
                             num_particles=num_particles,
                             pos_cut_sigmas=coll_manager.colldb.gap[coll],
                             dr_sigmas=0.002,
                             side='+-')

# Other plane: generate gaussian distribution in normalized coordinates
y_norm = np.random.normal(scale=0.01, size=num_particles)
py_norm = np.random.normal(scale=0.01, size=num_particles)

part = xp.build_particles(
            tracker=coll_manager.tracker,
            x_norm=x_norm, px_norm=px_norm,
            y_norm=y_norm, py_norm=py_norm,
            scale_with_transverse_norm_emitt=coll_manager.colldb.emittance,
            at_element=coll,
            match_at_s=coll_manager.s_match[coll])



# Track
tracker.track(part, num_turns=10)


collimator_losses = part.s[part.state==-333]
aperture_losses = part.s[part.state==-1]
other_losses = part.s[(part.state<=0) & (part.state!=-333)  & (part.state!=-1)]
if len(other_losses) > 0:
    raise Exception("Unexpected losses type. Please investigate.")

## Plot histogram of losses along the acceleratores
## ------------------------------------------------------------------

wdth=0.1;
S, count = np.unique(np.floor(aperture_losses/wdth)*wdth, return_counts = True);

fig, ax=plt.subplots(1,1,figsize=(30,7));
plt.bar(S+wdth*.5, count, width = wdth);
ax.set_xlim(-10, 650);
ax.set_xlabel('S [m]');
ax.set_ylabel('No. of lost particles');
plt.show()

