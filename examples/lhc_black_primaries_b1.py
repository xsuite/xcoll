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
    colldb=xc.load_SixTrack_colldb('colldb/lhc_tcp_b1.dat', emit=3.5e-6),
    _buffer=buffer
    )

# Install collimators in line as black absorbers
coll_manager.install_black_absorbers(verbose=True)

# Build the tracker
tracker = coll_manager.build_tracker()

# Align the collimators
coll_manager.align_collimators_to('front')

# Set the collimator openings based on the colldb,
# or manually override with the option gaps={collname: gap}
coll_manager.set_openings()

# Create initial particles, starting at the vertical primary
n_sigmas = 10
n_part = 50000
x_norm = np.random.uniform(-n_sigmas, n_sigmas, n_part)
y_norm = np.random.uniform(-n_sigmas, n_sigmas, n_part)
part = xp.build_particles(tracker=tracker, x_norm=x_norm, y_norm=y_norm,
                          scale_with_transverse_norm_emitt=(3.5e-6, 3.5e-6),
                          at_element='tcp.d6l7.b1',
                          match_at_s=coll_manager.s_match['tcp.d6l7.b1']
                         )

# Track
tracker.track(part, num_turns=10)

# The survival flags are sorted as surviving particles first,
# hence we need to 'unsort' them using their IDs
surv = part.state.copy()
surv[part.particle_id] = part.state

# Plot the surviving particles as green
plt.figure(1,figsize=(12,12))
plt.plot(x_norm, y_norm, '.', color='red')
plt.plot(x_norm[surv>0], y_norm[surv>0], '.', color='green')
plt.axis('equal')
plt.axis([n_sigmas, -n_sigmas, -n_sigmas, n_sigmas])
plt.show()

