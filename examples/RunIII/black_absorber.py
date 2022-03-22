import numpy as np

import matplotlib.pyplot as plt

import xtrack as xt
import xpart as xp
import xcoll as xc

import sixtracktools as st




# Import Run III lattice
line = xt.Line.from_sixinput(st.SixInput('./'))

# Attach reference particle (a proton at 6.8 TeV)
line.particle_ref = xp.Particles(mass0 = xp.PROTON_MASS_EV, p0c=6.8e12)

# Switch on RF (needed to twiss)
line['acsca.d5l4.b1'].voltage = 16e6
line['acsca.d5l4.b1'].frequency = 1e6

# Initialise collmanager
coll_manager = xc.CollimatorManager(
    line=line,
    colldb=xc.load_SixTrack_colldb('CollDB-RunIII_B1.dat', 3.5e-6)
    )

# Install collimators in line as black absorbers
coll_manager.install_black_absorbers()

# Build the tracker
tracker = line.build_tracker()

# Set the collimator openings based on the colldb,
# or manually override with the option gaps={collname: gap}
coll_manager.set_openings()


# Create initial particles
n_sigmas = 10
n_part = 50000
x_norm = np.random.uniform(-n_sigmas, n_sigmas, n_part)
y_norm = np.random.uniform(-n_sigmas, n_sigmas, n_part)
part = xp.build_particles(tracker=tracker, x_norm=x_norm, y_norm=y_norm,
                          scale_with_transverse_norm_emitt=(3.5e-6, 3.5e-6),
                          at_element='tcp.d6l7.b1',
                          match_at_s=coll_manager.s_start_active['tcp.d6l7.b1']
                         )

# Track
tracker.track(part, num_turns=5)


# The survival flags are sorted as surviving particles first,
# hence we need to 'unsort' them using their IDs
surv = part.state.copy()
surv[part.particle_id] = part.state

# Plot the surviving particles as green
plt.figure(1,figsize=(12,12))
plt.plot(x_norm, y_norm, '.', color='red')
plt.plot(x_norm[surv>0], y_norm[surv>0], '.', color='green')
plt.axis('equal')
plt.show()
