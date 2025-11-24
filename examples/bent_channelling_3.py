import numpy as np
import matplotlib.pyplot as plt
import time

import xtrack as xt
import xcoll as xc


L_tot   = 0.0004
# Maybe there is a maximum number of steps due to numerical instability? NO, just checked
n_steps = 1000
ds      = L_tot / n_steps
x0      = -0.95e-10
theta0  = 12e-6
#s_axis  = np.linspace(0.0, L_tot, n_steps+1)

# These are the classes and the respective methods, orders, versions:
CLASS_MAP = {
    # ---- Method 2, Version 1 ----
    ('M2','V1', 2): xc.BentChannellingDevM2V1o02,
    ('M2','V1', 4): xc.BentChannellingDevM2V1o04,
    ('M2','V1', 6): xc.BentChannellingDevM2V1o06,
    ('M2','V1', 8): xc.BentChannellingDevM2V1o08,
    ('M2','V1',10): xc.BentChannellingDevM2V1o10,
    ('M2','V1',12): xc.BentChannellingDevM2V1o12,

    # ---- Method 2, Version 2 ----
    ('M2','V2', 2): xc.BentChannellingDevM2V2o02,
    ('M2','V2', 4): xc.BentChannellingDevM2V2o04,
    ('M2','V2', 6): xc.BentChannellingDevM2V2o06,
    ('M2','V2', 8): xc.BentChannellingDevM2V2o08,
    ('M2','V2',10): xc.BentChannellingDevM2V2o10,
    ('M2','V2',12): xc.BentChannellingDevM2V2o12,

    # ---- Method 3, Version 1 ----
    ('M3','V1', 2): xc.BentChannellingDevM3V1o02,
    ('M3','V1', 4): xc.BentChannellingDevM3V1o04,
    ('M3','V1', 6): xc.BentChannellingDevM3V1o06,
    ('M3','V1', 8): xc.BentChannellingDevM3V1o08,
    ('M3','V1',10): xc.BentChannellingDevM3V1o10,
    ('M3','V1',12): xc.BentChannellingDevM3V1o12,

    # ---- Method 3, Version 2 ----
    ('M3','V2', 2): xc.BentChannellingDevM3V2o02,
    ('M3','V2', 4): xc.BentChannellingDevM3V2o04,
    ('M3','V2', 6): xc.BentChannellingDevM3V2o06,
    ('M3','V2', 8): xc.BentChannellingDevM3V2o08,
    ('M3','V2',10): xc.BentChannellingDevM3V2o10,
    ('M3','V2',12): xc.BentChannellingDevM3V2o12,

    # ---- Method 4, Version 1 ----
    ('M4','V1', 2): xc.BentChannellingDevM4V1o02,
    ('M4','V1', 4): xc.BentChannellingDevM4V1o04,
    ('M4','V1', 6): xc.BentChannellingDevM4V1o06,
    ('M4','V1', 8): xc.BentChannellingDevM4V1o08,
    ('M4','V1',10): xc.BentChannellingDevM4V1o10,
    ('M4','V1',12): xc.BentChannellingDevM4V1o12,

    # ---- Method 4, Version 2 ----
    ('M4','V2', 2): xc.BentChannellingDevM4V2o02,
    ('M4','V2', 4): xc.BentChannellingDevM4V2o04,
    ('M4','V2', 6): xc.BentChannellingDevM4V2o06,
    ('M4','V2', 8): xc.BentChannellingDevM4V2o08,
    ('M4','V2',10): xc.BentChannellingDevM4V2o10,
    ('M4','V2',12): xc.BentChannellingDevM4V2o12,
}



# picking the class 

method  = 'M3'
version = 'V2'
order   = 4
ChClass = CLASS_MAP[(method, version, order)]

# create element
 
ch = ChClass(length=ds)


# particle 
part = xt.Particles(x=x0, px=theta0)

# initialisation
xs  = np.zeros(n_steps + 1)
xps = np.zeros(n_steps + 1)
xs[0], xps[0] = x0, theta0

t0 = time.perf_counter()
for i in range(n_steps):
    ch.track(part)
    xs[i+1]  = part.x[0]
    xps[i+1] = part.px[0]     
elapsed = time.perf_counter() - t0
print(f"Execution time: {1000*elapsed:.3f} ms with {method}{version} order {order:02d}")

# results

s_axis = np.linspace(0, L_tot, n_steps + 1)*1e3
fig, ax = plt.subplots(1, 2, figsize=(10, 4))

ax[0].plot(s_axis, xs, label='x(s)')
ax[0].set_xlabel('s [mm]')
ax[0].set_ylabel('x [m]')
ax[0].set_title('Transverse position')

ax[1].plot(s_axis, xps, label="θ(s)")
ax[1].set_xlabel('s [mm]')
ax[1].set_ylabel("x' [rad]")
ax[1].set_title('Angle evolution')

for a in ax:
    a.grid(True)
    a.legend(loc='upper right')

plt.tight_layout()
plt.show()
