import numpy as np

import xcoll as xc
import xpart as xp

coll = xc.Collimator(inactive_length_at_start=3.,
                     active_length = 2, inactive_length_at_end=1,
                     n_slices=10, angle=20,
                     a_min=-0.5e-2, a_max=0.5e-2, b_min=-2e-2, b_max=2e-2)

p = xp.Particles(p0c=7e12, x=[1e-3, 2e-3], y=[3e-3, 4e-3])
coll.track(p)
assert np.allclose(p.s, 6, atol=1e-14, rtol=0)
assert np.allclose(p.x, [1e-3, 2e-3], atol=1e-14, rtol=0)
assert np.allclose(p.y, [3e-3, 4e-3], atol=1e-14, rtol=0)

coll = xc.Collimator(inactive_length_at_start=1.,
                     active_length = 3, inactive_length_at_end=1,
                     n_slices=100, angle=90,
                     a_min=-0.5e-2, a_max=0.5e-2, b_min=-2e-2, b_max=2e-2)
p = xp.Particles(p0c=7e12, py=[-0.5e-2/3, 0.5e-2/3])
coll.track(p)
assert np.all(p.state==-333)
assert np.allclose(p.s, 3, atol=3./100, rtol=0)
assert np.allclose(p.y[p.particle_id], [-5e-3, 5e-3], atol=1e-4, rtol=0)
assert np.allclose(p.x[p.particle_id], 0, atol=1e-10, rtol=0)


coll = xc.Collimator(inactive_length_at_start=3.,
                     active_length = 2, inactive_length_at_end=1,
                     n_slices=10, angle=20, dx=2e-2, dy=1e-2,
                     a_min=-0.5e-2, a_max=0.5e-2, b_min=-2e-2, b_max=2e-2)
part_gen_range = 3e-2
n_part = 10000
part = xp.Particles(
        p0c=6500e9,
        x=np.random.uniform(-part_gen_range, part_gen_range, n_part),
        px = np.zeros(n_part),
        y=np.random.uniform(-part_gen_range, part_gen_range, n_part),
        py = np.zeros(n_part),
        )
coll.track(part)

import matplotlib.pyplot as plt
plt.close('all')
plt.figure(1)
plt.plot(part.x, part.y, '.', color='red')
plt.plot(part.x[part.state>0], part.y[part.state>0], '.', color='green')
plt.plot(coll.dx, coll.dy, 'ok')
plt.axis('equal')
plt.show()

