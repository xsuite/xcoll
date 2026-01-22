import numpy as np
from pathlib import Path

import xobjects as xo
import xtrack as xt
import xpart as xp
import xcoll as xc


path = Path(__file__).parent


# ============================================
# From line
# ============================================

# Get line and collimators
env = xt.load(path / 'machines' / 'lhc_run3_b1.json')
line = env['lhcb1']

colldb = xc.CollimatorDatabase.from_yaml(path / 'colldbs' / 'lhc_run3.yaml', beam=1)
colldb.install_everest_collimators(verbose=True, line=line)
df_with_coll = line.check_aperture()
assert not np.any(df_with_coll.has_aperture_problem)

# Start interaction record
impacts = xc.InteractionRecord.start(line=line)

# Build tracker, assign optics and generate particles 
line.build_tracker()
line.collimators.assign_optics()
part = line['tcp.d6l7.b1'].generate_pencil(5000)

# This is not needed, but is done here so that we can track with 12 treads.
line.discard_tracker()
line.build_tracker(_context=xo.ContextCpu(omp_num_threads=12))

# Track
line.scattering.enable()
line.track(part, num_turns=20, time=True, with_progress=1)
line.scattering.disable()
line.discard_tracker()
impacts.stop()

df = impacts.to_pandas()
df.to_csv('results/impacts_line.csv', index=False)

# ============================================
# With collimator
# ============================================
coll = xc.EverestCollimator(length=0.6, jaw=0.0013, material=xc.materials.MolybdenumGraphite, emittance=3.5e-6)

num_part = int(5000)
x_init   = np.random.normal(loc=1.5e-3, scale=75.e-6, size=num_part)
px_init  = np.random.uniform(low=-50.e-6, high=250.e-6, size=num_part)
y_init   = np.random.normal(loc=0., scale=1e-3, size=num_part)
py_init  = np.random.normal(loc=0., scale=5.e-6, size=num_part)
part = xp.Particles(x=x_init, px=px_init, y=y_init, py=py_init, delta=0, p0c=4e11)

impacts_coll = xc.InteractionRecord.start(elements=[coll], names='TPCH')

coll.track(part)
part.sort(interleave_lost_particles=True)

df = impacts_coll.to_pandas()
df[df.interaction_type == 'Enter Jaw L'].to_csv('results/impacts_coll_enter_jaw_L.csv', index=False)

# ============================================
# With crystal
# ============================================
coll_cry = xc.EverestCrystal(length=0.002, material=xc.materials.Silicon, bending_angle=149e-6,
                         width=0.002, height=0.05, side='+', lattice='strip', jaw=0.001)

num_part = int(5000)
x_init   = np.random.normal(loc=1.5e-3, scale=75.e-6, size=num_part)
px_init  = np.random.uniform(low=-50.e-6, high=250.e-6, size=num_part)
y_init   = np.random.normal(loc=0., scale=1e-3, size=num_part)
py_init  = np.random.normal(loc=0., scale=5.e-6, size=num_part)
part = xp.Particles(x=x_init, px=px_init, y=y_init, py=py_init, delta=0, p0c=4e11)

impacts_crystal = xc.InteractionRecord.start(elements=[coll_cry], names='TPCH')
coll_cry.track(part)
part.sort(interleave_lost_particles=True)

impacts_crystal.to_pandas()
df_crystal = impacts_crystal.interactions_per_collimator()
df_crystal.to_csv('results/impacts_crystal_interactions.csv', index=False)
