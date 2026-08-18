import numpy as np
from pathlib import Path

import xobjects as xo
import xtrack as xt
import xpart as xp
import xcoll as xc

path = Path(__file__).parent


context = xo.ContextCpu(omp_num_threads='auto')  # For CPU
# context = xo.ContextCupy()                     # For CUDA GPUs
# context = xo.ContextPyopencl()                 # For OpenCL GPUs


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
impacts = xc.InteractionRecord(line=line)

# Build tracker, assign optics and generate particles 
line.build_tracker()
line.xcoll.collimators.assign_optics()
part = line['tcp.d6l7.b1'].generate_pencil(5000)

# Move to a more efficient context for tracking
line.discard_tracker()
line.build_tracker(_context=context)

# Track
line.xcoll.scattering.enable()
line.track(part, num_turns=20, time=True, with_progress=1)
line.xcoll.scattering.disable()
line.discard_tracker()
impacts.stop()

df = impacts.to_pandas()
df.to_csv('results/impacts_line.csv', index=False)

# ============================================
# With collimator
# ============================================
coll = xc.EverestCollimator(length=0.6,
                            jaw=0.0013,
                            material=xc.materials.MolybdenumGraphite,
                            emittance=3.5e-6,
                            _context=context)

num_part = int(5000)
x_init   = np.random.normal(loc=1.5e-3, scale=75.e-6, size=num_part)
px_init  = np.random.uniform(low=-50.e-6, high=250.e-6, size=num_part)
y_init   = np.random.normal(loc=0., scale=1e-3, size=num_part)
py_init  = np.random.normal(loc=0., scale=5.e-6, size=num_part)
part = xp.Particles(x=x_init, px=px_init, y=y_init, py=py_init, delta=0,
                    p0c=4e11, _context=context)

impacts_coll = xc.InteractionRecord(elements=[coll], names='TPCH')

coll.track(part)
if not isinstance(context, xo.ContextCpu):
    part.move(_context=xo.ContextCpu())    # Not super fast
part.sort(interleave_lost_particles=True)  # Not super fast

df = impacts_coll.to_pandas()
df[df.interaction_type == 'Enter Jaw L'].to_csv('results/impacts_coll_enter_jaw_L.csv', index=False)

# ============================================
# With crystal
# ============================================
coll_cry = xc.EverestCrystal(length=0.002,
                             material=xc.materials.Silicon,
                             bending_angle=149e-6,
                             width=0.002,
                             height=0.05,
                             side='+',
                             lattice='strip',
                             jaw=0.001,
                             _context=context)

num_part = int(5000)
x_init   = np.random.normal(loc=1.5e-3, scale=75.e-6, size=num_part)
px_init  = np.random.uniform(low=-50.e-6, high=250.e-6, size=num_part)
y_init   = np.random.normal(loc=0., scale=1e-3, size=num_part)
py_init  = np.random.normal(loc=0., scale=5.e-6, size=num_part)
part = xp.Particles(x=x_init, px=px_init, y=y_init, py=py_init, delta=0, p0c=4e11, _context=context)

impacts_crystal = xc.InteractionRecord(elements=[coll_cry], names='TPCH')
coll_cry.track(part)
if not isinstance(context, xo.ContextCpu):
    part.move(_context=xo.ContextCpu())    # Not super fast
part.sort(interleave_lost_particles=True)  # Not super fast

impacts_crystal.to_pandas()
df_crystal = impacts_crystal.interactions_per_collimator()
df_crystal.to_csv('results/impacts_crystal_interactions.csv', index=False)
