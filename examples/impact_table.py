# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

import xobjects as xo
import xtrack as xt
import xpart as xp
import xcoll as xc

# ============================================
# With collimators
# ============================================

# Get line and collimators
line = xt.Line.from_json(xc._pkg_root / '..' / 'examples' / 'machines' / 'lhc_run3_b1.json')

coll_manager = xc.CollimatorDatabase.from_yaml(xc._pkg_root / '..' / 'examples' / 'colldb' / 'lhc_run3.yaml', beam=1)
coll_manager.install_everest_collimators(verbose=True, line=line)
df_with_coll = line.check_aperture()
assert not np.any(df_with_coll.has_aperture_problem)

# Start interaction record
impacts = xc.InteractionRecord.start(line)

# Build tracker, assign optics and generate particles 
line.build_tracker()
xc.assign_optics_to_collimators(line=line)
part = xc.generate_pencil_on_collimator(line, 'tcp.d6l7.b1', 50000)

# This is not needed, but is done here so that we can track with 12 treads.
line.discard_tracker()
line.build_tracker(_context=xo.ContextCpu(omp_num_threads=12))

# Track
xc.enable_scattering(line)
line.track(part, num_turns=20, time=True, with_progress=1)
xc.disable_scattering(line)
line.discard_tracker()

df = impacts.to_pandas()
df[(df.interaction_type == 'Enter Jaw L') & (df.ds == 0.0)]
df.to_csv('impacts.csv', index=False)

# ============================================
# With crystal
# ============================================
coll = xc.EverestCrystal(length=0.002, material=xc.materials.SiliconCrystal, bending_angle=149e-6,
                         width=0.002, height=0.05, side='+', lattice='strip', jaw=0.001)

num_part = int(50000)
x_init   = np.random.normal(loc=1.5e-3, scale=75.e-6, size=num_part)
px_init  = np.random.uniform(low=-50.e-6, high=250.e-6, size=num_part)
y_init   = np.random.normal(loc=0., scale=1e-3, size=num_part)
py_init  = np.random.normal(loc=0., scale=5.e-6, size=num_part)
part = xp.Particles(x=x_init, px=px_init, y=y_init, py=py_init, delta=0, p0c=4e11)

io_buffer = xt.new_io_buffer(capacity=int(2e7)) # 4-5 GB of memory
coll.record_scatterings = True
impacts_crystal = xt.start_internal_logging(elements=[coll], io_buffer=io_buffer, capacity=io_buffer.capacity)

coll.track(part)
part.sort(interleave_lost_particles=True)

impacts_crystal.to_pandas()
df_crystal = impacts_crystal.interactions_per_collimator()
df_crystal.to_csv('interactions_per_crystal.csv', index=False)