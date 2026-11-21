# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2026.                 #
# ######################################### #
"""CPU-built Everest elements moved to CUDA by Line.build_tracker()."""
import subprocess
import sys
import pytest
from xobjects.test_helpers import skip_if_forbid_compile


# Run each element in a fresh Python process: a previous crystal must not have
# raised the CUDA stack limit, and a CUDA fault must not poison the test worker.
SCRIPT = r'''
import sys
import numpy as np
import cupy as cp
import xobjects as xo
import xtrack as xt
import xcoll as xc
from xcoll.scattering_routines.everest import CRYSTAL_STACK_LIMIT_BYTES

case = sys.argv[1]
cpu = xo.ContextCpu()
gpu = xo.ContextCupy()
cp.cuda.runtime.deviceSetLimit(0x00, 1024)
assert cp.cuda.runtime.deviceGetLimit(0x00) == 1024
kwargs = dict(length=0.002, _context=cpu)
if case == 'block':
    el = xc.EverestBlock(material=xc.materials.CarbonFibreComposite, **kwargs)
elif case == 'collimator':
    el = xc.EverestCollimator(material=xc.materials.CarbonFibreComposite,
                              jaw=[0.001, -0.001], **kwargs)
else:
    el = xc.EverestCrystal(material=xc.materials.Silicon, jaw=0.001,
                          side='+', bending_angle=149e-6, width=0.002,
                          height=0.05, **kwargs)
line = xt.Line(elements=[xt.Drift(length=0.1), el, xt.Drift(length=0.1)],
               element_names=['before','everest','after'],
               particle_ref=xt.Particles(p0c=400e9))
assert line.get('everest')._context is cpu
assert cp.cuda.runtime.deviceGetLimit(0x00) == 1024
line.build_tracker(_context=gpu)
assert line.get('everest')._context is gpu
if case == 'crystal':
    assert cp.cuda.runtime.deviceGetLimit(0x00) >= CRYSTAL_STACK_LIMIT_BYTES
n = 128
p = line.build_particles(x=np.full(n, 0.0015), px=np.zeros(n),
                         y=np.zeros(n), py=np.zeros(n))
assert p._context is gpu
p._init_random_number_generator(seeds=13579+np.arange(n, dtype=np.uint32))
line.track(p, num_turns=1)
cp.cuda.runtime.deviceSynchronize()
state = gpu.nparray_from_context_array(p.state)
assert np.count_nonzero(state > 0) > 0
assert np.all(np.isfinite(gpu.nparray_from_context_array(p.x)[state > 0]))
assert np.all(gpu.nparray_from_context_array(p.at_turn)[state > 0] == 1)
'''


@pytest.mark.gpu
@pytest.mark.everest
@pytest.mark.parametrize('case', ['block','collimator','crystal'])
def test_cpu_line_build_gpu_tracker(case):
    skip_if_forbid_compile()
    cp = pytest.importorskip('cupy')
    if cp.cuda.runtime.getDeviceCount() == 0:
        pytest.skip('No CUDA device available')
    result = subprocess.run([sys.executable,'-c',SCRIPT,case],
                            capture_output=True, text=True, timeout=600)
    assert result.returncode == 0, result.stdout+'\n'+result.stderr
