# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2026.                 #
# ######################################### #

import numpy as np
import pytest
import xobjects as xo
import xtrack as xt
import xcoll as xc
from xobjects.test_helpers import for_all_test_contexts, skip_if_forbid_compile


class DechannellingLengthTest(xt.BeamElement):
    # Test-only kernel: depend on the real Everest implementation rather than
    # duplicating its formula in a new C function.
    _xofields = {'material': xc.Material}
    _depends_on = [xc.EverestCrystal]
    allow_no_prebuilt_kernel = True
    allow_track = True
    _extra_c_sources = ["""
        #include "xcoll/scattering_routines/everest/crystal_parameters.h"
        GPUFUN
        void DechannellingLengthTest_track_local_particle(
                DechannellingLengthTestData el, LocalParticle* part0) {
            START_PER_PARTICLE_BLOCK(part0, part);
                EverestData_ everest;
                double pc = LocalParticle_get_p0c(part) * 1.e-9;
                double coefficient = calculate_dechannelling_length(
                    &everest, DechannellingLengthTestData_getp_material(el), pc);
                LocalParticle_set_x(part, coefficient);
            END_PER_PARTICLE_BLOCK;
        }
    """]


@pytest.mark.gpu
@pytest.mark.everest
@for_all_test_contexts(excluding=('ContextPyopencl',))
def test_dechannelling_length_units(test_context):
    skip_if_forbid_compile()
    # Set I explicitly so this regression test is independent of the material
    # database's choice of excitation energy for silicon.
    material = xc.Material(Z=14, A=28.085, density=2.3296, state='solid',
                           excitation_energy=173., _context=test_context)
    probe = DechannellingLengthTest(material=material, _context=test_context)
    # pc is the momentum argument in GeV, not the energy in MeV.
    pc = np.array([1., 10., 450., 7000.])
    # Each particle has its own reference momentum. No scattering or random
    # sampling is involved; x is simply a test output buffer.
    part = xt.Particles(p0c=pc*1e9, x=np.zeros(len(pc)), _context=test_context)
    probe.track(part)
    actual = test_context.nparray_from_context_array(part.x)

    # Independent calculation with masses and I consistently in GeV.
    # The length prefactor is therefore directly in m/GeV; there is no
    # additional conversion in the logarithm (nor in the prefactor).
    electron_mass = 0.51099890e-3
    proton_mass = 0.9382720813
    excitation_energy = 173e-9
    gamma = np.hypot(pc, proton_mass)/proton_mass
    screening = 0.194e-10
    plane_distance = 1.920e-10
    classical_radius = 2.817940285e-15
    prefactor = (256/(9*np.pi**2) * screening*plane_distance
                 / (classical_radius*electron_mass))
    expected = prefactor/(np.log(2*electron_mass*gamma/excitation_energy)-1)
    np.testing.assert_allclose(actual, expected, rtol=2e-7, atol=0.)

    # Sanity check: the old /1000 error in the logarithm is large at every
    # tested momentum. In particular it cannot pass by tolerance at 7 TeV.
    old = prefactor/(np.log(2*electron_mass*gamma/excitation_energy/1000)-1)
    assert np.all(np.abs(actual-old) > 0.1*np.abs(expected))
