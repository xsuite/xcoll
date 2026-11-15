# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2026.                 #
# ######################################### #
from types import SimpleNamespace
import numpy as np
import pytest
import xobjects as xo
import xtrack as xt
import xcoll as xc


@pytest.mark.xcother
@pytest.mark.parametrize('enabled', [True, False])
@pytest.mark.parametrize('failure', ['constructor', 'refine', None])
def test_refinement_restores_time_dependent_vars(monkeypatch, enabled, failure):
    line = SimpleNamespace(_context=xo.ContextCpu(), enable_time_dependent_vars=enabled)
    part = SimpleNamespace(s=np.array([1.]), state=np.array([0]))
    lm = xc.LossMap()
    lm.interpolation = 0.2

    class Refinement:
        def __init__(self, refinement_line, **kwargs):
            assert refinement_line is line
            assert not line.enable_time_dependent_vars
            if failure == 'constructor':
                raise RuntimeError('refinement failed')

        def refine_loss_location(self, particles):
            assert not line.enable_time_dependent_vars
            if failure == 'refine':
                raise RuntimeError('refinement failed')

    monkeypatch.setattr(xt, 'LossLocationRefinement', Refinement)
    if failure is None:
        lm._interpolate(part, line, verbose=False)
    else:
        with pytest.raises(RuntimeError, match='refinement failed'):
            lm._interpolate(part, line, verbose=False)
    assert line.enable_time_dependent_vars is enabled
