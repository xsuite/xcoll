# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2026.                 #
# ######################################### #

import pytest

import xtrack as xt
import xcoll as xc


@pytest.mark.xcother
def test_environment_api_facade():
    env = xt.Environment()
    assert isinstance(env.xcoll, xc.XcollEnvironmentAPI)
    assert env.xcoll._env is env
