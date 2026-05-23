# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2026.                 #
# ######################################### #

import pytest
from _common_api import check_skip


def pytest_runtest_setup(item):
    # If any engine marker applies to this test item, run the skip check.
    if item.get_closest_marker("fluka") is not None:
        check_skip("fluka")

    if item.get_closest_marker("geant4") is not None:
        check_skip("geant4")


def pytest_collection_modifyitems(config, items):
    running_xdist = hasattr(config, "workerinput") or config.getoption("-n") not in (None, 0)
    if not running_xdist:
        return

    for item in items:
        if item.get_closest_marker("serial"):
            item.add_marker(pytest.mark.skip("Serial test cannot run under xdist"))
