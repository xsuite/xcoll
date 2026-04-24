# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2026.                 #
# ######################################### #

import pytest
import shutil
from pathlib import Path
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

@pytest.fixture
def register_cleanup():
    paths = []

    def add(path):
        paths.append(path)

    yield add

    for p in paths:
        if Path(p).is_dir():
            shutil.rmtree(p, ignore_errors=True)
        elif Path(p).is_file():
            try:
                Path(p).unlink()
            except Exception:
                pass
