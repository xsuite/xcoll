# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2026.                 #
# ######################################### #

import os
from pathlib import Path
from contextlib import contextmanager


@contextmanager
def temporary_cwd(path):
    """Temporarily change the process working directory."""
    old_cwd = Path.cwd()

    if path is None:
        yield old_cwd
        return

    os.chdir(path)
    try:
        yield Path.cwd()
    finally:
        os.chdir(old_cwd)
