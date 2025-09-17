# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2025                  #
# ######################################### #

import os
import sys
from contextlib import contextmanager


@contextmanager
def pin_python_stdio():
    ORIG_STDOUT_FD = os.dup(1)
    ORIG_STDERR_FD = os.dup(2)
    try:
        yield

    finally:
        # flush the duplicates
        try: sys.stdout.flush()
        except Exception: pass
        try: sys.stderr.flush()
        except Exception: pass

        # restore the originals *before* closing the duplicates
        sys.stdout = os.fdopen(os.dup(ORIG_STDOUT_FD), "w", buffering=1)
        sys.stderr = os.fdopen(os.dup(ORIG_STDERR_FD), "w", buffering=1)

        # now it's safe to close the duplicates
        try: ORIG_STDOUT_FD.close()
        except Exception: pass
        try:
            if ORIG_STDERR_FD is not ORIG_STDOUT_FD:
                ORIG_STDERR_FD.close()
        except Exception: pass
