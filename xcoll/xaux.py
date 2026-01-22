# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

# Tenporary file that defines xaux tools - to be used until xaux is in main Xsuite release cycle

import os
import base64
import shutil
import inspect
import functools
import numpy as np
from pathlib import PosixPath


class FsPath(PosixPath):
    def copy_to(self, other, **kwargs):
        if self.is_dir():
            shutil.copytree(self, other / self.name, dirs_exist_ok=True)
        else:
            shutil.copy(self, other)
    def move_to(self, other, **kwargs):
        shutil.move(self, other)
    def rmtree(self, *args, **kwargs):
        shutil.rmtree(self, *args, **kwargs)
    def __eq__(self, other):
        try:
            other = FsPath(other).expanduser().resolve()
        except:
            return False
        self = self.expanduser().resolve()
        return self.as_posix() == other.as_posix()


def ranID(*, length=12, size=1, only_alphanumeric=False):
    """Base64 encoded random ID.
    Args:
        length (int): Length of the ID string, rounded up to
            the closest multiple of 4. Default 12.
        size (int): Number of random IDs to generate.
            Default 1.
        only_alphanumeric (bool): If True, only alphanumeric
            characters are used. Default False.
    Returns:
        str: Random ID string.
    """
    if length < 1:
        raise ValueError("Length must be greater than 0!")
    if size < 1:
        raise ValueError("Size must be greater than 0!")
    if size > 1:
        return [ranID(length=length, only_alphanumeric=only_alphanumeric)
                for _ in range(size)]
    length = int(np.ceil(length/4))
    if only_alphanumeric:
        ran = ''
        for _ in range(length):
            while True:
                this_ran = ranID(length=4, only_alphanumeric=False)
                if this_ran.isalnum():
                    break
            ran += this_ran
        return ran
    else:
        random_bytes = os.urandom(3*length)
        return base64.urlsafe_b64encode(random_bytes).decode('utf-8')



def count_required_arguments(func):
    i = 0
    sig = inspect.signature(func)
    for param in sig.parameters.values():
        if (param.kind == inspect.Parameter.POSITIONAL_ONLY \
        or param.kind == inspect.Parameter.POSITIONAL_OR_KEYWORD) \
        and param.default == inspect.Parameter.empty:
            i += 1
    return i
