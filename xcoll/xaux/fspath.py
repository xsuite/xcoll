# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

import shutil
from pathlib import PosixPath, Path


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
    def __hash__(self):
        return hash(Path(self))
