# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2026.                 #
# ######################################### #

from xcoll.xaux import FsPath
from xcoll.package_env import BaseInterface


class DummyInterface(BaseInterface):
    _paths = {"required": 0}
    _optional_paths = {"optional": 0}
    _read_only_paths = {"readonly": 0}

    @property
    def compiled(self):
        return True

@pytest.fixture
def interface(tmp_path, monkeypatch):
    monkeypatch.setattr(
        DummyInterface, "_config_dir", FsPath(tmp_path / "config")
    )
    monkeypatch.setattr(
        DummyInterface, "_data_dir", FsPath(tmp_path / "data")
    )
    monkeypatch.setattr(
        DummyInterface, "_lib_dir", FsPath(tmp_path / "lib")
    )

    return DummyInterface()
