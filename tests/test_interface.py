# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2026.                 #
# ######################################### #

import os
import sys
import json
import pytest
from subprocess import CompletedProcess

from xcoll.xaux import FsPath, track_construction
from xcoll.package_env import BaseInterface
import xcoll.package_env as package_env


class DummyInterface(BaseInterface):
    _paths = {
        "required": 0,
    }
    _optional_paths = {
        "optional": 0,
    }
    _read_only_paths = {
        "readonly": 0,
    }
    @property
    def compiled(self):
        return True


class DummyUncompiledInterface(DummyInterface):
    @property
    def compiled(self):
        return False


@track_construction
class DummyPostInitInterface(BaseInterface):
    _paths = {}
    _optional_paths = {}
    _read_only_paths = {
        "detected": 0,
    }

    def __init__(self, detected):
        super().__init__()

        # This assignment happens after BaseInterface.__init__(),
        # but is still inside the subclass constructor.
        self._detected = detected

        # One deliberate save at the end.
        self.save()

    @property
    def compiled(self):
        return True


@pytest.fixture
def interface_dirs(tmp_path, monkeypatch):
    config_dir = FsPath(tmp_path / "config")
    data_dir = FsPath(tmp_path / "data")
    lib_dir = FsPath(tmp_path / "lib")
    monkeypatch.setattr(DummyInterface, "_config_dir", config_dir)
    monkeypatch.setattr(DummyInterface, "_data_dir", data_dir)
    monkeypatch.setattr(DummyInterface, "_lib_dir", lib_dir)
    # The uncompiled subclass inherits these, but set them explicitly
    # to make the test independent of inheritance details.
    monkeypatch.setattr(DummyUncompiledInterface, "_config_dir", config_dir)
    monkeypatch.setattr(DummyUncompiledInterface, "_data_dir", data_dir)
    monkeypatch.setattr(DummyUncompiledInterface, "_lib_dir", lib_dir)

    monkeypatch.setattr(DummyPostInitInterface, "_config_dir", config_dir)
    monkeypatch.setattr(DummyPostInitInterface, "_data_dir", data_dir)
    monkeypatch.setattr(DummyPostInitInterface, "_lib_dir", lib_dir)
    old_sys_path = sys.path.copy()
    old_environ = os.environ.copy()
    yield {
        "config": config_dir,
        "data": data_dir,
        "lib": lib_dir,
    }
    sys.path[:] = old_sys_path
    os.environ.clear()
    os.environ.update(old_environ)


def test_interface_initialisation(interface_dirs):
    interface = DummyInterface()
    assert interface.config_dir == interface_dirs["config"]
    assert interface.data_dir == interface_dirs["data"]
    assert interface.lib_dir == interface_dirs["lib"]
    assert interface.config_dir.exists()
    assert interface.data_dir.exists()
    assert interface.lib_dir.exists()
    assert interface.required is None
    assert interface.optional is None
    assert interface.readonly is None
    assert interface.config_file.exists()


def test_constructor_does_not_save_or_bruteforce(interface_dirs, monkeypatch):
    required = FsPath(interface_dirs["data"] / "required")
    optional = FsPath(interface_dirs["data"] / "optional")
    readonly = FsPath(interface_dirs["data"] / "readonly")
    required.mkdir(parents=True)
    optional.mkdir(parents=True)
    readonly.mkdir(parents=True)
    interface_dirs["config"].mkdir(parents=True, exist_ok=True)
    config_file = interface_dirs["config"] / "dummy.config.json"
    with open(config_file, "w") as fid:
        json.dump(
            {
                "paths": {"required": required.as_posix()},
                "optional_paths": {"optional": optional.as_posix()},
                "read_only_paths": {"readonly": readonly.as_posix()},
            },
            fid,
        )
    calls = []
    def save(self):
        calls.append("save")
    def brute_force_path(self, path):
        calls.append(("brute_force", path))
    monkeypatch.setattr(DummyInterface, "save", save)
    monkeypatch.setattr(DummyInterface, "brute_force_path", brute_force_path)
    interface = DummyInterface()
    assert calls == []
    assert interface.required == required
    assert interface.optional == optional
    assert interface.readonly == readonly


def test_setting_path_validates_and_saves(interface_dirs, monkeypatch):
    interface = DummyInterface()
    path = FsPath(interface_dirs["data"] / "required")
    path.mkdir()
    calls = []
    def brute_force_path(value):
        calls.append(value)
    interface.brute_force_path = brute_force_path
    interface.required = path
    assert interface.required == path
    assert calls == [path]
    with open(interface.config_file) as fid:
        config = json.load(fid)
    assert config["paths"]["required"] == path.as_posix()


def test_optional_path(interface_dirs):
    interface = DummyInterface()
    path = FsPath(interface_dirs["data"] / "optional")
    path.mkdir()
    interface.brute_force_path = lambda value: None
    interface.optional = path
    assert interface.optional == path
    with open(interface.config_file) as fid:
        config = json.load(fid)
    assert config["optional_paths"]["optional"] == path.as_posix()


def test_read_only_path(interface_dirs):
    interface = DummyInterface()
    path = FsPath(interface_dirs["data"] / "readonly")
    path.mkdir()
    with pytest.raises(AttributeError, match="read-only"):
        interface.readonly = path
    # Internal assignment is allowed.
    interface._readonly = path
    assert interface.readonly == path
    with open(interface.config_file) as fid:
        config = json.load(fid)
    assert config["read_only_paths"]["readonly"] == path.as_posix()


def test_delete_path(interface_dirs):
    interface = DummyInterface()
    path = FsPath(interface_dirs["data"] / "required")
    path.mkdir()
    interface.brute_force_path = lambda value: None
    interface.required = path
    assert interface.required == path
    del interface.required
    assert interface.required is None
    with open(interface.config_file) as fid:
        config = json.load(fid)
    assert config["paths"]["required"] is None


def test_save_load_roundtrip(interface_dirs):
    required = FsPath(interface_dirs["data"] / "required")
    optional = FsPath(interface_dirs["data"] / "optional")
    readonly = FsPath(interface_dirs["data"] / "readonly")
    required.mkdir(parents=True)
    optional.mkdir(parents=True)
    readonly.mkdir(parents=True)
    first = DummyInterface()
    first.brute_force_path = lambda value: None
    first.required = required
    first.optional = optional
    first._readonly = readonly
    second = DummyInterface()
    assert second.required == required
    assert second.optional == optional
    assert second.readonly == readonly


def test_initialised_and_ready(interface_dirs):
    interface = DummyInterface()
    assert not interface.initialised
    assert not interface.ready
    required = FsPath(interface_dirs["data"] / "required")
    required.mkdir()
    interface.brute_force_path = lambda value: None
    interface.required = required
    assert interface.initialised
    assert interface.ready


def test_assert_environment_ready(interface_dirs):
    interface = DummyInterface()
    with pytest.raises(RuntimeError, match="not initialised"):
        interface.assert_environment_ready()
    required = FsPath(interface_dirs["data"] / "required")
    required.mkdir()
    interface.brute_force_path = lambda value: None
    interface.required = required
    interface.assert_environment_ready()


def test_assert_environment_not_compiled(interface_dirs):
    interface = DummyUncompiledInterface()
    required = FsPath(interface_dirs["data"] / "required")
    required.mkdir()
    interface.brute_force_path = lambda value: None
    interface.required = required
    assert interface.initialised
    assert not interface.ready
    with pytest.raises(RuntimeError, match="not compiled"):
        interface.assert_environment_ready()


def test_temp_dir(interface_dirs, tmp_path):
    interface = DummyInterface()
    first = interface.temp_dir
    assert first.exists()
    parent = FsPath(tmp_path / "custom_tmp")
    parent.mkdir()
    interface.temp_dir = parent
    assert not first.exists()
    second = interface.temp_dir
    assert second.exists()
    assert second.parent == parent
    del interface.temp_dir
    assert not second.exists()
    assert interface._temp_dir is None


def test_preserve_environment(interface_dirs):
    interface = DummyInterface()
    sys_path_object = sys.path
    environ_object = os.environ
    old_sys_path = sys.path.copy()
    old_environ = os.environ.copy()
    with interface.preserve_environment():
        sys.path.append("/this/should/disappear")
        os.environ["XCOLL_INTERFACE_TEST"] = "temporary"

    # The special global objects themselves should not be replaced.
    assert sys.path is sys_path_object
    assert os.environ is environ_object
    assert sys.path == old_sys_path
    assert dict(os.environ) == old_environ


def test_lib_path_is_not_duplicated(interface_dirs):
    lib = interface_dirs["lib"].as_posix()
    before = sys.path.count(lib)
    DummyInterface()
    DummyInterface()
    assert sys.path.count(lib) == before + 1


def test_post_super_constructor_assignment_is_shielded(
    interface_dirs,
    monkeypatch,
):
    detected = interface_dirs["data"] / "detected"
    detected.mkdir(parents=True)
    calls = []
    def save(self):
        calls.append(self._detected)
    monkeypatch.setattr(DummyPostInitInterface, "save", save)
    interface = DummyPostInitInterface(detected)
    assert interface.detected == detected
    # _detected assignment itself must not trigger BaseInterface.save().
    # The only call is the explicit save() at the end of __init__.
    assert calls == [detected]


def _completed_process(
    args,
    *,
    returncode=0,
    stdout="",
    stderr="",
):
    return CompletedProcess(
        args=args,
        returncode=returncode,
        stdout=stdout.encode(),
        stderr=stderr.encode(),
    )


def test_assert_installed_success_and_cache(interface_dirs, monkeypatch):
    interface = DummyInterface()
    calls = []
    def fake_run(args, stdout=None, stderr=None):
        calls.append(args)
        if args == ["which", "dummy-program"]:
            return _completed_process(args, stdout="/usr/bin/dummy-program\n")
        if args == ["dummy-program", "--version"]:
            return _completed_process(args, stdout="1.2.3\n")
        raise AssertionError(f"Unexpected command: {args}")
    monkeypatch.setattr(package_env, "run", fake_run)
    assert interface.assert_installed("dummy-program") == ("/usr/bin/dummy-program", "1.2.3")
    # Cached: no new subprocess calls.
    assert interface.assert_installed("dummy-program") == ("/usr/bin/dummy-program", "1.2.3")
    assert calls == [["which", "dummy-program"], ["dummy-program", "--version"]]


def test_assert_installed_which_failure(interface_dirs, monkeypatch):
    interface = DummyInterface()
    def fake_run(args, stdout=None, stderr=None):
        return _completed_process(args, returncode=1, stderr="not found")
    monkeypatch.setattr(package_env, "run", fake_run)
    with pytest.raises(RuntimeError, match=r"which missing-program"):
        interface.assert_installed("missing-program")


def test_assert_installed_which_exception(interface_dirs, monkeypatch):
    interface = DummyInterface()
    def fake_run(*args, **kwargs):
        raise OSError("boom")
    monkeypatch.setattr(package_env, "run", fake_run)
    with pytest.raises(RuntimeError, match=r"which dummy-program") as exc_info:
        interface.assert_installed("dummy-program")
    assert isinstance(exc_info.value.__cause__, OSError)


def test_assert_installed_version_failure(interface_dirs, monkeypatch):
    interface = DummyInterface()
    def fake_run(args, stdout=None, stderr=None):
        if args == ["which", "dummy-program"]:
            return _completed_process(args, stdout="/usr/bin/dummy-program\n")
        return _completed_process(args, returncode=2, stderr="bad version command")
    monkeypatch.setattr(package_env, "run", fake_run)
    with pytest.raises(RuntimeError, match=r"dummy-program --version"):
        interface.assert_installed("dummy-program")


@pytest.mark.parametrize(
    "method_name,env_var,program",
    [
        (
            "assert_gcc_installed",
            "CC",
            "test-gcc",
        ),
        (
            "assert_gxx_installed",
            "CXX",
            "test-g++",
        ),
        (
            "assert_gfortran_installed",
            "FC",
            "test-gfortran",
        ),
    ],
)
def test_old_compiler_check_is_repeatable(interface_dirs, monkeypatch,
                                          method_name, env_var, program):
    interface = DummyInterface()
    monkeypatch.setenv(env_var, program)
    calls = []
    def fake_run(args, stdout=None, stderr=None):
        calls.append(args)
        if args == ["which", program]:
            return _completed_process(args, stdout=f"/usr/bin/{program}\n")
        if args == [program, "-dumpversion"]:
            return _completed_process(args, stdout="8.5.0\n")
        raise AssertionError(f"Unexpected command: {args}")
    monkeypatch.setattr(package_env, "run", fake_run)
    method = getattr(interface, method_name)
    for _ in range(2):
        with pytest.raises(RuntimeError, match="9 or higher"):
            method(minimum_version=9)
    # assert_installed() result remains cached and valid.
    assert calls == [["which", program], [program, "-dumpversion"]]
