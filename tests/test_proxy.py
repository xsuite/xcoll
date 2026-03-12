# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2026.                 #
# ######################################### #

import pytest
from xcoll import Proxy


class Target:
    def __init__(self):
        self.value = 3
        self._prop = 8
        self._secret = 10

    @property
    def prop(self):
        return self._prop

    @prop.setter
    def prop(self, value):
        self._prop = value

    @prop.deleter
    def prop(self):
        self._prop = 8

    def multiply(self, x):
        return self.value * x

    def __len__(self):
        return 7


class LoggingProxy(Proxy[Target]):
    __slots__ = ("_enabled",)

    def __init__(self, target, enabled=True):
        super().__init__(target)
        self._enabled = enabled

    @property
    def enabled(self):
        return self._enabled

    @enabled.setter
    def enabled(self, value):
        self._enabled = bool(value)

    def _on_get(self, name, value):
        if self.enabled:
            print(f"GET {name} -> {value!r}")
        return value

    def _on_set(self, name, value):
        if self.enabled:
            print(f"SET {name} = {value!r}")
        return value

    def _on_del(self, name):
        if self.enabled:
            print(f"DEL {name}")
        return None

    def multiply(self, x):
        if self.enabled:
            print(f"CALL multiply({x!r})")
        return self._target.multiply(x)

    def __len__(self):
        if self.enabled:
            print("CALL __len__()")
        return len(self._target)

def test_proxy_delegates_reads_writes_and_deletes():
    t = Target()
    p = Proxy(t)

    assert p.value == 3
    assert p.prop == 8
    assert p._secret == 10
    assert p.multiply(5) == 15
    assert len(p) == 7

    p.value = 4
    assert p.value == 4
    assert t.value == 4

    p.prop = 11
    assert p.prop == 11
    assert t.prop == 11

    del p.prop
    assert p.prop == 8
    assert t.prop == 8


def test_custom_proxy(capsys):
    t = Target()
    p = LoggingProxy(t)

    assert p.enabled is True
    p.enabled = False
    assert p.enabled is False

    # Target should not receive this attribute
    assert not hasattr(t, "enabled")

    p.enabled = True

    # Proxy should be transparent
    assert p._secret == 10
    p._secret = 99
    assert p._secret == 99
    assert t._secret == 99
    capsys.readouterr()

    # Verify side effects of the proxy
    assert p.value == 3
    assert capsys.readouterr().out == "GET value -> 3\n"

    p.value = 4
    assert capsys.readouterr().out == "SET value = 4\n"

    assert p.prop == 8
    assert capsys.readouterr().out == "GET prop -> 8\n"

    p.prop = -99
    assert capsys.readouterr().out == "SET prop = -99\n"

    del p.prop
    assert capsys.readouterr().out == "DEL prop\n"

    assert p.prop == 8
    assert capsys.readouterr().out == "GET prop -> 8\n"

    assert p.multiply(5) == 20
    assert capsys.readouterr().out == "CALL multiply(5)\n"

    assert len(p) == 7
    assert capsys.readouterr().out == "CALL __len__()\n"

    p.enabled = False

    assert p.value == 4
    assert capsys.readouterr().out == ""

    p.value = 78
    assert capsys.readouterr().out == ""

    assert p.multiply(5) == 390
    assert capsys.readouterr().out == ""

    assert len(p) == 7
    assert capsys.readouterr().out == ""


def test_target_cannot_be_overwritten():
    p = Proxy(Target())

    with pytest.raises(AttributeError, match="_target"):
        p._target = Target()

    with pytest.raises(AttributeError, match="_target"):
        del p._target
