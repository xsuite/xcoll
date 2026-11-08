# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2026.                #
# ######################################### #

import pytest

from xcoll.xaux import (
    track_construction,
    skip_if_being_constructed,
    super_if_being_constructed,
)


def test_track_construction():
    observations = []
    @track_construction
    class Example:
        def __init__(self):
            observations.append(self._being_constructed())
    obj = Example()
    assert observations == [True]
    assert not obj._being_constructed()


def test_track_construction_is_nesting_safe():
    observations = []

    @track_construction
    class Parent:
        def __init__(self):
            observations.append(("parent", self._being_constructed()))

    @track_construction
    class Child(Parent):
        def __init__(self):
            observations.append(("child-before", self._being_constructed()))
            super().__init__()
            observations.append(("child-after", self._being_constructed()))

    obj = Child()
    assert observations == [
        ("child-before", True),
        ("parent", True),
        ("child-after", True),
    ]
    assert not obj._being_constructed()


def test_track_construction_resets_after_exception():
    instances = []

    @track_construction
    class Broken:
        def __init__(self):
            instances.append(self)
            assert self._being_constructed()
            raise RuntimeError("constructor failed")

    with pytest.raises(RuntimeError, match="constructor failed"):
        Broken()
    assert len(instances) == 1
    assert not instances[0]._being_constructed()


def test_skip_if_being_constructed():
    @track_construction
    class Example:
        def __init__(self):
            self._value = 0
            # Must be ignored during construction.
            self.value = 123

        @property
        def value(self):
            return self._value

        @value.setter
        @skip_if_being_constructed
        def value(self, value):
            self._value = value

    obj = Example()
    assert obj.value == 0
    obj.value = 456
    assert obj.value == 456


def test_super_if_being_constructed():
    class Parent:
        def action(self, value):
            self.calls.append(("parent", value))
            return "parent"

    @track_construction
    class Child(Parent):
        def __init__(self):
            self.calls = []
            self.constructor_result = self.action("constructor")

        @super_if_being_constructed
        def action(self, value):
            self.calls.append(("child", value))
            return "child"

    obj = Child()
    assert obj.constructor_result == "parent"
    assert obj.calls == [("parent", "constructor")]
    assert obj.action("runtime") == "child"
    assert obj.calls == [("parent", "constructor"), ("child", "runtime")]


def test_super_if_being_constructed_with_decorated_subclass():
    class Parent:
        def action(self, value):
            self.calls.append(("parent", value))
            return "parent"

    @track_construction
    class Child(Parent):
        def __init__(self):
            self.calls = []
            self.child_result = self.action("child")

        @super_if_being_constructed
        def action(self, value):
            self.calls.append(("child", value))
            return "child"

    @track_construction
    class GrandChild(Child):
        def __init__(self):
            super().__init__()
            # Child's construction wrapper has exited here, but the
            # GrandChild construction wrapper is still active.
            self.grandchild_result = (self.action("grandchild"))

    obj = GrandChild()
    assert obj.child_result == "parent"
    assert obj.grandchild_result == "parent"
    assert obj.calls == [("parent", "child"), ("parent", "grandchild")]
    assert not obj._being_constructed()
    assert obj.action("runtime") == "child"
    assert obj.calls[-1] == ("child", "runtime")


def test_super_if_being_constructed_class_call():
    class Parent:
        def action(self, value):
            return ("parent", value)

    @track_construction
    class Child(Parent):
        def __init__(self):
            pass

        @super_if_being_constructed
        def action(self, value):
            return ("child", value)

    obj = Child()
    # __call__ on the descriptor supports Class.method(instance, ...).
    assert Child.action(obj, 123) == ("child", 123)
