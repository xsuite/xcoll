# copyright ############################### #
# This file is part of the Xcoll Package.  #
# Copyright (c) CERN, 2026.                #
# ######################################### #

from contextvars import ContextVar
from functools import wraps, update_wrapper
from types import MethodType


# Tuple rather than a set, because the same object can enter nested
# constructors (e.g. subclass __init__ -> superclass __init__).
_objects_under_construction = ContextVar(
    "_objects_under_construction",
    default=(),
)


def _is_being_constructed(obj):
    """Return whether `obj` is currently inside a tracked constructor."""
    return any(current is obj for current in _objects_under_construction.get())


def track_construction(cls):
    """Class decorator tracking execution of the class' constructor.

    The decorated class gets a ``_being_constructed()`` method that can be
    queried from methods and property setters.

    The mechanism is nesting-safe, so both a subclass and its parent can be
    decorated.
    """
    original_init = cls.__init__

    @wraps(original_init)
    def wrapped_init(self, *args, **kwargs):
        current = _objects_under_construction.get()
        token = _objects_under_construction.set(current + (self,))
        try:
            return original_init(self, *args, **kwargs)
        finally:
            _objects_under_construction.reset(token)

    cls.__init__ = wrapped_init

    def _being_constructed(self):
        return _is_being_constructed(self)

    cls._being_constructed = _being_constructed

    return cls


def skip_if_being_constructed(method):
    """Skip a method call while the instance is being constructed."""

    @wraps(method)
    def wrapper(self, *args, **kwargs):
        if self._being_constructed():
            return None
        return method(self, *args, **kwargs)

    return wrapper


class super_if_being_constructed:
    """Method decorator calling the parent implementation during construction.

    Outside construction, the decorated implementation is called normally.

    This is implemented as a descriptor so that the class in which the
    decorated method is defined is known. This makes

        super(<defining class>, self)

    correct also when the method is inherited by subclasses.
    """

    def __init__(self, method):
        self._method = method
        self._owner = None
        self._name = None
        update_wrapper(self, method)

    def __set_name__(self, owner, name):
        self._owner = owner
        self._name = name

    def _invoke(self, instance, *args, **kwargs):
        if instance._being_constructed():
            parent = super(self._owner, instance)
            method = getattr(parent, self._name)
            return method(*args, **kwargs)

        return self._method(instance, *args, **kwargs)

    def __get__(self, instance, owner=None):
        if instance is None:
            return self
        return MethodType(self._invoke, instance)

    def __call__(self, instance, *args, **kwargs):
        # Allows Class.method(instance, ...) as well as instance.method(...).
        return self._invoke(instance, *args, **kwargs)
