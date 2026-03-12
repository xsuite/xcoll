# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2026.                 #
# ######################################### #

from __future__ import annotations

from typing import Any, Callable, TypeVar, Generic


T = TypeVar("T")

class Proxy(Generic[T]):
    """Delegating proxy with interception hooks.

    This class forwards attribute access to an underlying target instance
    whenever the requested attribute is not defined on the proxy itself.
    It is intended as a lightweight base for wrappers that should behave
    almost transparently like another object, while still allowing custom
    behaviour on delegated reads, writes, and deletions.

    The proxy follows a "proxy-first, target-second" lookup policy:

    - If an attribute is defined on the proxy class or one of its bases,
      it is handled by the proxy itself.
    - Otherwise, the operation is forwarded to the target instance.
    - Delegated reads, writes, and deletions may be intercepted through
      the `_on_get`, `_on_set`, and `_on_del` hooks.

    Special methods such as `__len__` or `__getitem__` are not
    automatically delegated by Python's normal attribute lookup rules.
    To make operators, builtins, and protocols behave transparently,
    a collection of common special methods is injected onto the proxy
    class and forwarded explicitly to the target.

    Parameters
    ----------
    target : T
        The target instance to which non-proxy attributes and methods are
        delegated.

    Notes
    -----
    The target is stored in the internal attribute `_target` and is
    intended to be immutable after initialisation. Normal attribute
    assignment to `_target` is rejected after construction.

    This protection prevents accidental reassignment in subclasses, but
    it is not meant as a security boundary: deliberate use of
    `object.__setattr__` can still bypass it.

    The proxy does not wrap delegated callables. If a delegated attribute
    is a bound method, the original bound method of the target is
    returned. This preserves its normal behaviour and usually gives
    better introspection, signatures, and debugging output than wrapping
    it in an additional closure.

    Attributes
    ----------
    _target : T
        The wrapped target instance.

    Examples
    --------
    A minimal proxy delegates unknown attributes to the target:

    >>> class Target:
    ...     def __init__(self):
    ...         self.value = 3
    ...         self._prop = 8
    ...         self._secret = 10
    ...
    ...     @property
    ...     def prop(self):
    ...         return self._prop
    ...
    ...     @prop.setter
    ...     def prop(self, value):
    ...         self._prop = value
    ...
    ...     def multiply(self, x):
    ...         return self.value * x
    ...
    ...     def __len__(self):
    ...         return 7
    ...
    >>> p = Proxy(Target())
    >>> p.value
    3
    >>> p.prop
    8
    >>> p.prop = 11
    >>> p.prop
    11
    >>> p._secret
    10
    >>> p.multiply(5)
    15
    >>> len(p)
    7

    Proxy-defined attributes remain on the proxy itself. A subclass can
    add its own state and intercept delegated operations:

    >>> class LoggingProxy(Proxy[Target]):
    ...     __slots__ = ("_enabled",)
    ...
    ...     def __init__(self, target, enabled=True):
    ...         super().__init__(target)
    ...         self._enabled = enabled
    ...
    ...     @property
    ...     def enabled(self):
    ...         return self._enabled
    ...
    ...     @enabled.setter
    ...     def enabled(self, value):
    ...         self._enabled = bool(value)
    ...
    ...     def _on_get(self, name, value):
    ...         if self.enabled:
    ...             print(f"GET {name} -> {value!r}")
    ...         return value
    ...
    ...     def _on_set(self, name, value):
    ...         if self.enabled:
    ...             print(f"SET {name} = {value!r}")
    ...         return value
    ...
    ...     def multiply(self, x):
    ...         if self.enabled:
    ...             print(f"CALL multiply({x!r})")
    ...         return self._target.multiply(x)
    ...
    >>> t = Target()
    >>> p = LoggingProxy(t)
    >>> p.value
    GET value -> 3
    3
    >>> p.value = 4
    SET value = 4
    >>> t.value
    4
    >>> p.prop
    GET prop -> 8
    8
    >>> p.prop = -99
    SET prop = -99
    >>> t.prop
    -99
    >>> p._secret
    GET _secret -> 10
    10
    >>> p.multiply(5)
    CALL multiply(5)
    15
    >>> len(p)
    7
    >>> p.enabled = False
    >>> p.value
    4

    Attributes defined on the proxy class are handled by the proxy, while
    all other attributes, including underscore-prefixed ones, are forwarded
    to the target unless they are true proxy internals such as `_target`.
    """

    __slots__ = ("_target",)
    _PROXY_INTERNALS = frozenset({
        "_target",
        "_on_get",
        "_on_set",
        "_on_del",
        "__class__",
        "__dict__",
        "__slots__",
    })

    def __init__(self, target: T) -> None:
        object.__setattr__(self, "_target", target)

    # ----------------
    # Optional hooks
    # ----------------

    def _on_get(self, name: str, value: Any) -> Any:
        """Hook for delegated attribute reads."""
        return value

    def _on_set(self, name: str, value: Any) -> Any:
        """Hook for delegated attribute writes."""
        return value

    def _on_del(self, name: str) -> None:
        """Hook for delegated attribute deletions."""
        return None

    # -------------------------
    # Attribute access logic
    # -------------------------

    def __getattribute__(self, name: str) -> Any:
        # Never delegate true proxy internals.
        if name in object.__getattribute__(self, "_PROXY_INTERNALS"):
            return object.__getattribute__(self, name)

        # First try the proxy itself.
        try:
            return object.__getattribute__(self, name)
        except AttributeError:
            pass

        # Otherwise delegate to the target.
        target = object.__getattribute__(self, "_target")
        value = getattr(target, name)
        hook = object.__getattribute__(self, "_on_get")
        return hook(name, value)

    def __setattr__(self, name: str, value: Any) -> None:
        # _target is proxy-internal and write-once.
        if name == "_target":
            if hasattr(self, "_target"):
                raise AttributeError(
                    f"{type(self).__name__}._target is read-only after initialisation"
                )
            object.__setattr__(self, name, value)
            return

        tp = type(self)

        # If the proxy class defines this name, keep it on the proxy.
        descr = _get_proxy_attr(tp, name)
        if descr is not _MISSING:
            # Data descriptor, e.g. property with setter.
            if hasattr(descr, "__set__"):
                descr.__set__(self, value)
                return

            # Slot or plain proxy attribute.
            object.__setattr__(self, name, value)
            return

        # Otherwise forward to the target.
        value = object.__getattribute__(self, "_on_set")(name, value)
        setattr(object.__getattribute__(self, "_target"), name, value)

    def __delattr__(self, name: str) -> None:
        if name == "_target":
            raise AttributeError(
                f"{type(self).__name__}._target cannot be deleted"
            )

        tp = type(self)

        # If the proxy class defines this name, delete it on the proxy.
        descr = _get_proxy_attr(tp, name)
        if descr is not _MISSING:

            if hasattr(descr, "__delete__"):
                descr.__delete__(self)
                return

            object.__delattr__(self, name)
            return

        # Otherwise forward to the target.
        object.__getattribute__(self, "_on_del")(name)
        delattr(object.__getattribute__(self, "_target"), name)

    # -------------------------
    # Convenience
    # -------------------------

    def __dir__(self) -> list[str]:
        own = set(object.__dir__(self))
        target = set(dir(object.__getattribute__(self, "_target")))
        return sorted(own | target)

    def __repr__(self) -> str:
        return f"{type(self).__name__}({self._target!r})"


_MISSING = object()

def _get_proxy_attr(tp: type, name: str) -> Any:
    """Return the descriptor/class attribute for `name` on the proxy MRO."""
    for cls in tp.__mro__:
        if cls is object:
            break
        if name in cls.__dict__:
            return cls.__dict__[name]
    return _MISSING

def _make_forwarder(name: str) -> Callable[..., Any]:
    """Create a special-method forwarder to the underlying target."""
    def method(self: Proxy[Any], *args: Any, **kwargs: Any) -> Any:
        return getattr(self._target, name)(*args, **kwargs)
    method.__name__ = name
    return method


# A fairly comprehensive (but still finite) list of “common” special
# methods to forward.
_SPECIALS_TO_FORWARD: tuple[str, ...] = (
    # Basic protocol / representation
    "__bytes__",
    "__format__",
    "__str__",

    # Container / iteration
    "__len__",
    "__iter__",
    "__reversed__",
    "__contains__",
    "__getitem__",
    "__setitem__",
    "__delitem__",

    # Callable
    "__call__",

    # Context manager
    "__enter__",
    "__exit__",

    # Comparisons
    "__lt__",
    "__le__",
    "__eq__",
    "__ne__",
    "__gt__",
    "__ge__",
    "__hash__",

    # Numeric / unary
    "__neg__",
    "__pos__",
    "__abs__",
    "__invert__",

    # Numeric / binary
    "__add__",
    "__sub__",
    "__mul__",
    "__matmul__",
    "__truediv__",
    "__floordiv__",
    "__mod__",
    "__pow__",
    "__and__",
    "__or__",
    "__xor__",
    "__lshift__",
    "__rshift__",

    # Reflected numeric / binary
    "__radd__",
    "__rsub__",
    "__rmul__",
    "__rmatmul__",
    "__rtruediv__",
    "__rfloordiv__",
    "__rmod__",
    "__rpow__",
    "__rand__",
    "__ror__",
    "__rxor__",
    "__rlshift__",
    "__rrshift__",

    # In-place numeric / binary
    "__iadd__",
    "__isub__",
    "__imul__",
    "__imatmul__",
    "__itruediv__",
    "__ifloordiv__",
    "__imod__",
    "__ipow__",
    "__iand__",
    "__ior__",
    "__ixor__",
    "__ilshift__",
    "__irshift__",
)

# Inject special-method forwarders onto the Proxy class.
# (These must be on the class for Python to use them in operators/builtins.)
for _name in _SPECIALS_TO_FORWARD:
    if not hasattr(Proxy, _name):
        setattr(Proxy, _name, _make_forwarder(_name))
