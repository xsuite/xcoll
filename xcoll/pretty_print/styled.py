# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

from .style import style
from .graphemes import _take_graphemes, _slice_graphemes
from .helpers import _disp_len, _split_string_spec, _wrap_like


class Styled(str):
    """
    Wrap any value; format it first (respecting .2f, g, %, etc.),
    then style only the inner content, and finally do width/align padding
    using the visible length (ignoring ANSI codes).
    """
    __slots__ = ("_value", "_style_kwargs")
    _STYLE_DEFAULTS = {
        "color": None,
        "bg": None,
        "bold": False,
        "italic": False,
        "underline": False,
        "dim": False,
        "enabled": None,
        "reset": True,
    }

    # names we must never proxy (our own attributes/methods & dunders)
    _OWN_ATTRS = frozenset({
        "_value", "_style_kwargs", "_STYLE_DEFAULTS",
        "__class__", "__dict__", "__slots__", "__new__", "__init__",
        "__str__", "__repr__", "__len__", "__format__", "__getitem__",
        "__getattribute__", "__setattr__", "__reduce__", "__reduce_ex__",
        "__add__", "__radd__", "__iadd__", "join",
    })

    # computed once: public callable names on str to proxy
    _STR_CALLABLES = frozenset({
        name for name in dir(str)
        if not name.startswith("__") and callable(getattr(str, name, None))
    })

    def __new__(cls, value, **style_kwargs):
        if value is None:
            value = ''
        value = str(value)
        if value == '':
            # empty string: no need to style
            return value
        return super().__new__(cls, value)

    def __init__(self, value, **style_kwargs):
        object.__setattr__(self, '_value', value)
        if "colour" in style_kwargs:
            if "color" in style_kwargs:
                raise ValueError("Cannot specify both 'color' and 'colour'")
            style_kwargs["color"] = style_kwargs.pop("colour")
        bad = [k for k in style_kwargs if k not in self._STYLE_DEFAULTS]
        if bad:
            raise TypeError(f"Styled() unexpected kwargs: {', '.join(bad)}")
        kw = self._STYLE_DEFAULTS.copy()
        kw.update(style_kwargs)
        object.__setattr__(self, "_style_kwargs", kw)

    def __repr__(self) -> str:
        kw = [ f"{kk}={vv!r}" for kk, vv in self._style_kwargs.items()
              if vv != self._STYLE_DEFAULTS[kk]]
        if kw:
            return f"Styled(value={self._value!r}, {', '.join(kw)})"
        return f"Styled(value={self._value!r})"

    def __str__(self):
        return style(str(self._value), **self._style_kwargs)

    def __len__(self):
        return _disp_len(str(self._value))

    def __add__(self, other):
        from .multi import MultiStyled
        return MultiStyled(self, other)

    def __radd__(self, other):
        from .multi import MultiStyled
        return MultiStyled(other, self)

    def __iadd__(self, other):
        from .multi import MultiStyled
        return MultiStyled(self, other)

    def __getitem__(self, key):
        return Styled(_slice_graphemes(str(self._value), key), **self._style_kwargs)

    def __format__(self, spec: str) -> str:
        fill, align, width, prec, inner = _split_string_spec(spec)

        # Numeric formatting etc. (no width/align here)
        core_plain = format(self._value, inner)

        if prec is not None:
            # GRAPHEME truncation for precision:
            core_plain = _take_graphemes(core_plain, prec)

        # Style the truncated core
        core_styled = style(core_plain, **self._style_kwargs)
        if width is None:
            return core_styled
        pad = width - _disp_len(core_plain)
        if pad <= 0: return core_styled
        if align == '>':  return f"{fill*pad}{core_styled}"
        if align == '^':  return f"{fill*(pad//2)}{core_styled}{fill*(pad - pad//2)}"
        return f"{core_styled}{fill*pad}"

    def join(self, iterable):
        from .multi import MultiStyled
        out = []
        first = True
        for item in iterable:
            if not first:
                out.append(self)  # the separator (Styled)
            first = False
            if isinstance(item, MultiStyled):
                out.extend(item._values)
            elif isinstance(item, Styled):
                out.append(item)
            else:
                out.append(str(item))
        style = self._style_kwargs
        for v in out:
            if isinstance(v, Styled) and v._style_kwargs != style:
                # Different styles; must preserve as MultiStyled
                return MultiStyled(*out, enabled=self.enabled if hasattr(self, "enabled") else None)
        # Single style for all parts; can collapse to one Styled
        out = [v._value if isinstance(v, Styled) else str(v) for v in out]
        return Styled(''.join(out), **style)

    def __getattribute__(self, name):
        # 1) our own attributes always win
        if name in Styled._OWN_ATTRS or name in type(self).__dict__:
            return object.__getattribute__(self, name)

        # 2) style fields as attributes (colour, bold, etc.)
        try:
            skw = object.__getattribute__(self, "_style_kwargs")
            if name in skw:
                return skw[name]
        except Exception:
            pass  # during early init, fall through

        # 3) if it's a public str method, return a wrapper bound to the *plain* core
        if name in Styled._STR_CALLABLES:
            # bind to the underlying plain string (your stored core)
            val = object.__getattribute__(self, "_value")
            bound = getattr(str(val), name, None)
            if callable(bound):
                def _wrapped(*a, **k):
                    out = bound(*a, **k)  # operate on the plain core
                    return _wrap_like(self, out)
                return _wrapped
            else:
                return bound

        # 4) otherwise, default lookup (covers properties on str, etc.)
        return object.__getattribute__(self, name)

    def __setattr__(self, name, value):
        if name == "colour": name = "color"
        skw = object.__getattribute__(self, "_style_kwargs")
        if name in skw:
            skw[name] = value
            return
        raise AttributeError(f"'{type(self).__name__}' has no attribute "
                             f"'{name}' and no __dict__ for setting new "
                             f"attributes")
