# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

from typing import Iterable, Union

from .styled import Styled
from .helpers import (_disp_len, _split_string_spec, _unwrap, _as_plain_and_spans,
                      _slice_ms_by_plain_span, _split_spans_plain, _rsplit_spans_plain,
                      _splitlines_spans_plain, _partition_spans_plain, _rpartition_spans_plain)
from .graphemes import _truncate_parts_graphemes, _grapheme_spans_ms


class MultiStyled(str):
    __slots__ = ("_values", "_enabled")

    # names we must never proxy (our own attributes/methods & dunders)
    _OWN_ATTRS = frozenset({
        "_values", "_enabled", "enabled",
        "__class__", "__dict__", "__slots__", "__new__", "__init__",
        "__str__", "__repr__", "__len__", "__format__", "__getitem__",
        "__getattribute__", "__setattr__", "__reduce__", "__reduce_ex__",
        "__add__", "__radd__", "__iadd__","replace", "join",
    })

    # Registry of behaviours
    _PER_CHUNK = frozenset({
        "upper","lower","casefold","swapcase","title","capitalize",
        "translate","expandtabs",
    })
    _BOUNDARY = frozenset({
        "strip","lstrip","rstrip","removeprefix","removesuffix",
        "center","ljust","rjust","zfill",
    })
    _SCALAR = frozenset({
        "find","rfind","index","rindex","count","startswith","endswith",
        "isalnum","isalpha","isascii","isdecimal","isdigit","isidentifier",
        "islower","isnumeric","isprintable","isspace","istitle","isupper",
    })
    _SEGMENT = frozenset({"split","rsplit","splitlines","partition","rpartition"})
    _SPECIAL = frozenset({"replace","join"})  # handled explicitly

    def __new__(cls, *values, enabled=None):
        if all(isinstance(v, Styled) for v in values):
            if all(v._style_kwargs == values[0]._style_kwargs for v in values):
                return Styled(''.join(v._value for v in values), **values[0]._style_kwargs)
        return super().__new__(cls)

    def __init__(self, *values, enabled=None):
        vals = []
        flags = set()
        for v in values:
            if isinstance(v, MultiStyled):
                vals.extend(v._values)
                flags.add(v.enabled)
            elif isinstance(v, Styled):
                vals.append(v)
                flags.add(v.enabled)
            else:
                vals.append(v)
        self._values = vals
        if enabled is None:
            # No enabled given; infer from components
            flags.discard(None)  # Ignore None values
            if len(flags) > 1:
                raise ValueError("Cannot combine Styled with differing enabled states")
            enabled = next(iter(flags), None)
        self.enabled = enabled

    @property
    def enabled(self) -> bool | None:
        return self._enabled

    @enabled.setter
    def enabled(self, val: bool | None):
        for v in self._values:
            if isinstance(v, Styled):
                v.enabled = val
        object.__setattr__(self, "_enabled", val)

    @enabled.deleter
    def enabled(self):
        self.enabled = None

    def __repr__(self) -> str:
        if self.enabled is not None:
            return  f"MultiStyled({', '.join(repr(v) for v in self._values)}, enabled={self.enabled!r})"
        return  f"MultiStyled({', '.join(repr(v) for v in self._values)})"

    def __str__(self):         return "".join(str(v) for v in self._values)
    def __len__(self):         return sum(_disp_len(str(_unwrap(v))) for v in self._values)
    def __add__(self, other):  return MultiStyled(*self._values, other)
    def __radd__(self, other): return MultiStyled(other, *self._values)
    def __iadd__(self, other): return MultiStyled(*self._values, other)

    def __getitem__(self, key):
        _, gspans = _grapheme_spans_ms(self)
        ng = len(gspans)

        # normalise a single index
        if isinstance(key, int):
            idx = key + ng if key < 0 else key
            if idx < 0 or idx >= ng:
                raise IndexError("MultiStyled grapheme index out of range")
            a, b = gspans[idx]
            piece = _slice_ms_by_plain_span(self, a, b)
            # Collapse to a single Styled/str if it didn’t cross a boundary
            if len(piece._values) == 1:
                return piece._values[0]
            return piece

        # slicing by grapheme indices
        if isinstance(key, slice):
            # Build the list of grapheme indices following Python semantics
            indices = range(ng)[key]  # this handles start/stop/step/negatives
            # Fast-path: empty selection
            try:
                first = next(iter(indices))
            except StopIteration:
                return ''

            # Collect contiguous spans to minimise fragmentation
            spans = []
            prev = first
            start_a, start_b = gspans[prev]
            run_a, run_b = start_a, start_b

            def flush():
                if spans and run_a == spans[-1][1]:
                    # merge with previous (shouldn’t happen with arbitrary steps)
                    spans[-1] = (spans[-1][0], run_b)
                else:
                    spans.append((run_a, run_b))

            for idx in list(indices)[1:]:
                a, b = gspans[idx]
                # If we’re stepping by 1 and contiguous, extend; else flush/start new run
                if idx == prev + 1 and a == run_b:
                    run_b = b
                else:
                    flush()
                    run_a, run_b = a, b
                prev = idx
            flush()

            # Slice each span, then flatten into one MultiStyled
            parts = []
            for a, b in spans:
                seg = _slice_ms_by_plain_span(self, a, b)
                parts.extend(seg._values)
            return MultiStyled(*parts, enabled=self.enabled)

        raise TypeError(f"indices must be integers or slices, not {type(key).__name__}")


    # outer string formatting on the WHOLE (width/align/precision)
    def __format__(self, spec: str) -> str:
        fill, align, width, prec, tail = _split_string_spec(spec)
        if tail:
            # any exotic type → fallback
            return format(str(self), spec)

        # precision across visible width, preserving styles
        parts = self._values
        if prec is not None:
            # keep styles; slice last chunk
            parts = _truncate_parts_graphemes(parts, prec)

        core = MultiStyled(*parts, enabled=self.enabled)  # re-wrap so we can compute len again
        vis = len(core)
        if width is None or vis >= width:
            return str(core)
        pad = width - vis
        if align == '>': return f"{fill*pad}{core}"
        if align == '^': return f"{fill*(pad//2)}{core}{fill*(pad - pad//2)}"
        return f"{core}{fill*pad}"

    def __getattribute__(self, name):
        # 1) our own attributes always win
        if name in MultiStyled._OWN_ATTRS or name in type(self).__dict__:
            return object.__getattribute__(self, name)

        # 2) Special methods we handle explicitly
        if name in MultiStyled._SPECIAL:
            return object.__getattribute__(self, name)

        # 3) Known categories
        if name in MultiStyled._PER_CHUNK:
            def _op(*a, **k):
                new = []
                for v in object.__getattribute__(self, "_values"):
                    s = _unwrap(v)
                    res = getattr(s, name)(*a, **k)
                    new.append(Styled(res, **v._style_kwargs) if isinstance(v, Styled) else res)
                return MultiStyled(*new, enabled=object.__getattribute__(self, "_enabled"))
            return _op

        # 4) Boundary-based operations like strip, padding, etc.
        if name in MultiStyled._BOUNDARY:
            def _op(*a, **k):
                plain, _ = _as_plain_and_spans(self)
                res = getattr(plain, name)(*a, **k)
                # derive kept span(s)
                if name == "lstrip":
                    start = len(plain) - len(plain.lstrip(*a, **k)); end = len(plain)
                    return _slice_ms_by_plain_span(self, start, end)
                if name == "rstrip":
                    end = len(plain.rstrip(*a, **k)); return _slice_ms_by_plain_span(self, 0, end)
                if name == "strip":
                    start = len(plain) - len(plain.lstrip(*a, **k))
                    end = len(plain.rstrip(*a, **k))
                    return _slice_ms_by_plain_span(self, start, end)
                if name == "removeprefix":
                    pref = a[0] if a else k.get("prefix","")
                    start = len(pref) if plain.startswith(pref) else 0
                    return _slice_ms_by_plain_span(self, start, len(plain))
                if name == "removesuffix":
                    suff = a[0] if a else k.get("suffix","")
                    end = len(plain) - len(suff) if plain.endswith(suff) else len(plain)
                    return _slice_ms_by_plain_span(self, 0, end)
                # padding: add unstyled margins from the actual result
                if name in {"center","ljust","rjust","zfill"}:
                    out = res
                    # compute left/right pad lengths by comparing positions
                    # simplest: locate 'plain' inside 'out' (first occurrence)
                    j = out.find(plain) if plain else (len(out) // 2)
                    if j < 0:  # fallback
                        return MultiStyled(out)
                    left, right = j, len(out) - j - len(plain)
                    return MultiStyled(out[:left], _slice_ms_by_plain_span(self, 0, len(plain)), out[len(out)-right:])
            return _op

        # 5) Scalar operations that return a single value
        if name in MultiStyled._SCALAR:
            def _op(*a, **k):
                plain, _ = _as_plain_and_spans(self)
                return getattr(plain, name)(*a, **k)
            return _op

        # 6) Segmentation operations that return multiple parts
        if name in MultiStyled._SEGMENT:
            def _op(*a, **k):
                plain, _ = _as_plain_and_spans(self)
                if name == "split":
                    sep = a[0] if a else k.get("sep", None)
                    maxsplit = a[1] if len(a) > 1 else k.get("maxsplit", -1)
                    spans = _split_spans_plain(plain, sep, maxsplit)
                    return [_slice_ms_by_plain_span(self, x, y) for x, y in spans]
                if name == "rsplit":
                    sep = a[0] if a else k.get("sep", None)
                    maxsplit = a[1] if len(a) > 1 else k.get("maxsplit", -1)
                    spans = _rsplit_spans_plain(plain, sep, maxsplit)
                    return [_slice_ms_by_plain_span(self, x, y) for x, y in spans]
                if name == "splitlines":
                    keepends = a[0] if a else k.get("keepends", False)
                    spans = _splitlines_spans_plain(plain, keepends)
                    return [_slice_ms_by_plain_span(self, x, y) for x, y in spans]
                if name == "partition":
                    sep = a[0] if a else k["sep"]
                    left, mid, right = _partition_spans_plain(plain, sep)
                    if mid is None:
                        return (self, "", "")
                    return (_slice_ms_by_plain_span(self, *left),
                            _slice_ms_by_plain_span(self, *mid),
                            _slice_ms_by_plain_span(self, *right))
                if name == "rpartition":
                    sep = a[0] if a else k["sep"]
                    left, mid, right = _rpartition_spans_plain(plain, sep)
                    if mid is None:
                        return ("", "", self)
                    return (_slice_ms_by_plain_span(self, *left),
                            _slice_ms_by_plain_span(self, *mid),
                            _slice_ms_by_plain_span(self, *right))
            return _op

        # 7) Fall back to our own dict / default behaviour
        return object.__getattribute__(self, name)

    # --- explicit specials: replace / join (need bespoke logic) ---
    def replace(self, old: str, new: Union[str,"Styled","MultiStyled"], count: int = -1):
        plain, _ = _as_plain_and_spans(self)
        if old == "":
            # match str semantics for empty old: insert 'new' between graphemes; we’ll keep it simple and defer
            parts = [plain[0:0]]
            for ch in plain:
                parts.extend([new, ch])
            parts.append(new)
            segs = []
            idx = 0
            for p in parts:
                if p is new:
                    segs.append(new)
                else:
                    n = len(p)
                    segs.append(_slice_ms_by_plain_span(self, idx, idx+n))
                    idx += n
        else:
            segs, i, L, done = [], 0, len(old), 0
            while True:
                j = plain.find(old, i)
                if j < 0 or (count >= 0 and done >= count):
                    segs.append(_slice_ms_by_plain_span(self, i, len(plain))); break
                segs.append(_slice_ms_by_plain_span(self, i, j))
                segs.append(new)
                i = j + L; done += 1
        flat = []
        for seg in segs:
            if isinstance(seg, MultiStyled):
                flat.extend(seg._values)
            else:
                flat.append(seg)
        return MultiStyled(*flat, enabled=self.enabled)

    def join(self, iterable: Iterable[Union[str,"Styled","MultiStyled"]]):
        out = []
        first = True
        for item in iterable:
            if not first:
                out.extend(self._values)   # add separator (may be styled)
            first = False
            if isinstance(item, MultiStyled):
                out.extend(item._values)
            elif isinstance(item, Styled):
                out.append(item)
            else:
                out.append(str(item))
        return MultiStyled(*out, enabled=self.enabled)
