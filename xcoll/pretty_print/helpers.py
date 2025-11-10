# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

import re
from typing import List, Tuple, Iterable, Iterator, Optional
try:
    from wcwidth import wcswidth as _visible_len
except Exception:
    _visible_len = None


def _disp_len(s: str) -> int:
    if _visible_len is None:
        return len(s)
    n = _visible_len(s)
    return 0 if n < 0 else n

def _unwrap(v):
    # plain view used for numeric formatting
    from .styled import Styled
    return v._value if isinstance(v, Styled) else v

def _wrap_like(self, out):
    from .styled import Styled
    if isinstance(out, str):
        return Styled(out, **self._style_kwargs)
    elif isinstance(out, list):
        return [Styled(x, **self._style_kwargs) if isinstance(x, str) else x for x in out]
    elif isinstance(out, tuple):
        return tuple(Styled(x, **self._style_kwargs) if isinstance(x, str) else x for x in out)
    elif isinstance(out, Iterator) or isinstance(out, Iterable) and not isinstance(out, (str, bytes)):
        # lazy wrapper
        def _gen():
            for x in out:
                yield Styled(x, **self._style_kwargs) if isinstance(x, str) else x
        return _gen()
    return out

def _split_string_spec(spec: str):
    # returns (fill, align, width, precision, rest)
    fill, align, width, prec = ' ', '<', None, None
    j = 0
    if spec:
        if len(spec) >= 2 and spec[1] in '<^>':
            fill, align = spec[0], spec[1]
            j = 2
        elif spec[0:1] in '<^>':
            align = spec[0]
            j = 1
        k = j
        while k < len(spec) and spec[k].isdigit():
            k += 1
        if k > j: width = int(spec[j:k])
        j = k
        if j < len(spec) and spec[j] == '.':
            j += 1; k = j
            while k < len(spec) and spec[k].isdigit():
                k += 1
            if k > j: prec = int(spec[j:k])
            j = k
        # tail = any exotic presentation type; we return it to signal fallback
        tail = spec[j:]
    else:
        tail = ""
    return fill, align, width, prec, tail

def _slice_visible(s: str, n: int) -> str:
    # naive codepoint walk; swap in wcwidth/grapheme if you need perfection
    if n <= 0: return ""
    acc = 0
    for i, ch in enumerate(s):
        w = _disp_len(ch)
        if acc + w > n: return s[:i]
        acc += w
    return s

def _truncate_parts(parts, prec: int):
    from .styled import Styled
    out = []
    remaining = prec
    for v in parts:
        core = str(_unwrap(v))
        w = _disp_len(core)
        if w <= remaining:
            out.append(v)
            remaining -= w
        else:
            cut = _slice_visible(core, remaining)
            if isinstance(v, Styled):
                out.append(Styled(cut, **v._style_kwargs))
            else:
                out.append(cut)
            break
    return out

def _as_plain_and_spans(ms: "MultiStyled") -> Tuple[str, List[Tuple[int, int, object]]]:
    """Plain string + list of (start,end,chunk) in plain indices."""
    from .styled import Styled
    parts, spans, pos = [], [], 0
    for v in ms._values:
        s = _unwrap(v)
        parts.append(s)
        end = pos + len(s)
        spans.append((pos, end, v))
        pos = end
    return "".join(parts), spans

def _slice_ms_by_plain_span(ms: "MultiStyled", start: int, end: int) -> "MultiStyled":
    """Build a MultiStyled covering [start,end) in plain indices, preserving per-chunk style."""
    from .styled import Styled
    from .multi import MultiStyled
    if start >= end:
        return ''  # empty piece
    _, spans = _as_plain_and_spans(ms)
    out = []
    for c0, c1, v in spans:
        if c1 <= start or c0 >= end:
            continue
        sub0 = max(start, c0) - c0
        sub1 = min(end, c1) - c0
        if isinstance(v, Styled):
            core = str(v._value)
            out.append(Styled(core[sub0:sub1], **v._style_kwargs))
        else:
            out.append(str(v)[sub0:sub1])
    return MultiStyled(*out, enabled=ms.enabled)

def _split_spans_plain(text: str, sep: Optional[str], maxsplit: int):
    if sep is None:
        spans = [(m.start(), m.end()) for m in re.finditer(r"\S+", text)]
        if maxsplit >= 0 and len(spans) > maxsplit + 1:
            spans = spans[:maxsplit+1]
        return spans
    if sep == "":  # str.split raises
        raise ValueError("empty separator")
    spans, n, i, L, cnt = [], len(text), 0, len(sep), 0
    while (maxsplit < 0 or cnt < maxsplit):
        j = text.find(sep, i)
        if j < 0:
            break
        spans.append((i, j))
        i = j + L
        cnt += 1
    spans.append((i, n))
    return spans

def _rsplit_spans_plain(text: str, sep: Optional[str], maxsplit: int):
    if sep is None:
        spans = _split_spans_plain(text, None, -1)
        if maxsplit >= 0 and len(spans) > maxsplit + 1:
            spans = spans[-(maxsplit+1):]
        return spans
    if sep == "":
        raise ValueError("empty separator")
    spans, i, L, cnt = [], len(text), len(sep), 0
    while (maxsplit < 0 or cnt < maxsplit):
        j = text.rfind(sep, 0, i)
        if j < 0:
            break
        spans.append((j+L, i))
        i = j
        cnt += 1
    spans.append((0, i))
    spans.reverse()
    return spans

def _splitlines_spans_plain(text: str, keepends: bool):
    spans, i, n = [], 0, len(text)
    while i < n:
        j = i
        while j < n and text[j] not in "\r\n":
            j += 1
        line_end = j
        if j < n:
            if text[j] == "\r" and j+1 < n and text[j+1] == "\n":
                j += 2
                end = j if keepends else line_end
            else:
                j += 1
                end = j if keepends else line_end
            spans.append((i, end))
            i = j
        else:
            spans.append((i, line_end))
            i = j
    return spans

def _partition_spans_plain(text: str, sep: str):
    j = text.find(sep)
    if j < 0:
        return (0, 0), None, (0, len(text))
    return (0, j), (j, j+len(sep)), (j+len(sep), len(text))

def _rpartition_spans_plain(text: str, sep: str):
    j = text.rfind(sep)
    if j < 0:
        return (0, 0), None, (0, len(text))
    return (0, j), (j, j+len(sep)), (j+len(sep), len(text))
