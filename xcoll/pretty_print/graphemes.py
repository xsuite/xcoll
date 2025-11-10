# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

try:
    import regex
except Exception:
    regex = None

from .helpers import _slice_visible, _truncate_parts, _as_plain_and_spans


def _grapheme_spans(text: str):
    """Return list of (start, end) spans for grapheme clusters in text.
       Falls back to per-codepoint spans if regex is unavailable."""
    if regex is None:
        n = len(text)
        return [(i, i+1) for i in range(n)]
    return [m.span() for m in regex.finditer(r"\X", text)]

def _grapheme_spans_ms(ms: "MultiStyled"):
    """Plain concatenation + grapheme spans over it."""
    plain, _ = _as_plain_and_spans(ms)
    return plain, _grapheme_spans(plain)

def _take_graphemes(s: str, n: int) -> str:
    """Return the first n graphemes of s."""
    if regex is None:
        return _slice_visible(s, n)

    if n <= 0:
        return ""
    count = 0
    end = 0
    for (a, b) in _grapheme_spans(s):
        end = b
        count += 1
        if count >= n:
            break
    return s[:end]

def _slice_graphemes(s: str, sl: slice | int) -> str:
    if regex is None:
        return s[sl]
    if isinstance(sl, int):
        # one grapheme at position key
        sl = slice(sl, sl+1, None)

    spans = list(_grapheme_spans(s))
    # translate grapheme indices to codepoint indices
    start, stop, step = sl.start, sl.stop, sl.step
    if step not in (None, 1, -1):
        # Non-unit steps with graphemes are messy; fall back to codepoints or raise
        return s[sl]
    n = len(spans)
    i0 = 0 if start is None else (start + n if start < 0 else start)
    i1 = n if stop  is None else (stop  + n if stop  < 0 else stop)
    i0 = max(0, min(n, i0)); i1 = max(0, min(n, i1))
    if (step or 1) == 1:
        a = spans[i0][0] if i0 < n else len(s)
        b = spans[i1][0] if i1 < n else len(s)
        return s[a:b]
    else:  # step == -1 (reverse)
        # Cheap fallback: reverse by grapheme clusters, then slice
        clusters = [s[a:b] for (a,b) in spans]
        return "".join(clusters[::-1][i0:i1])

def _truncate_parts_graphemes(parts, n_graphemes):
    """
    parts: iterable of chunks (str or Styled).
    Return a list of chunks truncated to the first n graphemes, preserving styles.
    """
    from .styled import Styled
    if regex is None:
        return _truncate_parts(parts, n_graphemes)

    out = []
    remaining = n_graphemes
    for v in parts:
        if remaining <= 0:
            break
        if isinstance(v, Styled):
            core = str(v._value)
            # count graphemes in this chunk
            count = sum(1 for _ in _grapheme_spans(core))
            if count <= remaining:
                out.append(v)
                remaining -= count
            else:
                sliced = _take_graphemes(core, remaining)
                out.append(Styled(sliced, **v._style_kwargs))
                remaining = 0
                break
        else:
            core = str(v)
            count = sum(1 for _ in _grapheme_spans(core))
            if count <= remaining:
                out.append(core)
                remaining -= count
            else:
                out.append(_take_graphemes(core, remaining))
                remaining = 0
                break
    return out
