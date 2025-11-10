# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

import sys
import types
from pathlib import Path

import pytest

from xcoll.pretty_print import Styled, MultiStyled

try:
    import regex  # noqa: F401
except ImportError:
    regex = None

# --- Utilities ---------------------------------------------------------------

@pytest.fixture(autouse=True)
def fake_style(monkeypatch):
    """
    Patch the module-level `style(text, **kwargs)` used by Styled.__str__/__format__
    to a deterministic marker-based formatter so we can assert reliably.

    Output form: <S|sorted_kwargs>text</S>
    Example: <S|bold=True;color=green>hello</S>
    """
    mod = sys.modules[Styled.__module__]

    def format_kwargs(kwargs):
        items = [f"{k}={repr(v)}" for k, v in sorted(kwargs.items()) if v is not None and v is not False]
        return ";".join(items)

    def _fake_style(text, **kwargs):
        # Keep 'enabled=None' invisible in tags, but if enabled=False then show it
        if kwargs.get("enabled", None) is None:
            kwargs = {k: v for k, v in kwargs.items() if k != "enabled"}
        return f"<S|{format_kwargs(kwargs)}>{text}</S>"

    monkeypatch.setattr(mod, "style", _fake_style, raising=True)
    yield

def strip_tags(rendered: str) -> str:
    # Helper: remove our fake tags
    if rendered.startswith("<S|") and rendered.endswith("</S>"):
        # drop leading <S|...> (find the first > after 3)
        i = rendered.find(">")
        return rendered[i+1:-4]
    return rendered


def plain_text(val):
    if isinstance(val, (Styled, MultiStyled)):
        return strip_tags(str(val))
    return val


_STRING_TRANSLATE_TABLE = str.maketrans({"a": "1", "b": "2"})

STRING_METHOD_CASES = [
    ("capitalize", "hELLo world", (), {}),
    ("casefold", "Straße", (), {}),
    ("center", "core", (9, "*"), {}),
    ("count", "banana", ("na",), {}),
    ("encode", "café", ("utf-8",), {}),
    ("endswith", "foobar", ("bar",), {}),
    ("expandtabs", "col1\tcol2", (4,), {}),
    ("find", "needle haystack", ("hay",), {}),
    ("format", "{} + {}", (1, 2), {}),
    ("format_map", "{x}-{y}", ({"x": 1, "y": 2},), {}),
    ("index", "hello", ("ll",), {}),
    ("isalnum", "abc123", (), {}),
    ("isalpha", "abc", (), {}),
    ("isascii", "abc", (), {}),
    ("isdecimal", "123", (), {}),
    ("isdigit", "123", (), {}),
    ("isidentifier", "foo_bar", (), {}),
    ("islower", "abc", (), {}),
    ("isnumeric", "123", (), {}),
    ("isprintable", "abc123", (), {}),
    ("isspace", "   ", (), {}),
    ("istitle", "Hello World", (), {}),
    ("isupper", "ABC", (), {}),
    ("ljust", "foo", (6, "."), {}),
    ("lower", "Foo", (), {}),
    ("lstrip", "  foo", (), {}),
    ("maketrans", "abc", ("abc", "123"), {}),
    ("partition", "foo=bar", ("=",), {}),
    ("removeprefix", "prefix_value", ("prefix_",), {}),
    ("removesuffix", "value_suffix", ("_suffix",), {}),
    ("replace", "banana", ("na", "NA"), {}),
    ("rfind", "banana", ("na",), {}),
    ("rindex", "bananas", ("na",), {}),
    ("rjust", "foo", (6, "."), {}),
    ("rpartition", "foo=bar=baz", ("=",), {}),
    ("rsplit", "a,b,c", (",", 1), {}),
    ("rstrip", "foo   ", (), {}),
    ("split", "a,b,,c", (",",), {}),
    ("splitlines", "a\nb\r\n", (True,), {}),
    ("startswith", "foobar", ("foo",), {}),
    ("strip", "  foo  ", (), {}),
    ("swapcase", "AbC", (), {}),
    ("title", "hello world", (), {}),
    ("translate", "abca", (_STRING_TRANSLATE_TABLE,), {}),
    ("upper", "foo", (), {}),
    ("zfill", "42", (5,), {}),
]

STRING_METHOD_NAMES = {case[0] for case in STRING_METHOD_CASES}


def make_multistyled(value: str, colour: str = "magenta") -> MultiStyled:
    split = max(1, len(value) // 2)
    return MultiStyled(Styled(value[:split], colour=colour), value[split:])


_MULTI_TRANSLATE_TABLE = str.maketrans({"a": "1", "b": "2"})

MULTI_PER_CHUNK_CASES = [
    ("upper", "Abc!", (), {}),
    ("lower", "ABc!", (), {}),
    ("casefold", "Straße", (), {}),
    ("swapcase", "aBc", (), {}),
    ("title", "hello world", (), {}),
    ("capitalize", "hello world", (), {}),
    ("translate", "abba", (_MULTI_TRANSLATE_TABLE,), {}),
    ("expandtabs", "a\tb\t", (4,), {}),
]

MULTI_SCALAR_CASES = [
    ("find", "hello world", ("lo",), {}),
    ("rfind", "banana", ("na",), {}),
    ("index", "hello world", ("wor",), {}),
    ("rindex", "bananas", ("na",), {}),
    ("count", "banana", ("na",), {}),
    ("startswith", "foobar", ("foo",), {}),
    ("endswith", "foobar", ("bar",), {}),
    ("isalnum", "abc123", (), {}),
    ("isalpha", "abc", (), {}),
    ("isascii", "abc", (), {}),
    ("isdecimal", "123", (), {}),
    ("isdigit", "123", (), {}),
    ("isidentifier", "foo_bar", (), {}),
    ("islower", "abc", (), {}),
    ("isnumeric", "123", (), {}),
    ("isprintable", "abc123", (), {}),
    ("isspace", "   ", (), {}),
    ("istitle", "Hello World", (), {}),
    ("isupper", "ABC", (), {}),
]

MULTI_SEGMENT_CASES = [
    ("split", "one two  three", (), {}),
    ("split", "a,,b", (",",), {}),
    ("rsplit", "a,,b,,c", (",", 2), {}),
    ("splitlines", "a\nb\r\nc", (), {}),
    ("splitlines", "a\nb\r\nc", (True,), {}),
    ("partition", "foo=bar=baz", ("=",), {}),
    ("partition", "foobar", ("=",), {}),
    ("rpartition", "foo=bar=baz", ("=",), {}),
    ("rpartition", "foobar", ("=",), {}),
]

# --- Styled: basics ----------------------------------------------------------

def test_styled_repr_and_attrs():
    s = Styled("hello", colour="green", bold=True, reset=False)
    assert "Styled(value='hello'" in repr(s)
    assert s.color == "green"
    assert s.bold is True
    assert s.reset is False
    # Setting colour alias
    s.colour = "blue"
    assert s.color == "blue"

def test_styled_str_renders_with_style_tags():
    s = Styled("hello", colour="green", bold=True)
    out = str(s)
    assert out.startswith("<S|")
    assert "color='green'" in out or "colour='green'" in out
    assert "bold=True" in out
    assert strip_tags(out) == "hello"

def test_styled_len_counts_visible_chars_not_tags():
    s = Styled("hello", colour="green", bold=True)
    assert len(s) == 5

def test_styled_empty_value_returns_plain_string():
    empty = Styled("")
    assert empty == ""
    assert type(empty) is str
    none_input = Styled(None)
    assert none_input == ""
    assert type(none_input) is str

def test_styled_string_methods_are_wrapped_and_preserve_style_kwargs():
    s = Styled("Hello", colour="green")
    lo = s.lower()
    assert isinstance(lo, Styled)
    assert lo.color == "green"
    assert strip_tags(str(lo)) == "hello"
    up_list = s.split("e")
    assert all(isinstance(p, Styled) for p in up_list)
    assert [strip_tags(str(p)) for p in up_list] == ["H", "llo"]

def test_styled_join_returns_multistyled_and_keeps_styles():
    a = Styled("A", colour="green")
    b = Styled("B", colour="blue")
    sep = Styled("-", dim=True)
    out = sep.join([a, b])
    assert isinstance(out, MultiStyled)
    pieces = [strip_tags(str(x)) if isinstance(x, Styled) else strip_tags(str(x)) for x in out._values]
    assert pieces == ["A", "-", "B"]
    # The middle is the separator Styled
    assert isinstance(out._values[1], Styled) and out._values[1].dim is True

def test_styled_join_collapses_to_single_styled_when_styles_match():
    pieces = [Styled("a", colour="green"), Styled("b", colour="green")]
    sep = Styled("-", colour="green")
    out = sep.join(pieces)
    assert isinstance(out, Styled)
    assert plain_text(out) == "a-b"
    assert out.color == "green"

def test_styled_add_variants_create_multistyled_instances():
    left = Styled("L", colour="cyan")
    right = Styled("R", colour="magenta")
    combo = left + right
    assert isinstance(combo, MultiStyled)
    assert [plain_text(chunk) for chunk in combo._values] == ["L", "R"]
    rcombo = "!" + left
    assert isinstance(rcombo, MultiStyled)
    assert [plain_text(chunk) for chunk in rcombo._values] == ["!", "L"]
    left_aug = left
    left_aug += "?"
    assert isinstance(left_aug, MultiStyled)
    assert [plain_text(chunk) for chunk in left_aug._values] == ["L", "?"]

def test_styled_format_numeric_then_style_and_pad():
    x = Styled(3.14159, colour="yellow")
    # Precision on number + right align width
    out = f"|{x:>8.2f}|"
    # Inside bars, there must be 8 characters of visible text, spaces on the left
    inner = out[1:-1]
    assert len(strip_tags(inner).lstrip()) == 4   # "3.14"
    assert len(inner) == len(strip_tags(inner)) + (len(inner) - len(strip_tags(inner)))  # sanity
    assert strip_tags(inner).endswith("3.14")
    # Ensure left padding spaces exist
    assert strip_tags(inner).startswith(" " * (8 - 4))

def test_styled_string_method_dataset_matches_builtin_surface_area():
    expected = {name for name in dir(str) if not name.startswith("__") and callable(getattr(str, name))}
    expected.discard("join")  # custom implementation tested separately
    assert STRING_METHOD_NAMES == expected

@pytest.mark.parametrize("method,value,args,kwargs", STRING_METHOD_CASES)
def test_styled_string_methods_match_python_str(method, value, args, kwargs):
    styled = Styled(value, colour="magenta", underline=True)
    result = getattr(styled, method)(*args, **kwargs)
    expected = getattr(value, method)(*args, **kwargs)
    if isinstance(expected, str):
        assert isinstance(result, Styled)
        assert plain_text(result) == expected
        assert result.color == "magenta"
        assert result.underline is True
    elif isinstance(expected, list):
        assert isinstance(result, list)
        assert [plain_text(item) for item in result] == expected
    elif isinstance(expected, tuple):
        assert isinstance(result, tuple)
        assert tuple(plain_text(item) for item in result) == tuple(expected)
    else:
        assert result == expected

# --- MultiStyled: construction, str, len -------------------------------------

def test_multistyled_construct_and_str_len():
    a = Styled("Hello", colour="red")
    b = " "
    c = Styled("world", colour="green")
    ms = MultiStyled(a, b, c)
    rendered = str(ms)
    assert strip_tags(rendered.replace("<S|", "").replace("</S>", "")) == "Hello world"
    assert len(ms) == len("Hello world")

def test_multistyled_enabled_propagates():
    a = Styled("x", enabled=None)
    b = Styled("y", enabled=None)
    ms = MultiStyled(a, b, enabled=True)
    assert all(v.enabled is True for v in ms._values if isinstance(v, Styled))
    ms.enabled = False
    assert all(v.enabled is False for v in ms._values if isinstance(v, Styled))

def test_multistyled_enabled_conflict_raises_value_error():
    a = Styled("a", enabled=True)
    b = Styled("b", enabled=False)
    with pytest.raises(ValueError):
        MultiStyled(a, b)

# --- MultiStyled: split/rsplit/splitlines/partition/rpartition ---------------

def test_multistyled_split_whitespace_and_maxsplit():
    ms = MultiStyled(Styled("Hello", colour="red"), "  ", Styled("world", colour="green"), "  !  ")
    parts = ms.split()
    plain = [strip_tags(str(p)) for p in parts]
    assert plain == ["Hello", "world", "!"]
    parts2 = ms.split(None, 1)
    assert [strip_tags(str(p)) for p in parts2] == ["Hello", "world  !  "]

def test_multistyled_split_fixed_sep_and_rsplit():
    ms = MultiStyled(Styled("foo", colour="red"), "oXo", Styled("bar", colour="blue"))
    parts = ms.split("o")
    assert [strip_tags(str(p)) for p in parts] == ["f", "", "X", "bar"]
    rparts = ms.rsplit("o", 1)
    assert [strip_tags(str(p)) for p in rparts] == ["fooX", "bar"]

def test_multistyled_splitlines_keepends():
    ms = MultiStyled(Styled("a\r\n", bold=True), "b\nc")
    no_ends = [strip_tags(str(p)) for p in ms.splitlines()]
    with_ends = [strip_tags(str(p)) for p in ms.splitlines(True)]
    assert no_ends == ["a", "b", "c"]
    assert with_ends == ["a\r\n", "b\n", "c"]

def test_multistyled_partition_and_rpartition():
    ms = MultiStyled("ab--", Styled("cd", colour="green"), "--ef")
    left, sep, right = ms.partition("--")
    assert [strip_tags(str(x)) for x in (left, sep, right)] == ["ab", "--", "cd--ef"]
    left2, sep2, right2 = ms.rpartition("--")
    assert [strip_tags(str(x)) for x in (left2, sep2, right2)] == ["ab--cd", "--", "ef"]

# --- MultiStyled: replace/join ------------------------------------------------

def test_multistyled_replace_with_plain_and_styled_insertion():
    ms = MultiStyled(Styled("foo", colour="red"), "-", Styled("bar", colour="blue"))
    # plain insertion
    out1 = ms.replace("-", "X")
    assert [strip_tags(str(x)) for x in out1._values] == ["foo", "X", "bar"]
    # styled insertion
    ins = Styled("+", underline=True)
    out2 = ms.replace("-", ins)
    assert isinstance(out2._values[1], Styled) and out2._values[1].underline is True

def test_join_plain_separator_loses_object_type_but_keeps_ansi_in_str():
    a = Styled("hello", colour="green")
    b = Styled("worldie", colour="blue")
    out = "\n".join([a, b])  # Python calls str.join → returns plain str
    assert isinstance(out, str)
    assert "hello" in out and "worldie" in out
    # These contain our fake tags because Styled.__str__ was invoked
    assert "<S|" in out

def test_join_styled_separator_keeps_multistyled():
    a = Styled("hello", colour="green")
    b = Styled("worldie", colour="blue")
    sep = Styled("\n", dim=True)
    out = sep.join([a, b])
    assert isinstance(out, MultiStyled)
    assert [strip_tags(str(x)) for x in out._values] == ["hello", "\n", "worldie"]
    assert isinstance(out._values[1], Styled) and out._values[1].dim is True

def test_multistyled_replace_with_empty_old_inserts_between_graphemes():
    base = MultiStyled(Styled("ab", colour="red"))
    insert = Styled("+", bold=True)
    out = base.replace("", insert)
    assert isinstance(out, MultiStyled)
    assert plain_text(out) == "+a+b+"
    assert sum(1 for chunk in out._values if chunk is insert) == 3
    core_chunks = [chunk for chunk in out._values if isinstance(chunk, Styled) and chunk is not insert]
    assert core_chunks and core_chunks[0].color == "red"

def test_multistyled_join_with_mixed_iterables_preserves_plain_text():
    sep = MultiStyled("(", Styled("|", colour="red"), ")")
    parts = [Styled("A", colour="blue"), MultiStyled("b", Styled("C", colour="green")), "D"]
    out = sep.join(parts)
    expected_plain = plain_text(sep).join(plain_text(p) for p in parts)
    assert plain_text(out) == expected_plain
    assert any(isinstance(v, Styled) and v.color == "red" for v in out._values)

# --- MultiStyled: string-method parity --------------------------------------

@pytest.mark.parametrize("method,value,args,kwargs", MULTI_PER_CHUNK_CASES)
def test_multistyled_per_chunk_methods_match_plain_string(method, value, args, kwargs):
    ms = make_multistyled(value, colour="blue")
    result = getattr(ms, method)(*args, **kwargs)
    expected = getattr(value, method)(*args, **kwargs)
    assert isinstance(result, MultiStyled)
    assert plain_text(result) == expected
    orig_colors = [chunk.color for chunk in ms._values if isinstance(chunk, Styled)]
    new_colors = [chunk.color for chunk in result._values if isinstance(chunk, Styled)]
    assert new_colors == orig_colors

@pytest.mark.parametrize(
    "method, expected_plain",
    [("strip", "core"), ("lstrip", "core  "), ("rstrip", "   core")],
)
def test_multistyled_strip_family_preserves_chunk_styles(method, expected_plain):
    ms = MultiStyled("   ", Styled("core", colour="red"), "  ")
    result = getattr(ms, method)()
    assert plain_text(result) == expected_plain
    if expected_plain.strip():
        assert any(isinstance(chunk, Styled) for chunk in result._values)

def test_multistyled_strip_all_whitespace_returns_plain_empty_string():
    ms = MultiStyled("  ", Styled(" ", colour="yellow"))
    assert ms.strip() == ""

def test_multistyled_removeprefix_and_suffix_preserve_interior_styles():
    ms = MultiStyled("##", Styled("value", colour="cyan"), "--")
    no_prefix = ms.removeprefix("##")
    assert plain_text(no_prefix) == "value--"
    assert isinstance(no_prefix._values[0], Styled)
    no_suffix = ms.removesuffix("--")
    assert plain_text(no_suffix) == "##value"
    assert isinstance(no_suffix._values[-1], Styled)

@pytest.mark.parametrize("method,args", [
    ("center", (8, ".")),
    ("ljust", (7, ".")),
    ("rjust", (7, ".")),
    ("zfill", (6,)),
])
def test_multistyled_padding_methods_add_plain_margins(method, args):
    ms = MultiStyled(Styled("core", colour="magenta"))
    result = getattr(ms, method)(*args)
    expected = getattr("core", method)(*args)
    assert plain_text(result) == expected
    styled_chunks = [chunk for chunk in result._values if isinstance(chunk, Styled)]
    assert styled_chunks and styled_chunks[0].color == "magenta"
    assert isinstance(result._values[0], str)

@pytest.mark.parametrize("method,value,args,kwargs", MULTI_SCALAR_CASES)
def test_multistyled_scalar_methods_match_plain_values(method, value, args, kwargs):
    ms = make_multistyled(value, colour="green")
    assert getattr(ms, method)(*args, **kwargs) == getattr(value, method)(*args, **kwargs)

@pytest.mark.parametrize("method,value,args,kwargs", MULTI_SEGMENT_CASES)
def test_multistyled_segment_methods_match_plain_string(method, value, args, kwargs):
    ms = make_multistyled(value, colour="orange")
    result = getattr(ms, method)(*args, **kwargs)
    expected = getattr(value, method)(*args, **kwargs)
    if method in {"partition", "rpartition"}:
        assert tuple(plain_text(part) for part in result) == expected
    else:
        assert [plain_text(part) for part in result] == expected

# --- Alignment / precision behaviour -----------------------------------------

def test_alignment_and_precision_on_styled():
    s = Styled("abcde", colour="yellow")
    # precision (truncate) then centre align
    out = f"|{s:^7.3s}|"
    inner = out[1:-1]
    visible = strip_tags(inner)
    assert visible == "  abc  "[:7]  # 3 chars centred in width 7
    # ensure only the central 3 chars are styled (our fake tags enclose just them)
    assert inner.count("<S|") == 1 and inner.count("</S>") == 1

def test_alignment_and_precision_on_multistyled_across_chunks():
    ms = MultiStyled(Styled("ab", colour="red"), Styled("cd", colour="green"), "ef")
    out = f"|{ms:>6.4s}|"   # keep first 4 codepoints: "abcd"; right-align in width 6
    inner = out[1:-1]
    assert strip_tags(inner) == "  abcd"
    # confirm resulting pieces are still two Styled chunks inside
    # (the last "ef" must be dropped by precision)
    trimmed = MultiStyled(*ms._values)  # copy to avoid mutating
    trimmed_str = f"{trimmed:.4s}"
    assert strip_tags(trimmed_str) == "abcd"

# --- Grapheme-aware precision (optional) -------------------------------------

@pytest.mark.skipif(regex is None, reason="regex not installed")
def test_grapheme_precision_on_styled(monkeypatch):
    # Only meaningful if you implemented grapheme-aware `take_graphemes` in Styled.__format__
    s = Styled("👍🏽a", colour="cyan")  # thumbs up + skin tone is one grapheme in UAX#29
    out = f"{s:.1s}"
    assert strip_tags(out) in {"👍🏽", "👍"}  # depending on your grapheme splitter

@pytest.mark.skipif(regex is None, reason="regex not installed")
def test_grapheme_precision_on_multistyled_across_chunks():
    # Assuming grapheme-aware parts truncation in MultiStyled.__format__
    a = Styled("👍🏽", colour="red")
    b = Styled("x", colour="blue")
    ms = MultiStyled(a, b)
    out = f"{ms:.1s}"
    assert strip_tags(out) in {"👍🏽", "👍"}  # first grapheme only

@pytest.mark.skipif(regex is None, reason="regex not installed")
def test_grapheme_precision_on_styled_with_multiple_clusters():
    s = Styled("👍🏽🙂x", colour="green")
    out = f"{s:.2s}"
    assert strip_tags(out) in {"👍🏽🙂", "👍🙂"}

@pytest.mark.skipif(regex is None, reason="regex not installed")
def test_grapheme_precision_on_multistyled_with_multiple_chunks():
    ms = MultiStyled(Styled("👍🏽", colour="red"), Styled("🙂", colour="blue"), "x")
    out = f"{ms:.2s}"
    assert strip_tags(out) in {"👍🏽🙂", "👍🙂"}

# --- String predicates / searches on MultiStyled ------------------------------

def test_scalar_predicates_and_searches_on_multistyled():
    ms = MultiStyled(Styled("hel", colour="red"), Styled("lo", colour="green"))
    assert ms.find("l") == "hello".find("l")
    assert ms.count("l") == 2
    assert ms.startswith("he")
    assert ms.endswith("lo")
    assert ms.isalpha() is True

# --- Error parity -------------------------------------------------------------

def test_split_with_empty_separator_raises():
    ms = MultiStyled("abc")
    with pytest.raises(ValueError):
        ms.split("")

# --- Chaining transforms ------------------------------------------------------

def test_chain_per_chunk_transforms_preserve_styles():
    ms = MultiStyled(Styled("a", colour="red"), Styled("B", colour="blue"))
    ms2 = ms.upper().lower().swapcase()
    assert isinstance(ms2, MultiStyled)
    # Styles preserved per original chunk
    assert isinstance(ms2._values[0], Styled) and ms2._values[0].color == "red"
    assert isinstance(ms2._values[1], Styled) and ms2._values[1].color == "blue"
