# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

import os
import re
import sys
from typing import Tuple, Iterable, Union


# ---------- ANSI basics ----------
_RESET = "\033[0m"

def _sgr(*codes: Union[int, str]) -> str:
    return f"\033[{';'.join(str(c) for c in codes)}m"

# Text attributes (SGR)
_SGR_CODES = {
    "bold": 1,
    "dim": 2,
    "italic": 3,
    "underline": 4,
}

# ANSI named colours: normal + bright
# You can extend this dict as you wish.
_NAMED = {
    "black": 0, "red": 1, "green": 2, "yellow": 3, "blue": 4, "magenta": 5, "cyan": 6, "white": 7,
    "bright_black": 8, "bright_red": 9, "bright_green": 10, "bright_yellow": 11,
    "bright_blue": 12, "bright_magenta": 13, "bright_cyan": 14, "bright_white": 15,
}

# Add friendly aliases mapped to 256-colour indices (xterm palette)
# (These go through the 38;5;<index> pathway.)
_ALIAS_256 = {
    # ---- Greys (common names) ----
    "grey": 8, "gray": 8,
    "dark_grey": 240, "dark_gray": 240,
    "slate_grey": 243, "slate_gray": 243,
    "light_grey": 250, "light_gray": 250,
    "gainsboro": 252,
    "silver": 251,
    "whitesmoke": 255,

    # Fine greyscale steps (darker -> lighter)
    "grey_0": 232, "grey_1": 233, "grey_2": 234, "grey_3": 235,
    "grey_4": 236, "grey_5": 237, "grey_6": 238, "grey_7": 239,
    "grey_8": 240, "grey_9": 241, "grey_10": 242, "grey_11": 243,
    "grey_12": 244, "grey_13": 245, "grey_14": 246, "grey_15": 247,
    "grey_16": 248, "grey_17": 249, "grey_18": 250, "grey_19": 251,
    "grey_20": 252, "grey_21": 253, "grey_22": 254, "grey_23": 255,
    "gray_0": 232, "gray_1": 233, "gray_2": 234, "gray_3": 235,
    "gray_4": 236, "gray_5": 237, "gray_6": 238, "gray_7": 239,
    "gray_8": 240, "gray_9": 241, "gray_10": 242, "gray_11": 243,
    "gray_12": 244, "gray_13": 245, "gray_14": 246, "gray_15": 247,
    "gray_16": 248, "gray_17": 249, "gray_18": 250, "gray_19": 251,
    "gray_20": 252, "gray_21": 253, "gray_22": 254, "gray_23": 255,

    # ---- Reds ----
    "maroon": 88, "firebrick": 124, "crimson": 160, "red2": 196,
    "tomato": 203, "salmon": 210, "light_coral": 210, "indian_red": 167,
    "rosy_brown": 138, "brown": 94,

    # ---- Oranges / Browns ----
    "dark_orange": 166, "orange": 208, "coral": 209, "sandy_brown": 215,
    "peru": 173, "chocolate": 166, "sienna": 130, "tan": 180,
    "burlywood": 180, "peach": 216,

    # ---- Yellows / Golds ----
    "gold": 220, "goldenrod": 178, "dark_khaki": 143, "khaki": 186,
    "lemon": 226, "light_yellow": 229, "wheat": 229, "cornsilk": 230,

    # ---- Greens ----
    "dark_green": 22, "forest_green": 22, "olive": 58, "olive_drab": 64,
    "green4": 28, "sea_green": 29, "dark_sea_green": 108, "spring_green": 48,
    "lawn_green": 118, "chartreuse": 118, "pale_green": 120,
    "medium_spring_green": 49, "yellow_green": 190, "mint_green": 121,

    # ---- Teals / Cyans ----
    "teal": 30, "dark_cyan": 30, "cadet_blue": 109, "light_sea_green": 37,
    "turquoise": 45, "medium_turquoise": 44, "aquamarine": 86,
    "cyan2": 51, "light_cyan": 195,

    # ---- Blues ----
    "navy": 17, "midnight_blue": 17, "dark_blue": 18, "blue2": 21,
    "royal_blue": 63, "dodger_blue": 33, "deep_sky_blue": 39,
    "cornflower_blue": 69, "steel_blue": 67, "light_steel_blue": 147,
    "sky_blue": 117,

    # ---- Purples / Violets ----
    "indigo": 54, "dark_violet": 92, "blue_violet": 93,
    "slate_blue": 62, "medium_purple": 104,
    "orchid": 170, "plum": 176, "violet": 177, "thistle": 182,

    # ---- Pinks / Magentas ----
    "magenta2": 201, "fuchsia": 201, "deep_pink": 199,
    "hot_pink": 205, "pink": 211, "light_pink": 218, "pale_violet_red": 168,

    # ---- Earth / Neutrals ----
    "beige": 230, "antique_white": 230, "linen": 255, "ivory": 255,
}



# Regex for stripping ANSI
_ANSI_RE = re.compile(r"\x1B\[[0-?]*[ -/]*[@-~]")

def strip_ansi(s: str) -> str:
    """Remove ANSI escape codes from a string."""
    return _ANSI_RE.sub("", s)


# ---------- Colour encoders ----------
def _encode_colour(
    colour: Union[str, int, Tuple[int, int, int], None], *, bg: bool = False
) -> Iterable[int]:
    """
    Return SGR code sequence (as ints) for a colour.

    Accepts:
      - Named colour: "red", "bright_blue", ...
      - 256-colour index: int 0..255
      - True-colour: (r,g,b) 0..255
      - Hex string: "#RRGGBB" or "RRGGBB"
      - None: no colour
    """
    if colour is None:
        return ()

    if isinstance(colour, str):
        s = colour.lower().strip()
        if s.startswith("#"):
            s = s[1:]
        if len(s) == 6 and all(c in "0123456789abcdef" for c in s):
            r = int(s[0:2], 16); g = int(s[2:4], 16); b = int(s[4:6], 16)
            return (48 if bg else 38, 2, r, g, b)
        if s in _NAMED:
            idx = _NAMED[s]
            # map 0..7 to 30..37 (fg) / 40..47 (bg), 8..15 to 90..97 / 100..107
            if idx < 8:
                base = 40 if bg else 30
                return (base + idx,)
            else:
                base = 100 if bg else 90
                return (base + (idx - 8),)
        if s in _ALIAS_256:
            return (48 if bg else 38, 5, _ALIAS_256[s])
        raise ValueError(f"Unknown colour name: {colour!r}")

    elif isinstance(colour, int):
        if not (0 <= colour <= 255):
            raise ValueError("256-colour index must be in 0..255")
        return (48 if bg else 38, 5, colour)

    elif isinstance(colour, tuple) and len(colour) == 3:
        r, g, b = colour
        for c in (r, g, b):
            if not (0 <= int(c) <= 255):
                raise ValueError("RGB components must be in 0..255")
        return (48 if bg else 38, 2, int(r), int(g), int(b))

    raise TypeError("Colour must be str|int|tuple[int,int,int]|None")


# ---------- Style core ----------
def style(
    text: str,
    *,
    color: Union[str, int, Tuple[int, int, int], None] = None,
    colour: Union[str, int, Tuple[int, int, int], None] = None,
    bg: Union[str, int, Tuple[int, int, int], None] = None,
    bold: bool = False,
    italic: bool = False,
    underline: bool = False,
    dim: bool = False,
    enabled: bool | None = None,
    reset: bool = True,
) -> str:
    """
    Wrap `text` with ANSI styles.

    - colo(u)r/bg: named ("red", "bright_blue"), 0..255, (r,g,b), or "#RRGGBB"
    - toggles: bold, italic, underline, dim
    - enabled: override auto-detection. If None, enable only when stdout is a TTY and NO_COLOR not set.
    - reset: append a trailing reset (recommended when printing standalone fragments)
    """
    if enabled is None:
        enabled = sys.stdout.isatty() and "NO_COLOR" not in os.environ

    if not enabled:
        return text

    codes: list[int] = []
    if colour is not None:
        if color is not None:
            raise ValueError("Cannot specify both 'color' and 'colour'")
        color = colour
    codes += list(_encode_colour(color, bg=False))
    codes += list(_encode_colour(bg, bg=True))
    for flag, name in ((bold, "bold"), (dim, "dim"), (italic, "italic"), (underline, "underline")):
        if flag:
            codes.append(_SGR_CODES[name])

    if not codes:
        return text

    prefix = _sgr(*codes)
    return f"{prefix}{text}{_RESET if reset else ''}"
