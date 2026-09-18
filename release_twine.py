#!/usr/bin/env python
# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

"""
Poetry-free release script for xcoll (alternative to release.py, using
`python -m build` + `twine`). Leaves the repository in the state expected by
the poetry-based flow: version bumped in pyproject.toml, xcoll/general.py and
tests/test_version.py, one commit pushed, annotated tag vX.Y.Z pushed,
sdist + wheel uploaded to PyPI. GitHub release notes are done by hand.

Prerequisites: pip install build twine; twine credentials (~/.pypirc or
TWINE_USERNAME=__token__ / TWINE_PASSWORD=<token>).

Usage (from the xcoll root directory, on a clean and pushed branch):
    python release_twine.py           # prompts for the new version
    python release_twine.py 0.12.5    # or give it directly
"""

import re
import subprocess
import sys
import tempfile
import tomllib
from pathlib import Path

PACKAGE = "xcoll"
VERSION_FILES = {
    Path("pyproject.toml"): "version = ",
    Path(f"{PACKAGE}/general.py"): "__version__ = ",
    Path("tests/test_version.py"): "    assert __version__ == ",
}


def run(*cmd, capture=False):
    print(f"  $ {' '.join(cmd)}")
    res = subprocess.run(cmd, text=True, check=True,
                         stdout=subprocess.PIPE if capture else None)
    return res.stdout.strip() if capture else None


def fail(msg):
    sys.exit(f"\nERROR: {msg}")


def confirm(question):
    print(f"{question} Type y to continue (or anything else to cancel):")
    if input().strip() not in ("y", "Y"):
        sys.exit("Cancelled.")


def set_version(new_ver):
    for path, prefix in VERSION_FILES.items():
        quote = '"' if path.suffix == ".toml" else "'"
        lines = path.read_text().splitlines(keepends=True)
        hits = [i for i, l in enumerate(lines) if l.startswith(prefix)]
        if len(hits) != 1:
            fail(f"Expected exactly one line starting with {prefix!r} in {path}")
        lines[hits[0]] = f"{prefix}{quote}{new_ver}{quote}\n"
        path.write_text("".join(lines))


def main(argv):
    # Where are we, and is everything clean?
    if not Path("pyproject.toml").is_file() or not Path(PACKAGE).is_dir():
        fail("Run this script from the xcoll root directory.")
    with open("pyproject.toml", "rb") as fid:
        current = tomllib.load(fid)["project"]["version"]
    for mod in ("build", "twine"):
        if subprocess.run([sys.executable, "-m", mod, "--help"],
                          capture_output=True).returncode != 0:
            fail(f"`{mod}` is not installed: pip install build twine")
    if run("git", "status", "--porcelain", "--untracked-files=no", capture=True):
        fail("Working tree has uncommitted changes.")
    run("git", "fetch", "origin")
    branch = run("git", "branch", "--show-current", capture=True)
    if run("git", "rev-list", "--left-right", "--count", f"{branch}...@{{u}}",
           capture=True) != "0\t0":
        fail(f"Branch {branch} is not in sync with its upstream. Pull/push first.")

    # New version number
    new_ver = argv[0] if argv else None
    while True:
        if new_ver is None:
            print(f"Current version is {current}. Type the new version (X.Y.Z), "
                  "or leave empty to cancel:")
            new_ver = input().strip() or sys.exit("Cancelled.")
        if not re.fullmatch(r"0\.\d+\.\d+", new_ver):
            print(f"  Invalid version {new_ver!r}: expected 0.Y.Z")
        elif new_ver == current:
            print(f"  {new_ver} is already the current version.")
        elif (run("git", "tag", "-l", f"v{new_ver}", capture=True)
              or run("git", "ls-remote", "--tags", "origin", f"v{new_ver}", capture=True)):
            print(f"  Tag v{new_ver} already exists.")
        else:
            break
        new_ver = None

    # Bump, commit, push, tag
    confirm(f"Bumping {PACKAGE} from {current} to {new_ver} on branch {branch}.")
    set_version(new_ver)
    run("git", "--no-pager", "diff")
    confirm("Commit, push, tag and upload to PyPI?")
    run("git", "add", *map(str, VERSION_FILES))
    run("git", "commit", "--no-verify", "-m", f"Updated version number to v{new_ver}.")
    run("git", "push")
    run("git", "tag", "-a", f"v{new_ver}", "-m", f"{PACKAGE.capitalize()} release {new_ver}")
    run("git", "push", "origin", f"v{new_ver}")

    # Build and upload
    with tempfile.TemporaryDirectory(prefix=f"{PACKAGE}_release_") as outdir:
        run(sys.executable, "-m", "build", "--outdir", outdir)
        files = [str(f) for f in sorted(Path(outdir).iterdir())]
        run(sys.executable, "-m", "twine", "check", *files)
        run(sys.executable, "-m", "twine", "upload", *files)

    print(f"\nAll done! Now create the release notes on GitHub:\n"
          f"  https://github.com/xsuite/{PACKAGE}/releases/new?tag=v{new_ver}")


if __name__ == "__main__":
    try:
        main(sys.argv[1:])
    except subprocess.CalledProcessError as e:
        fail(f"command failed: {' '.join(e.cmd)}")
