#!/usr/bin/env python
# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

import sys
import platform
from gh import *
# sys.tracebacklimit = 0

class VersionError(OSError):
    pass


# Check necessary setup and installs
assert_git_repo()
assert_poetry_installed()


# Check the script arguments
if len(sys.argv) != 2:
    raise ValueError("This script needs exactly one argument: the new version number or a bump. "
                     "The latter can be: patch, minor, major.")
bump = sys.argv[1]


# Check our working directory is clean
git_assert_working_tree_clean()
branch = git_current_branch()
if branch != "main":
    raise GitError("This script needs to be ran on the main branch.")
git_push()   # Sync with the remote to be sure we don't delete an incomplete branch later
git_pull()
expected_ver = poetry_get_expected_version(bump)


# Check that we are not accidentally bumping a major version
major_ver = int(expected_ver.split('.')[0])
if major_ver != 0:
    raise VersionError("Bumping a major version! If this is really what you want, "
                       "then adapt make_release_branch.py manually.")


# Confirm version bump
branch = f"release/v{expected_ver}"
expected_ver = f"{expected_ver}rc0"
current_ver = poetry_get_version()
print(f"Bumping from {current_ver} to {expected_ver}.")
print("Type y to continue (or anything else to cancel):")
answer = input()
if answer not in ["y", "Y"]:
    print("Cancelled.")
    sys.exit(1)


# Create release branch
print("Creating release branch...")
git_switch(branch, create=True)


# Update version in the release branch
print("Poetry version bump...")
poetry_bump_version(expected_ver)
new_ver = poetry_get_version()
if new_ver != expected_ver:
    raise VersionError(f"Fatal error: `poetry --dry-run` expected {expected_ver}, but result is {new_ver}..."
                        "Need to recover manually!")


# Adapt version files
for file, pattern in zip(["xcoll/general.py", "tests/test_version.py"],
                         ["__version__ = ", "    assert __version__ == "]):
    file_adapted = False
    with Path(file).open("r") as fid:
        lines = fid.readlines()
    with Path(file).open("w") as fid:
        for line in lines:
            if line.startswith(pattern):
                fid.write(f"{pattern}'{new_ver}'\n")
                file_adapted = True
            else:
                fid.write(line)
    if not file_adapted:
        raise VersionError(f"Fatal error: could not adapt {file}...")


# Commit and push
git_add(['pyproject.toml', 'xcoll/general.py', 'tests/test_version.py'])
git_commit(f"Created release branch release/v{new_ver}.", no_verify=True)
git_push(set_upstream=True)

print("All done!")

