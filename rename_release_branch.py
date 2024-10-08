#!/usr/bin/env python
# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

import sys
from gh import *
sys.tracebacklimit = 0

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
git_push()   # Sync with the remote to be sure we don't delete an incomplete branch later
git_pull()
branch = git_current_branch()
current_ver = poetry_get_version()
if not branch == f"release/v{current_ver[:-3]}":
    raise GitError("This script needs to be ran from a release branch.")


# Check that there are no conflicting PRs open
prs = gh_pr_list(base=branch)
if prs:
    raise GitError("There are open PRs to the current release branch:\n" \
                 + "\n".join([f"PR#{pr} from {br}" for pr, br in prs.items()]) \
                 + "\nThese would be automatically closed by this script, " \
                 + "as the target branch disappears. Please close them manually, " \
                 + "or change the target branch.")
prs = gh_pr_list(head=branch)
if prs:
    raise GitError("There are open PRs from the current release branch:\n" \
                 + "\n".join([f"PR#{pr} from {br}" for pr, br in prs.items()]) \
                 + "\nThese would be automatically closed by this script, " \
                 + "as the head branch disappears. Please close or resolve them " \
                 + "manually.")


# Check that we are not accidentally bumping a major version
expected_ver = poetry_get_expected_version(bump)
major_ver = int(expected_ver.split('.')[0])
if major_ver != 0:
    raise VersionError("Bumping a major version! If this is really what you want, "
                       "then adapt make_release_branch.py manually.")


# Confirm version bump
new_branch = f"release/v{expected_ver}"
expected_ver = f"{expected_ver}rc0"
print(f"Bumping from {current_ver} to {expected_ver}.")
print("Type y to continue (or anything else to cancel):")
answer = input()
if answer not in ["y", "Y"]:
    print("Cancelled.")
    sys.exit(1)


# Rename release branch
git_rename_current_branch(new_branch, set_upstream=True)


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
git_commit(f"Renamed release branch {branch} into {new_branch}.", no_verify=True)
git_push(set_upstream=True)

print("All done!")

