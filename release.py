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
assert_gh_installed()


# Check the script arguments
if len(sys.argv) < 2 or len(sys.argv) > 3:
    raise ValueError("This script needs exactly one argument: the new version number or a bump. "
                     "The latter can be: patch, minor, major.")
force = False
bump = sys.argv[1]
if len(sys.argv) == 3:
    force = True
    if sys.argv[1] == "--force":
        bump = sys.argv[2]
    elif sys.argv[2] != "--force":
        raise ValueError("Only '--force' is allowed as an option.")


# Check our working directory is clean
print("Verifying repository status...")
git_assert_working_tree_clean()
branch = git_current_branch()
if branch == "main":
    raise GitError("\nThis script cannot be ran on the main branch."
                      "Make a release branch and make the new release from there."
                      "Make sure that the release branch has an upstream version (i.e."
                      "push at least once before running this script), or this script"
                      "will fail.")
expected_ver = poetry_get_expected_version(bump)
if not force:
    if branch != f"release/v{expected_ver}":
        raise VersionError(f"\nYou are bumping to {expected_ver} but this branch is {branch}. "
                                "If this is intentional, use --force.")
git_push()   # Sync with the remote to be sure we don't delete an incomplete branch later
git_pull()
print("Repository is clean.")


# Check that there are no conflicting PRs open
prs = gh_pr_list(base=branch)
if prs:
    raise GitError(f"There are open PRs to the release branch:\n" \
                   + "\n".join([f"PR#{pr} from {br}" for pr, br in prs.items()]) \
                   + "\nThese would be automatically closed by this script, " \
                   + "as the target branch disappears. Please close them manually, " \
                   + "or change the target branch.")
prs = gh_pr_list(base='main', head=branch)
if prs:
    raise GitError(f"There are open PRs from the release branch to main:\n" \
                   + "\n".join([f"PR#{pr} from {br}" for pr, br in prs.items()]) \
                   + "\nThese would conflict with the versioning script. "
                   + "Please close them manually.")


# Check that we are not accidentally bumping a major version
major_ver = int(expected_ver.split('.')[0])
if major_ver != 0:
    raise VersionError("Bumping a major version! If this is really what you want, "
                       "then adapt version.py manually.")


# Confirm version bump
current_ver = poetry_get_version()
print(f"Bumping from {current_ver} to {expected_ver}.")
print("Type y to continue (or anything else to cancel):")
answer = input()
if answer not in ["y", "Y"]:
    print("Cancelled.")
    sys.exit(1)


print("Updating version in the release branch...")
# Make the bump
poetry_bump_version(bump)
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
git_commit(f"Updated version number to v{new_ver}.", no_verify=True)
git_push()


print("Creating and merging pull request to main branch...")
gh_pr_create('main', f"Release {new_ver}")
git_switch('main')
git_pull()
prs = gh_pr_list(base='main', head=branch)
if len(prs) != 1:
    raise GitError(f"Expected one PR from {branch} to main, found {len(prs)}:\n"
                  + "\n".join([f"PR#{pr} from {br}" for pr, br in prs.items()]))
gh_pr_merge(list(prs.keys())[0], admin=True, delete_branch=True)
git_pull()
git_make_tag(f"v{new_ver}")


print("Creating draft release and publishing to PyPi...")
gh_release_create(f"v{new_ver}", f"Xcoll release {new_ver}", draft=True)
poetry_publish(build=True)

print("All done!")
