#!/usr/bin/env python
# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

# HOW TO USE THIS SCRIPT
# ======================
#
# Prerequisites (only needed once):
#     - Make sure you have release rights for the package on PyPI and GitHub
#     - Install xaux and poetry (pip install xaux poetry)
#     - Install gh (GitHub Command Line interface, see https://cli.github.com/)
#     - Get an API token from PYPI (https://pypi.org/manage/account/ -> "Add API token"). Make
#       sure the package is covere by the scope of the token.
#     - Set the PyPI token with poetry: `poetry config pypi-token.pypi <my-token>`
#     - Get an API token from GitHub (github.com -> profile -> "Settings" -> "Developer Settings"
#       -> "Personal access tokens" -> "Tokens (classic)"). You will need all scopes under "repo",
#       "admin:org", "admin:public_key", and "admin:repo_hook".
#     - Set the GitHub token with gh: `gh auth login --with-token` and then paste the token.
#       In principle, this should be done for every terminal session, but an easy way around is to
#       store the token in a private file (not exposed to the network) and add
#       `gh auth login --with-token <  <path_to_file>` to your .bashrc
#
# Steps:
#     - Version numbers are of the form major.minor.patch
#     - Run `make_release_branch.py` with argument major, minor, patch, or a version number.
#       Alternatively, manually make a release branch named release/v<version_number>
#     - Make sure all changes are committed and pushed (working directory clean)
#     - In this branch, run `release.py` (again with argument major, minor, patch, or a version number).
#       This will make the release, push to main, push to PyPI, and create draft release notes on GitHub.


from xaux.dev_tools import make_release
# sys.tracebacklimit = 0


make_release("xcoll")
