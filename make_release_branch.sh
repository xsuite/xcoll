#!/bin/bash
# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

# Check necessary setup and installs
poetry --version || exit 1

if [ $# -eq 0 ]
then
    echo "This script needs one argument: the new version number or a bump."
    echo "The latter can be: patch, minor, major."
    exit 1
fi
bump=$1

# Check our working directory is clean
git diff --quiet
if [ $? -ne 0 ]
then
    echo "Some changes are not staged."
    echo "Stage and commit before bumping the version."
    exit 1
fi
git diff --staged --quiet
if [ $? -ne 0 ]
then
    echo "Some changes are staged."
    echo "Commit before bumping the version."
    exit 1
fi

# Check we are on main
branch=$(git status | head -1 | awk '{print $NF}')
if [[ "$branch" != "main" ]]
then
    echo "Error: this script needs to be ran on the main branch."
    exit 1
fi

# Get the release version
expected_ver=$(poetry version $bump --dry-run | awk '{print $6;}')

# Check that we are not accidentally bumping a major version
major_ver=(${expected_ver//./ })
if [ ${major_ver[0]} -ne 0 ]
then
    echo "Error: bumping a major version! If this is really what you want, then adapt version.sh"
    exit 1
fi

# Confirm version bump
branch=release/v$expected_ver
expected_ver=${expected_ver}rc0
current_ver=$( poetry version | awk '{print $2;}' )
echo "Bumping from $current_ver to $expected_ver"
read -n 1 -p "Press y to continue (or any other key to cancel): " answer
case ${answer:0:1} in
    y|Y )
        echo
    ;;
    * )
        echo "Cancelled."
        exit 1
    ;;
esac

# Kill script on first error
set -e

# Create release branch
echo "Creating release branch..."
git switch -c $branch

# Update version in the release branch
echo "Poetry version bump..."
poetry version $expected_ver
new_ver=$( poetry version | awk '{print $2;}' )
if [[ "$new_ver" != "$expected_ver" ]]
then
    echo "Fatal error: poetry --dry-run expected $expected_ver, but result is $new_ver..."
    exit 1
fi
echo "Adapting version files..."
sed -i "s/\(__version__ =\).*/\1 '"${new_ver}"'/"         xcoll/general.py
sed -i "s/\(assert __version__ ==\).*/\1 '"${new_ver}"'/" tests/test_version.py
echo "Committing version change..."
git reset
git add pyproject.toml xcoll/general.py tests/test_version.py
git commit --no-verify -m "Created release branch release/v"${new_ver}"."
git push --set-upstream origin $branch

echo "All done!"

