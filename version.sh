#!/bin/bash
# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

# Check necessary setup and installs
poetry --version || exit 1
gh pr list || exit 1



if [ $# -eq 0 ]
then
    echo "This script needs one argument: the new version number or a bump."
    echo "The latter can be: patch, minor, major."
    exit 1
fi
if [ $# -gt 1 ] && [[ "$2" != "--force" ]]
then
    echo "The second argument can only be '--force'!"
    exit 1
fi
bump=$1
force=$2

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

# Check we are not on main
branch=$(git status | head -1 | awk '{print $NF}')
if [[ "$branch" == "main" ]]
then
    echo "Error: this script cannot be ran on the main branch."
    echo "Make a release branch and make the new release from there."
    echo "Make sure that the release branch has an upstream version (i.e."
    echo "push at least once before running this script), or this script"
    echo "will fail."
    exit 1
fi

# Check we are on the correct branch
expected_ver=$(poetry version $bump --dry-run | awk '{print $6;}')
if [[ "$force" != "--force" ]]
then
    if [[ "$branch" != "release/v$expected_ver" ]]
    then
        echo "Error: you are bumping to $expected_ver but this branch is $branch."
        echo "If this is intentional, use version.sh 'value' --force."
        exit 1
    fi
fi

# Check that we are not accidentally bumping a major version
major_ver=(${expected_ver//./ })
if [ ${major_ver[0]} -ne 0 ]
then
    echo "Error: bumping a major version! If this is really what you want, then adapt version.sh"
    exit 1
fi

# Confirm version bump
current_ver=$( poetry version | awk '{print $2;}' )
echo "Bumping from $current_ver to $expected_ver."
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

# Update version in the release branch
echo "Poetry version bump..."
poetry version $bump
new_ver=$( poetry version | awk '{print $2;}' )
if [[ "$new_ver" != "$expected_ver" ]]
then
    echo "Fatal error: poetry --dry-run expected $expected_ver, but result is $new_ver..."
    exit 1
fi
echo "Adapting version files..."
if [[ "$OSTYPE" == "darwin"* ]]
then
    sed -i '' "s/\(__version__ =\).*/\1 '"${new_ver}"'/"         xcoll/general.py
    sed -i '' "s/\(assert __version__ ==\).*/\1 '"${new_ver}"'/" tests/test_version.py
else
    sed -i'' "s/\(__version__ =\).*/\1 '"${new_ver}"'/"         xcoll/general.py
    sed -i'' "s/\(assert __version__ ==\).*/\1 '"${new_ver}"'/" tests/test_version.py
fi
echo "Committing version change..."
git reset
git add pyproject.toml xcoll/general.py tests/test_version.py
git commit --no-verify -m "Updated version number to v"${new_ver}"."
git push

# Make and accept pull request
echo "Creating and merging pull request..."
gh pr create --base main --title "Release "${new_ver} --fill
git switch main
git pull
PR=$( gh pr list --head $branch | tail -1 | awk '{print $1;}' )
gh pr merge ${PR} --merge --admin --delete-branch
git pull

# Make a tag
git tag v${new_ver}
git push origin v${new_ver}

# Draft release
gh release create v${new_ver} --draft --generate-notes --title "Xcoll release ${new_ver}" --verify-tag

# Build release and publish to PyPi
poetry publish --build

echo "All done!"

