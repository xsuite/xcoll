#!/bin/bash

if [ $# -ne 1 ]
then
    echo "This script needs exactly one argument: the new version number or a bump."
    echo "The latter can be: patch, minor, major, prepatch, preminor, premajor, prerelease."
    exit 1
fi

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

# Update version in temporary branch
poetry version $1
new_ver=$( poetry version | awk '{print $2;}' )
sed -i "s/\(__version__ =\).*/\1 '"${new_ver}"'/"         xcoll/__init__.py
sed -i "s/\(assert __version__ ==\).*/\1 '"${new_ver}"'/" tests/test_version.py
git reset
git add pyproject.toml xcoll/__init__.py tests/test_version.py
git commit -m "Updated version number to v"${new_ver}"."
git push

# Make and accept pull request
gh pr create --base main --title "Release "${new_ver} --fill
git switch main
git pull
PR=$( gh pr list --head $branch | tail -1 | awk '{print $1;}' )
gh pr merge ${PR} --merge --admin --delete-branch
git pull

# Make a tag
git tag v${new_ver}
git push origin v${new_ver}

# Draft release notes
curl \
  -X POST \
  -H "Accept: application/vnd.github+json" \
  -H "Authorization: Bearer "$(cat ../github_token) \
  https://api.github.com/repos/xsuite/xcoll/releases \
  -d '{"tag_name":"v'${new_ver}'","target_commitish":"main","name":"Xcoll release '${new_ver}'","body":"","draft":true,"prerelease":false,"generate_release_notes":true}'

# Build release and publish to PyPi
poetry publish --build
