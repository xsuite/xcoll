#!/bin/bash

if [ $# -ne 1 ]
then
    echo "This script needs exactly one argument: the new version number or a bump."
    echo "The latter can be: patch, minor, major, prepatch, preminor, premajor, prerelease."
    exit 1
fi

poetry version $1
new_ver=$( poetry version | awk '{print $2;}' )
sed -i "s/\(__version__ =\).*/\1 '"${new_ver}"'/"         xcoll/__init__.py
sed -i "s/\(assert __version__ ==\).*/\1 '"${new_ver}"'/" tests/test_version.py
git reset
git add pyproject.toml xcoll/__init__.py tests/test_version.py
git commit -m "Updated version number to v"${new_ver}"."
git push
git tag v${new_ver}
git push origin v${new_ver}

#poetry publish --build
