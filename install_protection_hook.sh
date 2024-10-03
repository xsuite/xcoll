#!/bin/bash

# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2023.                 #
# ######################################### #


cat << 'EOF' > pre-commit
#!/bin/sh

branch="$(git rev-parse --abbrev-ref HEAD)"

if [ "$branch" = "main" ]
then
  echo "Cannot commit to main!"
  echo "Use pull requests or, if package admin, the versioning script."
  exit 1
fi

for file in pyproject.toml gh.py version.py make_release_branch.py LICENSE install_protection_hook.sh xcoll/general.py
do
  git diff --name-only | grep '^'${file}'$' &> /dev/null
  if [ $? -eq 0 ]
  then
    thisfile=${file/.\*/\*}
    echo "File $thisfile is protected but has local changes."
    echo "Restore the file with 'git restore "${thisfile}"' before commiting anything,"
    echo "or force this commit with 'git commit --no-verify'."
    exit 1
  fi
  git diff --cached --name-only | grep '^'${file}'$' &> /dev/null
  if [ $? -eq 0 ]
  then
    thisfile=${file/.\*/\*}
    echo "File $thisfile is protected but has local changes that are staged."
    echo "First unstage the file with 'git restore --staged "${thisfile}"', then "
    echo "restore the file with 'git restore "${thisfile}"' before commiting anything."
    echo "Alternatively, force this commit with 'git commit --no-verify'."
    exit 1
  fi
done
EOF

chmod +x pre-commit
mv pre-commit .git/hooks/
