#!/bin/bash

path=xcoll/lib
collimasim_path=$path/collimasim

# git submodule update --init --recursive

echo "Preparing $collimasim_path"
python -m pip install --user --editable $collimasim_path || (echo; echo; echo "Failed compiling collimasim!"; echo)
