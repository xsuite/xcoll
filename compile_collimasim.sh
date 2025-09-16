#!/bin/bash

path=xcoll/lib
collimasim_path=$path/collimasim

# conda install -c conda-forge compilers clang clangxx clang-tools ninja make cmake bdsim-g4
# pip install rpyc

# If compilation fails, one can try to do (before running this script):
# export CMAKE_GENERATOR="Unix Makefiles"

echo "Preparing $collimasim_path"
python -m pip install --user --editable $collimasim_path || (echo; echo; echo "Failed compiling collimasim!"; echo)
