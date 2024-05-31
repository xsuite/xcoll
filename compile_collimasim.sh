# Requires BDSIM to be installed
# and available to CMake. For install instructions:
# https://www.pp.rhul.ac.uk/bdsim/manual/
#
# using bash (not zsh) on CentOS:
# source /cvmfs/beam-physics.cern.ch/bdsim/x86_64-centos7-gcc11-opt/bdsim-env-v1.7.6-g4v10.7.2.3-ftfp-boost.sh
# source /cvmfs/beam-physics.cern.ch/bdsim/x86_64-centos7-gcc11-opt/bdsim-v1.7.6-g4v10.7.2.3-ftfp-boost/env-for-build.sh
# source /cvmfs/sft.cern.ch/lcg/releases/LCG_104/ninja/1.10.0/x86_64-centos7-gcc11-opt/ninja-env.sh
# source /cvmfs/sft.cern.ch/lcg/releases/LCG_104/jsonmcpp/3.10.5/x86_64-centos7-gcc11-opt/jsonmcpp-env.sh
# export nlohmann_json_DIR=/cvmfs/sft.cern.ch/lcg/releases/LCG_104/jsonmcpp/3.10.5/x86_64-centos7-gcc11-opt/lib64/cmake/nlohmann_json/
# unset PYTHONPATH
# unset PYTHONHOME
# if [ ! -z "$CONDA_PREFIX" ]; then export PATH="$CONDA_PREFIX/bin":$PATH; fi
# GITCOMMAND=/usr/bin/git


path=xcoll/scattering_routines/geant4

collimasim_path=$path/collimasim

git submodule update --init --recursive

echo "Preparing $collimasim_path"
python -m pip install --user --editable $collimasim_path || (echo; echo; echo "Failed installing collimasim. Make sure BDSIM is installed correctly (see https://www.pp.rhul.ac.uk/bdsim/manual/)"; echo)
