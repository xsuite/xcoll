# Requires BDSIM to be installed
# and available to CMake. For install instructions: 
# https://www.pp.rhul.ac.uk/bdsim/manual/

path=xcoll/scattering_routines/geant4

collimasim_path=$path/collimasim

git submodule update --init --recursive

echo "Preparing $collimasim_path"
python -m pip install --user --editable $collimasim_path