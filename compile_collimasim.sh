#!/bin/bash

set -euo pipefail

python - <<'PY'
from xcoll.scattering_routines.geant4.environment import Geant4Environment

env = Geant4Environment()
print("Compiling Geant4 interface ...")
env.compile(force=True)
print("Compilation completed. Built module located in:", env.collimasim)
PY
