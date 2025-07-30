# Collimasim

Collimasim provides Python bindings to BDSIM (Geant4) in order to enable collimation studies in pure tracking codes

## Requiremends
 - BDSIM installed on the system
 - python3
 - CMake >3.17

## Installation

To isntall use:

```bash
git clone --recurse-submodules https://gitlab.cern.ch/anabramo/collimasim.git
python -m pip install --user --editable ./
```

## Usage

See tests/test.py for an example of how to use