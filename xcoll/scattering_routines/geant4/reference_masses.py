# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

import xtrack.particles.pdg as pdg

from ..xoconstants import Constants, constant

class XcollGeant4Masses(Constants):
    _category_ = "geant4_mass"
    _plural_   = "geant4_masses"
    _reverse_  = "unique"         # builds particle_state_names
    _c_prefix_ = "XC"

    ELECTRON_MASS_EV_GEANT4 = constant(510998.91, "Electron mass in eV as used in Geant4 (PDG ID 11).")  # best known value: 510998.95069(16) eV


def get_geant4_mass_ev(pdgid: int) -> float | None:
    vals = [vv['value'] for _, vv in XcollGeant4Masses.types_meta.items()
            if vv['info'].split('PDG ID ')[1].split(')')[0] == f"{abs(pdgid)}"]
    if len(vals) == 1:
        return vals[0]
