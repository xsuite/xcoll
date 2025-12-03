# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

import xtrack.particles.pdg as pdg

from ...xoconstants import Constants, constant


class Geant4Masses(Constants):
    _category_ = "geant4_mass"
    _plural_   = "geant4_masses"
    _c_prefix_ = "XC_GEANT4"

    ELECTRON_MASS_EV  = constant(0.51099891e6,   "Electron mass in eV as used in Geant4.",  pdg_id=11)    # Best known value: 510998.95069(16) eV
    MUON_MASS_EV      = constant(105.6583715e6,  "Muon mass in eV as used in Geant4.",      pdg_id=13)

    PION_MASS_EV      = constant(139.5701e6,     "Pi+ mass in eV as used in Geant4.",       pdg_id=211)
    KAON_MASS_EV      = constant(493.677e6,      "K+ mass in eV as used in Geant4.",        pdg_id=321)
    D_MASS_EV         = constant(1869.58e6,      "D+ mass in eV as used in Geant4.",        pdg_id=411)
    Ds_MASS_EV        = constant(1968.27e6,      "Ds+ mass in eV as used in Geant4.",       pdg_id=431)
    B_MASS_EV         = constant(5279.29e6,      "B+ mass in eV as used in Geant4.",        pdg_id=521)

    PROTON_MASS_EV    = constant(938.272013e6,   "Proton mass in eV as used in Geant4.",    pdg_id=2212)  # Best known value: 938272088.16 eV
    SIGMA_POS_MASS_EV = constant(1189.37e6,      "Sigma+ mass in eV as used in Geant4.",    pdg_id=3222)
    SIGMA_NEG_MASS_EV = constant(1197.449e6,     "Sigma- mass in eV as used in Geant4.",    pdg_id=3112)
    XI_MASS_EV        = constant(1321.71e6,      "Xi- mass in eV as used in Geant4.",       pdg_id=3312)
    OMEGA_MASS_EV     = constant(1672.45e6,      "Omega- mass in eV as used in Geant4.",    pdg_id=3334)
    LAMBDAc_MASS_EV   = constant(2286.46e6,      "Lambda_c+ mass in eV as used in Geant4.", pdg_id=4122)
    XIc_MASS_EV       = constant(2467.71e6,      "Xi_c+ mass in eV as used in Geant4.",     pdg_id=4232)
    XIb_MASS_EV       = constant(5794.5e6,       "Xi_c- mass in eV as used in Geant4.",     pdg_id=5132)

    DEUTERON_MASS_EV  = constant(1.875613e9,     "H2 mass in eV as used in Geant4.",        pdg_id=1000010020)
    TRITON_MASS_EV    = constant(2.808921e9,     "H3 mass in eV as used in Geant4.",        pdg_id=1000010030)
    HE3_MASS_EV       = constant(2.808391e9,     "He3 mass in eV as used in Geant4.",       pdg_id=1000020030)
    HE4_MASS_EV       = constant(3.727379e9,     "He4 mass in eV as used in Geant4.",       pdg_id=1000020040)


class Geant4MassesAccessor:
    """Accessor to get Geant4 reference masses."""

    def __getitem__(self, pdgid: int | str) -> float | None:
        if isinstance(pdgid, str):
            pdgid = pdg.get_pdg_id_from_name(pdgid)
        vals = [vv['value'] for _, vv in geant4_masses_meta.items()
                if vv['pdg_id'] == pdgid]
        if len(vals) == 1:
            return vals[0]
        vals = [vv['value'] for _, vv in geant4_masses_meta.items()
                if vv['pdg_id'] == -pdgid]
        if len(vals) == 1:
            return vals[0]

    def __contains__(self, pdgid: int | str) -> bool:
        return self[pdgid] is not None
