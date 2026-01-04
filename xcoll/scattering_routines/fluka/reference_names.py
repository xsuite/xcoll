# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

from ...xoconstants import Constants, constant


class FlukaNames(Constants):
    _category_ = "fluka_name"
    _c_prefix_ = "XC_FLUKA"

    PHOTON_NAME                  = constant("PHOTON",    "Photon name as used in FLUKA.",                 pdg_id=22)

    ELECTRON_NAME                = constant("ELECTRON",  "Electron name as used in FLUKA.",               pdg_id=11)
    ELECTRON_NEUTRINO_NAME       = constant("NEUTRIE",   "Electron-neutrino name as used in FLUKA.",      pdg_id=12)
    MUON_NAME                    = constant("MUON-",     "Muon name as used in FLUKA.",                   pdg_id=13)
    MUON_NEUTRINO_NAME           = constant("NEUTRIM",   "Muon-neutrino name as used in FLUKA.",          pdg_id=14)
    TAU_NAME                     = constant("TAU-",      "Tau name as used in FLUKA.",                    pdg_id=15)
    TAU_NEUTRINO_NAME            = constant("NEUTRIT",   "Tau-neutrino name as used in FLUKA.",           pdg_id=16)

    POSITRON_NAME                = constant("POSITRON",  "Positron name as used in FLUKA.",               pdg_id=-11)
    ANTI_ELECTRON_NEUTRINO_NAME  = constant("ANEUTRIE",  "Anti-electron-neutrino name as used in FLUKA.", pdg_id=-12)
    ANTI_MUON_NAME               = constant("MUON+",     "Anti-muon name as used in FLUKA.",              pdg_id=-13)
    ANTI_MUON_NEUTRINO_NAME      = constant("ANEUTRIM",  "Anti-muon-neutrino name as used in FLUKA.",     pdg_id=-14)
    ANTI_TAU_NAME                = constant("TAU+",      "Anti-tau name as used in FLUKA.",               pdg_id=-15)
    ANTI_TAU_NEUTRINO_NAME       = constant("ANEUTRIT",  "Anti-tau-neutrino name as used in FLUKA.",      pdg_id=-16)

    PROTON_NAME                  = constant("PROTON",    "Proton name as used in FLUKA.",                 pdg_id=2212)
    ANTI_PROTON_NAME             = constant("APROTON",   "Anti-proton name as used in FLUKA.",            pdg_id=-2212)
    NEUTRON_NAME                 = constant("NEUTRON",   "Neutron name as used in FLUKA.",                pdg_id=2112)
    ANTI_NEUTRON_NAME            = constant("ANEUTRON",  "Anti-neutron name as used in FLUKA.",           pdg_id=-2112)
    TRITON_NAME                  = constant("TRITON",    "Triton name as used in FLUKA.",                 pdg_id=1000010030)
    DEUTERON_NAME                = constant("DEUTERON",  "Deuteron name as used in FLUKA.",               pdg_id=1000010020)
    HELIUM3_NAME                 = constant("3-HELIUM",  "He3 name as used in FLUKA.",                    pdg_id=1000020030)
    HELIUM4_NAME                 = constant("4-HELIUM",  "He4 name as used in FLUKA.",                    pdg_id=1000020040)

    PION_NAME                    = constant("PIZERO",    "π0 name as used in FLUKA.",                     pdg_id=111)
    PION_POS_NAME                = constant("PION+",     "π+ name as used in FLUKA.",                     pdg_id=211)
    PION_NEG_NAME                = constant("PION-",     "π- name as used in FLUKA.",                     pdg_id=-211)
    KAON_NAME                    = constant("KAONZERO",  "K0 name as used in FLUKA.",                     pdg_id=311)
    ANTI_KAON_NAME               = constant("AKAONZER",  "Anti-K0 name as used in FLUKA.",                pdg_id=-311)
    KAON_POS_NAME                = constant("KAON+",     "K+ name as used in FLUKA.",                     pdg_id=321)
    KAON_NEG_NAME                = constant("KAON-",     "K- name as used in FLUKA.",                     pdg_id=-321)
    KAON_LONG_NAME               = constant("KAONLONG",  "K long name as used in FLUKA.",                 pdg_id=130)
    KAON_SHORT_NAME              = constant("KAONSHRT",  "K short name as used in FLUKA.",                pdg_id=310)
    D_NAME                       = constant("D0",        "D0 name as used in FLUKA.",                     pdg_id=421)
    ANTI_D_NAME                  = constant("D0BAR",     "Anti-D0 name as used in FLUKA.",                pdg_id=-421)
    D_POS_NAME                   = constant("D+",        "D+ name as used in FLUKA.",                     pdg_id=411)
    D_NEG_NAME                   = constant("D-",        "D- name as used in FLUKA.",                     pdg_id=-411)
    D_S_POS_NAME                 = constant("DS+",       "Ds+ name as used in FLUKA.",                    pdg_id=431)
    D_S_NEG_NAME                 = constant("DS-",       "Ds- name as used in FLUKA.",                    pdg_id=-431)

    LAMBDA_NAME                  = constant("LAMBDA",    "Λ0 name as used in FLUKA.",                     pdg_id=3122)
    LAMBDA_C_POS_NAME            = constant("LAMBDAC+",  "Λc+ name as used in FLUKA.",                    pdg_id=4122)
    ANTI_LAMBDA_NAME             = constant("ALAMBDA",   "Anti-Λ0 name as used in FLUKA.",                pdg_id=-3122)
    ANTI_LAMBDA_C_NEG_NAME       = constant("ALAMBDC-",  "Anti-Λc- name as used in FLUKA.",               pdg_id=-4122)
    SIGMA_NAME                   = constant("SIGMAZER",  "Σ0 name as used in FLUKA.",                     pdg_id=3212)
    SIGMA_POS_NAME               = constant("SIGMA+",    "Σ+ name as used in FLUKA.",                     pdg_id=3222)
    SIGMA_NEG_NAME               = constant("SIGMA-",    "Σ- name as used in FLUKA.",                     pdg_id=3112)
    ANTI_SIGMA_NAME              = constant("ASIGMAZE",  "Anti-Σ0 name as used in FLUKA.",                pdg_id=-3212)
    ANTI_SIGMA_NEG_NAME          = constant("ASIGMA-",   "Anti-Σ- name as used in FLUKA.",                pdg_id=-3222)
    ANTI_SIGMA_POS_NAME          = constant("ASIGMA+",   "Anti-Σ+ name as used in FLUKA.",                pdg_id=-3112)
    XI_NAME                      = constant("XSIZERO",   "Ξ0 name as used in FLUKA.",                     pdg_id=3322)
    XI_NEG_NAME                  = constant("XSI-",      "Ξ- name as used in FLUKA.",                     pdg_id=3312)
    XI_C_NAME                    = constant("XSIC0",     "Ξc0 name as used in FLUKA.",                    pdg_id=4132)
    XI_C_POS_NAME                = constant("XSIC+",     "Ξc+ name as used in FLUKA.",                    pdg_id=4232)
    XI_C_PRIME_NAME              = constant("XSIPC0",    "Ξ'c0 name as used in FLUKA.",                   pdg_id=4312)
    XI_C_PRIME_POS_NAME          = constant("XSIPC+",    "Ξ'c+ name as used in FLUKA.",                   pdg_id=4322)
    ANTI_XI_NAME                 = constant("AXSIZERO",  "Anti-Ξ0 name as used in FLUKA.",                pdg_id=-3322)
    ANTI_XI_POS_NAME             = constant("AXSI+",     "Anti-Ξ+ name as used in FLUKA.",                pdg_id=-3312)
    ANTI_XI_C_NAME               = constant("AXSIC0",    "Anti-Ξc0 name as used in FLUKA.",               pdg_id=-4132)
    ANTI_XI_C_NEG_NAME           = constant("AXSIC-",    "Anti-Ξc- name as used in FLUKA.",               pdg_id=-4232)
    ANTI_XI_C_PRIME_NAME         = constant("AXSIPC0",   "Anti-Ξ'c0 name as used in FLUKA.",              pdg_id=-4312)
    ANTI_XI_C_PRIME_NEG_NAME     = constant("AXSIPC-",   "Anti-Ξ'c- name as used in FLUKA.",              pdg_id=-4322)
    OMEGA_NEG_NAME               = constant("OMEGA-",    "Ω- name as used in FLUKA.",                     pdg_id=3334)
    OMEGA_C_NAME                 = constant("OMEGAC0",   "Ωc0 name as used in FLUKA.",                    pdg_id=4332)
    ANTI_OMEGA_POS_NAME          = constant("AOMEGA+",   "Anti-Ω+ name as used in FLUKA.",                pdg_id=-3334)
    ANTI_OMEGA_C_NAME            = constant("AOMEGAC0",  "Anti-Ωc0 name as used in FLUKA.",               pdg_id=-4332)
