# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

from math import sqrt

import xtrack as xt
from xpart.pdg import get_name_from_pdg_id, get_properties_from_pdg_id
try:
    from xaux import singleton, FsPath  # TODO: once xaux is in Xsuite keep only this
except (ImportError, ModuleNotFoundError):
    from ...xaux import FsPath

from .reference_names import fluka_names
from .environment import format_fluka_float
from ...general import _pkg_root


def _is_proton(pdg_id):
    if isinstance(pdg_id, xt.Particles):
        pdg_id = pdg_id.pdg_id[0]
    return pdg_id >= 2212


def _is_ion(pdg_id):
    if isinstance(pdg_id, xt.Particles):
        pdg_id = pdg_id.pdg_id[0]
    return pdg_id >= 1000000000


def _is_lepton(pdg_id):
    if isinstance(pdg_id, xt.Particles):
        pdg_id = pdg_id.pdg_id[0]
    return 11 <= abs(pdg_id) <= 19


def get_include_files(particle_ref, include_files=[], verbose=True, lower_momentum_cut=None,
                      photon_lower_momentum_cut=None, electron_lower_momentum_cut=None,
                      include_showers=None, return_electrons=None, return_leptons=None,
                      return_neutrinos=None, return_protons=None, return_neutrons=None,
                      return_ions=None, return_exotics=None, return_all=None,
                      return_neutral=None, use_crystals=False):
    this_include_files = include_files.copy()
    # Required default include files
    if 'include_settings_beam.inp' not in [file.name for file in this_include_files]:
        this_include_files.append(_beam_include_file(particle_ref))
    if 'include_settings_physics.inp' in [file.name for file in this_include_files]:
        if photon_lower_momentum_cut is not None:
            raise ValueError("Physics include file already provided. Cannot change "
                            + "photon lower momentum cut.")
        if electron_lower_momentum_cut is not None:
            raise ValueError("Physics include file already provided. Cannot change "
                            + "electron lower momentum cut.")
        if lower_momentum_cut is not None:
            raise ValueError("Physics include file already provided. Cannot change "
                            + "hadron lower momentum cut.")
        if include_showers is not None:
            raise ValueError("Physics include file already provided. Cannot change "
                            + "shower settings.")
    else:
        if photon_lower_momentum_cut is None:
            photon_lower_momentum_cut = 1.e5
        elif photon_lower_momentum_cut < 1.e3:
            if verbose:
                print(f"Warning: Photon lower momentum cut of {photon_lower_momentum_cut/1.e3}keV "
                     + "is very low and will result in very long computation times.")
        if electron_lower_momentum_cut is None:
            if _is_lepton(particle_ref):
                electron_lower_momentum_cut = particle_ref.p0c[0] / 10
            else:
                electron_lower_momentum_cut = 1.e6
        elif electron_lower_momentum_cut < 1.e6:
            if verbose:
                print(f"Warning: Electron lower momentum cut of {electron_lower_momentum_cut/1.e6}MeV "
                     + "is very low and will result in very long computation times.")
        if lower_momentum_cut is None:
            lower_momentum_cut = particle_ref.p0c[0] / 10
            if _is_ion(particle_ref):
                _, A, _, _ = get_properties_from_pdg_id(particle_ref.pdg_id[0])
                lower_momentum_cut /= A
        elif lower_momentum_cut < 1.e9:
            if verbose:
                print(f"Warning: Hadron lower momentum cut of {lower_momentum_cut/1.e9}GeV "
                     + "is very low and will result in very long computation times.")
        if include_showers is None:
            include_showers = True if _is_lepton(particle_ref) else False
        physics_file = _physics_include_file(verbose=verbose, lower_momentum_cut=lower_momentum_cut,
                                             photon_lower_momentum_cut=photon_lower_momentum_cut,
                                             electron_lower_momentum_cut=electron_lower_momentum_cut,
                                             include_showers=include_showers)
        this_include_files.append(physics_file)
    if 'include_custom_scoring.inp' in [file.name for file in this_include_files]:
        if return_electrons is not None:
            raise ValueError("Custom scoring include file already provided. Cannot change "
                           + "`return_electrons`.")
        if return_leptons is not None:
            raise ValueError("Custom scoring include file already provided. Cannot change "
                           + "`return_leptons`.")
        if return_neutrinos is not None:
            raise ValueError("Custom scoring include file already provided. Cannot change "
                           + "`return_neutrinos`.")
        if return_protons is not None:
            raise ValueError("Custom scoring include file already provided. Cannot change "
                           + "`return_protons`.")
        if return_neutrons is not None:
            raise ValueError("Custom scoring include file already provided. Cannot change "
                           + "`return_neutrons`.")
        if return_ions is not None:
            raise ValueError("Custom scoring include file already provided. Cannot change "
                           + "`return_ions`.")
        if return_exotics is not None:
            raise ValueError("Custom scoring include file already provided. Cannot change "
                           + "`return_exotics`.")
        if return_all is not None:
            raise ValueError("Custom scoring include file already provided. Cannot change "
                           + "`return_all`.")
        if return_neutral is not None:
            raise ValueError("Custom scoring include file already provided. Cannot change "
                           + "`return_neutral`.")
    else:
        if return_electrons is None:
            return_electrons = _is_lepton(particle_ref)
        if return_leptons is None:
            return_leptons = _is_lepton(particle_ref)
        if return_neutrinos is None:
            return_neutrinos = False
        if return_protons is None:
            return_protons = _is_proton(particle_ref)
        if return_neutrons is None:
            return_neutrons = False
        if return_ions is None:
            return_ions = _is_ion(particle_ref)
        if return_exotics is None:
            return_exotics = False
        if return_all is None:
            return_all = False
        if return_neutral is None:
            return_neutral = False
        scoring_file = _scoring_include_file(return_electrons=return_electrons, return_leptons=return_leptons,
                                             return_neutrinos=return_neutrinos, return_protons=return_protons,
                                             return_neutrons=return_neutrons, return_ions=return_ions,
                                             return_exotics=return_exotics, return_all=return_all,
                                             return_neutral=return_neutral, use_crystals=use_crystals)
        this_include_files.append(scoring_file)
    # Add any additional include files
    for ff in (_pkg_root / 'scattering_routines' / 'fluka' / 'data').glob('include_*'):
        if ff.name not in [file.name for file in this_include_files]:
            this_include_files.append(ff)
    this_include_files = [FsPath(ff).resolve() for ff in this_include_files]
    for ff in this_include_files:
        if not ff.exists():
            raise FileNotFoundError(f"Include file not found: {ff}.")
        elif ff.parent != FsPath.cwd():
            ff.copy_to(FsPath.cwd())
    return this_include_files


def _beam_include_file(particle_ref):
    filename = FsPath("include_settings_beam.inp").resolve()
    pdg_id = particle_ref.pdg_id[0]
    momentum_cut = particle_ref.p0c[0] / 1e9 * 1.05
    hi_prope = "*"
    if pdg_id in fluka_names:
        name = fluka_names[pdg_id]
    elif _is_ion(pdg_id):
        _, A, Z, _ = get_properties_from_pdg_id(pdg_id)
        momentum_cut *= 3.2 / A # Upper limit (scaling is slightly arbitrary)
        name = "HEAVYION"
        hi_prope = "HI-PROPE  {Z:8}.0{A:8}.0"
    else:
        raise ValueError(f"Reference particle {get_name_from_pdg_id(pdg_id)} not "
                        + "supported by FLUKA.")
    beam = f"BEAM      {format_fluka_float(momentum_cut)}{50*' '}{name}"

    template = f"""\
******************************************************************************
*                          BEAM SETTINGS                                     *
******************************************************************************
* ..+....1....+....2....+....3....+....4....+....5....+....6....+....7..
* maximum momentum per nucleon (3000 for 3.5Z TeV, 6000 for 6.37Z TeV)
{beam}
{hi_prope}
*
BEAMPOS
*
* Only asking for loss map and touches map as in Sixtrack
* ..+....1....+....2....+....3....+....4....+....5....+....6....+....7..
SOURCE                                         87.       88.        1.
SOURCE           89.       90.       91.        0.       -1.       10.&
*SOURCE           0.0       0.0      97.0       1.0      96.0       1.0&&
"""
    with filename.open('w') as fp:
        fp.write(template)

    return filename


def _physics_include_file(*, verbose, lower_momentum_cut, photon_lower_momentum_cut,
                          electron_lower_momentum_cut, include_showers):
    filename = FsPath("include_settings_physics.inp").resolve()
    emf = "*EMF" if include_showers else "EMF"
    emfcut = "EMFCUT" if include_showers else "*EMFCUT"
    photon_lower_momentum_cut = format_fluka_float(photon_lower_momentum_cut/1.e9)
    # TODO: FLUKA electron mass
    electron_lower_energy_cut = sqrt(electron_lower_momentum_cut**2 + 511e3**2)
    electron_lower_energy_cut = format_fluka_float(electron_lower_energy_cut/1.e9)
    lower_momentum_cut /= 1.e9
    if verbose:
        print(f"Physics include file created with:\n"
             + f"  - Hadron and muon lower momentum cut: {lower_momentum_cut} GeV")
        if include_showers:
            print(f"\n  - EM showers: ON\n"
                + f"  - Photon lower momentum cut: {photon_lower_momentum_cut} GeV\n"
                + f"  - Electron lower momentum cut: {electron_lower_momentum_cut} GeV")

    template = f"""\
******************************************************************************
*                          PHYSICS SETTINGS                                  *
******************************************************************************
* ..+....1....+....2....+....3....+....4....+....5....+....6....+....7..
DEFAULTS                                                              PRECISIO
*
* Thresholds for secondary particle (electron, positron, photon) production
* applied to all materials; electron, positron: 1.0 MeV, photons: 0.1 MeV
* ..+....1....+....2....+....3....+....4....+....5....+....6....+....7..
{emfcut}    {electron_lower_energy_cut}{photon_lower_momentum_cut}       1.0       1.0  @LASTMAT       1.0PROD-CUT
{emfcut}    {electron_lower_energy_cut}{photon_lower_momentum_cut}       0.0       1.0  @LASTREG       1.0
*
* Kill EM showers
{emf}                                                                   EMF-OFF
*
* All particle transport thresholds up to 1 TeV
* ..+....1....+....2....+....3....+....4....+....5....+....6....+....7..
PART-THR  {format_fluka_float(-lower_momentum_cut)}            @LASTPAR                 0.0
PART-THR  {format_fluka_float(-2*lower_momentum_cut)}  DEUTERON                           0.0
PART-THR  {format_fluka_float(-3*lower_momentum_cut)}    TRITON                           0.0
PART-THR  {format_fluka_float(-3*lower_momentum_cut)}  3-HELIUM                           0.0
PART-THR  {format_fluka_float(-4*lower_momentum_cut)}  4-HELIUM                           0.0
*
*
* Activate single scattering
MULSOPT                                        1.0       1.0       1.0GLOBAL
*
* ..+....1....+....2....+....3....+....4....+....5....+....6....+....7..
* un-comment in case you use FLUKA pro
PHYSICS           1.                                                  COALESCE
PHYSICS           3.                                                  EVAPORAT
PHYSICS        1.D+5     1.D+5     1.D+5     1.D+5     1.D+5     1.D+5PEATHRES
PHYSICS           1.     0.005      0.15       2.0       2.0       3.0IONSPLIT
PHYSICS           2.                                                  EM-DISSO
"""
    with filename.open('w') as fp:
        fp.write(template)
    print("TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO ")
    print("PART-THR does not work with positive numbers (momentum cut), only with negative numbers (kinetic energy cut). Whyyy?")
    print("TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO ")
    return filename


def _scoring_include_file(*, return_electrons, return_leptons, return_neutrinos, return_protons,
                          return_neutrons, return_ions, return_exotics, return_all, return_neutral,
                          use_crystals=False):
    filename = FsPath("include_custom_scoring.inp").resolve()
    electrons = 'USRBDX' if return_electrons or return_leptons else '*USRBDX'
    leptons = 'USRBDX' if return_leptons else '*USRBDX'
    neutrinos = 'USRBDX' if return_neutrinos else '*USRBDX'
    protons = 'USRBDX' if return_protons else '*USRBDX'
    neutrons = 'USRBDX' if return_neutrons or (return_protons and return_neutral) else '*USRBDX'
    ions = 'USRBDX' if return_ions else '*USRBDX'
    exotics = 'USRBDX' if return_exotics else '*USRBDX'
    neutral_exotics = 'USRBDX' if return_exotics and return_neutral else '*USRBDX'
    if return_all:
        all_charged = 'USRBDX' if not return_neutral else '*USRBDX'
        all_particles = 'USRBDX' if return_neutral else '*USRBDX'
        electrons = leptons = neutrinos = protons = neutrons = ions = exotics \
                  = neutral_exotics = '*USRBDX'
    else:
        all_charged = all_particles = '*USRBDX'
    crystal = 'USRICALL' if use_crystals else '*USRICALL'

    template = f"""\
* Original Pascal
* include_custom_scoring
* USERDUMP        101.                 2.0                              SEAMU
* SERWEIG                             3.0
* hisix: send back all particles
*                                          boundary crossing    special sdum
* USRBDX          99.0  ALL-PART     -42.0   VAROUND  BLKROUND          BACK2ICO
*       above max energy (576000 for 7Z TeV)
*USRBDX      576000.0   1.0E-04     210.0
*
* ..+....1....+....2....+....3....+....4....+....5....+....6....+....7..
*
* Place in this file all the custom scoring you require
*
*USRBIN           10.    PROTON     -41.0        2.        2.     6900.Check
*USRBIN           -2.       -2.        0.      100.      100.     1000.&
*
* New way to give back particles to Icosim via fluscw.f routine (no more usrmed)
*     through a fake USRBDX estimator
USERWEIG                             3.0
* ..+....1....+....2....+....3....+....4....+....5....+....6....+....7....+....8
{all_charged}          99.0  ALL-CHAR     -42.0   VAROUND  TRANSF_D          BACK2ICO
{all_particles}          99.0  ALL-PART     -42.0   VAROUND  TRANSF_D          BACK2ICO
{electrons}          99.0  ELECTRON     -42.0   VAROUND  TRANSF_D          BACK2ICO
{electrons}          99.0  POSITRON     -42.0   VAROUND  TRANSF_D          BACK2ICO
{leptons}          99.0     MUON+     -42.0   VAROUND  TRANSF_D          BACK2ICO
{leptons}          99.0     MUON-     -42.0   VAROUND  TRANSF_D          BACK2ICO
{leptons}          99.0      TAU+     -42.0   VAROUND  TRANSF_D          BACK2ICO
{leptons}          99.0      TAU-     -42.0   VAROUND  TRANSF_D          BACK2ICO
{neutrinos}          99.0   NEUTRIE     -42.0   VAROUND  TRANSF_D          BACK2ICO
{neutrinos}          99.0  ANEUTRIE     -42.0   VAROUND  TRANSF_D          BACK2ICO
{neutrinos}          99.0   NEUTRIM     -42.0   VAROUND  TRANSF_D          BACK2ICO
{neutrinos}          99.0  ANEUTRIM     -42.0   VAROUND  TRANSF_D          BACK2ICO
{neutrinos}          99.0   NEUTRIT     -42.0   VAROUND  TRANSF_D          BACK2ICO
{neutrinos}          99.0  ANEUTRIT     -42.0   VAROUND  TRANSF_D          BACK2ICO
{exotics}          99.0     PION+     -42.0   VAROUND  TRANSF_D          BACK2ICO
{exotics}          99.0     PION-     -42.0   VAROUND  TRANSF_D          BACK2ICO
{neutral_exotics}          99.0    PIZERO     -42.0   VAROUND  TRANSF_D          BACK2ICO
{exotics}          99.0     KAON+     -42.0   VAROUND  TRANSF_D          BACK2ICO
{exotics}          99.0     KAON-     -42.0   VAROUND  TRANSF_D          BACK2ICO
{neutral_exotics}          99.0  KAONZERO     -42.0   VAROUND  TRANSF_D          BACK2ICO
{neutral_exotics}          99.0  AKAONZER     -42.0   VAROUND  TRANSF_D          BACK2ICO
{neutral_exotics}          99.0  KAONLONG     -42.0   VAROUND  TRANSF_D          BACK2ICO
{neutral_exotics}          99.0  KAONSHRT     -42.0   VAROUND  TRANSF_D          BACK2ICO
{protons}          99.0    PROTON     -42.0   VAROUND  TRANSF_D          BACK2ICO
{protons}        8000.0   1.0E-04     210.0
{protons}          99.0   APROTON     -42.0   VAROUND  TRANSF_D          BACK2ICO
{neutrons}          99.0   NEUTRON     -42.0   VAROUND  TRANSF_D          BACK2ICO
{ions}          99.0  DEUTERON     -42.0   VAROUND  TRANSF_D          BACK2ICO
{ions}          99.0    TRITON     -42.0   VAROUND  TRANSF_D          BACK2ICO
{ions}          99.0  3-HELIUM     -42.0   VAROUND  TRANSF_D          BACK2ICO
{ions}          99.0  4-HELIUM     -42.0   VAROUND  TRANSF_D          BACK2ICO
{ions}          99.0  HEAVYION     -42.0   VAROUND  TRANSF_D          BACK2ICO
{neutral_exotics}          99.0    LAMBDA     -42.0   VAROUND  TRANSF_D          BACK2ICO
{neutral_exotics}          99.0   ALAMBDA     -42.0   VAROUND  TRANSF_D          BACK2ICO
{exotics}          99.0  LAMBDAC+     -42.0   VAROUND  TRANSF_D          BACK2ICO
{exotics}          99.0  ALAMBDC-     -42.0   VAROUND  TRANSF_D          BACK2ICO
{exotics}          99.0    SIGMA-     -42.0   VAROUND  TRANSF_D          BACK2ICO
{exotics}          99.0    SIGMA+     -42.0   VAROUND  TRANSF_D          BACK2ICO
{neutral_exotics}          99.0  SIGMAZER     -42.0   VAROUND  TRANSF_D          BACK2ICO
{neutral_exotics}          99.0  ASIGMAZE     -42.0   VAROUND  TRANSF_D          BACK2ICO
{exotics}          99.0   ASIGMA-     -42.0   VAROUND  TRANSF_D          BACK2ICO
{exotics}          99.0   ASIGMA+     -42.0   VAROUND  TRANSF_D          BACK2ICO
{exotics}          99.0      XSI-     -42.0   VAROUND  TRANSF_D          BACK2ICO
{exotics}          99.0     AXSI+     -42.0   VAROUND  TRANSF_D          BACK2ICO
{neutral_exotics}          99.0   XSIZERO     -42.0   VAROUND  TRANSF_D          BACK2ICO
{neutral_exotics}          99.0  AXSIZERO     -42.0   VAROUND  TRANSF_D          BACK2ICO
{exotics}          99.0     XSIC+     -42.0   VAROUND  TRANSF_D          BACK2ICO
{exotics}          99.0    AXSIC-     -42.0   VAROUND  TRANSF_D          BACK2ICO
{neutral_exotics}          99.0     XSIC0     -42.0   VAROUND  TRANSF_D          BACK2ICO
{neutral_exotics}          99.0    AXSIC0     -42.0   VAROUND  TRANSF_D          BACK2ICO
{exotics}          99.0    XSIPC+     -42.0   VAROUND  TRANSF_D          BACK2ICO
{exotics}          99.0   AXSIPC-     -42.0   VAROUND  TRANSF_D          BACK2ICO
{neutral_exotics}          99.0    XSIPC0     -42.0   VAROUND  TRANSF_D          BACK2ICO
{neutral_exotics}          99.0   AXSIPC0     -42.0   VAROUND  TRANSF_D          BACK2ICO
{exotics}          99.0    OMEGA-     -42.0   VAROUND  TRANSF_D          BACK2ICO
{exotics}          99.0   AOMEGA+     -42.0   VAROUND  TRANSF_D          BACK2ICO
{neutral_exotics}          99.0   OMEGAC0     -42.0   VAROUND  TRANSF_D          BACK2ICO
{neutral_exotics}          99.0  AOMEGAC0     -42.0   VAROUND  TRANSF_D          BACK2ICO
{exotics}          99.0        D+     -42.0   VAROUND  TRANSF_D          BACK2ICO
{exotics}          99.0        D-     -42.0   VAROUND  TRANSF_D          BACK2ICO
{neutral_exotics}          99.0        D0     -42.0   VAROUND  TRANSF_D          BACK2ICO
{neutral_exotics}          99.0     D0BAR     -42.0   VAROUND  TRANSF_D          BACK2ICO
{exotics}          99.0       DS+     -42.0   VAROUND  TRANSF_D          BACK2ICO
{exotics}          99.0       DS-     -42.0   VAROUND  TRANSF_D          BACK2ICO
*
* Get back touches
* ..+....1....+....2....+....3....+....4....+....5....+....6....+....7....+....8
USERDUMP       100.0
*
* crystal scoring
{crystal}        50.0                                                  CRYSTAL
"""
    with filename.open('w') as fp:
        fp.write(template)
    print("TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO ")
    print("What about the second USRBDX with energy cut and number of energy bins?")
    print("TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO ")
    return filename