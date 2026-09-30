# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

from math import sqrt
from warnings import warn

import xtrack as xt
from xtrack.particles.pdg import (
    get_name_from_pdg_id,
    get_properties_from_pdg_id,
    is_ion,
)

from .environment import format_fluka_float
from .prototype import FlukaAssembly
from .reference_names import fluka_names_meta
from ...xaux import FsPath


_FLUKA_PDG_IDS = {
    value["value"]: value["pdg_id"]
    for value in fluka_names_meta.values()
}


def _is_ion(pdg_id):
    if isinstance(pdg_id, xt.Particles):
        pdg_id = pdg_id.pdg_id[0]
    return is_ion(pdg_id)


def get_include_files(particle_ref, include_files=None, *, verbose=True,
                      assemblies=None, bb_int=False, **kwargs):

    import xcoll as xc
    phys = xc.fluka.engine._physics_settings
    if include_files is None:
        include_files = []
    this_include_files = include_files.copy()
    if assemblies is None:
        assemblies = []

    # Required default include files
    if 'include_settings_beam.inp' not in [file.name for file in this_include_files]:
        this_include_files.append(_beam_include_file(particle_ref, bb_int=bb_int))
    if 'include_settings_physics.inp' not in [file.name for file in this_include_files]:
        physics_file = _physics_include_file(
                            verbose=verbose,
                            particle_ref=particle_ref,
                            hadron_lower_momentum_cut=phys.hadron_lower_momentum_cut,
                            photon_lower_momentum_cut=phys.photon_lower_momentum_cut,
                            electron_lower_momentum_cut=phys.electron_lower_momentum_cut,
                            include_showers=phys.include_showers,
                            include_single_coulomb=phys.include_single_coulomb,
                            include_multiple_coulomb=phys.include_multiple_coulomb,
                            include_elastic=phys.include_elastic,
                            include_inelastic=phys.include_inelastic,
                            include_pair_production=phys.include_pair_production,
                            include_bremsstrahlung=phys.include_bremsstrahlung,
                            include_ionisation_fluctuations=phys.include_ionisation_fluctuations,
                            extra_physics_cards=kwargs.get('extra_physics_cards', [])
                        )
        this_include_files.append(physics_file)
    if 'include_custom_scoring.inp' not in [file.name for file in this_include_files]:
        scoring_file = _scoring_include_file(verbose=verbose, return_list=phys,
                                             get_touches=kwargs.get('touches', False),
                                             use_crystals=any([assm.is_crystal for assm in assemblies]))
        this_include_files.append(scoring_file)
    if 'include_custom_assignmat.inp' not in [file.name for file in this_include_files]:
        material_file = _assignmat_include_file(particle_ref, assemblies=assemblies)
        this_include_files.append(material_file)

    # Add any additional include files
    if 'include_custom_biasing.inp' not in [file.name for file in this_include_files]:
        this_include_files.append(_biasing_include_file())
    if 'include_define.inp' not in [file.name for file in this_include_files]:
        this_include_files.append(_define_include_file())
    this_include_files = [FsPath(ff).resolve() for ff in this_include_files]
    for ff in this_include_files:
        if not ff.exists():
            raise FileNotFoundError(f"Include file not found: {ff}.")
        elif ff.parent != FsPath.cwd():
            ff.copy_to(FsPath.cwd(), method='mount')
    return this_include_files, kwargs


def _assignmat_include_file(particle_ref, assemblies=[]):
    from xcoll.materials.database import db as mdb
    template =  ''.join([mat._generated_fluka_code for mat in mdb.fluka.values()
                         if mat._generated_fluka_code is not None])
    if any([assm.is_crystal for assm in assemblies]):
        template += f"""\
* ..+....1....+....2....+....3....+....4....+....5....+....6....+....7..
* Crystal Card
* what (1) = REGNUM ( mandatory )
* what (2) = Crystal bending angle [ mrad ]
* what (3) = A parameter = Crystal length [ cm ]
* what (4) = B parameter = Torsion coefficient [ mrad / cm ]
* what (5) = C parameter = Quasi - Mosaic factor [ mrad / cm ]
* what (6) = D parameter = Crystal temperature [K]
* sdum = crystal type
* what (7, 8, 9) = U/V/ W of nominal crystal normal to the curvature plane
* what (10, 11, 12) = U /V/W of nominal crystal channel direction
* what (13, 14, 15) = X /Y/Z of crystal nominal position
* ..+....1....+....2....+....3....+....4....+....5....+....6....+....7..
* CRYSTAL what (1) what (2) what (3) what (4) what (5) what (6) sdum
* CRYSTAL what (7) what (8) what (9) what (10) what (11) what (12) &
* CRYSTAL what (13) what (14) what (15) &&
*
"""
    for crystal in assemblies:
        if not crystal.is_crystal:
            continue
        name = crystal.name
        if isinstance(crystal, FlukaAssembly) and crystal.is_generic():
            name = crystal.prototypes[1].name   # region is tied to crystal body
        l = crystal.length * 100
        bang = round(l/(crystal.bending_radius*100) *1000, 6) # mrad
        if crystal.fluka_position is not None:
            pos = crystal.fluka_position
        else:
            if isinstance(crystal, FlukaAssembly):
                raise ValueError(f"Crystal {name} has no position in the prototypes file.")
            if len(crystal.dependant_assemblies) == 0:
                raise ValueError(f"Crystal {name} has no position in the prototypes file.")
            elif len(crystal.dependant_assemblies) > 1:
                raise ValueError(f"Crystal {name} has multiple positions in the prototypes file. "
                                + "(too many dependant asemblies).")
            pos = crystal.dependant_assemblies[0].fluka_position
        pos1   = format_fluka_float(pos[3]-50.0+1e-4) # XXX Does it always work with the 50cm shift?
        pos2   = format_fluka_float(pos[4])
        pos3   = format_fluka_float(pos[5])
        pos4   = format_fluka_float(particle_ref.p0c[0] / 1e9)
        template += f"""\
* ..+....1....+....2....+....3....+....4....+....5....+....6....+....7..
CRYSTAL     {name:>8}{format_fluka_float(bang)}{format_fluka_float(l)}       0.0       0.0     300.0 110
CRYSTAL          0.0      -1.0       0.0       0.0       0.0       1.0 &
CRYSTAL   {pos1}{pos2}{pos3}{pos4}                        &&
"""
    filename = FsPath("include_custom_assignmat.inp").resolve()

    with filename.open('w') as fp:
        fp.write(template)
    return filename


def _beam_include_file(particle_ref, bb_int=False):
    import xcoll as xc
    filename = FsPath("include_settings_beam.inp").resolve()
    pdg_id = particle_ref.pdg_id[0]
    momentum_cut = particle_ref.p0c[0] / 1e9 * 1.05

    hi_prope = "*"
    if pdg_id in xc.fluka.particle_names:
        name = xc.fluka.particle_names[pdg_id]
    elif _is_ion(pdg_id):
        _, A, Z, _ = get_properties_from_pdg_id(pdg_id)
        momentum_cut *= 3.2 / A # Upper limit (scaling is slightly arbitrary)
        name = "HEAVYION"
        hi_prope = f"HI-PROPE  {Z:8}.0{A:8}.0"
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
"""
    if bb_int:
        template += f"""\
* ..+....1....+....2....+....3....+....4....+....5....+....6....+....7..
* W(1): interaction type, W(2): Index IR, W(3): Length of IR[cm], W(4): SigmaZ
* Interaction type: 1.0 (Inelastic), 10.0 (Elastic), 100.0 (EMD)
* SigmaZ: sigma_z (cm) for the Gaussian sampling of the collision position around the center of the insertion
SOURCE    {format_fluka_float(bb_int['int_type'])}{format_fluka_float(1)}{format_fluka_float(20)}{format_fluka_float(bb_int['sigma_z'])}
SOURCE           89.       90.       91.        0.                    &
* ..+....1....+....2....+....3....+....4....+....5....+....6....+....7..
* W(1): Theta2-> Polar angle (rad) between b2 and -z direction.
* W(2): Azimuthal angle (deg) defining the crossing plane.
* W(3): Sigma_x for Gaussian sampling in b2, dir (x', y', 1) wrt to px
* W(4): Sigma_y for Gaussian sampling in b2, dir (x', y', 1) wrt to py
* W(5): Z b2
* W(6): A b2
USRICALL  {format_fluka_float(bb_int['theta2'])}{format_fluka_float(bb_int['xs'])}{format_fluka_float(bb_int['sigma_p_x2'])}{format_fluka_float(bb_int['sigma_p_y2'])}{format_fluka_float(bb_int['Z'])}{format_fluka_float(bb_int['A'])}BBEAMCOL
USRICALL  {format_fluka_float(bb_int['betx'])}{format_fluka_float(bb_int['alfx'])}{format_fluka_float(bb_int['dx'])}{format_fluka_float(bb_int['dpx'])}{format_fluka_float(bb_int['rms_emx'])}          BBCOL_H
USRICALL  {format_fluka_float(bb_int['bety'])}{format_fluka_float(bb_int['alfy'])}{format_fluka_float(bb_int['dy'])}{format_fluka_float(bb_int['dpy'])}{format_fluka_float(bb_int['rms_emy'])}          BBCOL_V
* * ..+....1....+....2....+....3....+....4....+....5....+....6....+....7..
* *USRICALL      #OFFSET                                                        BBCOL_O
* W(1,2,3,4,5,6): X, X', Y,  Y', dT, pc
* USRICALL         0.0       0.0       0.0       0.0       0.0       0.0BBCOL_O
USRICALL  {format_fluka_float(bb_int['offset'][0])}{format_fluka_float(bb_int['offset'][1])}{format_fluka_float(bb_int['offset'][2])}{format_fluka_float(bb_int['offset'][3])}       0.0       0.0BBCOL_O
* *USRICALL      #SIGDPP                                                        BBCOL_L
* USRICALL         0.0       0.0       0.0       0.0       0.0       0.0BBCOL_L
* ..+....1....+....2....+....3....+....4....+....5....+....6....+....7..
USRICALL  {format_fluka_float(bb_int['sigma_dpp'])}                                                  BBCOL_L
*
"""
    else:
        template += f"""\
* Only asking for loss map and touches map as in Xsuite
* ..+....1....+....2....+....3....+....4....+....5....+....6....+....7..
SOURCE                                         87.       88.        1.
SOURCE           89.       90.       91.        0.       -1.       10.&
SOURCE           0.0       0.0      97.0       1.0      96.0       1.0&&
"""


    with filename.open('w') as fp:
        fp.write(template)
    return filename


def _physics_include_file(*, verbose, particle_ref, hadron_lower_momentum_cut,
                          photon_lower_momentum_cut, electron_lower_momentum_cut,
                          include_showers, include_single_coulomb, include_multiple_coulomb,
                          include_elastic, include_inelastic, include_pair_production,
                          include_bremsstrahlung, include_ionisation_fluctuations,
                          extra_physics_cards=None):
    filename = FsPath("include_settings_physics.inp").resolve()
    if extra_physics_cards is None:
        extra_physics_cards = []
    # Showers
    emf = "*EMF" if include_showers else "EMF"
    deltaray = "DELTARAY" if not include_showers else "*DELTARAY"
    emfcut = "EMFCUT" if include_showers else "*EMFCUT"
    # Ionisation losses
    if include_ionisation_fluctuations:
        # Ionisation losses are on by default
        ionisation_losses = ""
    else:
        ionisation_losses  = "* Deactivate ionisation losses\n"
        ionisation_losses += "IONFLUCT        -1.0      -1.0            BLCKHOLE  @LASTMAT\n*"
    # Coulomb scattering
    if include_single_coulomb and include_multiple_coulomb:
        coulomb  = "* Activate Coulomb scattering (single and multiple)\n"
        coulomb += "MULSOPT                                        1.0       1.0       1.0GLOBAL\n*"
    elif include_single_coulomb:
        coulomb  = "* Activate single Coulomb scattering (disable multiple Coulomb scattering)\n"
        coulomb += "MULSOPT          0.0       0.0       0.0       1.0       1.0 99999999.GLOBAL\n*"
    elif include_multiple_coulomb:
        coulomb  = "* Activate multiple Coulomb scattering (disable single Coulomb scattering)\n"
        coulomb += "MULSOPT                                       -1.0      -1.0      -1.0GLOBAL\n*"
    else:
        # No multiple Coulomb either
        coulomb  = "* Deactivate all Coulomb scattering\n"
        coulomb += "MULSOPT                                       -1.0      -1.0      -1.0GLOBAL\n"
        coulomb += "MULSOPT                    3.0       3.0  BLCKHOLE  @LASTMAT\n*"
    # Pair production and bremsstrahlung
    if include_pair_production and include_bremsstrahlung:
        pairbrem  = "* Activate pair production and bremsstrahlung by muons/hadrons\n"
        pairbrem += "PAIRBREM         3.0                      BLCKHOLE  @LASTMAT\n*"
    elif include_pair_production:
        pairbrem  = "* Activate pair production by muons/hadrons (disable bremsstrahlung)\n"
        pairbrem += "PAIRBREM         1.0                      BLCKHOLE  @LASTMAT\n*"
    elif include_bremsstrahlung:
        pairbrem  = "* Activate bremsstrahlung by muons/hadrons (disable pair production)\n"
        pairbrem += "PAIRBREM         2.0                      BLCKHOLE  @LASTMAT\n*"
    else:
        pairbrem  = "* Disable pair production and bremsstrahlung by muons/hadrons\n"
        pairbrem += "PAIRBREM        -3.0                      BLCKHOLE  @LASTMAT\n*"
    # Cuts
    photon_lower_momentum_cut = format_fluka_float(photon_lower_momentum_cut/1.e9)
    # TODO: FLUKA electron mass
    electron_lower_energy_cut = sqrt(electron_lower_momentum_cut**2 + 511e3**2)
    electron_lower_energy_cut = format_fluka_float(electron_lower_energy_cut/1.e9)
    hadron_lower_momentum_cut /= 1.e9
    # Hadronic interactions
    thresh = format_fluka_float(5*particle_ref.energy0[0] / 1.e9)
    if not include_elastic and not include_inelastic:
        hadron = "* Deactivate hadronic interactions\n"
        hadron += f"THRESHOLd                     {thresh}{thresh}\n*"
    elif not include_elastic:
            hadron = "* Deactivate elastic hadronic interactions\n"
            hadron += f"THRESHOLd                     {thresh}\n*"
    elif not include_inelastic:
            hadron = "* Deactivate inelastic hadronic interactions\n"
            hadron += f"THRESHOLd                               {thresh}\n*"
    else:
        # Hadronic interactions are on by default
        hadron = ''
    # EM dissociation
    emdisso = "PHYSICS" if _is_ion(particle_ref.pdg_id[0]) else "*PHYSICS"
    # Extra physics cards
    if extra_physics_cards:
        extra_physics_cards = ["* Extra physics cards", *extra_physics_cards, '*']
    extra_physics_cards = "\n".join(extra_physics_cards)

    if verbose:
        print(f"Physics include file created with:\n"
             + f"  - Hadron and muon lower momentum cut: {hadron_lower_momentum_cut} GeV")
        if include_showers:
            print(f"  - EM showers: ON\n"
                + f"  - Photon lower momentum cut: {photon_lower_momentum_cut} GeV\n"
                + f"  - Electron lower momentum cut: {electron_lower_energy_cut} GeV")
        else:
            print(f"  - EM showers: OFF")
        print(f"  - Single Coulomb scattering: {'ON' if include_single_coulomb else 'OFF'}")
        print(f"  - Multiple Coulomb scattering: {'ON' if include_multiple_coulomb else 'OFF'}")
        print(f"  - Pair production: {'ON' if include_pair_production else 'OFF'}")
        print(f"  - Bremsstrahlung: {'ON' if include_bremsstrahlung else 'OFF'}")
        print(f"  - Ionisation losses: {'ON' if include_ionisation_fluctuations else 'OFF'}")
        print(f"  - Elastic hadronic interactions: {'ON' if include_elastic else 'OFF'}")
        print(f"  - Inelastic hadronic interactions: {'ON' if include_inelastic else 'OFF'}")
        if extra_physics_cards:
            print(f"  - Extra physics cards:\n{extra_physics_cards}")

    template = f"""\
******************************************************************************
*                          PHYSICS SETTINGS                                  *
******************************************************************************
* ..+....1....+....2....+....3....+....4....+....5....+....6....+....7..
DEFAULTS                                                              PRECISIO
*
* Thresholds for secondary particle (electron, positron, photon) production
* applied to all materials; electron, positron: 1.0 MeV, photons: 0.1 MeV
{emfcut}    {electron_lower_energy_cut}{photon_lower_momentum_cut}       1.0       1.0  @LASTMAT       1.0PROD-CUT
{emfcut}    {electron_lower_energy_cut}{photon_lower_momentum_cut}       0.0       1.0  @LASTREG       1.0
*
* Kill EM showers
{emf}                                                                   EMF-OFF
{deltaray}          -1                      BLCKHOLE  @LASTMAT
*
{coulomb}
{pairbrem}
{ionisation_losses}
* Particle transport thresholds (this sets physics process thresholds as well)
PART-THR  {format_fluka_float(hadron_lower_momentum_cut)}            @LASTPAR                 0.0
PART-THR  {format_fluka_float(2*hadron_lower_momentum_cut)}  DEUTERON                           0.0
PART-THR  {format_fluka_float(3*hadron_lower_momentum_cut)}    TRITON                           0.0
PART-THR  {format_fluka_float(3*hadron_lower_momentum_cut)}  3-HELIUM                           0.0
PART-THR  {format_fluka_float(4*hadron_lower_momentum_cut)}  4-HELIUM                           0.0
*
{hadron}
* Activate ion interactions:
* COALESCE and EVAPORAT if beam or target is ion, EM-DISSO if beam is ion
PHYSICS           1.                                                  COALESCE
PHYSICS           3.                                                  EVAPORAT
{emdisso}           2.                                                  EM-DISSO
*
* PHYSICS limits automatically set (default -1) by BEAM card 
*
* No low-energy neutron transport
LOW-PWXS          -1
*
{extra_physics_cards}
"""
    with filename.open('w') as fp:
        fp.write(template)
    return filename


def _scoring_include_file(*, verbose, return_list, get_touches=False, use_crystals=False):
    filename = FsPath("include_custom_scoring.inp").resolve()
    template = f"""\
* New way to give back particles to Icosim via fluscw.f routine (no more usrmed)
*     through a fake USRBDX estimator
* ..+....1....+....2....+....3....+....4....+....5....+....6....+....7....+....8
USERWEIG                             3.0
"""

    return_all = return_list.return_all
    return_all_charged = return_list.return_all_charged
    explicit_return = return_list._extra_pdg_ids_to_return
    covered_ids = set(_FLUKA_PDG_IDS.values())

    uncovered_explicit = [pid for pid in explicit_return
                          if pid not in covered_ids and not is_ion(pid)]
    if len(uncovered_explicit) > 0:
        print("WARNING: The following explicitly requested PDG IDs are not "
              "covered by the FLUKA scoring include file:\n"
              + ", ".join(str(pid) for pid in uncovered_explicit)
              + "Changed to return_all=True to ensure they are returned, if "
                "they are produced in the simulation."
        )
        return_all = True

    # Find explicitly requested particles for which there is no dedicated FLUKA
    # scorer. Ions are handled through HEAVYION. Charged particles are also
    # already covered when ALL-CHAR is active.
    uncovered_explicit = []
    for pdg_id in explicit_return:
        if is_ion(pdg_id):
            continue
        if pdg_id in covered_ids:
            continue
        if return_all_charged:
            charge = get_properties_from_pdg_id(pdg_id)[0]
            if abs(charge) > 1.e-12:
                continue
        uncovered_explicit.append(pdg_id)
    fallback_to_all = len(uncovered_explicit) > 0
    score_all = return_all or fallback_to_all

    if fallback_to_all:
        warn("The following explicitly requested PDG IDs do not have "
            "dedicated FLUKA scoring cards: "
            + ", ".join(str(pid) for pid in sorted(uncovered_explicit))
            + ". Falling back to ALL-PART scoring; the normal "
              "PhysicsSettingsHelper mask will filter the returned "
              "particles afterwards.", RuntimeWarning, stacklevel=2)

    if score_all:
        template += f"""\
USRBDX          99.0  ALL-PART     -42.0   VAROUND  TRANSF_D          BACK2ICO
"""

    else:
        if return_all_charged:
            template += f"""\
USRBDX          99.0  ALL-CHAR     -42.0   VAROUND  TRANSF_D          BACK2ICO
"""

        for fluka_name, pdg_id in _FLUKA_PDG_IDS.items():
            # Charged particles are already covered by ALL-CHAR.
            if return_all_charged:
                charge = get_properties_from_pdg_id(pdg_id)[0]
                if abs(charge) > 1.e-12:
                    continue
            # Do not score if not requested
            if not return_list.pdg_id_is_returned(pdg_id):
                continue
            # Score
            template += f"""\
USRBDX          99.0{fluka_name:>10}     -42.0   VAROUND  TRANSF_D          BACK2ICO
"""
            # FLUKA has two additional photon-like scorers:
            if fluka_name == "PHOTON":
                template += f"""\
USRBDX          99.0  OPTIPHOT     -42.0   VAROUND  TRANSF_D          BACK2ICO
USRBDX          99.0       RAY     -42.0   VAROUND  TRANSF_D          BACK2ICO
"""

        # HEAVYION is a FLUKA class rather than one concrete PDG ID.
        light_ion_pdg_ids = {
            _FLUKA_PDG_IDS["DEUTERON"],
            _FLUKA_PDG_IDS["TRITON"],
            _FLUKA_PDG_IDS["3-HELIUM"],
            _FLUKA_PDG_IDS["4-HELIUM"],
        }
        explicit_heavy_ion = any(
            is_ion(pdg_id) and pdg_id not in light_ion_pdg_ids
            for pdg_id in return_list._extra_pdg_ids_to_return
        )
        if not return_all_charged and (return_list.return_ions or explicit_heavy_ion):
            template += f"""\
USRBDX          99.0  HEAVYION     -42.0   VAROUND  TRANSF_D          BACK2ICO
"""

    if use_crystals:
        template += f"""*
* crystal scoring
* ..+....1....+....2....+....3....+....4....+....5....+....6....+....7....+....8
USRICALL        50.0                                                  CRYSTAL
"""
    if get_touches:
        template += f"""*
* Get back touches
* ..+....1....+....2....+....3....+....4....+....5....+....6....+....7....+....8
USERDUMP       100.0
"""

    if verbose:
        print("Scoring include file created with:")
        if score_all:
            if fallback_to_all and not return_all:
                print("  - Particle scoring: all particles "
                      "(fallback for explicitly requested PDG IDs)")
            else:
                print("  - Particle scoring: all particles")
        elif return_all_charged:
            print("  - Particle scoring: all charged particles")
            neutral_explicit = [
                pdg_id for pdg_id in sorted(explicit_return)
                if abs(get_properties_from_pdg_id(pdg_id)[0]) < 1.e-12
            ]
            if neutral_explicit:
                print("  - Additional explicitly requested neutral PDG IDs: "
                    + ", ".join(str(pdg_id) for pdg_id in neutral_explicit))
        else:
            enabled = [flag for flag in return_list._return_leaf_flags
                       if getattr(return_list, flag)]
            print("  - Particle scoring: selected particle types")
            if enabled:
                print("  - Enabled return types: "
                    + ", ".join(flag.removeprefix("return_")
                                for flag in enabled)
                )
            if explicit_return:
                print("  - Explicitly returned PDG IDs: "
                    + ", ".join(str(pdg_id)
                                for pdg_id in sorted(explicit_return))
                )
            if return_list._extra_pdg_ids_to_kill:
                print("  - Explicitly excluded PDG IDs: "
                    + ", ".join(
                        str(pdg_id)
                        for pdg_id
                        in sorted(
                            return_list._extra_pdg_ids_to_kill
                        )
                    )
                )
        print(f"  - Use crystals: {'ON' if use_crystals else 'OFF'}")
        print(f"  - Get touches:  {'ON' if get_touches else 'OFF'}")

    with filename.open('w') as fp:
        fp.write(template)
    return filename


def _biasing_include_file():
    filename = FsPath(f"include_custom_biasing.inp").resolve()
    template = f"""\
******************************************************************************
*                             BIASING SETTINGS                               *
******************************************************************************
* ..+....1....+....2....+....3....+....4....+....5....+....6....+....7..
*
* Place in this file all the custom biasing you require
*
* --------------------------------
* Standard Leading particle biasing applies everywhere
*________+_________+_________+_________+_________+_________+_________+
*EMF-BIAS      1022.0     10.00     10.00       2.0  @LASTREG          LPBEMF
*
"""
    with filename.open('w') as fp:
        fp.write(template)
    return filename


def _define_include_file():
    filename = FsPath(f"include_define.inp").resolve()
    template = f"""\
******************************************************************************
*                                  DEFINES                                   *
******************************************************************************
* ..+....1....+....2....+....3....+....4....+....5....+....6....+....7..
*
* Place in this file all the custom defines you require
*
"""
    with filename.open('w') as fp:
        fp.write(template)
    return filename
