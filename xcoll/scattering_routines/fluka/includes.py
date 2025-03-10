# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

import numpy as np

from xpart.pdg import get_name_from_pdg_id, get_properties_from_pdg_id
try:
    from xaux import singleton, FsPath  # TODO: once xaux is in Xsuite keep only this
except (ImportError, ModuleNotFoundError):
    from ...xaux import FsPath

from .reference_names import fluka_names
from ...general import _pkg_root


def get_include_files(particle_ref, include_files=[]):
    this_include_files = include_files.copy()
    # Required default include files
    required_includes = ['include_settings_beam.inp', 'include_settings_physics.inp',
                        'include_custom_scoring.inp']
    for ff in required_includes:
        if ff not in [file.name for file in include_files]:
            if ff == 'include_settings_beam.inp':
                this_include_files.append(_beam_include_file(particle_ref))
            else:
                this_include_files.append(_pkg_root / 'scattering_routines' / 'fluka' / 'data' / ff)
    # Add any additional include files
    for ff in (_pkg_root / 'scattering_routines' / 'fluka' / 'data').glob('include_*'):
        if ff.name not in [file.name for file in include_files]:
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
    momentum_cut = particle_ref.p0c[0] / 1e9 * 1.2
    hi_prope = "*"
    if pdg_id in fluka_names:
        name = fluka_names[pdg_id]
    elif pdg_id >= 1000000000: # heavy ion
        _, A, Z, _ = get_properties_from_pdg_id
        momentum_cut *= 3.2 / A # Upper limit (scaling is slightly arbitrary)
        name = "HEAVYION"
        hi_prope = "HI-PROPE  {Z:8}.0{A:8}.0"
    else:
        raise ValueError(f"Reference particle {get_name_from_pdg_id(pdg_id)} not "
                        + "supported by FLUKA.")
    if momentum_cut < 1:
        momentum_cut = f"{momentum_cut:.4E}"
    else:
        momentum_cut = f"{int(np.ceil(momentum_cut)):9}."
    beam = f"BEAM      {momentum_cut}{50*' '}{name}"

    template = f"""\
******************************************************************************
*                          BEAM SETTINGS                                     *
******************************************************************************
*
* ..+....1....+....2....+....3....+....4....+....5....+....6....+....7..
*
*============================================================================*
*                        common to all cases                                 *
*============================================================================*
*
*
* maximum momentum per nucleon (3000 for 3.5Z TeV, 6000 for 6.37Z TeV)
{beam}
{hi_prope}
*
BEAMPOS
*
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