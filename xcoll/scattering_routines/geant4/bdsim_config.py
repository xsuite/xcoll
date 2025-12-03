# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2025                  #
# ######################################### #

import json
import numpy as np

from xtrack.particles.pdg import get_properties_from_pdg_id, is_proton, is_ion

from ...materials import RefMaterial
try:
    from xaux import FsPath  # TODO: once xaux is in Xsuite keep only this
except (ImportError, ModuleNotFoundError):
    from ...xaux import FsPath


_header_start = "! ** XCOLL START  **"
_header_stop  = "! ** XCOLL END  **"


def create_bdsim_config_file(element_dict, particle_ref, physics_list=None, extra_opts=[],
                             extra_input=[], verbose=True, _all_black=False, **kwargs):
    momentum = int(np.ceil((particle_ref.p0c[0] + 1) / 1.e9))  # in GeV
    pdg_id = particle_ref.pdg_id[0]
    if pdg_id is None or pdg_id == 0:
        raise ValueError('particle_ref must have a valid pdg_id')
    particle = pdg_id
    if is_proton(pdg_id):
        particle = '"proton"'
    elif is_ion(pdg_id):
        _, A, Z, _ = get_properties_from_pdg_id(pdg_id)
        q0 = particle_ref.q0
        if not np.isclose(Z, q0):
            if verbose:
                print("Warning: particle_ref is a partially stripped ion.")
        particle = f'"ion {A} {Z} {q0}"'
    else:
        particle = f'"{pdg_id}"'   # Need quoted strings
    gmad = [f"beam, particle={particle}, momentum={momentum}*GeV;"]
    if physics_list is None:
        if is_ion(pdg_id):
            physics_list = "em ftfp_bert hadronic_elastic decay neutron_tracking_cut em_extra ion stopping ion_em_dissociation"
            # physics_list = "em decay hadronic_elastic ftfp_bert stopping neutron_tracking_cut ion ion_em_dissociation"
        else:
            # This is the same as g4FTFP_BERT (verified with option, physicsVerbose=1;)
            # physics_list = "em ftfp_bert hadronic_elastic decay neutron_tracking_cut em_extra ion"
            physics_list = 'g4FTFP_BERT'
    gmad.append(f'option, physicsList="{physics_list}";')
    # gmad.append(f'option, seed={seed};')
    if _all_black:
        gmad.append('option, collimatorsAreInfiniteAbsorbers=1;')
    gmad.append('')
    if extra_opts:
        if hasattr(extra_opts, '__iter__') and not isinstance(extra_opts, str):
            gmad.extend(extra_opts)
        else:
            gmad.append(extra_opts)
        gmad.append('')
    gmad.extend(generate_material_definitions(element_dict, verbose=verbose))
    if extra_input:
        if hasattr(extra_input, '__iter__') and not isinstance(extra_input, str):
            gmad.extend(extra_input)
        else:
            gmad.append(extra_input)
    input_file = FsPath(FsPath.cwd() / 'geant4_input.gmad').resolve()
    with open(input_file, 'w') as fp:
        fp.write('\n'.join([*_generate_xcoll_header(element_dict), *gmad]))
    return input_file, kwargs


def generate_material_definitions(element_dict, verbose=True):
    from ...beam_elements import Geant4CollimatorTip
    code = []
    for name, el in element_dict.items():
        all_mat = [el.material]
        if isinstance(el, Geant4CollimatorTip):
            all_mat.append(el.tip_material)
        for mat in all_mat:
            if mat is None:
                raise ValueError(f"Material not set for element {name}!")
            elif not isinstance(mat, RefMaterial):
                if mat.geant4_name is None:
                    mat._generate_geant4_code()
                if mat._generated_geant4_code is not None \
                and mat._generated_geant4_code not in code:
                    code.append(mat._generated_geant4_code)
            # Verify that material has a Geant4 name (after possible generation)
            if mat.geant4_name is None:
                raise ValueError(f"Material for element {name} has no Geant4 name!")
        if verbose:
            mats = [m.geant4_name for m in all_mat]
            if len(mats) == 1:
                print(f"Element {name} (id {el.geant4_id}) uses material {mats[0]}")
            else:
                print(f"Element {name} (id {el.geant4_id}) uses materials {mats}")
    code.append('')
    return code


def get_collimators_from_input_file(input_file):
    with open(input_file, 'r') as fp:
        data = fp.read()
    if _header_start not in data or _header_stop not in data:
        raise ValueError("No Xcoll header found in input file. Regenerate input file!")
    commented_dict = data.split(_header_start)[1].split(_header_stop)[0].split('\n')[1:-1]
    cleaned_dict = "".join([val[3:] for val in commented_dict])
    return json.loads(cleaned_dict)


def _generate_xcoll_header(element_dict):
    from ...beam_elements import Geant4CollimatorTip
    header = ["!  DO NOT CHANGE THIS HEADER", _header_start, "!  {"]
    for kk, el in element_dict.items():
        vv = {
            "geant4_id": np.array(el.geant4_id).tolist(),
            "length": np.array(el.length).tolist(),
            "material": np.array(el.material.geant4_name).tolist(),
            "angle": np.array(el.angle).tolist(),
            "jaw": [np.array(el.jaw_L).tolist(), np.array(el.jaw_R).tolist()],
            "tilt": np.array(el.tilt).tolist(),
        }
        if isinstance(el, Geant4CollimatorTip):
            vv["tip_material"] = np.array(el.tip_material.geant4_name).tolist()
            vv["tip_thickness"]   = np.array(el.tip_thickness).tolist()
        header.append(f'!  "{kk}": ' + json.dumps(vv) + ',')
    header[-1] = header[-1][:-1]  # remove last comma
    header.append("!  }")
    header.append(_header_stop)
    header.append('')
    return header