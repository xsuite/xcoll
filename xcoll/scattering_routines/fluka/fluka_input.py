# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

import json
import numpy as np

from ...beam_elements import FlukaCollimator


_header_start = "*  XCOLL START  **"
_header_stop  = "*  XCOLL END  **"

def fluka_builder(line):
    from .paths import fluka_builder as fluka_builder_path
    import sys
    import os
    from importlib.util import spec_from_file_location, module_from_spec

    sys.path.append(fluka_builder_path)
    # Convert to string and append to sys.path
    path_str = str(fluka_builder_path)
    print("Appending to sys.path:", path_str)
    sys.path.append(path_str)

    # Check if the file exists
    file_path = fluka_builder_path / "FLUKA_builder_with_main.py"
    print("Checking if file exists at:", file_path)
    if file_path.exists():
        print("File found. Attempting import.")
        try:
            import FLUKA_builder_with_main_ads as fb
            print("Module imported successfully.")
        except ImportError as e:
            print("ImportError:", e)
        except Exception as e:
            print("An unexpected error occurred:", e)
    else:
        print("File not found at:", file_path)
    # module_path = os.path.join(fluka_builder_path, 'FLUKA_builder_with_main_ads.py')

    # if os.path.exists(module_path):
    #     print("Confirmed FLUKA_builder_with_main_ads.py exists in the directory.")
    #     spec = spec_from_file_location("FLUKA_builder_with_main_ads", module_path)
    #     fb = module_from_spec(spec)
    #     spec.loader.exec_module(fb)
    #     print("Module imported successfully.")
    # else:
    #     print("FLUKA_builder_with_main_ads.py does not exist in the directory. Check the file name and path.")
    # print(f"Attempting to append '{fluka_builder_path}' to sys.path")

    import xtrack as xt
    import math
    elements, names = line.get_elements_of_type(FlukaCollimator)

    collimator_dict = {}

    for element, name in zip(elements, names):

        geometrical_emittance = xt.PROTON_MASS_EV * element.nemitt_x / (element.optics['upstream']['p0c'])
        betx = element.optics['upstream']['betx'][0]
        # betx = element.optics['downstream']['betx'] ???
        bety = element.optics['upstream']['bety'][0]
        # bety = element.optics['downstream']['bety'] ???
        material = None
        length = element.length
        angle = float(math.radians(element.angle))
        sigma_x = float(np.sqrt
       (geometrical_emittance * betx))
       
        sigma_y = float(np.sqrt(geometrical_emittance * bety))
        offset = 0.0
        tilt_1 = 0 if not isinstance(element.tilt_L, float) else math.radians(element.tilt_L) # if 'tilt_L' in line[name] else 0.0,
        tilt_2 = 0 if not isinstance(element.tilt_R, float) else math.radians(element.tilt_R) # if 'tilt_R' in line[name] else 0.0,

        # tilt_2 = math.radians(element.tilt_R) # if 'tilt_R' in line[name] else 0.0,

        # nsig = element.gap[0] if isinstance(element.gap, list) else element.gap
        nsig = element.gap
        nsig = 999.0 if not isinstance(nsig, float) else nsig
        half_gap = nsig * np.sqrt((sigma_x*np.cos(angle))**2 + (sigma_y*np.sin(angle))**2) if isinstance(nsig, float) else 999.0

        collimator_dict[name] = {
            'name': name,
            'betx': betx,
            'bety': bety,
            'material': material,
            'length': length,
            'angle': angle,
            'sigma_x': sigma_x, 
            'sigma_y': sigma_y,
            'offset': offset,
            'tilt_1': tilt_1,
            'tilt_2': tilt_2,
            'nsig': nsig,
            'half_gap': half_gap
        }

    collimatorList = fb.CollimatorList()
    collimatorList.acquireCollxsuite(collimator_dict)

    args_fb = fb.args_fluka_builder()
    args_fb.geometrical_emittance = geometrical_emittance
    args_fb.collimatorList = collimatorList
    args_fb.prototype_file = 'prototypes.lbp'

    coll_dict = {}
    input_file, coll_dict = fb.fluka_builder(args_fb)

    return input_file, coll_dict

    # _write_xcoll_header_to_fluka_input(input_file, coll_dict)

    # print(collimatorList)

    # # loop over elements and get the collimator settings
    # print("Acquiring collimator settings...")
    # # print(all atributes in line)
    # # tw = line.twiss()
    # collimator_dict = {}
    # for element, name in zip(elements, names):
    #     collimator_dict[name] = {
    #         'name': name,
    #         # 'betx': tw[name].betx[0],
    #         # 'bety': tw[name].bety[0],
    #         'material': line[name].material.name,
    #         'length': line[name].length,
    #         'angle': line[name].angle,
    #         # 'sigma_x': np.sqrt(geometric_emittance * tw[name].betx[0]),
    #         # 'sigma_y': np.sqrt(geometric_emittance * tw[name].bety[0]),
    #         'offset': 0.0, # ??
    #         # 'tilt_1': math.radians(line[name].tilt_L),
    #         'tilt_1': math.radians(line.get(name, {}).get('tilt_L', 0.0)),
    #         # 'tilt_2': math.radians(line[name].tilt_R),  
    #         'tilt_2': math.radians(line.get(name, {}).get('tilt_R', 0.0)),
    #         'nsig': element.gap,
    #         'half_gap': element.gap
    #     }

def _write_xcoll_header_to_fluka_input(input_file, collimator_dict):
    header = ["*  DO NOT CHANGE THIS HEADER", _header_start, "*  {"]
    for kk, vv in collimator_dict.items():
        header.append(f'*  "{kk}": ' + json.dumps(vv).replace('" jaw"', '\n*          "jaw"') + ',')
    header[-1] = header[-1][:-1]  # remove last comma
    header.append("*  }")
    header.append(_header_stop)

    with open(input_file, 'r') as fp:
        data = fp.read()
    with open(input_file, 'w') as fp:
        fp.write("\n".join(header) + "\n*\n" + data)


def create_fluka_input(line, cwd=None):
    elements, names = line.get_elements_of_type(FlukaCollimator)
    if len(elements) == 0:
        raise ValueError('No FlukaCollimator elements found in line')
    input_file = None
    collimator_dict = None
    input_file, collimator_dict = fluka_builder(line)

    # _read_fort3(file) # Need to add jaw as calculated in FLUKA builder
    # ....
    _write_xcoll_header_to_fluka_input(input_file, collimator_dict)
    return input_file


def get_collimators_from_input_file(input_file):
    with open(input_file, 'r') as fp:
        data = fp.read()
    if _header_start not in data or _header_stop not in data:
        raise ValueError("No XCOLL header found in input file. Regenerate input file!")
    commented_dict = data.split(_header_start)[1].split(_header_stop)[0].split('\n')[1:-1]
    cleaned_dict = "".join([val[3:] for val in commented_dict])
    return json.loads(cleaned_dict)


def verify_insertion_file(insertion_file, collimator_dict):
    all_fluka_ids = []
    with open(insertion_file, 'r') as fid:
        for line in fid.readlines():
            all_fluka_ids.append(int(line.split()[0]))
    for name, val in collimator_dict.items():
        if val['fluka_id'] not in all_fluka_ids:
            raise ValueError(f'FlukaCollimator {name} not found in insertion file!')


def _read_fort3(file):
    if not file.exists():
        raise ValueError(f"File {file} not found, cannot deduce FLUKA collimator ids!")
    with file.open('r') as fid:
        lines = fid.readlines()
    collimators = {}
    for line in lines:
        line = [l.strip() for l in line.split()]
        if line[0].startswith('/') or line[0].startswith('!'):
            continue
        collimators[line[0]] = {
            'fluka_id': int(line[2]),
            'length': float(line[3])
        }
    return collimators


    # def _read_gaps(self, input_file):
    #     # TODO: very hacky
    #     with input_file.open('r') as fid:
    #         lines = fid.readlines()
    #     for coll in self._collimators:
    #         try:
    #             idx = [i for i, line in enumerate(lines) if coll.upper() in line][0]
    #             hgap = float([line for line in lines[idx:] if line.startswith('* hGap')][0].split()[3])*1.e-3
    #         except:
    #             print(f"Warning: {coll} not found. Aperture set to 0!")
    #             hgap = 0.0
    #         self._collimators[coll]['jaw'] = hgap

# tmpAuxFile = open ('colldb_lhc_run3.yaml' ,'r'  )
# CDByamlflag = True

# twiss_fname = args.TFile[0]
# tmpOptFile = AuxFile( twiss_fname )
# NREMIT = float( args.EMIn[0] )    # normalised emmittance in [m]
# MOMENTUM = float( args.pc[0] )    # [GeV/c]
# GEOEMIT = NREMIT / (gamma*beta)   # gamma*beta = Pc / m

# # Collision IPs not done!

# BeamLine = TWISS_ELEMENTs()

#     # Prototype list (collimator names and positions in PARKING region)
#     lbp_fname  = "prototypes.lbp" #There is a better way...but enough for now!


#     # output files:
#     #     1- fort.3 insertion point list
#     sixt_fname = outFileName + ".fort3.list"
#     #     2- insertion.txt file
#     injec_fname = outFileName + ".insertion.txt"
#     #     3- FLUKA input file
#     out_fname  = outFileName + ".inp"
#     #     4- Log file
#     log_fname  = outFileName + ".log"


#                 element = TWISS_ELEMENT()
#                 COLLID += 1
#                 element.ID = COLLID
#                 element.NAME = name.upper()
#                 if type(settings['gap']) is str:
#                     if settings['gap'] not in data['Families']:
#                         error_message ( "The gap family setting '" + settings['gap'] + "' for collimator " + name + " is not defined in Families!", True )
#                     element.NSIG = NullGap( data['Families'][ settings['gap'] ]['gap'] )
#                 else:
#                     element.NSIG = NullGap( settings['gap'] )
#                 element.MATERIAL = settings['material']
#                 element.LENGTH = settings['length']
#                 element.ANGLE = math.radians(settings.get('angle',0))
#                 if float(settings.get('offset',0)) != 0:
#                     print ' collimator offset ignored for the time being...'
#                 element.TILT1 = math.radians(settings.get('tilt_right',0))
#                 element.TILT2 = math.radians(settings.get('tilt_left',0))
#                 element.BETX = 0.0
#                 element.BETY = 0.0
#                 element.SIGX = 0.0
#                 element.SIGY = 0.0
#                 element.HALFGAP = 0.0

# EXCLUDED_COLL = [ 'TCRYO' ]\
