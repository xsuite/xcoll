# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

import io
import re
import json
import numpy as np
import pandas as pd
from pathlib import Path

import xtrack as xt

from .beam_elements import BlackAbsorber, BlackCrystal, EverestCollimator, EverestCrystal, \
                           FlukaCollimator, BaseCollimator, BaseCrystal, _all_collimator_classes
from .install import install_elements
from .scattering_routines.everest.materials import SixTrack_to_xcoll
from .scattering_routines.fluka import FlukaEngine


def _initialise_None(dct):
    fields = {'gap': None, 'angle': 0, 'offset': 0, 'parking': 1, 'jaw': None, 'family': None}
    fields.update({'overwritten_keys': [], 'side': 'both', 'material': None, 'stage': None})
    fields.update({'length': 0, 'collimator_type': None, 'active': True, 'crystal': None, 'tilt': 0})
    fields.update({'bending_radius': None, 'bending_angle': None, 'width': 0, 'height': 0, 'miscut': 0})
    for f, val in fields.items():
        if f not in dct.keys():
            dct[f] = val
    for key in dct.keys():
        if key not in fields.keys():
            raise ValueError(f"Illegal setting {key} in collimator!")


def _dict_keys_to_lower(dct):
    if isinstance(dct, dict):
        return {k.lower(): _dict_keys_to_lower(v) for k,v in dct.items()}
    else:
        return dct


def _get_coll_dct_by_beam(coll, beam):
    # The dictionary can be a CollimatorDatabase for a single beam (beam=None)
    # or for both beams (beam='b1' or beam='b2)
    if beam is not None:
        if isinstance(beam, int) or len(beam) == 1:
            beam = f'b{beam}'
        beam = beam.lower()
    beam_in_db = list(coll.keys())

    if beam_in_db == ['b1','b2']:
        if beam is None:
            raise ValueError("Need to specify a beam, because the given dict is for both beams!")
        return coll[beam]

    elif len(beam_in_db) == 1:
        if beam is None:
            beam = beam_in_db[0].lower()
        elif beam != beam_in_db[0].lower():
            raise ValueError("Variable 'beam' does not match beam specified in CollimatorDatabase!")
        return coll[beam]

    elif beam is not None:
        print("Warning: Specified a beam, but the CollimatorDatabase is for a single beam only!")
    return coll


class CollimatorDatabase:

    def __init__(self, **kwargs):
        # Get all arguments
        for key in ['collimator_dict', 'nemitt_x', 'nemitt_y']:
            if key not in kwargs.keys():
                raise ValueError(f"CollimatorDatabase is missing required argument '{key}'!")

        self._parse_dict(kwargs['collimator_dict'],
                         kwargs.get('family_dict', {}),
                         kwargs.get('beam', None),
                         kwargs.get('_yaml_merged', False),
                         kwargs.get('ignore_crystals', True))
        self.nemitt_x = kwargs['nemitt_x']
        self.nemitt_y = kwargs['nemitt_y']
        self._elements = {}


    def _parse_dict(self, coll, fam, beam=None, _yaml_merged=False, ignore_crystals=True):
        # We make all keys case-insensitive to avoid confusion between different conventions
        coll = _dict_keys_to_lower(coll)
        fam  = _dict_keys_to_lower(fam)

        # The dictionary can be a CollimatorDatabase for a single beam (beam=None)
        # or for both beams (beam='b1' or beam='b2)
        coll = _get_coll_dct_by_beam(coll, beam)

        # Apply family settings
        crystals = []
        for thiscoll, settings in coll.items():
            settings = {k.lower(): v for k,v in settings.items()}
            if 'family' in settings.keys() and settings['family'] is not None:
                settings['family'] = settings['family'].lower()
                thisfam = settings['family']
                if thisfam not in fam.keys():
                    raise ValueError(f"Collimator {thiscoll} depends on family {thisfam}, "
                                   + f"but the latter is not defined!")

                # Check if some family settings are overwritten for this collimator
                # Only do this check if we didn't do a YAML merge earlier (because then it
                # is already taken care of)
                if not _yaml_merged:
                    overwritten_keys = [key.lower() for key in settings.keys() if key in fam[thisfam]]
                    if len(overwritten_keys) > 0:
                        settings['overwritten_keys'] = overwritten_keys

                # Load family settings, potentially overwriting settings for this collimator
                settings = {**fam[thisfam], **settings}

            else:
                settings['family'] = None
            coll[thiscoll] = settings

            # Save list of crystals
            if 'crystal' in settings:
                if settings['crystal'] != 0.0:
                    crystals += [thiscoll]
                else:
                    settings['crystal'] = None

        # Remove crystals from colldb
        if ignore_crystals:
            for thiscoll in crystals:
                del coll[thiscoll]

        # Check that all collimators have gap settings
        if not np.all(['gap' in val.keys() or 'jaw' in val.keys() for val in coll.values()]):
            raise ValueError("Ill-defined CollimatorDatabase: Not all collimators have a gap or "
                           + "jaw setting, (or the keys / structure of the dictionary is wrong)!")

        # Update collimators with default values for missing keys
        for name, collimator in coll.items():
            # Change all values to lower case
            for key, val in collimator.items():
                collimator[key] = val.lower() if isinstance(val, str) else val
            _initialise_None(collimator)

        self._collimator_dict = coll
        self._family_dict = fam

    # =======================================
    # ====== Loading/dumping functions ======
    # =======================================

    @classmethod
    def from_yaml(cls, file, **kwargs):
        # Only do the import here, as to not force people to install
        # ruamel if they don't load CollimatorDatabase yaml's
        from ruamel.yaml import YAML
        yaml = YAML(typ='safe')
        if isinstance(file, io.IOBase):
            dct = yaml.load(file)
        else:
            file = Path(file).resolve()
            with file.open('r') as fid:
                dct = yaml.load(fid)
        dct = _dict_keys_to_lower(dct)

        # If the CollimatorDatabase uses YAML merging, we need a bit of hackery to get the
        # family names from the tags (anchors/aliases)
        _yaml_merged = False
        if 'families' in dct.keys() and not isinstance(dct['families'], dict):
            _yaml_merged = True

            # First we load a round-trip yaml
            yaml = YAML()
            if isinstance(file, io.IOBase):
                full_dct = yaml.load(file)
            else:
                with open(file, 'r') as fid:
                    full_dct = yaml.load(fid)
            famkey = [key for key in full_dct.keys() if key.lower() == 'families'][0]
            collkey = [key for key in full_dct.keys() if key.lower() == 'collimators'][0]
            families = {}

            # We loop over family names to get the name tag ('anchor') of each family
            for fam, full_fam in zip(dct['families'], full_dct[famkey]):
                if full_fam.anchor.value is None:
                    raise ValueError("Missing name tag / anchor in "
                                   + "CollimatorDatabase['families']!")
                # We get the anchor from the rt yaml, and use it as key in the families dict
                families[full_fam.anchor.value.lower()] = {f.lower(): v for f, v in fam.items()}
            dct['families'] = families

            # Now we need to loop over each collimator, and verify which family was used
            beam = kwargs.get('beam', None)
            coll_dct = _get_coll_dct_by_beam(dct['collimators'], beam)
            full_coll_dct = _get_coll_dct_by_beam(full_dct[collkey], beam)
            for coll, full_coll in zip(coll_dct.values(), full_coll_dct.values()):
                if not isinstance(coll['gap'], (int,float)):
                    coll['gap'] = None
                if 'family' in coll.keys():
                    raise ValueError(f"Error in {coll}: Cannot use merging for families "
                                    + "and manually specify family as well!")
                elif len(full_coll.merge) > 0:
                    coll['family'] = full_coll.merge[0][1].anchor.value.lower()
                    # Check if some family settings are overwritten for this collimator
                    overwritten_keys = [key.lower() for key in full_coll.keys()
                                        if full_coll._unmerged_contains(key)
                                        and key.lower() in families[coll['family']].keys()]
                    if len(overwritten_keys) > 0:
                        coll['overwritten_keys'] = overwritten_keys
                else:
                    coll['family'] = None

        return cls.from_dict(dct, _yaml_merged=_yaml_merged, **kwargs)


    @classmethod
    def from_json(cls, file, **kwargs):
        if isinstance(file, io.IOBase):
            dct = json.load(file)
        else:
            file = Path(file).resolve()
            with file.open('r') as fid:
                dct = json.load(fid)
        dct = _dict_keys_to_lower(dct)
        return cls.from_dict(dct, **kwargs)


    @classmethod
    def from_dict(cls, dct, beam=None, _yaml_merged=False, nemitt_x=None, nemitt_y=None, ignore_crystals=True):
        # We make all keys case-insensitive to avoid confusion between different conventions
        # The families are optional
        fam = {}
        dct = _dict_keys_to_lower(dct)

        # Get the emittance
        if nemitt_x is None and nemitt_y is None:
            if 'emittance' not in dct.keys():
                raise ValueError("Missing emittance info! Add 'emittance' as a key to "
                               + "the CollimatorDatabase file, or specify it as 'nemitt_x' "
                               + "and 'nemitt_y' to the loader!")
            nemitt_x = dct['emittance']['x']
            nemitt_y = dct['emittance']['y']
        elif nemitt_x is None or nemitt_y is None:
            raise ValueError("Need to provide both 'nemitt_x' and 'nemitt_y'!")
        elif 'emittance' in dct.keys():
            if dct['emittance']['x'] != nemitt_x or dct['emittance']['y'] != nemitt_y:
                raise ValueError("Emittance in CollimatorDatabase file different from "
                               + "'nemitt_x' and 'nemitt_y'!")

        # Get family and collimator dicts
        if 'families' in dct.keys():
            if not 'collimators' in dct.keys():
                raise ValueError("Could not find 'collimators' dict in CollimatorDatabase!")
            fam  = dct['families']
            coll = dct['collimators']
        elif 'collimators' in dct.keys():
            coll = dct['collimators']
        else:
            coll = dct

        return cls(collimator_dict=coll, family_dict=fam, nemitt_x=nemitt_x, nemitt_y=nemitt_y,
                   beam=beam, _yaml_merged=_yaml_merged, ignore_crystals=ignore_crystals)


    @classmethod
    def from_SixTrack(cls, file, ignore_crystals=True, **kwargs):
        file = Path(file).resolve()
        with file.open('r') as fp:
            coll_data_string = ''
            family_settings = {}
            family_types = {}
            side = {}
            cry_fields = ['bending_radius', 'width', 'height', 'thick', 'tilt', 'miscut', 'crystal']
            cry = {}

            for line in fp:
                if line.startswith('#'):
                    continue # Comment
                sline = line.split()
                if len(sline) > 0:
                    if sline[0].lower() == 'nsig_fam':
                        family_settings[sline[1]] = float(sline[2])
                        family_types[sline[1]] = sline[3]
                    elif sline[0].lower() == 'onesided':
                        side[sline[1]] = int(sline[2])
                    elif sline[0].lower() == "crystal":
                        cry[sline[1]] = {}
                        for i, key in enumerate(cry_fields):
                            idx = i+2
                            if i < 6:
                                cry[sline[1]][key] = float(sline[idx])
                            else:
                                cry[sline[1]][key] = int(sline[idx])
                    elif sline[0].lower() == 'settings':
                        pass # Acknowledge and ignore this line
                    elif len(sline) == 6:
                        # Standard collimator definition
                        coll_data_string += line
                    else:
                        print(f"Unknown setting {line}")

        defaults = {}
        _initialise_None(defaults)
        defaults['thick'] = 0

        famdct = {key: {'gap': None if family_settings[key] > 900 else family_settings[key],
                        'stage': family_types[key]} for key in family_settings}
        names = ['name', 'gap', 'material', 'length', 'angle', 'offset']

        df = pd.read_csv(io.StringIO(coll_data_string), sep=r'\s+', index_col=False, names=names)
        df['family'] = df['gap'].copy()
        df['family'] = df['family'].apply(lambda s: None if re.match(r'^-?\d+(\.\d+)?$', str(s)) else s)
        df.insert(5,'stage', df['gap'].apply(lambda s: None if s in family_types else 'UNKNOWN'))

        df['gap'] = df['gap'].apply(lambda s: None if not isinstance(s, str) and s > 900 else s)
        df['gap'] = df['gap'].apply(lambda s: None if isinstance(s, str) else s)

        # TODO this breaks code if a key has upper case, e.g. gap_L
        df['name'] = df['name'].str.lower() # Make the names lowercase for easy processing
        df['parking'] = 0.025
        if ignore_crystals:
            df = df[~df.name.isin(list(cry.keys()))]
        else:
            for key in cry_fields:
                df[key] = [cry[name][key] if name in cry else defaults[key]
                           for name in df['name']]
            if not np.allclose(np.unique(df.thick.values), 0):
                print("Warning: Keyword 'thick' is currently not supported in xcoll! Ignoring.")
            df = df.drop('thick', axis=1)
            df['crystal'] = ['strip' if s==1 else s for s in df['crystal']]
            df['crystal'] = ['quasi-mosaic' if s==2 else s for s in df['crystal']]
        df['side'] = [side[name] if name in side else defaults['side']
                      for name in df['name']]
        df['side'] = ['both'  if s==0 else s for s in df['side']]
        df['side'] = ['left'  if s==1 else s for s in df['side']]
        df['side'] = ['right' if s==2 else s for s in df['side']]
        if not np.allclose(np.unique(df.offset.values), 0):
            print("Warning: Keyword 'offset' is currently not supported in xcoll! Ignoring.")
        df = df.drop('offset', axis=1)
        df = df.set_index('name')
        df = df.replace(np.nan, None)

        colldict = df.transpose().to_dict()
        # Remove Nonetype families
        colldict = {coll: {kk: vv for kk, vv in coll_settings.items()
                           if kk != 'family' or vv is not None}
                    for coll, coll_settings in colldict.items()}
        # Remove None gaps and stages if the collimator is assigned a family (they will be set in _parse_dict)
        colldict = {coll: {kk: vv for kk, vv in coll_settings.items()
                           if kk != 'gap' or vv is not None or 'family' not in coll_settings}
                    for coll, coll_settings in colldict.items()}
        colldict = {coll: {kk: vv for kk, vv in coll_settings.items()
                           if kk != 'stage' or vv is not None or 'family' not in coll_settings}
                    for coll, coll_settings in colldict.items()}

        return cls.from_dict({'collimators': colldict, 'families': famdct}, \
                             ignore_crystals=ignore_crystals, **kwargs)

    def to_pandas(self):
        return pd.DataFrame(self._collimator_dict).transpose()

    def to_yaml(self, out, lhc_style=True):
        """
        Writes a colldb in memory to disk in the yaml format.

        > colldb_object.write_to_yaml(<path+name>, lhc_style=Bool)

        if lhc_style == True, it will add comments assuming that the collimators are named
        as in the lhc.

        The function can dump b1, b2 and a general bx, however multi-beam functionality is not yet
        added to the collmanager. TODO

        If any of the dumped keys contains capital letters (e.g. gap_L), it will not be possible
        to load it back into xcoll, since all keys are set to lowercase when importing TODO
        """
        # Dumps collimator database to a YAML file with optional LHC style formatting
        import re

        # Local helper functions
        def _format_dict_entry(key, value, spacing='', mapping=False, key_width=15):
            # Formats a dictionary entry into a string for YAML output
            formatted_values = ',    '.join(f"{k}: {v}" for k, v in value.items())
            formatted_values = re.sub(r'none', 'null', formatted_values, flags=re.IGNORECASE)
            # Ensure key has a fixed width for alignment
            if mapping:
                formatted_key = f'{key}'.ljust(key_width)
            else:
                formatted_key = f'{key}:'.ljust(key_width)
            #formatted_values = formatted_values.ljust(key_width)
            return f"{spacing}{formatted_key} {{ {formatted_values} }}\n"

        def _print_values(keys, dct, file, spacing='', mapping=False):
            # Writes formatted dictionary entries to a file
            for key in keys:
                file.write(_format_dict_entry(key, dct[key], spacing=spacing, mapping=mapping))

        def _print_colls(colls, dcts, beam, file):
            # Filters and formats collimator data, then writes to a file
            coll_items_to_print = ['<<','gap','angle','material','active','length','side']
            file.write(f'  {beam}:\n')
            for coll in colls:
                coll_dict = dcts.to_pandas().transpose().to_dict()[coll]
                fam = coll_dict['family']
                fam_keys = []
                if fam is not None:
                    fam_keys = dcts._family_dict[fam].keys()
                    coll_dict = {**{'<<': '*'+fam}, **coll_dict}
                temp_items_to_print = []
                if coll_dict['crystal'] and str(coll_dict['crystal'])!='nan':
                    temp_items_to_print = ['bending_radius','width','height','miscut','crystal']
                # if 'angle_L' in coll_dict and coll_dict['angle_L'] == coll_dict['angle_R']:
                #     coll_dict.update({'angle': coll_dict['angle_L']})
                # else:
                #     temp_items_to_print = temp_items_to_print + ['angle_L','angle_R']
                # if coll_dict['gap_L'] == coll_dict['gap_R']:
                #     coll_dict.update({'gap': coll_dict['gap_L']})
                # elif coll_dict['gap_L'] is None and coll_dict['gap_R'] is not None:
                #     coll_dict.update({'gap': coll_dict['gap_R']})
                # elif coll_dict['gap_L'] is not None and coll_dict['gap_R'] is None:
                #     coll_dict.update({'gap': coll_dict['gap_L']})
                # else:
                #     temp_items_to_print = temp_items_to_print + ['gap_L','gap_R']
                value = {}
                overwritten_keys = coll_dict['overwritten_keys']
                for key, val in coll_dict.items():
                    if key == 'active_length':
                        key = 'length' 
                    if (key in coll_items_to_print+temp_items_to_print) and (key not in (set(fam_keys)-set(overwritten_keys))) and (val != 'both'):
                        value.update({key: val})
                file.write(_format_dict_entry(coll, value, spacing='    '))
            file.write('\n')

        LHC_families = ['tcp3', 'tcsg3', 'tcsm3', 'tcla3', 'tcp7', 'tcsg7', 'tcsm7', 'tcla7', 'tcli', 'tdi', 'tcdq', 'tcstcdq', 'tcth1', 'tcth2', 'tcth5', 'tcth8', 'tctv1', 'tctv2', 'tctv5', 'tctv8', 'tclp', 'tcxrp', 'tcryo', 'tcl4', 'tcl5', 'tcl6', 'tct15', 'tct2', 'tct8', 'tcsp', 'tcld']
        with open(f'{out}.yaml', 'w') as file:
            if '_family_dict' in self.__dict__.keys():
                file.write('families:\n')
                if lhc_style:
                    printed_families = []
                    fams_in_dict = self._family_dict.keys()

                    # Momentum cleaning
                    file.write('  # Momentum cleaning\n')
                    sel_fam = [fam for fam in LHC_families if re.match('.*3', fam) and (fam in fams_in_dict)]
                    printed_families += sel_fam
                    _print_values(sel_fam, self._family_dict, file, spacing='  - &', mapping=True)

                    # Betatron cleaning
                    file.write('  # Betatron cleaning\n')
                    sel_fam = [fam for fam in LHC_families if re.match('.*7', fam) and (fam in fams_in_dict)]
                    printed_families += sel_fam
                    _print_values(sel_fam, self._family_dict, file, spacing='  - &', mapping=True)

                    # Injection protection
                    file.write('  # Injection protection\n')
                    sel_fam = [fam for fam in LHC_families if (fam in ['tcli', 'tdi']) and (fam in fams_in_dict)]
                    printed_families += sel_fam
                    _print_values(sel_fam, self._family_dict, file, spacing='  - &', mapping=True)

                    # Dump protection
                    file.write('  # Dump protection\n')
                    sel_fam = [fam for fam in LHC_families if (fam in ['tcdq', 'tcsp', 'tcstcdq']) and (fam in fams_in_dict)]
                    printed_families += sel_fam
                    _print_values(sel_fam, self._family_dict, file, spacing='  - &', mapping=True)

                    # Physics background / debris
                    file.write('  # Physics background / debris\n')
                    sel_fam = [fam for fam in LHC_families if ((re.match('tc[lt][0-9dp].*', fam)) or (fam in ['tcryo', 'tcxrp'])) and (fam in fams_in_dict)]
                    printed_families += sel_fam
                    _print_values(sel_fam, self._family_dict, file, spacing='  - &', mapping=True)

                    # Other families
                    if set(printed_families) != set(fams_in_dict):
                        file.write('  # Other families\n')
                        _print_values(set(fams_in_dict) - set(printed_families), self._family_dict, file, spacing='  - &', mapping=True)
                else:
                    file.write('  # Families\n')
                    _print_values(self._family_dict.keys(), self._family_dict, file, spacing='  - &', mapping=True)

            # Emittance section
            ex = self.nemitt_x
            ey = self.nemitt_y
            file.write(f'\nemittance:\n  x: {ex}\n  y: {ey}\n')

            # Collimators section
            file.write('\ncollimators:\n')
            b1_colls, b2_colls, bx_colls = [], [], []
            for coll in self.to_pandas().index:
                if coll == 'tclia.4r2' or coll == 'tclia.4l8':    # TODO: hardcoded!!!
                    b1_colls.append(coll)
                    b2_colls.append(coll)
                elif coll[-2:] == 'b1':
                    b1_colls.append(coll)
                elif coll[-2:] == 'b2':
                    b2_colls.append(coll)
                else:
                    bx_colls.append(coll)

            # Handle special cases for collimators
            if (('tclia.4r2' in b1_colls) or ('tclia.4l8' in b1_colls)) and (len(b1_colls) <= 2):
                b1_colls = []
            if (('tclia.4r2' in b2_colls) or ('tclia.4l8' in b2_colls)) and (len(b2_colls) <= 2):
                b2_colls = []

            # Print collimators for each beam
            if len(b1_colls) > 0:
                _print_colls(b1_colls, self, 'b1', file)
            if len(b2_colls) > 0:
                _print_colls(b2_colls, self, 'b2', file)
            if len(bx_colls) > 0:
                _print_colls(bx_colls, self, 'bx', file)
                print('WARNING -- some collimators could not be assigned to b1 or b2. Tracking might not work with those collimators. Please manually change the output file if necessary.')


    # ====================================
    # ====== Installing collimators ======
    # ====================================

    def _get_names_from_line(self, line, names, families):
        if names is None and families is None:
            names = self.collimator_names
        elif names is None:
            names = self.get_collimators_from_family(families)
        elif families is not None:
            names.append(self.get_collimators_from_family(families))
        return list(set(names)) # Remove duplicates

    def _check_installed(self, line, name, collimator_class):
            # Check that collimator is not installed as different type
            # TODO: automatically replace collimator type and print warning
            if isinstance(line[name], _all_collimator_classes):
                raise ValueError(f"Trying to install {name} as {collimator_class.__name__}, "
                               + f"but it is already installed as {line[name].__class__.__name__}!\n"
                               + f"Please reconstruct the line.")
            # TODO: only allow Marker elements, no Drifts!!
            #       How to do this with importing a line for MAD-X or SixTrack...?
            #       Maybe we want a DriftCollimator type in Xtrack as a general placeholder
            elif not isinstance(line[name], (xt.Marker, xt.Drift)):
                raise ValueError(f"Trying to install {name} as {collimator_class.__name__}, "
                               + f"but the line element to replace is not an xtrack.Marker "
                               + f"(or xtrack.Drift)!\nPlease check the name, or correct the "
                               + f"element.")

    def _create_collimator(self, line, collimator_class, name, **kwargs):
        assert issubclass(collimator_class, BaseCollimator)
        self._check_installed(line, name, collimator_class)
        if kwargs.pop('verbose', False):
            print(f"Installing {name:20} as {collimator_class.__name__}")
        el = collimator_class(gap=self[name]['gap'], angle=self[name]['angle'],
                              length=self[name]['length'], side=self[name]['side'],
                              _tracking=False, **kwargs)
        el.emittance = [self.nemitt_x, self.nemitt_y]
        self._elements[name] = el

    def _create_crystal(self, line, crystal_class, name, **kwargs):
        assert issubclass(crystal_class, BaseCrystal)
        self._check_installed(line, name, crystal_class)
        if kwargs.pop('verbose', False):
            print(f"Installing {name:20} as {crystal_class.__name__}")
        el = crystal_class(gap=self[name]['gap'], angle=self[name]['angle'],
                           length=self[name]['length'], side=self[name]['side'],
                           bending_radius=self[name]['bending_radius'],
                           width=self[name]['width'], height=self[name]['height'],
                           _tracking=False, **kwargs)
        el.emittance = [self.nemitt_x, self.nemitt_y]
        self._elements[name] = el

    def install_black_absorbers(self, line, *, names=None, families=None, verbose=False, need_apertures=True):
        names = self._get_names_from_line(line, names, families)
        for name in names:
            if self[name]['bending_radius'] is None:
                self._create_collimator(line, BlackAbsorber, name, verbose=verbose)
            else:
                self._create_crystal(line, BlackCrystal, name, verbose=verbose)
        elements = [self._elements[name] for name in names]
        install_elements(line, names, elements, need_apertures=need_apertures)

    def install_everest_collimators(self, line, *, names=None, families=None, verbose=False, need_apertures=True):
        names = self._get_names_from_line(line, names, families)
        for name in names:
            mat = SixTrack_to_xcoll(self[name]['material'])
            if self[name]['bending_radius'] is None:
                self._create_collimator(line, EverestCollimator, name, material=mat[0],
                                        verbose=verbose)
            else:
                self._create_crystal(line, EverestCrystal, name, material=mat[1],
                                     lattice=self[name]['crystal'], verbose=verbose,
                                     miscut=self[name]['miscut'])
        elements = [self._elements[name] for name in names]
        install_elements(line, names, elements, need_apertures=need_apertures)

    def install_fluka_collimators(self, line, *, names=None, families=None, verbose=False, need_apertures=True,
                                  fluka_input_file=None, remove_missing=True):
        # Check server
        if FlukaEngine.is_running():
            print("Warning: FLUKA server is already running. Stopping server to install collimators.")
            FlukaEngine.stop_server(fluka_input_file)
        names = self._get_names_from_line(line, names, families)
        for name in names:
            self._create_collimator(line, FlukaCollimator, name, verbose=verbose)
        elements = [self._elements[name] for name in names]
        install_elements(line, names, elements, need_apertures=need_apertures)

    # ==================================
    # ====== Accessing attributes ======
    # ==================================

    @property
    def collimator_names(self):
        return list(self._collimator_dict.keys())

    @property
    def collimator_families(self):
        families = {fam: [] for fam in self._family_dict.keys()}
        families["no family"] = []
        for name in self.collimator_names:
            if 'family' not in self[name] or self[name]['family'].lower() == 'unknown':
                families["no family"].append(name)
            else:
                families[self[name]['family']].append(name)
        return families

    def get_collimators_from_family(self, family):
        if not hasattr(family, '__iter__') and not isinstance(family, str):
            family = [family]
        result = []
        for fam in family:
            if fam not in self.collimator_families:
                raise ValueError(f"Family '{fam}' not found in CollimatorDatabase!")
            result += self.collimator_families[fam]
        return result

    @property
    def properties(self):
        return {attr for d in self._collimator_dict.values() for attr in d.keys()}

    def __getattr__(self, attr):
        if attr in self.properties:
            # TODO: include families
            return {kk: vv.get(attr, None) for kk, vv in self._collimator_dict.items()}
        else:
            raise ValueError(f"Property `{attr}` not present in CollimatorDatabase!")

    def __getitem__(self, name):
        if name in self._family_dict:
            return self._family_dict[name]
        elif name in self._collimator_dict:
            return self._collimator_dict[name]
        else:
            raise ValueError(f"Family nor collimator `{name}` found in CollimatorDatabase!")

