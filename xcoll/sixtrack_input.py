# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

import numpy as np
import pathlib
import re

import xobjects as xo
import xcoll as xc
import xtrack as xt
from .beam_elements.k2 import _K2Collimator
from .beam_elements.k2 import EverestCrystal

def _get_spaces(str1,str2, total_spaces):
    length1 = len(str1)
    length2 = len(str2)
    middle_spaces = total_spaces - length2 - length1
    return middle_spaces

def create_dat_file(line, path_out, elements=None, names=None):
    # check if families
    # check which collimators we have in which family
    # write to file; file.write(f'  {beam}:\n')

    # get collimators
    # get opening (gap), material, length, angle, offset 
    # write to file

    # check for old collimators that are not in use 
    # write to file

    # check additional info: onesided or crystal (extra data)
    # write to tfile
    # Local helper functions
    with open(f'{path_out}.dat', 'w') as file:
        onesided = []
        crystal  = []
        for name in names:
            if line[name].side != 'both':
                onesided.append(name)
                continue
            if line[name].__class__ == EverestCrystal:
                crystal.append(name)
                continue
            
            gap = str(line[name].gap)
            if gap == str(None):
                gap = "null"
            mat = line[name].material
            length = line[name].length
            angle = line[name].angle
            offset = (line[name].jaw_L + line[name].jaw_R)/2 
            file.write(f"{name}" + " " * _get_spaces(name,gap, 65) + f" {mat}" \
                       + " " * _get_spaces(mat, "1.111", 15) + f"{length}" \
                       + " " * _get_spaces(length, angle, 14) + f"{angle}" \
                       + " " * _get_spaces(angle, offset, 15) + f"{offset} \n")
        if len(onesided) > 0:
            for name in onesided:
                file.write("SETTINGS \n")
                file.write(f"ONESIDED {name}" + " " * _get_spaces(name,"1", 50) \
                           + f"{line[name]._side} \n")
        if len(crystal) > 0:
            if len(onesided) == 0:
                file.write("SETTINGS \n")
            for name in crystal:
                file.write(f"CRYSTAL {name}" + f" {line[name].bending_radius}"\
                           + f" {line[name].width}" + f" {line[name].height}" \
                           + f" {line[name].miscut}" + f" {line[name].width}")
    return file

#  ['bending_radius', 'width', 'height', 'thick', 'miscut', 'crystal']
def get_collimators_from_input_file(input_file):
    with open(input_file,'r') as fid:
        data = fid.read()
    pass





    