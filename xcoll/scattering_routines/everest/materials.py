# copyright ############################### #
# This file is part of the Xcoll Package.  #
# Copyright (c) CERN, 2023.                 #
# ######################################### #

import xobjects as xo

# xo.Strings are complicated
# They come in predefined size of 7, 15, 23, 31, ... characters
# Once the string is assigned for the first time, the max size is defined,
# and extending the string is not possible.
# This means that, e.g. if we create a K2Collimator with Carbon and later change
# that to Aluminium, the name of the latter will be clipped to 'Alumini' (7 chars)


# CLASS INHERITANCE DOES NOT WORK WITH DYNAMIC XOFIELDS
# Because then the memory offsets between parent and child might no longer be the same,
# as xobjects enforces a strict order: static fields first, and then dynamic fields.
# See struct.py, in __new__ of MetaStruct

class GeneralMaterial(xo.HybridClass):
    _xofields = {
        'Z':                        xo.Float64,     # zatom
        'A':                        xo.Float64,     # anuc
        'density':                  xo.Float64,     # rho (g cm-3)
        'excitation_energy':        xo.Float64,     # exenergy
        'nuclear_radius':           xo.Float64,     # emr
        'nuclear_elastic_slope':    xo.Float64,     # bnref (g cm-2)
                # slope from https://journals.aps.org/prd/pdf/10.1103/PhysRevD.21.3010
        'cross_section':            xo.Float64[6],  # csref [barn]
                # Index 0:Total, 1:absorption, 2:nuclear elastic, 3:pp or pn elastic
                #       4:Single Diffractive pp or pn, 5:Coulomb for t above mcs
        'hcut':                     xo.Float64,
        'name':                     xo.String
    }

    def __init__(self, **kwargs):
        if '_xobject' not in kwargs:
            kwargs.setdefault('hcut', 0.02)
            kwargs.setdefault('cross_section', [0., 0., 0., 0., 0., 0.])
            kwargs.setdefault('name', "NO NAME")
            kwargs['name'] = kwargs['name'].ljust(55)  # Pre-allocate 64 byte using whitespace
        super().__init__(**kwargs)
        self.name = self.name.strip()


class Material(GeneralMaterial):
    _xofields = { **GeneralMaterial._xofields,
        'radiation_length':         xo.Float64      # radl
    }

    _depends_on = [GeneralMaterial]

    def __init__(self, **kwargs):
        super().__init__(**kwargs)



class CrystalMaterial(GeneralMaterial):
    _xofields = { **GeneralMaterial._xofields,
        'crystal_radiation_length': xo.Float64,     # dlri
        'crystal_nuclear_length':   xo.Float64,     # dlyi
        'crystal_plane_distance':   xo.Float64,     # ai  [mm]
        'crystal_potential':        xo.Float64,     # eum  [eV]
        'nuclear_collision_length': xo.Float64      # collnt [m]
    }

    _depends_on = [GeneralMaterial]

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    @classmethod
    def from_material(cls, material, **kwargs):
        thisdict = material.to_dict()
        thisdict.pop('radiation_length')    # not in crystals
        thisdict.update(kwargs)
        return cls(**thisdict)





# BE
Beryllium = Material(
        name = 'Beryllium',
        Z = 4.00,
        A = 9.01,
        density = 1.848,
        excitation_energy = 63.7e-9,
        radiation_length = 0.353,
        nuclear_radius = 0.22,
        nuclear_elastic_slope = 74.7,
        cross_section = [0.271, 0.192, 0, 0, 0, 0.0035e-2]
)

# AL
Aluminium = Material(
        name = 'Aluminium',
        Z = 13.00,
        A = 26.98,
        density = 2.70,
        excitation_energy = 166.0e-9,
        radiation_length = 0.089,
        nuclear_radius = 0.302,
        nuclear_elastic_slope = 120.3,
        cross_section = [0.643, 0.418, 0, 0, 0, 0.0340e-2]
)

# CU
Copper = Material(
        name = 'Copper',
        Z = 29.00,
        A = 63.55,
        density = 8.96,
        excitation_energy = 322.0e-9,
        hcut = 0.01,
        nuclear_radius = 0.366,
        radiation_length = 0.0143,
        nuclear_elastic_slope = 217.8,
        cross_section = [1.253, 0.769, 0, 0, 0, 0.1530e-2]
)

# W
Tungsten = Material(
        name = 'Tungsten',
        Z = 74.00,
        A = 183.85,
        density = 19.30,
        excitation_energy = 727.0e-9,
        hcut = 0.01,
        nuclear_radius = 0.520, 
        radiation_length = 0.0035,
        nuclear_elastic_slope = 440.3,
        cross_section = [2.765, 1.591, 0, 0, 0, 0.7680e-2]
)
TungstenCrystal = CrystalMaterial.from_material( Tungsten,
        crystal_radiation_length = 0.0035,
        crystal_nuclear_length = 0.096,
        crystal_plane_distance = 0.56e-7,
        crystal_potential = 21.0,
        nuclear_collision_length = 0
)

# PB
Lead = Material(
        name = 'Lead',
        Z = 82.00,
        A = 207.19,
        density = 11.35,
        excitation_energy = 823.0e-9,
        hcut = 0.01,
        nuclear_radius = 0.542,
        radiation_length = 0.0056,
        nuclear_elastic_slope = 455.3,
        cross_section = [3.016, 1.724, 0, 0, 0, 0.9070e-2]
)

# C
Carbon = Material(
        name = 'Carbon',
        Z = 6.00,
        A = 12.01,
        density = 1.67,
        excitation_energy = 78.0e-9,
        radiation_length = 0.2557,
        nuclear_radius = 0.25,
        nuclear_elastic_slope = 70.0,
        cross_section = [0.337, 0.232, 0, 0, 0, 0.0076e-2]
)
CarbonCrystal = CrystalMaterial.from_material( Carbon,
        crystal_radiation_length = 0.188,
        crystal_nuclear_length = 0.400,
        crystal_plane_distance = 0.63e-7,
        crystal_potential = 21.0,
        nuclear_collision_length = 0
)

# C2
Carbon2 = Material(
        name = 'Carbon2',
        Z = 6.00,
        A = 12.01,
        density = 4.52,
        excitation_energy = 78.0e-9,
        radiation_length = 0.094,
        nuclear_radius = 0.25,
        nuclear_elastic_slope = 70.0,
        cross_section = [0.337, 0.232, 0, 0, 0, 0.0076e-2]
)

# Si
Silicon = Material(
        name = 'Silicon',
        Z = 14.00,
        A = 28.08,
        density = 2.33,
        excitation_energy = 173.0e-9,
        radiation_length = 1,
        nuclear_radius = 0.441,
        nuclear_elastic_slope = 120.14,
        cross_section = [0.664, 0.430, 0, 0, 0, 0.0390e-2]
)
SiliconCrystal = CrystalMaterial.from_material( Silicon,
        crystal_radiation_length = 0.0937,
        crystal_nuclear_length = 0.4652,
        crystal_plane_distance = 0.96e-7,
        crystal_potential = 21.34,
        nuclear_collision_length = 0.3016
)

# Ge
Germanium = Material(
        name = 'Germanium',
        Z = 32.00,
        A = 72.63,
        density = 5.323,
        excitation_energy = 350.0e-9,
        radiation_length = 1,
        nuclear_radius = 0.605,
        nuclear_elastic_slope = 226.35,
        cross_section = [1.388, 0.844, 0, 0, 0, 0.1860e-2]
)
GermaniumCrystal = CrystalMaterial.from_material( Germanium,
        crystal_radiation_length = 0.02302,
        crystal_nuclear_length = 0.2686,
        crystal_plane_distance = 1.0e-7,
        crystal_potential = 40.0,
        nuclear_collision_length = 0.1632
)

# MoGR
MolybdenumGraphite = Material(
        name = 'MolybdenumGraphite',
        Z = 6.65,
        A = 13.53,
        density = 2.500,
        excitation_energy = 87.1e-9,
        radiation_length = 0.1193,
        nuclear_radius = 0.25,
        nuclear_elastic_slope = 76.7,
        cross_section = [0.362, 0.247, 0, 0, 0, 0.0094e-2]
)

# CuCD
CopperDiamond = Material(
        name = 'CopperDiamond',
        Z = 11.90,
        A = 25.24,
        density = 5.40,
        excitation_energy = 152.9e-9,
        radiation_length = 0.0316,
        nuclear_radius = 0.308,
        nuclear_elastic_slope = 115.0,
        cross_section = [0.572, 0.370, 0, 0, 0, 0.0279e-2]
)

# Mo
Molybdenum = Material(
        name = 'Molybdenum',
        Z = 42.00,
        A = 95.96,
        density = 10.22,
        excitation_energy = 424.0e-9,
        radiation_length = 0.0096,
        nuclear_radius = 0.481,
        nuclear_elastic_slope = 273.9,
        cross_section = [1.713, 1.023, 0, 0, 0, 0.2650e-2]
)

# Glid
Glidcop = Material(
        name = 'Glidcop',
        Z = 28.80,
        A = 63.15,
        density = 8.93,
        excitation_energy = 320.8e-9,
        radiation_length = 0.0144,
        nuclear_radius = 0.418,
        nuclear_elastic_slope = 208.7,
        cross_section = [1.246, 0.765, 0, 0, 0, 0.1390e-2]
)

# Iner
Inermet = Material(
        name = 'Inermet',
        Z = 67.70,
        A = 166.70,
        density = 18.00,
        excitation_energy = 682.2e-9,
        radiation_length = 0.00385,
        nuclear_radius = 0.578,
        nuclear_elastic_slope = 392.1,
        cross_section = [2.548, 1.473, 0, 0, 0, 0.5740e-2]
)


SixTrack_to_xcoll = {
    "be":   [Beryllium],
    "al":   [Aluminium],
    "cu":   [Copper],
    "w":    [Tungsten, TungstenCrystal],
    "pb":   [Lead],
    "c":    [Carbon, CarbonCrystal],
    "c2":   [Carbon2],
    "si":   [Silicon, SiliconCrystal],
    "ge":   [Germanium, GermaniumCrystal],
    "mogr": [MolybdenumGraphite],
    "cucd": [CopperDiamond],
    "mo":   [Molybdenum],
    "glid": [Glidcop],
    "iner": [Inermet]
}
