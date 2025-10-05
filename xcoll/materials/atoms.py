# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

from xtrack.particles.pdg import get_element_name_from_Z, get_element_full_name_from_Z

from .material import Material
from .database import db


# Atomic mass from https://iupac.qmul.ac.uk/AtWt/ rounded to 5 digits
# Densities (in g/cm3), temperature (in Kelvin), and pressure (in atm = 101.325 kPa) from:
#     https://en.wikipedia.org/wiki/List_of_chemical_elements
#     https://en.wikipedia.org/wiki/Densities_of_the_elements_(data_page)
Hydrogen      = Material(Z=  1,   A=1.0080,   state='gas',   density=89.88e-6,  excitation_energy=19.2, temperature=273.15, pressure=1)
Helium        = Material(Z=  2,   A=4.0026,   state='gas',   density=178.5e-6,  excitation_energy=41.8, temperature=273.15, pressure=1)
Lithium       = Material(Z=  3,   A=6.94,     state='solid', density=0.534,     excitation_energy=40,   temperature=293.15)
Beryllium     = Material(Z=  4,   A=9.0122,   state='solid', density=1.848,     excitation_energy=63.7, temperature=293.15)
Boron         = Material(Z=  5,   A=10.81,    state='solid', density=2.34,      excitation_energy=76,   temperature=293.15)
Carbon        = Material(Z=  6,   A=12.011,   state='solid', density=2.265,     excitation_energy=81,   temperature=293.15)  # Amorphous
Nitrogen      = Material(Z=  7,   A=14.007,   state='gas',   density=0.0012506, excitation_energy=82,   temperature=273.15, pressure=1)
Oxygen        = Material(Z=  8,   A=15.999,   state='gas',   density=0.001429,  excitation_energy=95,   temperature=273.15, pressure=1)
Fluorine      = Material(Z=  9,   A=18.998,   state='gas',   density=0.001696,  excitation_energy=115,  temperature=273.15, pressure=1)
Neon          = Material(Z= 10,   A=20.180,   state='gas',   density=0.0009002, excitation_energy=137,  temperature=273.15, pressure=1)
Sodium        = Material(Z= 11,   A=22.990,   state='solid', density=0.968,     excitation_energy=149,  temperature=293.15)
Magnesium     = Material(Z= 12,   A=24.305,   state='solid', density=1.738,     excitation_energy=156,  temperature=293.15)
Aluminium     = Material(Z= 13,   A=26.982,   state='solid', density=2.70,      excitation_energy=166,  temperature=293.15)
Silicon       = Material(Z= 14,   A=28.085,   state='solid', density=2.3296,    excitation_energy=173,  temperature=293.15)
Phosphorus    = Material(Z= 15,   A=30.974,   state='solid', density=2.2,       excitation_energy=173,  temperature=293.15)  # Red phosphorus
Sulfur        = Material(Z= 16,   A=32.06,    state='solid', density=2.067,     excitation_energy=180,  temperature=293.15)
Chlorine      = Material(Z= 17,   A=35.45,    state='gas',   density=0.003214,  excitation_energy=174,  temperature=273.15, pressure=1)
Argon         = Material(Z= 18,   A=39.948,   state='gas',   density=0.001784,  excitation_energy=188,  temperature=273.15, pressure=1)
Potassium     = Material(Z= 19,   A=39.098,   state='solid', density=0.862,     excitation_energy=190,  temperature=293.15)
Calcium       = Material(Z= 20,   A=40.078,   state='solid', density=1.54,      excitation_energy=191,  temperature=293.15)
Scandium      = Material(Z= 21,   A=44.956,   state='solid', density=2.99,      excitation_energy=216,  temperature=293.15)
Titanium      = Material(Z= 22,   A=47.867,   state='solid', density=4.506,     excitation_energy=233,  temperature=293.15)
Vanadium      = Material(Z= 23,   A=50.942,   state='solid', density=6.11,      excitation_energy=245,  temperature=293.15)
Chromium      = Material(Z= 24,   A=51.996,   state='solid', density=7.19,      excitation_energy=257,  temperature=293.15)
Manganese     = Material(Z= 25,   A=54.938,   state='solid', density=7.26,      excitation_energy=272,  temperature=293.15)
Iron          = Material(Z= 26,   A=55.845,   state='solid', density=7.874,     excitation_energy=286,  temperature=293.15)
Cobalt        = Material(Z= 27,   A=58.933,   state='solid', density=8.90,      excitation_energy=297,  temperature=293.15)
Nickel        = Material(Z= 28,   A=58.693,   state='solid', density=8.908,     excitation_energy=311,  temperature=293.15)
Copper        = Material(Z= 29,   A=63.546,   state='solid', density=8.96,      excitation_energy=322,  temperature=293.15)
Zinc          = Material(Z= 30,   A=65.38,    state='solid', density=7.134,     excitation_energy=330,  temperature=293.15)
Gallium       = Material(Z= 31,   A=69.723,   state='solid', density=5.91,      excitation_energy=334,  temperature=293.15)
Germanium     = Material(Z= 32,   A=72.630,   state='solid', density=5.323,     excitation_energy=350,  temperature=293.15)
Arsenic       = Material(Z= 33,   A=74.922,   state='solid', density=5.727,     excitation_energy=347,  temperature=293.15)
Selenium      = Material(Z= 34,   A=78.971,   state='solid', density=4.809,     excitation_energy=348,  temperature=293.15)
Bromine       = Material(Z= 35,   A=79.904,   state='liquid',density=3.1028,    excitation_energy=343,  temperature=293.15)
Krypton       = Material(Z= 36,   A=83.798,   state='gas',   density=0.003749,  excitation_energy=352,  temperature=273.15, pressure=1)
Rubidium      = Material(Z= 37,   A=85.468,   state='solid', density=1.532,     excitation_energy=363,  temperature=293.15)
Strontium     = Material(Z= 38,   A=87.62,    state='solid', density=2.64,      excitation_energy=366,  temperature=293.15)
Yttrium       = Material(Z= 39,   A=88.906,   state='solid', density=4.469,     excitation_energy=379,  temperature=293.15)
Zirconium     = Material(Z= 40,   A=91.224,   state='solid', density=6.506,     excitation_energy=393,  temperature=293.15)
Niobium       = Material(Z= 41,   A=92.906,   state='solid', density=8.57,      excitation_energy=417,  temperature=293.15)
Molybdenum    = Material(Z= 42,   A=95.95,    state='solid', density=10.223,    excitation_energy=424,  temperature=293.15)
Technetium    = Material(Z= 43,   A=97,       state='solid', density=11.5,      excitation_energy=428,  temperature=293.15)   # Only unstable isotopes
Ruthenium     = Material(Z= 44,   A=101.07,   state='solid', density=12.37,     excitation_energy=441,  temperature=293.15)
Rhodium       = Material(Z= 45,   A=102.91,   state='solid', density=12.41,     excitation_energy=449,  temperature=293.15)
Palladium     = Material(Z= 46,   A=106.42,   state='solid', density=12.02,     excitation_energy=470,  temperature=293.15)
Silver        = Material(Z= 47,   A=107.87,   state='solid', density=10.49,     excitation_energy=470,  temperature=293.15)
Cadmium       = Material(Z= 48,   A=112.41,   state='solid', density=8.65,      excitation_energy=469,  temperature=293.15)
Indium        = Material(Z= 49,   A=114.82,   state='solid', density=7.31,      excitation_energy=488,  temperature=293.15)
Tin           = Material(Z= 50,   A=118.71,   state='solid', density=7.365,     excitation_energy=488,  temperature=293.15)
Antimony      = Material(Z= 51,   A=121.76,   state='solid', density=6.697,     excitation_energy=487,  temperature=293.15)
Tellurium     = Material(Z= 52,   A=127.60,   state='solid', density=6.24,      excitation_energy=485,  temperature=293.15)
Iodine        = Material(Z= 53,   A=126.90,   state='solid', density=4.933,     excitation_energy=491,  temperature=293.15)
Xenon         = Material(Z= 54,   A=131.29,   state='gas',   density=0.005894,  excitation_energy=482,  temperature=273.15, pressure=1)
Cesium        = Material(Z= 55,   A=132.91,   state='solid', density=1.93,      excitation_energy=488,  temperature=293.15)
Barium        = Material(Z= 56,   A=137.33,   state='solid', density=3.62,      excitation_energy=491,  temperature=293.15)
Lanthanum     = Material(Z= 57,   A=138.91,   state='solid', density=6.15,      excitation_energy=501,  temperature=293.15)
Cerium        = Material(Z= 58,   A=140.12,   state='solid', density=6.77,      excitation_energy=523,  temperature=293.15)
Praseodymium  = Material(Z= 59,   A=140.91,   state='solid', density=6.77,      excitation_energy=535,  temperature=293.15)
Neodymium     = Material(Z= 60,   A=144.24,   state='solid', density=7.01,      excitation_energy=546,  temperature=293.15)
Promethium    = Material(Z= 61,   A=145,      state='solid', density=7.26,      excitation_energy=560,  temperature=293.15)   # Only unstable isotopes
Samarium      = Material(Z= 62,   A=150.36,   state='solid', density=7.52,      excitation_energy=574,  temperature=293.15)
Europium      = Material(Z= 63,   A=151.96,   state='solid', density=5.24,      excitation_energy=580,  temperature=293.15)
Gadolinium    = Material(Z= 64,   A=157.25,   state='solid', density=7.90,      excitation_energy=591,  temperature=293.15)
Terbium       = Material(Z= 65,   A=158.93,   state='solid', density=8.23,      excitation_energy=614,  temperature=293.15)
Dysprosium    = Material(Z= 66,   A=162.50,   state='solid', density=8.55,      excitation_energy=628,  temperature=293.15)
Holmium       = Material(Z= 67,   A=164.93,   state='solid', density=8.80,      excitation_energy=650,  temperature=293.15)
Erbium        = Material(Z= 68,   A=167.26,   state='solid', density=9.07,      excitation_energy=658,  temperature=293.15)
Thulium       = Material(Z= 69,   A=168.93,   state='solid', density=9.32,      excitation_energy=674,  temperature=293.15)
Ytterbium     = Material(Z= 70,   A=173.05,   state='solid', density=6.97,      excitation_energy=684,  temperature=293.15)
Lutetium      = Material(Z= 71,   A=174.97,   state='solid', density=9.84,      excitation_energy=694,  temperature=293.15)
Hafnium       = Material(Z= 72,   A=178.49,   state='solid', density=13.31,     excitation_energy=705,  temperature=293.15)
Tantalum      = Material(Z= 73,   A=180.95,   state='solid', density=16.65,     excitation_energy=718,  temperature=293.15)
Tungsten      = Material(Z= 74,   A=183.84,   state='solid', density=19.254,    excitation_energy=727,  temperature=293.15)
Rhenium       = Material(Z= 75,   A=186.21,   state='solid', density=21.02,     excitation_energy=736,  temperature=293.15)
Osmium        = Material(Z= 76,   A=190.23,   state='solid', density=22.59,     excitation_energy=746,  temperature=293.15)
Iridium       = Material(Z= 77,   A=192.22,   state='solid', density=22.56,     excitation_energy=757,  temperature=293.15)
Platinum      = Material(Z= 78,   A=195.08,   state='solid', density=21.45,     excitation_energy=790,  temperature=293.15)
Gold          = Material(Z= 79,   A=196.97,   state='solid', density=19.32,     excitation_energy=790,  temperature=293.15)
Mercury       = Material(Z= 80,   A=200.59,   state='liquid',density=13.534,    excitation_energy=800,  temperature=293.15)
Thallium      = Material(Z= 81,   A=204.38,   state='solid', density=11.85,     excitation_energy=810,  temperature=293.15)
Lead          = Material(Z= 82,   A=207.19,   state='solid', density=11.348,    excitation_energy=823,  temperature=293.15)
Bismuth       = Material(Z= 83,   A=208.98,   state='solid', density=9.78,      excitation_energy=823,  temperature=293.15)
Polonium      = Material(Z= 84,   A=209,      state='solid', density=9.32,      excitation_energy=830,  temperature=293.15)   # Only unstable isotopes
Astatine      = Material(Z= 85,   A=210,      state='solid', density=7.0,       excitation_energy=825,  temperature=293.15)   # Only unstable isotopes
Radon         = Material(Z= 86,   A=222,      state='gas',   density=0.00973,   excitation_energy=794,  temperature=273.15, pressure=1)   # Only unstable isotopes
Francium      = Material(Z= 87,   A=223,      state='solid', density=2.48,      excitation_energy=827,  temperature=293.15)   # Only unstable isotopes
Radium        = Material(Z= 88,   A=226,      state='solid', density=5.50,      excitation_energy=826,  temperature=293.15)   # Only unstable isotopes
Actinium      = Material(Z= 89,   A=227,      state='solid', density=10.07,     excitation_energy=841,  temperature=293.15)   # Only unstable isotopes
Thorium       = Material(Z= 90,   A=232.04,   state='solid', density=11.72,     excitation_energy=847,  temperature=293.15)
Protactinium  = Material(Z= 91,   A=231.04,   state='solid', density=15.37,     excitation_energy=878,  temperature=293.15)
Uranium       = Material(Z= 92,   A=238.03,   state='solid', density=18.95,     excitation_energy=890,  temperature=293.15)
Neptunium     = Material(Z= 93,   A=237,      state='solid', density=20.45,     excitation_energy=902,  temperature=293.15)   # Only unstable isotopes
Plutonium     = Material(Z= 94,   A=244,      state='solid', density=19.84,     excitation_energy=921,  temperature=293.15)   # Only unstable isotopes
Americium     = Material(Z= 95,   A=243,      state='solid', density=13.69,     excitation_energy=934,  temperature=293.15)   # Only unstable isotopes
Curium        = Material(Z= 96,   A=247,      state='solid', density=13.51,     excitation_energy=939,  temperature=293.15)   # Only unstable isotopes
Berkelium     = Material(Z= 97,   A=247,      state='solid', density=14.78,     excitation_energy=952,  temperature=293.15)   # Only unstable isotopes
Californium   = Material(Z= 98,   A=251,      state='solid', density=15.1,      excitation_energy=966,  temperature=293.15)   # Only unstable isotopes
Einsteinium   = Material(Z= 99,   A=252,      state='solid', density=8.84,      temperature=293.15)   # Only unstable isotopes
Fermium       = Material(Z=100,   A=257,      state='solid', density=9.7,       excitation_energy=994,  temperature=293.15)   # Only unstable isotopes
Mendelevium   = Material(Z=101,   A=258,      state='solid', density=10.3,      temperature=293.15)   # Only unstable isotopes
Nobelium      = Material(Z=102,   A=259,      state='solid', density=9.9,       temperature=293.15)   # Only unstable isotopes
Lawrencium    = Material(Z=103,   A=266,      state='solid', density=15.6,      temperature=293.15)   # Only unstable isotopes
Rutherfordium = Material(Z=104,   A=267,      state='solid', density=23.2,      temperature=293.15)   # Only unstable isotopes
Dubnium       = Material(Z=105,   A=268,      state='solid', density=29.3,      temperature=293.15)   # Only unstable isotopes
Seaborgium    = Material(Z=106,   A=269,      state='solid', density=35.0,      temperature=293.15)   # Only unstable isotopes
Bohrium       = Material(Z=107,   A=270,      state='solid', density=37.1,      temperature=293.15)   # Only unstable isotopes
Hassium       = Material(Z=108,   A=277,      state='solid', density=41.0,      temperature=293.15)   # Only unstable isotopes
Meitnerium    = Material(Z=109,   A=278,      state='solid', density=37.4,      temperature=293.15)   # Only unstable isotopes
Darmstadtium  = Material(Z=110,   A=281,      state='solid', density=34.8,      temperature=293.15)   # Only unstable isotopes
Roentgenium   = Material(Z=111,   A=282,      state='solid', density=28.7,      temperature=293.15)   # Only unstable isotopes
Copernicium   = Material(Z=112,   A=285,      state='solid', density=14.0,      temperature=293.15)   # Only unstable isotopes
Nihonium      = Material(Z=113,   A=286,      state='solid', density=16.0,      temperature=293.15)   # Only unstable isotopes
Flerovium     = Material(Z=114,   A=289,      state='solid', density=14.0,      temperature=293.15)   # Only unstable isotopes
Moscovium     = Material(Z=115,   A=290,      state='solid', density=13.5,      temperature=293.15)   # Only unstable isotopes
Livermorium   = Material(Z=116,   A=293,      state='solid', density=12.9,      excitation_energy=1213, temperature=293.15)   # Only unstable isotopes
Tennessine    = Material(Z=117,   A=294,      state='solid', density=7.2,       temperature=293.15)   # Only unstable isotopes
Oganesson     = Material(Z=118,   A=294,      state='solid', density=5.0,       temperature=293.15)   # Only unstable isotopes


# Extra parameters for Everest
Beryllium.adapt(inplace=True,  nuclear_radius=0.22, nuclear_elastic_slope=74.7,
                               cross_section=[0.271, 0.192, 0, 0, 0, 0.0035e-2])
Aluminium.adapt(inplace=True,  nuclear_radius=0.302, nuclear_elastic_slope=120.3,
                               cross_section=[0.643, 0.418, 0, 0, 0, 0.0340e-2])
Silicon.adapt(inplace=True,    nuclear_radius=0.441, nuclear_elastic_slope=120.14,
                               cross_section=[0.664, 0.430, 0, 0, 0, 0.0390e-2])
Copper.adapt(inplace=True,     nuclear_radius=0.366, nuclear_elastic_slope=217.8,
                               cross_section=[1.253, 0.769, 0, 0, 0, 0.1530e-2], hcut=0.01)
Germanium.adapt(inplace=True,  nuclear_radius=0.605, nuclear_elastic_slope=226.35,
                               cross_section=[1.388, 0.844, 0, 0, 0, 0.1860e-2])
Molybdenum.adapt(inplace=True, nuclear_radius=0.481, nuclear_elastic_slope=273.9,
                               cross_section=[1.713, 1.023, 0, 0, 0, 0.2650e-2])
Tungsten.adapt(inplace=True,   nuclear_radius=0.520, nuclear_elastic_slope=440.3,
                               cross_section=[2.765, 1.591, 0, 0, 0, 0.7680e-2], hcut=0.01)
Lead.adapt(inplace=True,       nuclear_radius=0.542, nuclear_elastic_slope=455.3,
                               cross_section=[3.016, 1.724, 0, 0, 0, 0.9070e-2], hcut=0.01)


# Give name and store in database
for name, obj in list(globals().items()):  # Have to wrap in list to take a snapshot (avoid updating in-place)
    if isinstance(obj, Material):
        obj.name = name
        if name == 'Aluminium':
            assert obj.Z == 13
        else:
            assert get_element_full_name_from_Z(obj.Z).lower() == name.lower(), \
            f'Inconsistency between material name {name} and Z={obj.Z} ({get_element_full_name_from_Z(obj.Z)})'
        obj.short_name = get_element_name_from_Z(obj.Z)
        # These elements are pre-defined in FLUKA
        if obj.Z in [1, 2, 4, 5, 6, 7, 8, 11, 12, 14, 15, 16, 18, 20, 22, 23, 24, 25, 26, 27, 28, 29, 30, 40,
                     41, 42, 47, 50, 72, 73, 74, 79, 80, 82]:
            obj.fluka_name = name[:8].upper().ljust(8) # FLUKA material name max length is 8
        obj.geant4_name = obj.short_name   # Not lowercase because case-sensitive in Geant4
        db[obj.short_name] = obj
        db.geant4[f'G4_{obj.geant4_name}'] = obj  # Additional Geant4 aliases

db['Aluminum'] = Aluminium  # Allow American spelling

Aluminium.fluka_name = 'ALUMINUM'
Phosphorus.fluka_name = 'PHOSPHO'
