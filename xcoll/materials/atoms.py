# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

from xtrack.particles.pdg import get_element_name_from_Z, get_element_full_name_from_Z

from .material import Material
from .database import db


# Atomic mass from https://iupac.qmul.ac.uk/AtWt/ rounded to 5 digits
# Densities (in g/cm3) from:
#     https://en.wikipedia.org/wiki/List_of_chemical_elements
#     https://en.wikipedia.org/wiki/Densities_of_the_elements_(data_page)
Hydrogen      = Material(Z=1,     A=1.0080,   state='gas',     density=89.88e-6)
Helium        = Material(Z=2,     A=4.0026,   state='gas',     density=178.5e-6)
Lithium       = Material(Z=3,     A=6.94,     state='solid',   density=0.534)
Beryllium     = Material(Z=4,     A=9.0122,   state='solid',   density=1.848)
Boron         = Material(Z=5,     A=10.81,    state='solid',   density=2.34)
Carbon        = Material(Z=6,     A=12.011,   state='solid',   density=2.265)   # Amorphous
Nitrogen      = Material(Z=7,     A=14.007,   state='gas',     density=0.0012506)
Oxygen        = Material(Z=8,     A=15.999,   state='gas',     density=0.001429)
Fluorine      = Material(Z=9,     A=18.998,   state='gas',     density=0.001696)
Neon          = Material(Z=10,    A=20.180,   state='gas',     density=0.0009002)
Sodium        = Material(Z=11,    A=22.990,   state='solid',   density=0.968)
Magnesium     = Material(Z=12,    A=24.305,   state='solid',   density=1.738)
Aluminium     = Material(Z=13,    A=26.982,   state='solid',   density=2.70)
Silicon       = Material(Z=14,    A=28.085,   state='solid',   density=2.3296)
Phosphorus    = Material(Z=15,    A=30.974,   state='solid',   density=2.2)     # Red phosphorus
Sulfur        = Material(Z=16,    A=32.06,    state='solid',   density=2.067)
Chlorine      = Material(Z=17,    A=35.45,    state='gas',     density=0.003214)
Argon         = Material(Z=18,    A=39.948,   state='gas',     density=0.001784)
Potassium     = Material(Z=19,    A=39.098,   state='solid',   density=0.862)
Calcium       = Material(Z=20,    A=40.078,   state='solid',   density=1.54)
Scandium      = Material(Z=21,    A=44.956,   state='solid',   density=2.99)
Titanium      = Material(Z=22,    A=47.867,   state='solid',   density=4.506)
Vanadium      = Material(Z=23,    A=50.942,   state='solid',   density=6.11)
Chromium      = Material(Z=24,    A=51.996,   state='solid',   density=7.19)
Manganese     = Material(Z=25,    A=54.938,   state='solid',   density=7.26)
Iron          = Material(Z=26,    A=55.845,   state='solid',   density=7.874)
Cobalt        = Material(Z=27,    A=58.933,   state='solid',   density=8.90)
Nickel        = Material(Z=28,    A=58.693,   state='solid',   density=8.908)
Copper        = Material(Z=29,    A=63.546,   state='solid',   density=8.96)
Zinc          = Material(Z=30,    A=65.38,    state='solid',   density=7.134)
Gallium       = Material(Z=31,    A=69.723,   state='solid',   density=5.91)
Germanium     = Material(Z=32,    A=72.63,    state='solid',   density=5.323)
Arsenic       = Material(Z=33,    A=74.922,   state='solid',   density=5.727)
Selenium      = Material(Z=34,    A=78.971,   state='solid',   density=4.809)
Bromine       = Material(Z=35,    A=79.904,   state='liquid',  density=3.1028)
Krypton       = Material(Z=36,    A=83.798,   state='gas',     density=0.003749)
Rubidium      = Material(Z=37,    A=85.468,   state='solid',   density=1.532)
Strontium     = Material(Z=38,    A=87.62,    state='solid',   density=2.64)
Yttrium       = Material(Z=39,    A=88.906,   state='solid',   density=4.469)
Zirconium     = Material(Z=40,    A=91.224,   state='solid',   density=6.506)
Niobium       = Material(Z=41,    A=92.906,   state='solid',   density=8.57)
Molybdenum    = Material(Z=42,    A=95.95,    state='solid',   density=10.223)
Technetium    = Material(Z=43,    A=97,       state='solid',   density=11.5)    # Only unstable isotopes
Ruthenium     = Material(Z=44,    A=101.07,   state='solid',   density=12.37)
Rhodium       = Material(Z=45,    A=102.91,   state='solid',   density=12.41)
Palladium     = Material(Z=46,    A=106.42,   state='solid',   density=12.02)
Silver        = Material(Z=47,    A=107.87,   state='solid',   density=10.49)
Cadmium       = Material(Z=48,    A=112.41,   state='solid',   density=8.65)
Indium        = Material(Z=49,    A=114.82,   state='solid',   density=7.31)
Tin           = Material(Z=50,    A=118.71,   state='solid',   density=7.365)
Antimony      = Material(Z=51,    A=121.76,   state='solid',   density=6.697)
Tellurium     = Material(Z=52,    A=127.60,   state='solid',   density=6.24)
Iodine        = Material(Z=53,    A=126.90,   state='solid',   density=4.933)
Xenon         = Material(Z=54,    A=131.29,   state='gas',     density=0.005894)
Cesium        = Material(Z=55,    A=132.91,   state='solid',   density=1.93)
Barium        = Material(Z=56,    A=137.33,   state='solid',   density=3.62)
Lanthanum     = Material(Z=57,    A=138.91,   state='solid',   density=6.15)
Cerium        = Material(Z=58,    A=140.12,   state='solid',   density=6.77)
Praseodymium  = Material(Z=59,    A=140.91,   state='solid',   density=6.77)
Neodymium     = Material(Z=60,    A=144.24,   state='solid',   density=7.01)
Promethium    = Material(Z=61,    A=145,      state='solid',   density=7.26)    # Only unstable isotopes
Samarium      = Material(Z=62,    A=150.36,   state='solid',   density=7.52)
Europium      = Material(Z=63,    A=151.96,   state='solid',   density=5.24)
Gadolinium    = Material(Z=64,    A=157.25,   state='solid',   density=7.90)
Terbium       = Material(Z=65,    A=158.93,   state='solid',   density=8.23)
Dysprosium    = Material(Z=66,    A=162.50,   state='solid',   density=8.55)
Holmium       = Material(Z=67,    A=164.93,   state='solid',   density=8.80)
Erbium        = Material(Z=68,    A=167.26,   state='solid',   density=9.07)
Thulium       = Material(Z=69,    A=168.93,   state='solid',   density=9.32)
Ytterbium     = Material(Z=70,    A=173.05,   state='solid',   density=6.97)
Lutetium      = Material(Z=71,    A=174.97,   state='solid',   density=9.84)
Hafnium       = Material(Z=72,    A=178.49,   state='solid',   density=13.31)
Tantalum      = Material(Z=73,    A=180.95,   state='solid',   density=16.65)
Tungsten      = Material(Z=74,    A=183.84,   state='solid',   density=19.3)
Rhenium       = Material(Z=75,    A=186.21,   state='solid',   density=21.02)
Osmium        = Material(Z=76,    A=190.23,   state='solid',   density=22.59)
Iridium       = Material(Z=77,    A=192.22,   state='solid',   density=22.56)
Platinum      = Material(Z=78,    A=195.08,   state='solid',   density=21.45)
Gold          = Material(Z=79,    A=196.97,   state='solid',   density=19.32)
Mercury       = Material(Z=80,    A=200.59,   state='liquid',  density=13.534)
Thallium      = Material(Z=81,    A=204.38,   state='solid',   density=11.85)
Lead          = Material(Z=82,    A=207.19,   state='solid',   density=11.348)
Bismuth       = Material(Z=83,    A=208.98,   state='solid',   density=9.78)
Polonium      = Material(Z=84,    A=209,      state='solid',   density=9.32)    # Only unstable isotopes
Astatine      = Material(Z=85,    A=210,      state='solid',   density=7.0)     # Only unstable isotopes
Radon         = Material(Z=86,    A=222,      state='gas',     density=0.00973) # Only unstable isotopes
Francium      = Material(Z=87,    A=223,      state='solid',   density=2.48)    # Only unstable isotopes
Radium        = Material(Z=88,    A=226,      state='solid',   density=5.50)    # Only unstable isotopes
Actinium      = Material(Z=89,    A=227,      state='solid',   density=10.07)   # Only unstable isotopes
Thorium       = Material(Z=90,    A=232.04,   state='solid',   density=11.72)
Protactinium  = Material(Z=91,    A=231.04,   state='solid',   density=15.37)
Uranium       = Material(Z=92,    A=238.03,   state='solid',   density=18.95)
Neptunium     = Material(Z=93,    A=237,      state='solid',   density=20.45)   # Only unstable isotopes
Plutonium     = Material(Z=94,    A=244,      state='solid',   density=19.84)   # Only unstable isotopes
Americium     = Material(Z=95,    A=243,      state='solid',   density=13.69)   # Only unstable isotopes
Curium        = Material(Z=96,    A=247,      state='solid',   density=13.51)   # Only unstable isotopes
Berkelium     = Material(Z=97,    A=247,      state='solid',   density=14.78)   # Only unstable isotopes
Californium   = Material(Z=98,    A=251,      state='solid',   density=15.1)    # Only unstable isotopes
Einsteinium   = Material(Z=99,    A=252,      state='solid',   density=8.84)    # Only unstable isotopes
Fermium       = Material(Z=100,   A=257,      state='solid',   density=9.7)     # Only unstable isotopes
Mendelevium   = Material(Z=101,   A=258,      state='solid',   density=10.3)    # Only unstable isotopes
Nobelium      = Material(Z=102,   A=259,      state='solid',   density=9.9)     # Only unstable isotopes
Lawrencium    = Material(Z=103,   A=266,      state='solid',   density=15.6)    # Only unstable isotopes
Rutherfordium = Material(Z=104,   A=267,      state='solid',   density=23.2)    # Only unstable isotopes
Dubnium       = Material(Z=105,   A=268,      state='solid',   density=29.3)    # Only unstable isotopes
Seaborgium    = Material(Z=106,   A=269,      state='solid',   density=35.0)    # Only unstable isotopes
Bohrium       = Material(Z=107,   A=270,      state='solid',   density=37.1)    # Only unstable isotopes
Hassium       = Material(Z=108,   A=277,      state='solid',   density=41.0)    # Only unstable isotopes
Meitnerium    = Material(Z=109,   A=278,      state='solid',   density=37.4)    # Only unstable isotopes
Darmstadtium  = Material(Z=110,   A=281,      state='solid',   density=34.8)    # Only unstable isotopes
Roentgenium   = Material(Z=111,   A=282,      state='solid',   density=28.7)    # Only unstable isotopes
Copernicium   = Material(Z=112,   A=285,      state='solid',   density=14.0)    # Only unstable isotopes
Nihonium      = Material(Z=113,   A=286,      state='solid',   density=16.0)    # Only unstable isotopes
Flerovium     = Material(Z=114,   A=289,      state='solid',   density=14.0)    # Only unstable isotopes
Moscovium     = Material(Z=115,   A=290,      state='solid',   density=13.5)    # Only unstable isotopes
Livermorium   = Material(Z=116,   A=293,      state='solid',   density=12.9)    # Only unstable isotopes
Tennessine    = Material(Z=117,   A=294,      state='solid',   density=7.2)     # Only unstable isotopes
Oganesson     = Material(Z=118,   A=294,      state='solid',   density=5.0)     # Only unstable isotopes


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
        if obj.state == 'gas':
            obj.temperature = 273.15
            obj.pressure = 1
        elif obj.state == 'liquid' or obj.state == 'solid':
            obj.temperature = 293.15
        # These elements are pre-defined in FLUKA
        if obj.Z in [1, 2, 4, 5, 6, 7, 8, 11, 12, 14, 15, 16, 18, 20, 22, 23, 24, 25, 26, 27, 28, 29, 30, 40,
                     41, 42, 47, 50, 72, 73, 74, 79, 80, 82]:
            obj.fluka_name = name[:8].upper() # FLUKA material name max length is 8
        obj.geant4_name = obj.short_name   # Not lowercase because case-sensitive in Geant4
        db[obj.short_name] = obj
        db.geant4[f'G4_{obj.geant4_name}'] = obj  # Additional Geant4 aliases

db['Aluminum'] = Aluminium  # Allow American spelling

Aluminium.fluka_name = 'ALUMINUM'
Phosphorus.fluka_name = 'PHOSPHO'
