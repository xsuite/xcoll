# copyright ############################### #
# This file is part of the Xcoll Package.  #
# Copyright (c) CERN, 2023.                 #
# ######################################### #

import numpy as np

# TODO: very much non-robust!!
def make_ion_from_properties(q, m):
    amu = 931.494061e6
    A = np.array(np.round(m/amu), dtype=np.int32)
    Z = np.array(np.round(q), dtype=np.int32)
    # A, Z, PDG
    return A, Z, np.array([ 1000000000 + z*10000 + a*10 for z, a in zip(Z,A)], dtype=np.int32)
    
# TODO: mass info ?
def get_particle_info_from_pdgid(pdgid):
    # q, A, Z, name
    if pdgid > 1000000000:
        # Ion
        tmpid = pdgid - 1000000000
        L = int(tmpid/1e7)
        tmpid -= L*1e7
        Z = int(tmpid /1e4)
        tmpid -= Z*1e4
        A = int(tmpid /10)
        tmpid -= A*10
        isomer_level = int(tmpid)
        return Z, A, Z, f'{get_element_name_from_Z(Z)}{A}'#, L, isomer_level, get_element_full_name_from_Z(Z)

    elif pdgid == 11:
        return -1, 0, 0, 'electron'
    elif pdgid == -11:
        return 1, 0, 0, 'positron'
    elif pdgid == 13:
        return -1, 0, 0, 'muon'
    elif pdgid == 15:
        return -1, 0, 0, 'tau'
    elif pdgid == 111:
        return 0, 0, 0, 'pi0'
    elif pdgid == 211:
        return 1, 0, 0, 'pi+'
    elif pdgid == -211:
        return -1, 0, 0, 'pi-'
    elif pdgid == 311:
        return 0, 0, 0, 'K0'
    elif pdgid == 321:
        return 1, 0, 0, 'K+'
    elif pdgid == -321:
        return -1, 0, 0, 'K-'
    elif pdgid == 2212:
        return 1, 1, 1, 'proton'
    elif pdgid == 2112:
        return 0, 1, 0, 'neutron'
    elif pdgid == 2224:
        return 2, 0, 0, 'Delta++'
    elif pdgid == 2214:
        return 1, 0, 0, 'Delta+'
    elif pdgid == 2114:
        return 0, 0, 0, 'Delta0'
    elif pdgid == 1114:
        return -1, 0, 0, 'Delta-'
    elif pdgid == 3122:
        return 0, 0, 0, 'Lambda'
    elif pdgid == 3222:
        return 1, 0, 0, 'Sigma+'
    elif pdgid == 3212:
        return 0, 0, 0, 'Sigma0'
    elif pdgid == 3112:
        return -1, 0, 0, 'Sigma-'
    elif pdgid == 3322:
        return 0, 0, 0, 'Xi'
    elif pdgid == 3312:
        return -1, 0, 0, 'Xi-'
    elif pdgid < 0:
        antipart = get_particle_info_from_pdgid(pdgid)
        return antipart[0], antipart[1], antipart[2], f'anti-{antipart[3]}'
    else:
        raise ValueError(f"PDG ID {pdgid} not recognised!")


def get_element_name_from_Z(z):
    ele = {1:   "H",  2:   "He", 3:   "Li", 4:   "Be", 5:   "B",  6:   "C",  7:   "N",  8:   "O",  9:   "F",  10:  "Ne",
           11:  "Na", 12:  "Mg", 13:  "Al", 14:  "Si", 15:  "P",  16:  "S",  17:  "Cl", 18:  "Ar", 19:  "K",  20:  "Ca",
           21:  "Sc", 22:  "Ti", 23:  "V",  24:  "Cr", 25:  "Mn", 26:  "Fe", 27:  "Co", 28:  "Ni", 29:  "Cu", 30:  "Zn",
           31:  "Ga", 32:  "Ge", 33:  "As", 34:  "Se", 35:  "Br", 36:  "Kr", 37:  "Rb", 38:  "Sr", 39:  "Y",  40:  "Zr",
           41:  "Nb", 42:  "Mo", 43:  "Tc", 44:  "Ru", 45:  "Rh", 46:  "Pd", 47:  "Ag", 48:  "Cd", 49:  "In", 50:  "Sn",
           51:  "Sb", 52:  "Te", 53:  "I",  54:  "Xe", 55:  "Cs", 56:  "Ba", 57:  "La", 58:  "Ce", 59:  "Pr", 60:  "Nd",
           61:  "Pm", 62:  "Sm", 63:  "Eu", 64:  "Gd", 65:  "Tb", 66:  "Dy", 67:  "Ho", 68:  "Er", 69:  "Tm", 70:  "Yb",
           71:  "Lu", 72:  "Hf", 73:  "Ta", 74:  "W",  75:  "Re", 76:  "Os", 77:  "Ir", 78:  "Pt", 79:  "Au", 80:  "Hg",
           81:  "Tl", 82:  "Pb", 83:  "Bi", 84:  "Po", 85:  "At", 86:  "Rn", 87:  "Fr", 88:  "Ra", 89:  "Ac", 90:  "Th",
           91:  "Pa", 92:  "U",  93:  "Np", 94:  "Pu", 95:  "Am", 96:  "Cm", 97:  "Bk", 98:  "Cf", 99:  "Es", 100: "Fm",
           101: "Md", 102: "No", 103: "Lr", 104: "Rf", 105: "Db", 106: "Sg", 107: "Bh", 108: "Hs", 109: "Mt", 110: "Ds",
           111: "Rg", 112: "Cn", 113: "Nh", 114: "Fl", 115: "Mc", 116: "Lv", 117: "Ts", 118: "Og"
          }
    if z not in ele:
        raise ValueError(f"Element with {z} protons not known.")
    return ele[z]

def get_element_full_name_from_Z(z):
    ele_long = {  1: "Hydrogen",        2: "Helium",          3: "Lithium",         4: "Beryllium",       5: "Boron",         
                  6: "Carbon",          7: "Nitrogen",        8: "Oxygen",          9: "Fluorine",       10: "Neon",          
                 11: "Sodium",         12: "Magnesium",      13: "Aluminum",       14: "Silicon",        15: "Phosphorus",    
                 16: "Sulfur",         17: "Chlorine",       18: "Argon",          19: "Potassium",      20: "Calcium",       
                 21: "Scandium",       22: "Titanium",       23: "Vanadium",       24: "Chromium",       25: "Manganese",     
                 26: "Iron",           27: "Cobalt",         28: "Nickel",         29: "Copper",         30: "Zinc",          
                 31: "Gallium",        32: "Germanium",      33: "Arsenic",        34: "Selenium",       35: "Bromine",       
                 36: "Krypton",        37: "Rubidium",       38: "Strontium",      39: "Yttrium",        40: "Zirconium",     
                 41: "Niobium",        42: "Molybdenum",     43: "Technetium",     44: "Ruthenium",      45: "Rhodium",       
                 46: "Palladium",      47: "Silver",         48: "Cadmium",        49: "Indium",         50: "Tin",           
                 51: "Antimony",       52: "Tellurium",      53: "Iodine",         54: "Xenon",          55: "Cesium",        
                 56: "Barium",         57: "Lanthanum",      58: "Cerium",         59: "Praseodymium",   60: "Neodymium",     
                 61: "Promethium",     62: "Samarium",       63: "Europium",       64: "Gadolinium",     65: "Terbium",       
                 66: "Dysprosium",     67: "Holmium",        68: "Erbium",         69: "Thulium",        70: "Ytterbium",     
                 71: "Lutetium",       72: "Hafnium",        73: "Tantalum",       74: "Tungsten",       75: "Rhenium",       
                 76: "Osmium",         77: "Iridium",        78: "Platinum",       79: "Gold",           80: "Mercury",       
                 81: "Thallium",       82: "Lead",           83: "Bismuth",        84: "Polonium",       85: "Astatine",      
                 86: "Radon",          87: "Francium",       88: "Radium",         89: "Actinium",       90: "Thorium",       
                 91: "Protactinium",   92: "Uranium",        93: "Neptunium",      94: "Plutonium",      95: "Americium",     
                 96: "Curium",         97: "Berkelium",      98: "Californium",    99: "Einsteinium",   100: "Fermium",       
                101: "Mendelevium",   102: "Nobelium",      103: "Lawrencium",    104: "Rutherfordium", 105: "Dubnium",       
                106: "Seaborgium",    107: "Bohrium",       108: "Hassium",       109: "Meitnerium",    110: "Darmstadtium",  
                111: "Roentgenium",   112: "Copernicium",   113: "Nihonium",      114: "Flerovium",     115: "Moscovium",     
                116: "Livermorium",   117: "Tennessine",    118: "Oganesson"
               }
    if z not in ele_long:
        raise ValueError(f"Element with {z} protons not known.")
    return ele_long[z]


