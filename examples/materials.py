import xcoll as xc


# Examples of material definitions
# ================================

# Example compound material: Ethanol (C2H6O)
mat = xc.Material(components=['C', 'H', 'O'], n_atoms=[2, 6, 1], density=0.78945, name='Ethanol',
                  state='liquid', temperature=293.15)
# # or equivalently (CH3CH2OH):
# mat = xc.Material(components=['C', 'H', 'C', 'H', 'O', 'H'], n_atoms=[1, 3, 1, 2, 1, 1], density=0.78945,
#                   name='Ethanol', state='liquid', temperature=293.15)
print(mat)

mat = xc.Material(components=['H', 'C', 'O', 'Na', 'Mg', 'Al', 'Si', 'K', 'Ca', 'Fe'],
                mass_fractions=[0.01, 0.001, 0.529107, 0.016, 0.002, 0.033872, 0.337021, 0.013, 0.044, 0.014],
                name='Concrete', density=2.35, state='solid')
print(mat)
