# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

import numpy as np
import matplotlib.pyplot as plt

import xtrack as xt
import xcoll as xc


# Comparing K2 material properties to new material properties
# ===========================================================

# Check where properties differ:
for mat in xc.materials.db:
    if mat.name is not None and mat.name.startswith('K2'):
        old_dct = mat._xobject._to_dict()
        if mat.name[2:] not in xc.materials.db:
            print(f"Material {mat.name}: no base material found for comparison")
            continue
        new_dct = xc.Material.from_dict(xc.materials.db[mat.name[2:]].to_dict())._xobject._to_dict()
        for key in ['_density', '_ZA_mean', '_Z2_eff', '_atoms_per_volume', '_num_nucleons_eff',
                    '_radiation_length', '_excitation_energy', '_nuclear_radius', '_nuclear_elastic_slope',
                    '_cross_section', '_hcut']:
            if key in old_dct and key in new_dct:
                old_val = old_dct[key]
                new_val = new_dct[key]
                if old_val != new_val:
                    mess = f"Material {mat.name}: field {key} differs: {old_val} != {new_val}"
                    print(f"{mess:112}(diff: {100*abs(old_val-new_val)/new_val}%)")

# Only compounds are clearly different, so let's investigate those:


# Plotting scattering differences between old and new material definitions
# ========================================================================

def plot_material_scattering(mat1, mat2, length, p0c, n_points=100000, savefig=None):
    coll1 = xc.EverestBlock(length=length, material=mat1)
    coll2 = xc.EverestBlock(length=length, material=mat2)
    part_init = xt.Particles(x=np.zeros(n_points), p0c=p0c)
    part_init._init_random_number_generator()

    part1 = part_init.copy()
    part2 = part_init.copy()
    coll1.track(part1)
    coll2.track(part2)

    mask1 = part1.state > 0
    mask2 = part2.state > 0

    fig, ax = plt.subplots(3, 1, figsize=(12, 5))
    r1 = np.sqrt(part1.x[mask1]**2 + part1.y[mask1]**2)
    r2 = np.sqrt(part2.x[mask2]**2 + part2.y[mask2]**2)
    ax[0].hist(r1, bins=np.logspace(-8, -3.4, 1000), density=True, histtype='step', label=mat1.name)
    ax[0].hist(r2, bins=np.logspace(-8, -3.4, 1000), density=True, histtype='step', label=mat2.name)
    ax[0].set_xlabel(r'$\Delta r$ [$\mu$m]')
    ax[0].set_ylabel('Density')
    ax[0].set_xscale('log')
    ax[0].legend()
    phi1 = np.arctan2(part1.y[mask1], part1.x[mask1])
    phi2 = np.arctan2(part2.y[mask2], part2.x[mask2])
    ax[1].hist(phi1, bins=1000, range=(-np.pi, np.pi), density=True, histtype='step', label=mat1.name)
    ax[1].hist(phi2, bins=1000, range=(-np.pi, np.pi), density=True, histtype='step', label=mat2.name)
    ax[1].set_xlabel(r'$\Delta \phi$ [rad]')
    ax[1].set_ylabel('Density')
    ax[1].legend()
    ax[2].hist(-part1.delta[mask1], bins=np.logspace(-4, -2, 1000), density=True, histtype='step', label=mat1.name)
    ax[2].hist(-part2.delta[mask2], bins=np.logspace(-4, -2, 1000), density=True, histtype='step', label=mat2.name)
    ax[2].set_xlabel(r'$-\delta$ [-]')
    ax[2].set_ylabel('Density')
    ax[2].set_xscale('log')
    ax[2].legend()
    fig.suptitle(f"{mat1.name} ({mask1.sum()/n_points:.2%} survived) vs {mat2.name} ({mask2.sum()/n_points:.2%} survived)",
                 fontsize=11)
    plt.tight_layout()
    if savefig is not None:
        plt.savefig(savefig, dpi=300)
    plt.show()


for length in [0.1, 0.25, 1.2]:
    for energy in [20e9, 450e9, 7e12]:
        plot_material_scattering(xc.materials.MolybdenumGraphite, xc.materials.K2MolybdenumGraphite,
                                 length, energy, 1_000_000,
                                 f'scattering_mogr_{length}m_{int(energy/1e9)}GeV.png')
        plot_material_scattering(xc.materials.CopperDiamond, xc.materials.K2CopperDiamond,
                                 length, energy, 1_000_000,
                                 f'scattering_cucd_{length}m_{int(energy/1e9)}GeV.png')
        plot_material_scattering(xc.materials.Glidcop15, xc.materials.K2Glidcop15,
                                 length, energy, 1_000_000,
                                 f'scattering_glid_{length}m_{int(energy/1e9)}GeV.png')
        plot_material_scattering(xc.materials.Inermet180, xc.materials.K2Inermet180,
                                 length, energy, 1_000_000,
                                 f'scattering_iner_{length}m_{int(energy/1e9)}GeV.png')
