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

def plot_material_scattering(materials, length, p0c, n_points=100000, savefig=None):
    colls = {}
    for mat in materials:
        colls[mat.name] = xc.EverestBlock(length=length, material=mat)
    part_init = xt.Particles(x=np.zeros(n_points), p0c=p0c)
    part_init._init_random_number_generator()

    parts = {}
    masks = {}
    for mat in materials:
        parts[mat.name] = part_init.copy()
        colls[mat.name].track(parts[mat.name])
        masks[mat.name] = parts[mat.name].state > 0

    fig, ax = plt.subplots(3, 1, figsize=(12, 5))
    title = []
    for mat, part in parts.items():
        mask = masks[mat]
        # Plot radial deflection
        r = np.sqrt(part.x[mask]**2 + part.y[mask]**2)
        r_sorted = np.sort(r)
        low = np.log10(r_sorted[int(0.005*len(r_sorted))])
        high = np.log10(r_sorted[int(0.985*len(r_sorted))])
        low -= 0.3*(high - low)
        ax[0].hist(r, bins=np.logspace(low, high, 1000), density=True, histtype='step', label=mat)
        ax[0].set_xlabel(r'$\Delta r$ [$\mu$m]')
        ax[0].set_ylabel('Density')
        ax[0].set_xscale('log')
        ax[0].legend()
        # Plot angular deflection
        phi = np.arctan2(part.y[mask], part.x[mask])
        ax[1].hist(phi, bins=1000, range=(-np.pi, np.pi), density=True, histtype='step', label=mat)
        ax[1].set_xlabel(r'$\Delta \phi$ [rad]')
        ax[1].set_ylabel('Density')
        ax[1].legend()
        # Plot energy loss
        delta_sorted = np.sort(np.array(-part.delta[mask]))
        low  = np.log10(delta_sorted[int(0.005*len(delta_sorted))])
        high = np.log10(delta_sorted[int(0.985*len(delta_sorted))])
        low -= 0.3*(high - low)
        ax[2].hist(-part.delta[mask], bins=np.logspace(low, high, 1000), density=True, histtype='step', label=mat)
        ax[2].set_xlabel(r'$-\delta$ [-]')
        ax[2].set_ylabel('Density')
        ax[2].set_xscale('log')
        ax[2].legend()
        title.append(f"{mat} ({mask.sum()/n_points:.2%} surv.)")
    fig.suptitle(f"{length}m {int(energy/1e9)}GeV " + ' vs '.join(title), fontsize=11)
    plt.tight_layout()
    if savefig is not None:
        plt.savefig(savefig, dpi=300)
    plt.show()


for length in [0.1, 0.25, 1.2]:
    for energy in [20e9, 450e9, 7e12]:
        plot_material_scattering([xc.materials.MolybdenumGraphite, xc.materials.MolybdenumGraphite6400, xc.materials.K2MolybdenumGraphite],
                                 length, energy, 10_000_000,
                                 f'scattering_mogr_{length}m_{int(energy/1e9)}GeV.png')
        plot_material_scattering([xc.materials.CopperDiamond, xc.materials.K2CopperDiamond],
                                 length, energy, 10_000_000,
                                 f'scattering_cucd_{length}m_{int(energy/1e9)}GeV.png')
        plot_material_scattering([xc.materials.Glidcop15, xc.materials.K2Glidcop15],
                                 length, energy, 10_000_000,
                                 f'scattering_glid_{length}m_{int(energy/1e9)}GeV.png')
        plot_material_scattering([xc.materials.Inermet180, xc.materials.K2Inermet180],
                                 length, energy, 10_000_000,
                                 f'scattering_iner_{length}m_{int(energy/1e9)}GeV.png')
