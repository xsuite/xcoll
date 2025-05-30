# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

import numpy as np
from numbers import Number
import matplotlib.pyplot as plt


_NORMS = ["total", "coll_max", "max", "none", "raw"]


def _plot_lossmap_base(lossmap: dict, *, norm="total", ax=None, xlim=None, ylim=None,
                       legend=True, grid=True, energy=False, show=True, savefig=None):
    if norm in _NORMS:
        norm = norm.lower()
    elif not isinstance(norm, Number):
        raise ValueError(f"Norm must be one of {_NORMS} or a number, not {norm}.")

    coll_s = lossmap['collimator']['s']
    coll_val = lossmap['collimator']['e'] if energy else lossmap['collimator']['n']
    if lossmap['interpolation'] is None:
        cold_s = lossmap['aperture']['s']
        cold_val = lossmap['aperture']['e'] if energy else lossmap['aperture']['n']
    else:
        cold_s = lossmap['aperture']['s_bins']
        cold_val = lossmap['aperture']['e_bins'] if energy else lossmap['aperture']['n_bins']
    warm_s = np.array([])  # Placeholder for warm losses, if needed
    warm_val = np.array([])  # Placeholder for warm losses, if needed
    if norm != "raw":
        if lossmap['interpolation']:
            cold_val = cold_val / lossmap['interpolation']
            # warm_val = warm_val / lossmap['interpolation']
        coll_val = coll_val / lossmap['collimator']['length']
        if norm == "total":
            scale = coll_val.sum() + cold_val.sum() + warm_val.sum()
        elif norm == "coll_max":
            scale = coll_val.max()
        elif norm == "max":
            scale = max(coll_val.max(), cold_val.max(), warm_val.max())
        elif norm == "none":
            scale = 1
        else:
            scale = norm
        coll_val = coll_val / scale
        cold_val = cold_val / scale
        warm_val = warm_val / scale
        if lossmap['interpolation'] is not None:
            cold_val = cold_val / lossmap['aperture']['length_bins']

    L = lossmap['machine_length']
    xlim = xlim if xlim else [-0.01*L, 1.01*L]
    ylim = ylim if ylim else [1.e-7, 1.e1]

    font = {'family': 'serif', 'size': 17}
    format_dict = {f"font.{prop}": font[prop] for prop in font}
    with plt.rc_context(format_dict):
        if ax is None:
            _, ax = plt.subplots(figsize=(16, 4))

        if grid:
            ax.grid(axis="y", which="major", zorder=0)

        bar_common_kwargs = dict(width = 0.8, lw = 1, bottom = 1.e-9)
        ax.bar(coll_s, coll_val, color="k", edgecolor="k", label="Collimator", zorder=10, **bar_common_kwargs)
        ax.bar(cold_s, cold_val, color="b", edgecolor="b", label="Cold",  zorder=11, **bar_common_kwargs)
        # ax.bar(warm_s, warm_val, color="r", edgecolor="r", label="Warm",  zorder=12, **bar_common_kwargs)

        ax.set_yscale("log")
        ax.set_ylim(ylim)
        ax.set_xlim(xlim)

        ax.set_yticks([10**i for i in range(int(np.log10(ylim[0])), 1 + int(np.log10(ylim[1])))])
        ax.set_xlabel("s [m]")
        ax.set_ylabel(_label_from_norm(norm, energy))

        if legend:
            ax.legend(loc='best', fancybox=True, framealpha=0.8, fontsize="small")

        plt.tight_layout()
        if savefig:
            plt.savefig(savefig, dpi=300, bbox_inches='tight')
        if show:
            plt.show()
        else:
            plt.close()
        return ax.figure, ax


def _label_from_norm(norm, energy):
    if norm == "none":
        if energy:
            return "Deposited energy [eV/m]"
        else:
            return "Particles absorbed [1/m]"
    elif norm == "raw":
        if energy:
            return "Deposited energy [eV]"
        else:
            return "Particles absorbed [-]"
    elif norm == "total":
        return "Norm. inefficiency [1/m]"
    elif norm == "max":
        return "Inefficiency (norm. by max) [-]"
    elif norm == "coll_max":
        return "Inefficiency (norm. by coll max) [-]"
    else:
        return f"Inefficiency (norm. by {norm}) [1/m]"
