# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

import numpy as np
from numbers import Number


_NORMS = ["total", "coll_max", "max", "none", "raw"]


def plot_lossmap(lossmap: dict, *, norm="total", ax=None, xlim=None, ylim=None,
                       legend=True, grid=True, energy=False, show=True, zoom=None,
                       cold_regions=None, warm_regions=None, savefig=None):
    import matplotlib.pyplot as plt

    font = {'family': 'serif', 'size': 17}
    format_dict = {f"font.{prop}": font[prop] for prop in font}
    with plt.rc_context(format_dict):
        if zoom is not None:
            if not hasattr(zoom, "__iter__") or len(zoom) != 2 or isinstance(zoom, str):
                raise ValueError("Zoom must be a list/array of 2 elements [s_min, s_max].")
            if ax is None:
                _, ax = plt.subplots(2, 1, figsize=(16, 8))
            fig = ax[0].figure
            if not isinstance(ax, np.ndarray) or ax.ndim != 1 or ax.shape[0] != 2:
                raise ValueError("When zoom is used, ax must be a list/array of 2 axes.")
            _plot_lossmap_base(lossmap, norm=norm, ax=ax[0], xlim=xlim, ylim=ylim,
                            legend=legend, grid=grid, energy=energy,
                            cold_regions=cold_regions, warm_regions=warm_regions)
            _plot_lossmap_base(lossmap, norm=norm, ax=ax[1], xlim=zoom, ylim=ylim,
                            legend=legend, grid=grid, energy=energy,
                            cold_regions=cold_regions, warm_regions=warm_regions)
        else:
            if ax is None:
                _, ax = plt.subplots(figsize=(16, 4))
            fig = ax.figure
            _plot_lossmap_base(lossmap, norm=norm, ax=ax, xlim=xlim, ylim=ylim,
                            legend=legend, grid=grid, energy=energy,
                            cold_regions=cold_regions, warm_regions=warm_regions)
        plt.tight_layout()
        if savefig:
            plt.savefig(savefig, dpi=300, bbox_inches='tight')
        if show:
            plt.show()
        else:
            plt.close()
        return fig, ax


def _plot_lossmap_base(lossmap: dict, *, norm="total", ax=None, xlim=None, ylim=None,
                       legend=True, grid=True, energy=False, cold_regions=None,
                       warm_regions=None):
    import matplotlib.pyplot as plt
    if not isinstance(ax, plt.Axes):
        raise ValueError("ax must be a matplotlib Axes instance.")
    if norm in _NORMS:
        norm = norm.lower()
    elif not isinstance(norm, Number):
        raise ValueError(f"Norm must be one of {_NORMS} or a number, not {norm}.")

    coll_s = lossmap['collimator']['s']
    coll_val = lossmap['collimator']['e'] if energy else lossmap['collimator']['n']
    if lossmap['interpolation']:
        aper_s = lossmap['aperture']['s_bins']
        aper_val = lossmap['aperture']['e_bins'] if energy else lossmap['aperture']['n_bins']
        aper_length = lossmap['aperture']['length_bins']
    else:
        aper_s = lossmap['aperture']['s']
        aper_val = lossmap['aperture']['e'] if energy else lossmap['aperture']['n']
        aper_length = 1 #lossmap['aperture']['length']    TODO
    if len(coll_s) == 0 and len(aper_s) == 0:
        raise ValueError("Empty loss map.")

    cold_s = np.array([], dtype=aper_s.dtype)
    cold_val = np.array([], dtype=aper_val.dtype)
    cold_length = 1
    warm_s = np.array([], dtype=aper_s.dtype)
    warm_val = np.array([], dtype=aper_val.dtype)
    warm_length = 1
    if cold_regions is not None:
        if warm_regions is not None:
            raise ValueError("Provide only one of cold_regions or warm_regions.")
        starts = np.sort(cold_regions[:, 0])
        ends   = np.sort(cold_regions[:, 1])
        num_starts = np.searchsorted(starts, aper_s, side='right') # For each s, count how many intervals start before or at s
        num_ends   = np.searchsorted(ends, aper_s, side='left')    # Count how many intervals end before s
        mask = num_starts > num_ends
        cold_s = aper_s[mask]
        cold_val = aper_val[mask]
        cold_length = aper_length if len(aper_length)==1 else aper_length[mask]
        warm_s = aper_s[~mask]
        warm_val = aper_val[~mask]
        warm_length = aper_length if len(aper_length)==1 else aper_length[~mask]
        aper_s = np.array([], dtype=aper_s.dtype)
        aper_val = np.array([], dtype=aper_val.dtype)
        aper_length = 1
    elif warm_regions is not None:
        starts = np.sort(warm_regions[:, 0])
        ends   = np.sort(warm_regions[:, 1])
        num_starts = np.searchsorted(starts, aper_s, side='right') # For each s, count how many intervals start before or at s
        num_ends   = np.searchsorted(ends, aper_s, side='left')    # Count how many intervals end before s
        mask = num_starts > num_ends
        warm_s = aper_s[mask]
        warm_val = aper_val[mask]
        warm_length = aper_length if len(aper_length)==1 else aper_length[mask]
        cold_s = aper_s[~mask]
        cold_val = aper_val[~mask]
        cold_length = aper_length if len(aper_length)==1 else aper_length[~mask]
        aper_s = np.array([], dtype=aper_s.dtype)
        aper_val = np.array([], dtype=aper_val.dtype)
        aper_length = 1

    if norm != "raw":
        if norm == "total":
            scale = coll_val.sum() + cold_val.sum() + warm_val.sum() + aper_val.sum()
        elif norm == "coll_max":
            scale = coll_val.max()
        elif norm == "max":
            scale = np.concatenate([coll_val, aper_val, cold_val, warm_val]).max()
        elif norm == "none":
            scale = 1
        else:
            scale = norm
        cold_val = cold_val / scale / cold_length
        warm_val = warm_val / scale / warm_length
        aper_val = aper_val / scale / aper_length
        coll_val = coll_val / scale / lossmap['collimator']['length']

    L = lossmap['machine_length']
    if xlim is None:
        xlim = [-0.01*L, 1.01*L]
    if ylim is None:
        minimum = np.concatenate([coll_val, aper_val, cold_val, warm_val]).min()
        maximum = np.concatenate([coll_val, aper_val, cold_val, warm_val]).max()
        # ylim = [1.e-7, 1.e1]
        ylim = [minimum * 0.9 if minimum>0 else 1.e-9, maximum * 1.5]

    if grid:
        ax.grid(axis="y", which="major", zorder=0)

    bar_common_kwargs = dict(width = 0.8, lw = 1, bottom = 1.e-9)
    if len(coll_s) > 0:
        ax.bar(coll_s, coll_val, color="k", edgecolor="k", label="Collimator", zorder=9, **bar_common_kwargs)
    if len(aper_s) > 0:
        ax.bar(aper_s, aper_val, color="tab:orange", edgecolor="tab:orange", label="Aperture", zorder=10, **bar_common_kwargs)
    if len(warm_s) > 0:
        ax.bar(warm_s, warm_val, color="tab:red", edgecolor="tab:red", label="Warm",  zorder=11, **bar_common_kwargs)
    if len(cold_s) > 0:
        ax.bar(cold_s, cold_val, color="blue", edgecolor="blue", label="Cold",  zorder=12, **bar_common_kwargs)

    ax.set_yscale("log")
    ax.set_ylim(ylim)
    ax.set_xlim(xlim)

    ax.set_yticks([10**i for i in range(int(np.log10(ylim[0])), 1 + int(np.log10(ylim[1])))])
    ax.set_xlabel("s [m]")
    ax.set_ylabel(_label_from_norm(norm, energy))

    if legend:
        ax.legend(loc='best', fancybox=True, framealpha=0.8, fontsize="small")


def _label_from_norm(norm, energy):
    if norm == "raw":
        if energy:
            return "Deposited energy [eV]"
        else:
            return "Particles absorbed [-]"
    elif norm == "none":
        if energy:
            return "Deposited energy [eV/m]"
        else:
            return "Particles absorbed [1/m]"
    elif norm == "total":
        if energy:
            return "Norm. energy [1/m]"
        else:
            return "Norm. inefficiency [1/m]"
    elif norm == "max":
        if energy:
            return "Energy (max norm) [1/m]"
        else:
            return "Inefficiency (max norm) [1/m]"
    elif norm == "coll_max":
        if energy:
            return "Energy (coll max norm) [1/m]"
        else:
            return "Inefficiency (coll max norm) [1/m]"
    else:
        if energy:
            return f"Energy (norm. by {norm}) [1/m]"
        else:
            return f"Inefficiency (norm. by {norm}) [1/m]"


# def _plot_arrow_beam_direction(lossmap, axes, beam_direction, maximum):
#     if type(beam_direction) is not dict:
#         beam_direction = {}
#     # arrow style
#     if 'arrow_opts' in beam_direction.keys():
#         arrow_opts = beam_direction['arrow_opts']
#     else:
#         arrow_opts = dict(color           = 'black',
#                           arrowstyle      = 'simple, head_width=.75, head_length=.75',
#                           connectionstyle = 'arc3,rad=0')
#     annotate_kwargs = dict(xycoords   = 'data',
#                            textcoords = 'data',
#                            arrowprops = arrow_opts,
#                            size       = 15)
#     arrow_hpos = 3000/27000*lossmap['machine_length']
#     arrow_vpos = 2.e-1 * maximum

#     # draw arrow
#     if lossmap.beam.name == "BEAM1":
#         text_alignment = 'left'
#         annotate_kwargs['xy']     = (arrow_hpos+1500, arrow_vpos)
#         annotate_kwargs['xytext'] = (arrow_hpos,      arrow_vpos)
#         xtext = arrow_hpos

#     elif lossmap.beam.name == "BEAM2":
#         annotate_kwargs['xy']     = (arrow_hpos     , arrow_vpos)
#         annotate_kwargs['xytext'] = (arrow_hpos+1500, arrow_vpos)
#         text_alignment = 'right'
#         xtext = arrow_hpos + 1500

#     axes[0].text(s="Beam", x=xtext, y=1.5*arrow_vpos, horizontalalignment=text_alignment, fontsize=17)
#     axes[0].annotate("", **annotate_kwargs)
