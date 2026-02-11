# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

import numpy as np
from numbers import Number


_NORMS = ["total", "coll_max", "max", "none", "raw", "deposited_energy", "deposited_energy_per_length"]


def plot_lossmap(lossmap: dict, *,
                 ax=None,
                 show=True,
                 zoom=None,
                 savefig=None,
                 titles=None,
                 norm="total",
                 energy=False,
                 rescaled_beam_intensity=None,
                 cold_regions=None,
                 warm_regions=None,
                 **kwargs):
    import matplotlib.pyplot as plt

    # Resolve zoom and titles parameters
    zoom = _resolve_zoom(zoom)
    if titles is None:
        titles = []
    if not hasattr(titles, "__iter__") or isinstance(titles, str):
        titles = [titles]
    if zoom and titles:
        if len(titles) != len(zoom) + 1:
            raise ValueError("When zoom is used, titles must be a list of as many titles as zoom intervals plus one.")

    font = {'family': 'serif', 'size': 17}
    format_dict = {f"font.{prop}": font[prop] for prop in font}
    with plt.rc_context(format_dict):
        if ax is None:
            _, ax = plt.subplots(len(zoom) + 1, 1, figsize=(16, (len(zoom) + 1)*4))
        if not zoom:
            ax = np.array([ax])
        fig = ax[0].figure
        if not isinstance(ax, np.ndarray) or ax.ndim != 1 or ax.shape[0] != len(zoom) + 1:
            raise ValueError("When zoom is used, ax must be a list of as many axes as zoom intervals.")
        xlim = kwargs.pop("xlim", None)
        _plot_lossmap_base(lossmap, ax=ax[0], norm=norm, energy=energy, xlim=xlim,
                           rescaled_beam_intensity=rescaled_beam_intensity, cold_regions=cold_regions,
                           warm_regions=warm_regions, **kwargs)
        if titles:
            ax[0].set_title(titles[0])
        for i, zz in enumerate(zoom):
            _plot_lossmap_base(lossmap, ax=ax[i+1], norm=norm, energy=energy, xlim=zz,
                               rescaled_beam_intensity=rescaled_beam_intensity, cold_regions=cold_regions,
                               warm_regions=warm_regions, **kwargs)
            if titles:
                ax[i+1].set_title(titles[i+1])

        plt.tight_layout()
        if savefig:
            plt.savefig(savefig, dpi=300, bbox_inches='tight')
        if show:
            plt.show()
        else:
            plt.close()
    if zoom:
        return fig, ax
    else:
        return fig, ax[0]


def _plot_lossmap_base(lossmap: dict, *,
                       ax=None,
                       grid=True,
                       norm="total",
                       normalise_by_length=None,
                       energy=False,
                       rescaled_beam_intensity=None,
                       cold_regions=None,
                       warm_regions=None,
                       xlim=None,
                       ylim=None,
                       x_label=None,
                       y_label=None,
                       xticks=None,
                       yticks=None,
                       legend=True):
    import matplotlib.pyplot as plt
    if not isinstance(ax, plt.Axes):
        raise ValueError("ax must be a matplotlib Axes instance.")

    # Resolve normalisation
    if isinstance(norm, str) and norm.lower() in _NORMS:
        norm = norm.lower()
    elif not isinstance(norm, Number):
        raise ValueError(f"Norm must be one of {_NORMS} or a number, not {norm}.")
    if norm == "deposited_energy" or norm == "deposited_energy_per_length":
        energy = True
        if rescaled_beam_intensity is None:
            raise ValueError("When norm is 'deposited_energy', beam_intensity " \
                           + "or rescaled_beam_intensity must be provided.")
    if normalise_by_length is None:
        if norm == "raw" or norm == "deposited_energy":
            normalise_by_length = False
        else:
            normalise_by_length = True
    elif norm == 'raw' and normalise_by_length:
        normalise_by_length = False
        print("Warning: norm is 'raw', normalise_by_length set to False.")
    elif norm == 'deposited_energy' and normalise_by_length:
        normalise_by_length = False
        print("Warning: norm is 'deposited_energy', normalise_by_length set to False.")
    elif norm == 'deposited_energy_per_length' and not normalise_by_length:
        normalise_by_length = True
        print("Warning: norm is 'deposited_energy_per_length', normalise_by_length set to True.")

    # Create arrays for plotting
    coll_s = lossmap['collimator']['s']
    coll_val = lossmap['collimator']['e'] if energy else lossmap['collimator']['n']
    if lossmap['interpolation']:
        aper_s = lossmap['aperture']['s_bins']
        aper_val = lossmap['aperture']['e_bins'] if energy else lossmap['aperture']['n_bins']
        aper_length = lossmap['aperture']['length_bins']
    else:
        aper_s = lossmap['aperture']['s']
        aper_val = lossmap['aperture']['e'] if energy else lossmap['aperture']['n']
        aper_length = lossmap['aperture']['length']
    if len(coll_s) == 0 and len(aper_s) == 0:
        raise ValueError("Empty loss map.")

    L = lossmap['machine_length']

    cold_s = np.array([], dtype=aper_s.dtype)
    cold_val = np.array([], dtype=aper_val.dtype)
    cold_length = 1
    warm_s = np.array([], dtype=aper_s.dtype)
    warm_val = np.array([], dtype=aper_val.dtype)
    warm_length = 1
    if cold_regions is not None:
        if warm_regions is not None:
            raise ValueError("Provide only one of cold_regions or warm_regions.")
        cold_regions = _unwrap_and_check(cold_regions, L)
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
        warm_regions = _unwrap_and_check(warm_regions, L)
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

    if norm == "total":
        scale = coll_val.sum() + cold_val.sum() + warm_val.sum() + aper_val.sum()
    elif norm == "coll_max":
        scale = coll_val.max()
    elif norm == "max":
        scale = np.concatenate([coll_val, aper_val, cold_val, warm_val]).max()
    elif norm == "deposited_energy" or norm == "deposited_energy_per_length":
        scale = 1 / (rescaled_beam_intensity * 1.60218e-19)
    elif norm == "none" or norm == 'raw':
        scale = 1
    else:  # numeric
        scale = norm
    cold_val = cold_val / scale
    warm_val = warm_val / scale
    aper_val = aper_val / scale
    coll_val = coll_val / scale
    if normalise_by_length:
        cold_val = cold_val / cold_length
        warm_val = warm_val / warm_length
        aper_val = aper_val / aper_length
        coll_val = coll_val / lossmap['collimator']['length']

    if xlim is None:
        xlim = [-0.01*L, 1.01*L]
    if xlim[1] < xlim[0]:
        xlim[0] -= L
        cold_s = np.concatenate([cold_s - L, cold_s])
        warm_s = np.concatenate([warm_s - L, warm_s])
        aper_s = np.concatenate([aper_s - L, aper_s])
        coll_s = np.concatenate([coll_s - L, coll_s])
        cold_val = np.concatenate([cold_val, cold_val])
        warm_val = np.concatenate([warm_val, warm_val])
        aper_val = np.concatenate([aper_val, aper_val])
        coll_val = np.concatenate([coll_val, coll_val])
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

    if xticks is not None:
        ax.set_xticks(xticks)
    if yticks is None:
        yticks = [10**i for i in range(int(np.log10(ylim[0])), 1 + int(np.log10(ylim[1])))]
    ax.set_yticks(yticks)
    if x_label is None:
        x_label = "s [m]"
    ax.set_xlabel(x_label)
    if y_label is None:
        y_label = _label_from_norm(norm, energy)
    ax.set_ylabel(y_label)

    if legend:
        ax.legend(loc='upper right', fancybox=True, framealpha=0.8, fontsize="small")



def _valid_single_zoom(zz, allow_str=False):
    if isinstance(zz, str):
        return allow_str
    return hasattr(zz, "__iter__")  and len(zz) == 2 and not hasattr(zz[0], "__iter__") \
           and not hasattr(zz[1], "__iter__")

def _resolve_zoom(zoom, allow_str=False):
    if _valid_single_zoom(zoom, allow_str=allow_str):
        zoom = [zoom]
    if zoom is None:
        zoom = []
    if not hasattr(zoom, "__iter__") or isinstance(zoom, str) \
    or np.any([not _valid_single_zoom(zz, allow_str=allow_str) for zz in zoom]):
        raise ValueError("Zoom must be a list of 2 elements [s_min, s_max], or a list of such lists.")
    return zoom


def _label_from_norm(norm, energy):
    if norm == "raw":
        if energy:
            return "Deposited energy [eV]"
        else:
            return "Particles absorbed [-]"
    elif norm == "deposited_energy":
        return "Deposited energy [J]"
    elif norm == "deposited_energy_per_length":
        return "Deposited energy [J/m]"
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
            return "Energy (norm. by max) [1/m]"
        else:
            return "Inefficiency (norm. by max) [1/m]"
    elif norm == "coll_max":
        if energy:
            return "Energy (norm. by coll max) [1/m]"
        else:
            return "Inefficiency (norm. by coll max) [1/m]"
    else:
        if energy:
            return f"Deposited energy [eV/m]"
        else:
            return f"Inefficiency [1/m]"


def _unwrap_and_check(pairs, L):
    # 1) unwrap to linear intervals in [0, L]
    linear = []
    for a, b in pairs:
        if a <= b:
            linear.append([a, b])
        else:
            linear.append([a, L])
            linear.append([0, b])

    if not linear:
        return []

    # 2) sort + merge on [0, L]
    linear.sort(key=lambda x: (x[0], x[1]))
    merged = [linear[0]]

    for a, b in linear[1:]:
        last = merged[-1]
        overlaps = a <= last[1]
        if overlaps:
            last[1] = max(last[1], b)
        else:
            merged.append([a, b])

    return np.array(merged)


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
