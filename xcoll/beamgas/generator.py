# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2026.                 #
# ######################################### #

from functools import lru_cache
from warnings import warn

import numpy as np
from scipy import constants as sc

import xtrack as xt
from xtrack.particles import pdg

from ..beam_elements.beamgas import (BeamGasMarker, _resolve_process,
                                     _resolve_n_scattering_events)
from .analysis import (BeamGasResult, _cutoff_scan, _event_losses, _lost,
                       _sample_sums, _tracked)
from .cross_sections import (ELECTRON_MASS_EV, BremsstrahlungCalculator,
                             CoulombScatteringCalculator,
                             _energy_from_momentum, atomic_number_from_symbol)

PDG_ID_ELECTRON = 11
PDG_ID_POSITRON = -11

# Tolerance on the position of the last marker, at the end of the line [m]
_S_TOL = 1e-6

# Resolution of BeamGasResult.cutoff_scan
_CUTOFF_SCAN_BINS_PER_DECADE = 10
# Warn when more than this fraction of the loss rate comes from events
# within a factor two of the lower cut
_LOWER_CUT_WARN_FRACTION = 1e-3
# Warn when the Coulomb interaction rate above theta_max exceeds this
# fraction of the loss rate
_UPPER_CUT_WARN_FRACTION = 1e-2


def _lepton_charge(particle_ref):
    """
    Charge of the beam particles, checked to be electrons or positrons.

    Parameters
    ----------
    particle_ref : xtrack.Particles
        Reference particle of the line.

    Returns
    -------
    q0 : float
        ``-1`` for electrons, ``+1`` for positrons.
    """
    pdg_id = int(np.atleast_1d(particle_ref.pdg_id)[0])
    if pdg_id == 0:
        raise ValueError(
            "The reference particle has no PDG id set. The beam-gas models "
            "need to know whether the beam is made of electrons or "
            "positrons.")
    if pdg_id not in (PDG_ID_ELECTRON, PDG_ID_POSITRON):
        raise ValueError(
            "The beam-gas models implemented in Xcoll are only valid for "
            "electron and positron beams, but the reference particle is a "
            f"{pdg.get_name_from_pdg_id(pdg_id)} (PDG id {pdg_id}).")

    # q0 sets the sign of the McKinley-Feshbach interference term, so it
    # must agree with the PDG id rather than be trusted blindly
    q0 = -1.0 if pdg_id == PDG_ID_ELECTRON else 1.0
    if not np.isclose(float(particle_ref.q0), q0):
        raise ValueError(
            f"The reference particle is a {pdg.get_name_from_pdg_id(pdg_id)} "
            f"(PDG id {pdg_id}), which must have q0={q0:+.0f}, but it has "
            f"q0={float(particle_ref.q0):+g}.")
    return q0


@lru_cache(maxsize=None)
def _cross_section(process, Z, p0c, q0, window):
    """
    Cross-section model of a gas species, and its total cross section.

    The models depend only on these parameters, so they are built once and
    shared by all the markers.

    Parameters
    ----------
    process : {'brems', 'coulomb'}
        Beam-gas process.
    Z : int
        Atomic number of the species.
    p0c : float
        Reference momentum [eV].
    q0 : float
        Charge of the beam particles.
    window : tuple of float
        ``(brems_energy_cut,)`` [eV] for bremsstrahlung, or
        ``(theta_min, theta_max)`` [rad] for Coulomb scattering.

    Returns
    -------
    calculator : BremsstrahlungCalculator or CoulombScatteringCalculator
        Cross-section model, with its samplers.
    xsec : float
        Total cross section of the generated part of the process [m^2].
    """
    if process == 'brems':
        calculator = BremsstrahlungCalculator(Z, p0c, energy_cut=window[0])
    else:
        calculator = CoulombScatteringCalculator(
            Z, p0c, q0=q0, theta_lim=window)
    return calculator, calculator.compute_xsec()


def _stream(seed, name, replica):
    """
    Random-number generator of a marker and replica.

    The stream is derived from the configured ``seed`` with
    :class:`numpy.random.SeedSequence`, keyed by the name of the marker and
    by the replica index. The stream of a marker therefore does not depend on
    which other markers or replicas are generated, nor on the order of the
    generation: independent jobs reproduce the particles of a single one.

    Parameters
    ----------
    seed : int
        Seed configured on the marker.
    name : str
        Name of the marker.
    replica : int
        Index of the replica.

    Returns
    -------
    rng : numpy.random.Generator
        Random-number generator of the stream.
    """
    return np.random.default_rng(np.random.SeedSequence(
        entropy=int(seed), spawn_key=(int(replica), *name.encode())))


def _check_configured(name, marker, particle_ref):
    """
    Check that a marker is configured for the reference particle of the line.
    """
    if not marker.p0c > 0:
        raise ValueError(
            f"BeamGasMarker '{name}' is not configured: call "
            "`line.xcoll.beamgas.configure(...)` first.")
    p0c = particle_ref.p0c[0]
    if abs(marker.p0c - p0c) > 1e-10 * p0c:
        raise ValueError(
            f"BeamGasMarker '{name}' is configured for p0c = "
            f"{marker.p0c:.10g} eV, but the reference particle of the line "
            f"has p0c = {p0c:.10g} eV: configure the markers again.")


def _generate_marker(name, marker, *, particle_ref, q0, at_element, s, rng,
                     weight_scale=1.0):
    """
    Generate the beam-gas-scattered particles of a marker.

    Parameters
    ----------
    name : str
        Name of the marker.
    marker : BeamGasMarker
        Configured marker.
    particle_ref : xtrack.Particles
        Reference particle, defining the species of the particles.
    q0 : float
        Charge of the beam particles.
    at_element : int
        Index of the marker in the line.
    s : float
        Position of the marker in the line [m].
    rng : numpy.random.Generator
        Random-number generator.
    weight_scale : float, optional
        Factor applied to the weights of the particles.

    Returns
    -------
    particles : xtrack.Particles
        Particles that have undergone one beam-gas interaction, with a
        capacity of four times their number to leave room for secondaries.
    log : xtrack.Table
        Record of the interactions, one row per particle: ``name``, ``s``,
        ``particle_id``, ``gas``, ``theta``, ``photon_energy`` and
        ``weight``.
    """
    n = int(marker.n_scattering_events)
    if n <= 0:
        raise ValueError(
            f"BeamGasMarker '{name}': `n_scattering_events` must be positive "
            "to generate beam-gas-scattered particles.")

    process = marker.process
    if process == 'brems':
        window = (float(marker.brems_energy_cut),)
    else:
        window = (float(marker.coulomb_theta_min),
                  float(marker.coulomb_theta_max))
    densities = marker.atomic_densities
    species = list(densities)
    models = {kk: _cross_section(process, atomic_number_from_symbol(kk),
                                 float(marker.p0c), q0, window)
              for kk in species}
    contributions = np.array([densities[kk] * models[kk][1]
                              for kk in species])
    if contributions.sum() <= 0:
        raise ValueError(
            f"No beam-gas interactions possible at marker '{name}': all the "
            "gas densities are zero.")

    x, px, y, py, zeta, delta = marker._generate_primaries(n, rng)

    # Assign a gas species to each event, proportionally to its
    # contribution to the local interaction rate
    i_species = rng.choice(len(species), size=n,
                           p=contributions / contributions.sum())

    weight = np.empty(n)
    theta = np.empty(n)
    photon_energy = np.empty(n)
    # NOTE: an object array (and not a fixed-width unicode array) is used
    # on purpose: `xtrack.view.View`, which wraps the elements accessed
    # through `line[name]`, corrupts numpy unicode arrays on conversion.
    gas = np.empty(n, dtype=object)

    for ii, kk in enumerate(species):
        mask = i_species == ii
        if not mask.any():
            continue
        sample = models[kk][0].sample_deflections(
            px[mask], py[mask], delta[mask], rng)
        px[mask] = sample.px
        py[mask] = sample.py
        delta[mask] = sample.delta
        weight[mask] = sample.weight
        theta[mask] = sample.theta
        photon_energy[mask] = sample.photon_energy
        gas[mask] = kk

    # Normalise the weights to the interaction rate of the section
    weight *= marker.interaction_rate / n * weight_scale

    # Leave room for the secondaries that collimators may produce
    particles = xt.Particles(
        _capacity=4*n,
        p0c=marker.p0c,
        mass0=particle_ref.mass0,
        q0=particle_ref.q0,
        pdg_id=particle_ref.pdg_id,
        x=x, px=px, y=y, py=py, zeta=zeta, delta=delta,
        weight=weight,
        s=s,
    )
    particles.at_element = at_element

    log = xt.Table({
        'name': np.full(n, name, dtype=object),
        's': np.full(n, float(s)),
        'particle_id': particles.particle_id[:n].copy(),
        'gas': gas,
        'theta': theta,
        'photon_energy': photon_energy,
        'weight': weight,
    })
    return particles, log


def _log_table(name, replica, log):
    """
    Columns of the record of the interactions of a sample.
    """
    particle_id = np.asarray(log['particle_id'], dtype=np.int64)
    n = len(particle_id)
    return {
        'name': np.full(n, name, dtype=object),
        'replica': np.full(n, replica),
        's': np.asarray(log['s'], dtype=float),
        'particle_id': particle_id,
        'gas': np.asarray(log['gas'], dtype=object),
        'theta': np.asarray(log['theta'], dtype=float),
        'photon_energy': np.asarray(log['photon_energy'], dtype=float),
        'weight': np.asarray(log['weight'], dtype=float),
    }


def _warn_on_truncation(process, window, x, contribution, cutoff_scan,
                        rate_above_theta_max):
    """
    Warn when the generation window visibly truncates the loss rate.

    Parameters
    ----------
    process : {'brems', 'coulomb'}
        Beam-gas process.
    window : tuple of float
        Range of the generation variable.
    x : ndarray
        Generation variable of each event.
    contribution : ndarray
        Loss-rate contribution of each event, i.e. its weight if lost and
        zero otherwise [1/s].
    cutoff_scan : xtrack.Table
        Cutoff scan of the events.
    rate_above_theta_max : float or None
        Coulomb interaction rate above ``theta_max`` [1/s].

    Returns
    -------
    None
    """
    rate = float(np.sum(contribution))
    if rate <= 0:
        return

    lo = window[0]
    if process == 'coulomb':
        name = '`coulomb_theta[0]`'
        what = 'scattering angles'
        hint = (" The log-uniform angular sampling makes a lower cut "
                "cheap: the statistics only degrade as "
                "log(theta_max/theta_min).")
    else:
        name = '`brems_energy_cut`'
        what = 'photon energies'
        hint = ""

    near_lower = float(np.sum(contribution[x < 2*lo])) / rate
    if near_lower > _LOWER_CUT_WARN_FRACTION:
        warn(f"{100*near_lower:.2g}% of the tracked loss rate comes from "
             f"events generated within a factor 2 of the lower cut "
             f"{name}={lo:.3g}. The cut lies inside the range of {what} that "
             f"cause losses, so the losses from below it are missing and the "
             f"lifetime is overestimated. Lower {name}, and check that "
             f"`BeamGasResult.cutoff_scan.lifetime_tracking` is flat in its "
             f"first rows.{hint}", stacklevel=3)

    if rate_above_theta_max is not None:
        above_upper = rate_above_theta_max / rate
        if above_upper > _UPPER_CUT_WARN_FRACTION:
            warn(f"The Coulomb interaction rate above `coulomb_theta[1]`="
                 f"{window[1]:.3g} rad, which is not generated, amounts to "
                 f"{100*above_upper:.2g}% of the tracked loss rate. Up to that "
                 f"much loss rate is missing, and all of it if these events "
                 f"are always lost (loss probability in the last "
                 f"`cutoff_scan` row: {cutoff_scan.loss_probability[-1]:.2f}). "
                 f"Raise `coulomb_theta[1]`, or add `rate_above_theta_max` to "
                 f"the loss rate.", stacklevel=3)


def _resolve_seed(seed):
    """
    Seed of the beam-gas random streams: drawn if ``None``.
    """
    if seed is None:
        return int(np.random.randint(1, 2**31 - 1))
    if isinstance(seed, (bool, np.bool_)) or not isinstance(
            seed, (int, np.integer)):
        raise TypeError(
            f"`seed` must be an integer (got {type(seed).__name__}).")
    if seed < 0:
        raise ValueError("`seed` must be non-negative.")
    return int(seed)


def _check_gas_density(gas_density):
    """
    Validate a gas-density table.

    Parameters
    ----------
    gas_density : xtrack.Table
        Table with a column ``s`` and one column per gas species.

    Returns
    -------
    species : list of str
        Names of the gas-species columns.
    """
    if not isinstance(gas_density, xt.Table):
        raise TypeError("`gas_density` must be an `xt.Table` object.")

    col_names = list(gas_density._col_names)
    if 's' not in col_names:
        raise ValueError("`gas_density` must contain an `s` column.")

    species = [cc for cc in col_names if cc not in ('name', 's')]
    if len(species) == 0:
        raise ValueError(
            "`gas_density` must contain at least one gas-species column, "
            "named after the chemical symbol of the species "
            "(e.g. 'H', 'N', 'Ar').")

    s = np.asarray(gas_density['s'], dtype=float)
    if np.any(np.diff(s) < 0):
        raise ValueError("`gas_density.s` must be non-decreasing.")

    for cc in ['s'] + species:
        try:
            vals = np.asarray(gas_density[cc], dtype=float)
        except (TypeError, ValueError):
            raise TypeError(
                f"`{cc}` column must be numeric (cannot coerce to float).")
        if np.isnan(vals).any():
            bad = list(np.where(np.isnan(vals))[0][:5])
            raise ValueError(f"`{cc}` contains NaN at indices {bad}.")
        if np.isinf(vals).any():
            bad = list(np.where(np.isinf(vals))[0][:5])
            raise ValueError(f"`{cc}` contains inf at indices {bad}.")
        if cc != 's' and (vals < 0).any():
            bad = list(np.where(vals < 0)[0][:5])
            raise ValueError(
                f"`{cc}` contains negative densities at indices {bad}.")

    for cc in species:
        # Raises a clear error for an unknown chemical symbol
        atomic_number_from_symbol(cc)

    return species


def _integrated_densities(gas_density, species, s_start, s_end):
    """
    Integrate the atomic density of each species over a section.

    The tabulated density profile is treated as piecewise linear in ``s``,
    so the trapezoidal rule evaluated on the union of the tabulated positions
    inside the section and the two section edges is exact.

    Parameters
    ----------
    gas_density : xtrack.Table
        Gas-density table.
    species : list of str
        Gas-species columns.
    s_start, s_end : float
        Limits of the section [m].

    Returns
    -------
    integrals : dict
        Mapping ``{species: integral}`` [atoms/m^2].
    """
    if s_end <= s_start:
        return {kk: 0.0 for kk in species}

    s_table = np.asarray(gas_density['s'], dtype=float)
    s_inner = s_table[(s_table > s_start) & (s_table < s_end)]
    s_nodes = np.concatenate(([s_start], s_inner, [s_end]))

    return {kk: float(np.trapezoid(
                np.interp(s_nodes, s_table,
                          np.asarray(gas_density[kk], dtype=float)),
                s_nodes))
            for kk in species}


def _is_sequence(value):
    return isinstance(value, (list, tuple, np.ndarray))


def _per_marker(value, n_markers, name):
    """
    One value per marker: a sequence is checked, a scalar is repeated.
    """
    if _is_sequence(value):
        value = list(value)
        if len(value) != n_markers:
            raise ValueError(
                f"`{name}` has {len(value)} entries, but {n_markers} "
                "markers are installed.")
        return value
    return [value] * n_markers


def _is_free(env, name):
    """
    Whether `name` can be given to a new element of the environment.
    """
    if name in env.elements:
        return False
    try:
        env._check_name_clashes(name)
    except ValueError:
        return False
    return True


def _free_names(env, prefix, n_names):
    """
    The `n_names` names ``f'{prefix}.{i}'`` with the lowest free indices.
    """
    names = []
    ii = 0
    while len(names) < n_names:
        name = f'{prefix}.{ii}'
        if _is_free(env, name):
            names.append(name)
        ii += 1
    return names


class BeamGasGenerator:
    """
    Beam-gas scattering simulation of a line, based on the
    :class:`xcoll.BeamGasMarker` elements installed in it.

    An instance bound to a line is available as ``line.xcoll.beamgas``. The
    generator holds no state other than the line: the configuration is stored
    in the markers, so that a line saved to JSON after it has been configured
    generates the same particles once loaded again.

    Parameters
    ----------
    line : xtrack.Line
        Line containing, or receiving, the beam-gas markers.
    """

    def __init__(self, line):
        if not isinstance(line, xt.Line):
            raise TypeError(
                f"`line` must be an xtrack.Line (got {type(line).__name__}).")
        self._line = line

    @property
    def line(self):
        """
        Line of the beam-gas simulation.

        Returns
        -------
        line : xtrack.Line
            Line bound to this generator.
        """
        return self._line

    @property
    def marker_names(self):
        """
        Names of the beam-gas markers of the line.

        Returns
        -------
        names : list of str
            Names of the :class:`xcoll.BeamGasMarker` elements, in line order.
        """
        element_dict = self.line.element_dict
        names = [nn for nn in self.line.element_names
                 if isinstance(element_dict[nn], BeamGasMarker)]
        if len(set(names)) != len(names):
            repeated = sorted({nn for nn in names if names.count(nn) > 1})
            raise ValueError(
                f"Beam-gas markers {repeated} appear more than once in the "
                "line: each marker must be a distinct element.")
        return names

    def install_markers(self, at, *, from_=None, anchor=None,
                        from_anchor=None, names=None, prefix='beamgas',
                        close_ring=True):
        """
        Install beam-gas markers in the line.

        Each marker is a new :class:`xcoll.BeamGasMarker` element, placed with
        ``env.place`` and inserted with ``line.insert``: ``at``, ``from_``,
        ``anchor`` and ``from_anchor`` have the same meaning as there.
        Installing the markers by hand in the same way is equivalent.

        Each marker represents the section of the line from the previous
        marker (from the start of the line for the first one) to itself, so
        markers are best placed densely where the gas density or the optics
        change.

        Parameters
        ----------
        at : float, str, xdeps reference, or sequence of them
            Position of each marker, evaluated at installation: an ``s``
            position [m], an expression of the environment variables, the
            name of an element (to place the marker at its centre), or a
            reference ``'name@anchor'`` with ``anchor`` one of ``'start'``,
            ``'center'`` and ``'end'``. With ``from_``, the offset from the
            reference element.
        from_ : str or sequence of str, optional
            Element from which ``at`` is measured. ``'name@anchor'`` selects
            the anchor of the element.
        anchor : {'start', 'center', 'end'} or sequence of them, optional
            Anchor of the marker placed at ``at``. Markers are thin, so this
            only matters to follow the conventions of ``env.place``.
        from_anchor : {'start', 'center', 'end'} or sequence of them, optional
            Anchor of ``from_`` from which ``at`` is measured.
        names : str or sequence of str, optional
            Names of the new markers. By default ``f'{prefix}.{i}'``, with the
            lowest indices whose names are free in the environment.
        prefix : str, optional
            Prefix of the default names.
        close_ring : bool, optional
            If ``True`` (default), a marker is added at the end of the line
            when, after the insertion, the last marker is not there: otherwise
            the part of the line after it would belong to no section and
            :meth:`configure` refuses the line.

        Returns
        -------
        names : list of str
            Names of the installed markers, the closing one included.
        """
        line = self.line
        env = line.env

        at = list(at) if _is_sequence(at) else [at]
        n_markers = len(at)
        if n_markers == 0:
            raise ValueError("`at` must give at least one position.")
        from_ = _per_marker(from_, n_markers, 'from_')
        anchor = _per_marker(anchor, n_markers, 'anchor')
        from_anchor = _per_marker(from_anchor, n_markers, 'from_anchor')

        if names is None:
            names = _free_names(env, prefix, n_markers)
        else:
            names = _per_marker(
                [names] if isinstance(names, str) else names,
                n_markers, 'names')
            if len(set(names)) != len(names):
                raise ValueError("`names` contains repeated names.")
            for nn in names:
                if not _is_free(env, nn):
                    raise ValueError(
                        f"The name '{nn}' is already used in the environment.")

        places = []
        for nn, aa, ff, an, fa in zip(names, at, from_, anchor, from_anchor):
            env.elements[nn] = BeamGasMarker()
            places.append(env.place(
                nn, at=aa, from_=ff, anchor=an, from_anchor=fa))
        line.insert(places)

        installed = list(names)
        if close_ring and not self._last_marker_at_end():
            end_name, = _free_names(env, prefix, 1)
            env.elements[end_name] = BeamGasMarker()
            line.insert(end_name, at=line.get_length())
            installed.append(end_name)

        return installed

    def configure(self, *, gas_density=None, process='brems', sigma_z=None,
                  sigma_delta=None, n_scattering_events=None,
                  nemitt_x=None, nemitt_y=None, gemitt_x=None, gemitt_y=None,
                  bunch_intensity=1.0, brems_energy_cut=10e3,
                  coulomb_theta=(1e-7, 50e-3), seed=None, twiss=None,
                  method='6d', at=None, verbose=True, **kwargs):
        """
        Configure the beam-gas markers of the line.

        Integrates the atomic density of each gas species over the section of
        each marker, and stores on the marker the resulting interaction rate,
        the local gas composition, optics and closed orbit, the beam parameters
        and the Monte Carlo settings. The section of a marker goes from the
        previous marker to it, and the section of the first one from
        ``s = 0``: the last marker must sit at the end of the line, which
        :meth:`install_markers` ensures.

        The whole configuration is stored on the markers: a line saved to JSON
        after this call generates the same particles once loaded again.

        Parameters
        ----------
        gas_density : xtrack.Table
            Residual-gas density profile. It must contain a column ``s`` with
            the longitudinal positions [m] and one column per gas species,
            named after the chemical symbol of the species (``'H'``, ``'N'``,
            ``'Ar'``, ...) and holding the *atomic* density [atoms/m^3]. For a
            molecular gas, use the atomic density of each constituent, e.g.
            ``N = 2 * n_molecules`` for N2.
        process : {'brems', 'coulomb'}, optional
            Beam-gas process to simulate: bremsstrahlung on the gas nuclei, or
            elastic Coulomb scattering off them. Default ``'brems'``.
        sigma_z : float
            RMS bunch length [m].
        sigma_delta : float
            RMS relative momentum spread.
        n_scattering_events : int
            Number of Monte Carlo interactions generated at each marker.
        nemitt_x, nemitt_y : float, optional
            Normalized emittances. Mutually exclusive with ``gemitt_x`` and
            ``gemitt_y``.
        gemitt_x, gemitt_y : float, optional
            Geometric emittances. Mutually exclusive with ``nemitt_x`` and
            ``nemitt_y``.
        bunch_intensity : float, optional
            Number of particles in the bunch. The interaction and loss rates
            are proportional to it; the lifetime is not. Default 1, which
            makes the rates per stored particle.
        brems_energy_cut : float, optional
            Lower limit of the generated bremsstrahlung photon energy [eV].
            Only used when ``process='brems'``. Default 10 keV.
        coulomb_theta : tuple of float, optional
            Range of the generated Coulomb scattering angle [rad]. Only used
            when ``process='coulomb'``. Losses from angles outside this range
            are missing from the result, so it must extend well below the
            smallest angle that can cause a loss and well above the angle
            beyond which every particle is lost. Default ``(1e-7, 50e-3)``.
        seed : int or None, optional
            Seed of the random-number streams of the markers. If ``None``, a
            seed is drawn now and stored, so that the generation stays
            reproducible.
        twiss : xtrack.TwissTable or None, optional
            Twiss table of the line. If ``None``, it is computed with
            ``line.twiss(method=method, **kwargs)``.
        method : {'6d', '4d'}, optional
            Twiss method, used if ``twiss`` is ``None``.
        at : str, sequence of str, or None, optional
            Markers to configure. If ``None`` (default), all the markers. The
            sections are always defined by all the markers of the line.
        verbose : bool, optional
            Whether to print the progress of the configuration.
        **kwargs
            Keyword arguments forwarded to ``line.twiss`` if ``twiss`` is
            ``None``.

        Returns
        -------
        None
        """
        line = self.line
        if getattr(line, 'particle_ref', None) is None:
            raise ValueError("`line` must have a `particle_ref`.")
        for name, value in (('gas_density', gas_density),
                            ('sigma_z', sigma_z),
                            ('sigma_delta', sigma_delta),
                            ('n_scattering_events', n_scattering_events)):
            if value is None:
                raise ValueError(f"`{name}` is required.")
        particle_ref = line.particle_ref
        q0 = _lepton_charge(particle_ref)
        process = _resolve_process(process)
        species = _check_gas_density(gas_density)
        n_scattering_events = _resolve_n_scattering_events(
            n_scattering_events=n_scattering_events, default=None)
        bunch_intensity = float(bunch_intensity)
        if bunch_intensity < 0:
            raise ValueError("`bunch_intensity` must be non-negative.")
        seed = _resolve_seed(seed)

        nemitt_given = nemitt_x is not None and nemitt_y is not None
        gemitt_given = gemitt_x is not None and gemitt_y is not None
        if nemitt_given and gemitt_given:
            raise ValueError(
                "Provide either normalised emittances (nemitt_x, nemitt_y) "
                "OR geometric emittances (gemitt_x, gemitt_y), not both.")
        if not (nemitt_given or gemitt_given):
            raise ValueError(
                "You must provide either both normalised emittances "
                "(nemitt_x, nemitt_y) OR both geometric emittances "
                "(gemitt_x, gemitt_y).")
        if nemitt_given:
            beta_gamma = particle_ref.beta0[0] * particle_ref.gamma0[0]
            gemitt_x = float(nemitt_x / beta_gamma)
            gemitt_y = float(nemitt_y / beta_gamma)

        # Check the markers before the (long) twiss
        markers = self._markers()
        self._check_sections_cover_line()
        selected = markers if at is None else self._markers(at)

        if twiss is None:
            twiss = line.twiss(method=method, **kwargs)

        p0c = float(particle_ref.p0c[0])
        window = ((float(brems_energy_cut),) if process == 'brems'
                  else (float(coulomb_theta[0]), float(coulomb_theta[1])))
        atomic_numbers = {kk: atomic_number_from_symbol(kk) for kk in species}
        xsecs = {kk: _cross_section(process, Z, p0c, q0, window)[1]
                 for kk, Z in atomic_numbers.items()}
        if process == 'coulomb' and window[1] < np.pi:
            xsecs_above = {
                kk: _cross_section(process, Z, p0c, q0,
                                   (window[1], np.pi))[1]
                for kk, Z in atomic_numbers.items()}
        else:
            xsecs_above = {kk: 0.0 for kk in species}

        f_rev = particle_ref.beta0[0] * sc.c / line.get_length()
        s_table = np.asarray(gas_density['s'], dtype=float)
        tab = line.get_table()
        index = {nn: ii for ii, nn in enumerate(markers)}

        for nn in selected:
            s = float(tab['s', nn])
            ii = index[nn]
            s_start = float(tab['s', markers[ii - 1]]) if ii > 0 else 0.0

            integrated = _integrated_densities(
                gas_density, species, s_start, s)
            atomic_densities = {
                kk: float(np.interp(s, s_table,
                                    np.asarray(gas_density[kk], dtype=float)))
                for kk in species}
            inverse_mfp = sum(atomic_densities[kk] * xsecs[kk]
                              for kk in species)

            if verbose:
                print(f'Configuring BeamGasMarker {nn}')

            line.element_dict[nn]._configure(
                p0c=p0c,
                bunch_intensity=bunch_intensity,
                ds=s - s_start,
                alfx=twiss['alfx', nn], betx=twiss['betx', nn],
                alfy=twiss['alfy', nn], bety=twiss['bety', nn],
                dx=twiss['dx', nn], dpx=twiss['dpx', nn],
                dy=twiss['dy', nn], dpy=twiss['dpy', nn],
                x_co=twiss['x', nn], px_co=twiss['px', nn],
                y_co=twiss['y', nn], py_co=twiss['py', nn],
                zeta_co=twiss['zeta', nn], delta_co=twiss['delta', nn],
                gemitt_x=gemitt_x, gemitt_y=gemitt_y,
                sigma_z=float(sigma_z), sigma_delta=float(sigma_delta),
                n_scattering_events=n_scattering_events,
                interaction_rate=float(bunch_intensity * f_rev * sum(
                    integrated[kk] * xsecs[kk] for kk in species)),
                rate_above_theta_max=float(bunch_intensity * f_rev * sum(
                    integrated[kk] * xsecs_above[kk] for kk in species)),
                mean_free_path=(np.inf if inverse_mfp <= 0
                                else 1.0 / inverse_mfp),
                process=process,
                brems_energy_cut=float(brems_energy_cut),
                coulomb_theta_min=float(coulomb_theta[0]),
                coulomb_theta_max=float(coulomb_theta[1]),
                atomic_densities=atomic_densities,
                seed=seed,
            )

    def generate(self, at=None, *, replica=0, n_replicas=1, seed=None,
                 return_log=False):
        """
        Generate beam-gas-scattered particles at the markers.

        At each marker, ``n_scattering_events`` particles are drawn from the
        local matched Gaussian beam and each of them interacts once with the
        residual gas, on a species drawn according to its contribution to the
        local interaction rate. The weights of the particles of a marker sum
        to the interaction rate of its section [1/s], divided by
        ``n_replicas``.

        The particles of a marker are tracked from the marker itself, e.g.
        ``line.track(part, ele_start=name, ele_stop=name, num_turns=...)``,
        and the samples of different markers or replicas, generated in the
        same or in different jobs, can then be merged with
        :meth:`xtrack.Particles.merge` and given to :meth:`analyse`.

        Parameters
        ----------
        at : str, sequence of str, or None, optional
            Markers at which particles are generated. If ``None`` (default),
            all the markers of the line.
        replica : int, optional
            Index of the replica: replicas of a marker are independent
            samples of the same section, e.g. generated in different jobs.
            Default ``0``.
        n_replicas : int, optional
            Number of replicas of each marker that are merged. The weights are
            divided by it, so that the merged replicas represent the
            interaction rate once. Default ``1``.
        seed : int or None, optional
            Seed of the random-number streams. If ``None`` (default), the seed
            stored on each marker by :meth:`configure`. Each marker and
            replica has its own stream, derived from the seed, the name of the
            marker and the replica: the particles of a marker do not depend on
            the other markers that are generated, nor on the order.
        return_log : bool, optional
            If ``True``, also return the record of the interactions, needed by
            :meth:`analyse` for the per-event diagnostics. The record of the
            last generation is also kept, in memory, as the ``scatter_log`` of
            each marker.

        Returns
        -------
        particles : xtrack.Particles or dict
            If ``at`` is a string, the particles generated at that marker.
            Otherwise, a mapping ``{name: xtrack.Particles}`` in line order.
        log : xtrack.Table or dict
            Only if ``return_log`` is ``True``: the record of the
            interactions, one row per particle, with columns ``name``, ``s``,
            ``particle_id``, ``gas``, ``theta`` [rad], ``photon_energy`` [eV]
            and ``weight`` [1/s]; or a mapping ``{name: xtrack.Table}``.
            ``particle_id`` identifies the particle in the sample of its
            marker, and survives tracking.
        """
        line = self.line
        names = self._markers(at)

        if (isinstance(n_replicas, (bool, np.bool_))
                or not isinstance(n_replicas, (int, np.integer))
                or n_replicas < 1):
            raise ValueError("`n_replicas` must be a positive integer.")
        if (isinstance(replica, (bool, np.bool_))
                or not isinstance(replica, (int, np.integer))
                or not 0 <= replica < n_replicas):
            raise ValueError(
                f"`replica` must be an integer in [0, n_replicas) "
                f"(got {replica} with n_replicas={n_replicas}).")
        if seed is not None:
            seed = _resolve_seed(seed)

        particle_ref = line.particle_ref
        q0 = _lepton_charge(particle_ref)
        markers = {nn: line.element_dict[nn] for nn in names}
        for nn, marker in markers.items():
            _check_configured(nn, marker, particle_ref)

        tab = line.get_table()
        element_index = {nn: ii for ii, nn in enumerate(line.element_names)}

        particles, logs = {}, {}
        for nn in names:
            marker = markers[nn]
            particles[nn], logs[nn] = _generate_marker(
                nn, marker,
                particle_ref=particle_ref,
                q0=q0,
                at_element=element_index[nn],
                s=tab['s', nn],
                rng=_stream(marker.seed if seed is None else seed,
                            nn, replica),
                weight_scale=1.0 / n_replicas,
            )
            marker.scatter_log = logs[nn]

        if isinstance(at, str):
            particles, logs = particles[at], logs[at]
        if return_log:
            return particles, logs
        return particles

    def analyse(self, particles=None, *, at=None, logs=None):
        """
        Analyse the beam-gas interaction rates and the tracked particles.

        Without particles, the result contains the interaction rates
        configured on the markers. With the particles generated at the
        markers and tracked, it also contains the loss rate: the total weight
        of the scattered particles that are lost (any ``state <= 0``, e.g. on
        an aperture or in a collimator). Secondaries produced while tracking
        carry the weight of their parent and are not counted.

        ``particles`` can be the sample of one marker, the particles of
        several markers merged with :meth:`xtrack.Particles.merge`, possibly
        coming from separate jobs, or the samples of each marker separately.
        The particles do not record the marker at which they were generated,
        so the per-marker counts and the Monte Carlo error of the loss rate
        are only available for the samples given separately, or for the
        particles of a single marker selected with ``at``. The per-event
        diagnostics (``cutoff_scan``, ``interaction_log`` and the warnings on
        the generation window) also need the records of the interactions,
        returned by ``generate(..., return_log=True)``, given as ``logs`` with
        the same structure as ``particles``.

        Parameters
        ----------
        particles : xtrack.Particles, dict, or None, optional
            The particles: ``None`` (only the configured rates), an
            ``xtrack.Particles`` object, or a mapping ``{name: particles}``
            from marker names to the ``xtrack.Particles`` of the marker or a
            list of its replicas.
        at : str, sequence of str, or None, optional
            Markers whose particles are given, when ``particles`` is not a
            mapping. If ``None`` (default), all the markers of the line.
        logs : dict or None, optional
            Records of the interactions, as a mapping ``{name: log}`` or
            ``{name: [logs of the replicas]}`` matching ``particles``. Each
            record can be an ``xtrack.Table``, a ``pandas.DataFrame`` or a
            dict of arrays.

        Returns
        -------
        result : BeamGasResult
            Interaction rates and, with tracked particles, the loss rate and
            the beam-gas lifetime.
        """
        line = self.line

        samples = None
        if isinstance(particles, dict):
            if at is not None:
                raise ValueError(
                    "`at` cannot be given with a mapping of particles: its "
                    "keys select the markers.")
            names = self._markers(list(particles))
            samples = {}
            for nn in names:
                replicas = particles[nn]
                if not isinstance(replicas, (list, tuple)):
                    replicas = [replicas]
                if len(replicas) == 0:
                    raise ValueError(f"No particles given for marker '{nn}'.")
                samples[nn] = list(replicas)
            all_particles = [pp for vv in samples.values() for pp in vv]
        else:
            names = self._markers(at)
            all_particles = [] if particles is None else [particles]
        for pp in all_particles:
            if not isinstance(pp, xt.Particles):
                raise TypeError(
                    "The particles must be `xtrack.Particles` objects "
                    f"(got {type(pp).__name__}).")

        sample_logs = None
        if logs is not None:
            if samples is None:
                raise ValueError(
                    "`logs` need the particles given marker by marker, as a "
                    "mapping `{name: particles}`: the merged particles cannot "
                    "be matched to the records.")
            if set(logs) != set(samples):
                raise ValueError(
                    "`logs` must have the same markers as `particles`.")
            sample_logs = {}
            for nn in names:
                replica_logs = logs[nn]
                if isinstance(replica_logs, (list, tuple)):
                    replica_logs = list(replica_logs)
                else:
                    replica_logs = [replica_logs]
                if len(replica_logs) != len(samples[nn]):
                    raise ValueError(
                        f"Marker '{nn}': {len(samples[nn])} samples but "
                        f"{len(replica_logs)} records of interactions.")
                sample_logs[nn] = replica_logs

        markers = {nn: line.element_dict[nn] for nn in names}
        for nn, marker in markers.items():
            _check_configured(nn, marker, line.particle_ref)
        for label, values in (
                ('processes', {mm.process for mm in markers.values()}),
                ('bunch intensities',
                 {float(mm.bunch_intensity) for mm in markers.values()})):
            if len(values) > 1:
                raise ValueError(
                    f"The markers are configured with different {label}: "
                    f"{sorted(values)}.")
        first = markers[names[0]]
        process = first.process
        bunch_intensity = float(first.bunch_intensity)

        f_rev = (line.particle_ref.beta0[0] * sc.c / line.get_length())
        tab = line.get_table()
        mean_free_path = np.array([mm.mean_free_path
                                   for mm in markers.values()])
        with np.errstate(divide='ignore'):
            local_rate = bunch_intensity * f_rev / mean_free_path
        data = {
            'name': np.array(names),
            's': np.array([float(tab['s', nn]) for nn in names]),
            'ds': np.array([mm.ds for mm in markers.values()]),
            'mean_free_path': mean_free_path,
            'local_rate': local_rate,
            'interaction_rate': np.array([mm.interaction_rate
                                          for mm in markers.values()]),
        }
        rate_scattering = float(np.sum(data['interaction_rate']))
        rate_above_theta_max = None
        if process == 'coulomb':
            rate_above_theta_max = float(sum(
                mm.rate_above_theta_max for mm in markers.values()))

        tracked = any(_tracked(pp) for pp in all_particles)
        rate_generated = None
        rate_tracking = None
        rate_tracking_error = None
        cutoff_scan = None
        interaction_log = None

        if samples is not None:
            sums = {nn: [_sample_sums(pp) for pp in samples[nn]]
                    for nn in names}

            def per_marker(key):
                return np.array([sum(ss[key] for ss in sums[nn])
                                 for nn in names])

            data['n_replicas'] = np.array([len(samples[nn]) for nn in names])
            data['num_particles'] = per_marker('num_particles')
            data['sum_weight'] = per_marker('sum_weight')
            rate_generated = float(np.sum(data['sum_weight']))
            if tracked:
                data['num_lost_particles'] = per_marker('num_lost_particles')
                data['sum_lost_weight'] = per_marker('sum_lost_weight')
                variance = per_marker('lost_weight_variance')
                data['sum_lost_weight_error'] = np.sqrt(variance)
                rate_tracking = float(np.sum(data['sum_lost_weight']))
                rate_tracking_error = float(np.sqrt(np.sum(variance)))
        elif particles is not None:
            sums = _sample_sums(particles)
            rate_generated = sums['sum_weight']
            if tracked:
                rate_tracking = sums['sum_lost_weight']
            if len(names) == 1:
                # Events of a single marker are identically distributed,
                # also when replicas are merged
                data['num_particles'] = np.array([sums['num_particles']])
                data['sum_weight'] = np.array([sums['sum_weight']])
                if tracked:
                    data['num_lost_particles'] = np.array(
                        [sums['num_lost_particles']])
                    data['sum_lost_weight'] = np.array(
                        [sums['sum_lost_weight']])
                    rate_tracking_error = float(
                        np.sqrt(sums['lost_weight_variance']))

        if sample_logs is not None:
            column = 'theta' if process == 'coulomb' else 'photon_energy'
            if process == 'coulomb':
                window = (float(first.coulomb_theta_min),
                          float(first.coulomb_theta_max))
            else:
                # Up to the kinetic energy of the beam
                window = (float(first.brems_energy_cut), float(
                    _energy_from_momentum(first.p0c) - ELECTRON_MASS_EV))
            x, weight, lost, i_sample, tables = [], [], [], [], []
            for nn in names:
                for rr, (pp, log) in enumerate(zip(samples[nn],
                                                   sample_logs[nn])):
                    xx, ww, ll = _event_losses(pp, log, column)
                    x.append(xx)
                    weight.append(ww)
                    lost.append(ll)
                    i_sample.append(np.full(len(xx), len(tables)))
                    tables.append(_log_table(nn, rr, log))
            x, weight, lost, i_sample = (
                np.concatenate(vv) for vv in (x, weight, lost, i_sample))
            interaction_log = xt.Table({
                kk: np.concatenate([tt[kk] for tt in tables])
                for kk in tables[0]})
            if tracked:
                cutoff_scan = _cutoff_scan(
                    window, x, weight, lost, i_sample, bunch_intensity,
                    _CUTOFF_SCAN_BINS_PER_DECADE)
                _warn_on_truncation(process, window, x, weight * lost,
                                    cutoff_scan, rate_above_theta_max)

        def lifetime(rate):
            return np.inf if rate == 0 else float(bunch_intensity / rate)

        lifetime_tracking = None
        lifetime_tracking_error = None
        if rate_tracking is not None:
            lifetime_tracking = lifetime(rate_tracking)
            if rate_tracking_error is not None and rate_tracking > 0:
                lifetime_tracking_error = float(
                    lifetime_tracking * rate_tracking_error / rate_tracking)

        return BeamGasResult(
            element_names=names,
            local_rates=xt.Table(data),
            rate_scattering=rate_scattering,
            lifetime_scattering=lifetime(rate_scattering),
            rate_tracking=rate_tracking,
            lifetime_tracking=lifetime_tracking,
            tracked=tracked,
            rate_generated=rate_generated,
            rate_tracking_error=rate_tracking_error,
            lifetime_tracking_error=lifetime_tracking_error,
            cutoff_scan=cutoff_scan,
            rate_above_theta_max=rate_above_theta_max,
            interaction_log=interaction_log,
        )

    def run(self, n_turns, *, at=None, n_replicas=1, keep_particles=False,
            with_progress=False):
        """
        Generate, track and analyse the beam-gas particles of the markers.

        Convenience helper doing in one call the steps of a beam-gas
        simulation: the particles of each marker (and replica) are generated
        with :meth:`generate`, tracked from the marker back to it with
        ``line.track(part, ele_start=name, ele_stop=name, num_turns=n_turns)``
        and analysed with :meth:`analyse`, together with the records of the
        interactions. Doing these steps explicitly gives the same result, and
        allows splitting them over separate jobs or refining the loss
        locations before the analysis.

        Parameters
        ----------
        n_turns : int
            Number of turns to track.
        at : str, sequence of str, or None, optional
            Markers at which particles are generated. If ``None`` (default),
            all the markers of the line.
        n_replicas : int, optional
            Number of replicas of each marker.
        keep_particles : bool, optional
            If ``True``, the result keeps the particles of each marker (its
            replicas merged), all of them merged, and the lost ones (if any).
        with_progress : bool, optional
            Forwarded to ``line.track``.

        Returns
        -------
        result : BeamGasResult
            Result of :meth:`analyse` for the particles of each marker, with
            the records of the interactions.
        """
        if (isinstance(n_turns, (bool, np.bool_))
                or not isinstance(n_turns, (int, np.integer)) or n_turns < 1):
            raise ValueError("`n_turns` must be a positive integer.")

        line = self.line
        names = self._markers(at)
        samples = {nn: [] for nn in names}
        logs = {nn: [] for nn in names}
        for replica in range(n_replicas):
            generated, records = self.generate(
                at=names, replica=replica, n_replicas=n_replicas,
                return_log=True)
            for nn, part in generated.items():
                line.track(part, ele_start=nn, ele_stop=nn,
                           num_turns=n_turns, with_progress=with_progress)
                samples[nn].append(part)
                logs[nn].append(records[nn])

        result = self.analyse(samples, logs=logs)
        if keep_particles:
            result.particles_by_element = {
                nn: vv[0] if len(vv) == 1 else xt.Particles.merge(vv)
                for nn, vv in samples.items()}
            result.particles = xt.Particles.merge(
                list(result.particles_by_element.values()))
            lost = _lost(result.particles)
            if np.any(lost):
                result.lost_particles = result.particles.filter(lost)
        return result

    def _check_sections_cover_line(self):
        """
        Check that every part of the line belongs to the section of a marker.

        The section of a marker goes from the previous one to it, and the
        section of the first one from ``s = 0``. So the last marker must sit
        at the end of the line: otherwise the part of the line after it is not
        integrated, and its interaction rate is silently missing.

        Returns
        -------
        None
        """
        if self._last_marker_at_end():
            return
        last = self._markers()[-1]
        tab = self.line.get_table()
        raise ValueError(
            f"The last BeamGasMarker ('{last}') is at "
            f"s = {float(tab['s', last]):.6g} m, but the line ends at "
            f"s = {float(tab.s[-1]):.6g} m: the last BeamGasMarker must sit "
            "at the end of the line, otherwise the part of the line after it "
            "is not accounted for in the interaction rate. "
            "`line.xcoll.beamgas.install_markers(...)` adds it, or place a "
            "marker at `at=line.get_length()`.")

    def _last_marker_at_end(self):
        """
        Whether the last marker sits at the end of the line.

        Returns
        -------
        at_end : bool
            ``True`` if the line has markers and the last one is at the end
            of the line.
        """
        names = self.marker_names
        if len(names) == 0:
            return False
        tab = self.line.get_table()
        return abs(float(tab['s', names[-1]]) - float(tab.s[-1])) <= _S_TOL

    def _markers(self, at=None):
        """
        Validated names of the selected markers, in line order.

        Parameters
        ----------
        at : str, sequence of str, or None
            Selected markers. ``None`` selects all the markers of the line.

        Returns
        -------
        names : list of str
            Names of the selected markers, in line order.
        """
        names = self.marker_names
        if at is None:
            if len(names) == 0:
                raise ValueError(
                    "The line does not contain any BeamGasMarker. Install "
                    "them with `line.xcoll.beamgas.install_markers(...)`.")
            return names

        if isinstance(at, str):
            selected = [at]
        elif hasattr(at, '__iter__'):
            selected = list(at)
        else:
            raise TypeError(
                "Beam-gas markers are selected by name "
                f"(got {type(at).__name__}).")
        if len(selected) == 0:
            raise ValueError("No beam-gas marker selected.")
        position = {nn: ii for ii, nn in enumerate(names)}
        for nn in selected:
            if not isinstance(nn, str):
                raise TypeError(
                    "Beam-gas markers are selected by name "
                    f"(got {type(nn).__name__}).")
            if nn in position:
                continue
            if nn in set(self.line.element_names):
                raise TypeError(
                    f"Element '{nn}' is not a BeamGasMarker "
                    f"(got {type(self.line.element_dict[nn]).__name__}).")
            raise ValueError(f"Element '{nn}' is not present in the line.")

        return sorted(set(selected), key=position.__getitem__)
