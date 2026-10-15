# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2026.                 #
# ######################################### #

from dataclasses import dataclass

import numpy as np

import xtrack as xt


def _allocated(particles):
    """
    Mask of the slots of the particle buffer that hold a particle.

    The generated particles leave room for four times their number; the
    unused slots have a negative particle id and state.
    """
    return particles.particle_id >= 0


def _primary(particles):
    """
    Mask of the beam-gas-scattered particles, excluding any secondaries
    produced while tracking (e.g. by collimator scattering codes).
    """
    return _allocated(particles) & (
        particles.parent_particle_id == particles.particle_id)


def _lost(particles):
    """
    Mask of the lost particles: any state <= 0, i.e. on an aperture (0) or
    in other ways (negative states, e.g. absorbed in a collimator).
    """
    return _allocated(particles) & (particles.state <= 0)


def _tracked(particles):
    """
    Whether the particles have been tracked: some are lost, or have
    completed turns.
    """
    allocated = _allocated(particles)
    return bool(np.any(particles.state[allocated] <= 0)
                or np.any(particles.at_turn[allocated] > 0))


def _sample_sums(particles):
    """
    Counts, weights and variance of the scattered particles of a sample.

    Only the scattered particles (primaries) are counted: secondaries carry
    the weight of their parent. The particles are independent events, so
    the variance of the lost weight is the number of particles times the
    variance of their contributions.

    Parameters
    ----------
    particles : xtrack.Particles
        Particle sample, of one marker or of several replicas of a marker.

    Returns
    -------
    sums : dict
        ``num_particles``, ``sum_weight``, ``num_lost_particles``,
        ``sum_lost_weight`` and ``lost_weight_variance``.
    """
    primary = _primary(particles)
    lost = _lost(particles) & primary
    contribution = particles.weight[primary] * lost[primary]
    n = len(contribution)
    if n > 1:
        variance = float(n / (n - 1) * (np.sum(contribution**2)
                                       - np.sum(contribution)**2 / n))
    else:
        variance = 0.0
    return {
        'num_particles': int(np.sum(primary)),
        'sum_weight': float(np.sum(particles.weight[primary])),
        'num_lost_particles': int(np.sum(lost)),
        'sum_lost_weight': float(np.sum(particles.weight[lost])),
        'lost_weight_variance': max(variance, 0.0),
    }


def _event_losses(particles, log, column):
    """
    Join the record of the interactions of a sample with its tracked states.

    Xtrack reorders the particles while tracking, so the events are matched
    to the particles by ``particle_id``.

    Parameters
    ----------
    particles : xtrack.Particles
        Tracked particles of one sample.
    log : xtrack.Table, pandas.DataFrame or dict
        Record of the interactions of the sample, as returned by
        ``generate(..., return_log=True)``.
    column : {'theta', 'photon_energy'}
        Generation variable.

    Returns
    -------
    x, weight : ndarray
        Generation variable and weight of each event.
    lost : ndarray of bool
        Whether the particle of each event is lost.
    """
    log_ids = np.asarray(log['particle_id'], dtype=np.int64)
    allocated = np.flatnonzero(_allocated(particles))
    ids = particles.particle_id[allocated]
    order = np.argsort(ids)
    position = np.searchsorted(ids, log_ids, sorter=order)
    position = np.minimum(position, len(ids) - 1)
    idx = allocated[order[position]]
    weight = np.asarray(log['weight'], dtype=float)
    # The weights are not changed by tracking: they identify the sample
    if (not np.array_equal(particles.particle_id[idx], log_ids)
            or not np.array_equal(particles.weight[idx], weight)):
        raise ValueError(
            "The record of the interactions does not match the particles: "
            "give each sample with the record returned with it by "
            "`generate(..., return_log=True)`.")
    return (np.asarray(log[column], dtype=float), weight,
            particles.state[idx] <= 0)


def _cutoff_scan(window, x, weight, lost, i_sample, bunch_intensity,
                 bins_per_decade):
    """
    Tabulate the loss rate as a function of the lower generation cut.

    The events above a cut are an unbiased importance sample of a simulation
    generated with its lower cut set there, so the loss rate of every such
    simulation follows from the single tracked sample. The loss rate of each
    sample is a sum of independent event contributions, so its variance is
    the number of events times their sample variance.

    Parameters
    ----------
    window : tuple of float
        Range of the generation variable.
    x, weight, lost, i_sample : ndarray
        Generation variable, weight, loss flag and sample index of each event.
    bunch_intensity : float
        Number of particles in the bunch.
    bins_per_decade : int
        Resolution of the scan.

    Returns
    -------
    table : xtrack.Table
        See :attr:`BeamGasResult.cutoff_scan`.
    """
    lo, hi = window
    n_bins = max(int(np.ceil(bins_per_decade * np.log10(hi / lo))), 1)
    edges = np.geomspace(lo, hi, n_bins + 1)
    i_bin = np.clip(np.searchsorted(edges, x, side='right') - 1,
                    0, n_bins - 1)
    contribution = weight * lost
    n_samples = int(np.max(i_sample)) + 1

    def per_bin(values, sample_wise=False):
        if not sample_wise:
            return np.bincount(i_bin, weights=values, minlength=n_bins)
        flat = np.bincount(i_sample * n_bins + i_bin, weights=values,
                           minlength=n_samples * n_bins)
        return flat.reshape(n_samples, n_bins)

    def above(binned):
        # Sum over all the bins from the current one upwards
        return np.flip(np.cumsum(np.flip(binned, axis=-1), axis=-1),
                       axis=-1)

    num_events = per_bin(np.ones_like(x))
    weight_in_bin = per_bin(weight)
    with np.errstate(invalid='ignore', divide='ignore'):
        loss_probability = per_bin(contribution) / weight_in_bin

    s1 = above(per_bin(contribution, sample_wise=True))
    s2 = above(per_bin(contribution**2, sample_wise=True))
    n = np.bincount(i_sample, minlength=n_samples)[:, None]
    with np.errstate(invalid='ignore', divide='ignore'):
        variance = np.where(n > 1, n * (s2 - s1**2 / n) / (n - 1), 0.0)
    rate = s1.sum(axis=0)
    with np.errstate(divide='ignore'):
        lifetime = np.where(rate > 0, bunch_intensity / rate, np.inf)

    return xt.Table({
        'cut': edges[:-1],
        'num_events': num_events.astype(int),
        'num_lost': per_bin(lost.astype(float)).astype(int),
        'loss_probability': loss_probability,
        'rate_tracking': rate,
        'rate_tracking_error': np.sqrt(np.maximum(variance, 0.0)
                                       .sum(axis=0)),
        'lifetime_tracking': lifetime,
    }, index='cut')


@dataclass
class BeamGasResult:
    """
    Result of a beam-gas simulation.

    It is returned by :meth:`xcoll.BeamGasGenerator.analyse` and
    :meth:`xcoll.BeamGasGenerator.run`. The result separates quantities
    inferred from the cross sections (the interaction rates configured on the
    markers) from quantities obtained from the particles (the generated and
    the lost weight).

    Parameters
    ----------
    element_names : list
        Names of the beam-gas markers included in the result.
    local_rates : xtrack.Table
        Per-marker diagnostics table. It always contains ``name``, ``s``,
        ``ds``, ``mean_free_path``, ``local_rate`` and ``interaction_rate``.
        When the particles of each marker are known, it also contains
        ``num_particles`` and ``sum_weight`` and, once tracked,
        ``num_lost_particles`` and ``sum_lost_weight``.
    rate_scattering : float
        Total beam-gas interaction rate represented by the result [1/s]. Only
        the part of the cross section that is actually generated is counted,
        i.e. photon energies above ``brems_energy_cut`` or scattering angles
        within ``coulomb_theta``.
    lifetime_scattering : float
        Beam lifetime that would result if every generated interaction led to
        a loss [s]. It is a lower limit on the beam-gas lifetime.
    rate_tracking : float or None
        Beam-gas loss rate: weighted rate of the scattered particles lost in
        tracking [1/s]. ``None`` if the particles were not tracked.
    lifetime_tracking : float or None
        Beam-gas lifetime, inferred from the loss rate ``rate_tracking`` [s].
        ``None`` if the particles were not tracked.
    tracked : bool
        Whether the particles were tracked.
    rate_generated : float or None
        Total weight of the scattered particles [1/s]. ``None`` without
        particles.
    rate_tracking_error : float or None
        Monte Carlo standard error of ``rate_tracking`` [1/s]. The scattered
        particles are independent events, so it follows from the spread of
        their contributions within each sample. ``None`` unless the particles
        of each marker are known.
    lifetime_tracking_error : float or None
        Monte Carlo standard error of ``lifetime_tracking`` [s], propagated
        from ``rate_tracking_error`` to first order.
    cutoff_scan : xtrack.Table or None
        Loss rate as a function of the lower generation cut, obtained from the
        single tracked sample: since every event carries its importance
        weight, the events above a cut ``c`` are an unbiased sample of a
        simulation generated with its lower cut set to ``c``. One row per
        log-spaced ``cut``, from the configured lower cut (``theta_min`` for
        Coulomb scattering, ``brems_energy_cut`` for bremsstrahlung) upwards,
        with columns ``cut``, ``num_events`` and ``num_lost`` (in ``[cut,
        next cut)``), ``loss_probability`` (weighted fraction of those events
        that is lost), and ``rate_tracking``, ``rate_tracking_error`` and
        ``lifetime_tracking`` that a simulation with its lower cut set to
        ``cut`` would give. If ``lifetime_tracking`` has not saturated in the
        first rows, the configured cut is too high. ``None`` unless the
        records of the interactions are given.
    rate_above_theta_max : float or None
        Interaction rate for Coulomb scattering above ``theta_max`` [1/s],
        which is not generated. It is an upper limit on the loss rate missing
        because of the upper cut, and equals it when every event near
        ``theta_max`` is lost. ``None`` for bremsstrahlung.
    particles_by_element : dict or None
        Mapping ``{element_name: xtrack.Particles}`` with the particles of
        each marker, when kept.
    particles : xtrack.Particles or None
        Merged particles, when kept.
    lost_particles : xtrack.Particles or None
        Subset of ``particles`` lost in tracking, when kept and if any is
        lost. It is a copy taken right after tracking.
    interaction_log : xtrack.Table or None
        Record of the interactions, with columns ``name``, ``replica``,
        ``s``, ``particle_id``, ``gas``, ``theta``, ``photon_energy`` and
        ``weight``. ``particle_id`` refers to the sample of the marker and
        replica of the same row. ``None`` unless the records are given.
    gas_density : xtrack.Table or None
        Gas-density profile of the deprecated :class:`xcoll.BeamGasStudy`.

    Attributes
    ----------
    element_names, local_rates, rate_scattering, lifetime_scattering,
    rate_tracking, lifetime_tracking, tracked, rate_generated,
    rate_tracking_error, lifetime_tracking_error, cutoff_scan,
    rate_above_theta_max, particles_by_element, particles, lost_particles,
    interaction_log, gas_density
        See above.
    """
    element_names: list
    local_rates: xt.Table
    rate_scattering: float
    lifetime_scattering: float
    rate_tracking: float | None
    lifetime_tracking: float | None
    tracked: bool
    rate_generated: float | None = None
    rate_tracking_error: float | None = None
    lifetime_tracking_error: float | None = None
    cutoff_scan: xt.Table | None = None
    rate_above_theta_max: float | None = None
    particles_by_element: dict | None = None
    particles: xt.Particles | None = None
    lost_particles: xt.Particles | None = None
    interaction_log: xt.Table | None = None
    gas_density: xt.Table | None = None
