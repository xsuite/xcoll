# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2026.                 #
# ######################################### #

from dataclasses import dataclass
from warnings import warn

import numpy as np
from scipy import constants as sc

import xtrack as xt
from xtrack.particles import pdg

from ..beam_elements.beamgas import (BeamGasScattering, _resolve_process,
                                     _resolve_n_scattering_events)
from .cross_sections import (BremsstrahlungCalculator,
                             CoulombScatteringCalculator,
                             atomic_number_from_symbol)

C_LIGHT = sc.c

PDG_ID_ELECTRON = 11
PDG_ID_POSITRON = -11


def _allocated_mask(particles):
    """
    Mask selecting the allocated slots of a particles object.

    An :class:`xtrack.Particles` object created with a spare capacity (to make
    room for the secondaries produced by collimators) carries unallocated
    slots, flagged with a particle id and a state of ``-999999999``. They must
    be excluded from any statistics.

    Parameters
    ----------
    particles : xtrack.Particles
        Particles object to inspect.

    Returns
    -------
    mask : ndarray of bool
        ``True`` for the allocated slots.
    """
    return particles.particle_id >= 0


def _lost_mask(particles):
    """
    Mask selecting the allocated particles that are not alive any more.

    Xcoll flags the different loss mechanisms with distinct non-positive
    states (see :data:`xcoll.particle_states`), so all of them are counted as
    losses.

    Parameters
    ----------
    particles : xtrack.Particles
        Particles object to inspect.

    Returns
    -------
    mask : ndarray of bool
        ``True`` for the allocated particles that are lost.
    """
    return _allocated_mask(particles) & (particles.state <= 0)


@dataclass
class BeamGasResult:
    """
    Result returned by :meth:`BeamGasStudy.run`.

    The result separates quantities inferred from the scattering model from
    quantities obtained after tracking the generated particles. The particle
    samples are optional and are stored only when :meth:`BeamGasStudy.run` is
    called with ``keep_particles=True``.

    Parameters
    ----------
    element_names : list
        Names of the beam-gas scattering elements included in the result.
    gas_density : xtrack.Table
        Gas-density profile used by the study.
    local_rates : xtrack.Table
        Per-element diagnostics table. It always contains ``name``, ``s``,
        ``ds``, ``mean_free_path``, ``local_rate`` and ``interaction_rate``.
        When particles are generated it also contains ``num_particles`` and
        ``sum_weight``. When tracking is enabled it additionally contains
        ``num_lost_particles`` and ``sum_lost_weight``.
    rate_scattering : float
        Total beam-gas interaction rate represented by the result [1/s]. Only
        the part of the cross section that is actually generated is counted,
        i.e. photon energies above ``brems_energy_cut`` or scattering angles
        within ``coulomb_theta``. ``local_rates.sum_weight`` reproduces it up
        to the Monte Carlo noise of the importance sampling.
    lifetime_scattering : float
        Beam lifetime that would result if every generated interaction led to
        a loss [s]. It is a lower bound on the beam-gas lifetime, and is a
        useful diagnostic only when the generated part of the cross section is
        restricted to the loss-relevant range.
    rate_tracking : float or None
        Weighted loss rate after tracking the generated particles [1/s].
        ``None`` when tracking is disabled.
    lifetime_tracking : float or None
        Beam-gas lifetime inferred from ``rate_tracking`` [s]. ``None`` when
        tracking is disabled. This is the physically meaningful result of the
        study.
    tracked : bool
        Whether the generated particles were tracked to determine losses.
    rate_tracking_error : float or None
        Monte Carlo standard error of ``rate_tracking`` [1/s]. ``None`` when
        tracking is disabled.
    lifetime_tracking_error : float or None
        Monte Carlo standard error of ``lifetime_tracking`` [s], propagated
        from ``rate_tracking_error`` to first order. ``None`` when tracking is
        disabled.
    cutoff_scan : xtrack.Table or None
        Loss rate as a function of the lower generation cut, obtained from
        the single tracked sample: since every event carries its importance
        weight, the events above a cut ``c`` are an unbiased sample of a
        study generated with its lower cut set to ``c``. One row per
        log-spaced ``cut``, from the configured lower cut (``theta_min`` for
        Coulomb scattering, ``brems_energy_cut`` for bremsstrahlung) upwards,
        with columns ``cut``, ``num_events`` and ``num_lost`` (in
        ``[cut, next cut)``), ``loss_probability`` (weighted fraction of
        those events that is lost), and ``rate_tracking``,
        ``rate_tracking_error`` and ``lifetime_tracking`` that a study with
        its lower cut set to ``cut`` would give. If ``lifetime_tracking``
        has not saturated in the first rows, the configured cut is too high.
        Only losses of the generated primaries are included. ``None`` when
        tracking is disabled.
    rate_above_theta_max : float or None
        Interaction rate for Coulomb scattering above ``theta_max`` [1/s],
        which is not generated. It is an upper bound on the loss rate missing
        because of the upper cut, and equals it when every event near
        ``theta_max`` is lost (``cutoff_scan.loss_probability`` close to one
        in the last rows). ``None`` unless ``process='coulomb'`` and
        tracking is enabled.
    particles_by_element : dict or None
        Mapping ``{element_name: xtrack.Particles}`` with the generated
        particles for each scattering element. ``None`` unless
        ``keep_particles=True``.
    particles : xtrack.Particles or None
        Merged generated particle sample. ``None`` unless
        ``keep_particles=True``.
    lost_particles : xtrack.Particles or None
        Subset of ``particles`` lost during tracking. ``None`` unless both
        ``track=True`` and ``keep_particles=True``. It is a copy taken right
        after tracking, so operations that modify ``particles`` in place (such
        as :class:`xtrack.LossLocationRefinement`) are not reflected in it.
    interaction_log : xtrack.Table or None
        Per-event record of the generated interactions, with columns
        ``name``, ``s``, ``particle_id``, ``gas``, ``theta``,
        ``photon_energy`` and ``weight``. ``None`` unless
        ``keep_particles=True``. ``particle_id`` refers to the particle sample
        of the element named in the same row, i.e. to
        ``particles_by_element[name]``; :meth:`xtrack.Particles.merge`
        renumbers the particles, so the log cannot be joined to ``particles``
        by id.

    Attributes
    ----------
    element_names, gas_density, local_rates, rate_scattering,
    lifetime_scattering, rate_tracking, lifetime_tracking, tracked,
    rate_tracking_error, lifetime_tracking_error, cutoff_scan,
    rate_above_theta_max, particles_by_element, particles, lost_particles,
    interaction_log
        See above.
    """
    element_names: list
    gas_density: xt.Table
    local_rates: xt.Table
    rate_scattering: float
    lifetime_scattering: float
    rate_tracking: float | None
    lifetime_tracking: float | None
    tracked: bool
    rate_tracking_error: float | None = None
    lifetime_tracking_error: float | None = None
    cutoff_scan: xt.Table | None = None
    rate_above_theta_max: float | None = None
    particles_by_element: dict | None = None
    particles: xt.Particles | None = None
    lost_particles: xt.Particles | None = None
    interaction_log: xt.Table | None = None


class BeamGasStudy:
    """
    Monte Carlo study of beam-residual-gas scattering in a line.

    See :meth:`__init__` for the full description of the parameters.
    """

    def __init__(self, line=None, gas_density=None, process='brems',
                 elements=None, twiss=None,
                 nemitt_x=None, nemitt_y=None,
                 gemitt_x=None, gemitt_y=None,
                 sigma_z=None, sigma_delta=None,
                 bunch_intensity=1.0,
                 n_scattering_events=None,
                 brems_energy_cut=10e3,
                 coulomb_theta=(1e-7, 50e-3),
                 seed=None, **kwargs):
        """
        Build a beam-gas study and validate the line, optics, beam parameters
        and gas-density profile.

        The configured study computes the local beam-gas interaction rate at
        each :class:`xcoll.BeamGasScattering` element, integrates it over the
        lattice section represented by the element, and configures the element
        with the local optics, beam parameters and gas composition needed to
        generate weighted Monte Carlo particles.

        After construction, call :meth:`run` to obtain a
        :class:`BeamGasResult`. With tracking disabled, :meth:`run` returns
        the integrated interaction rate and a per-element diagnostics table.
        With tracking enabled, it also returns the tracked loss rate and the
        corresponding beam-gas lifetime. The generated and lost particle
        samples are returned only when ``keep_particles=True`` is passed to
        :meth:`run`.

        Parameters
        ----------
        line : xtrack.Line
            Line containing the :class:`xcoll.BeamGasScattering` elements to
            configure. It must have a ``particle_ref`` with its PDG id set to
            an electron (``11``) or a positron (``-11``).
        gas_density : xtrack.Table
            Residual-gas density profile. It must contain a column ``s`` with
            the longitudinal positions [m] and one column per gas species,
            named after the chemical symbol of the species (``'H'``, ``'N'``,
            ``'Ar'``, ...) and holding the *atomic* density [atoms/m^3]. For
            a molecular gas, use the atomic density of each constituent, e.g.
            ``N = 2 * n_molecules`` for N2.
        process : {'brems', 'coulomb'}, optional
            Beam-gas process to simulate: bremsstrahlung on the gas nuclei, or
            elastic Coulomb scattering off them. Default ``'brems'``.
        elements : str, sequence of str, or None, optional
            Beam-gas scattering elements included in the study. If ``None``,
            all :class:`xcoll.BeamGasScattering` elements in the line are used.
        twiss : xtrack.TwissTable or None, optional
            Twiss table used for the local optics. If ``None``, it is computed
            when :meth:`initialise_beamgas` is called.
        nemitt_x, nemitt_y : float, optional
            Normalised horizontal and vertical emittances. Mutually exclusive
            with ``gemitt_x`` and ``gemitt_y``.
        gemitt_x, gemitt_y : float, optional
            Geometric horizontal and vertical emittances. Mutually exclusive
            with ``nemitt_x`` and ``nemitt_y``.
        sigma_z : float
            RMS bunch length [m].
        sigma_delta : float
            RMS relative momentum spread.
        bunch_intensity : float, optional
            Number of real particles in the bunch. The interaction and loss
            rates are proportional to it; the lifetime is not, since it
            cancels. Default 1, which makes the reported rates per stored
            particle.
        n_scattering_events : int
            Number of Monte Carlo interactions generated at each scattering
            element.
        brems_energy_cut : float, optional
            Lower limit of the generated bremsstrahlung photon energy [eV].
            Only used when ``process='brems'``. Setting it close to the
            momentum acceptance of the machine concentrates the statistics on
            the events that can lead to a loss. Default 10 keV.
        coulomb_theta : tuple of float, optional
            Range of the generated Coulomb scattering angle [rad]. Only used
            when ``process='coulomb'``. Losses from angles outside this range
            are missing from the result, so ``theta_min`` must lie below the
            smallest angle that can cause a loss and ``theta_max`` high enough
            that the cross section above it is negligible. Thanks to the
            log-uniform importance sampling, widening the range only costs
            statistics as ``log(theta_max/theta_min)``, so err on the wide
            side. :meth:`run` reports both truncation effects (see
            ``BeamGasResult.cutoff_scan`` and ``rate_above_theta_max``).
            Default ``(1e-7, 50e-3)``.
        seed : int or None, optional
            Seed of the random-number generator. If ``None`` (default), a
            fresh unseeded generator is used.
        **kwargs
            Additional keyword arguments forwarded to ``line.twiss()`` when a
            Twiss table is computed internally.

        Returns
        -------
        None

        Notes
        -----
        Each :class:`xcoll.BeamGasScattering` element represents the lattice
        section between the previous scattering element (or the start of the
        line) and itself. To cover the whole ring, place one scattering
        element at the end of the line; a warning is issued otherwise.

        References
        ----------
        .. [1] Geant4 Collaboration, "Geant4 Physics Reference Manual",
           Rev. 11.x. https://geant4.web.cern.ch/docs/
        """
        # Input validation
        if line is None:
            raise ValueError("`line` is required.")
        if getattr(line, 'particle_ref', None) is None:
            raise ValueError("`line` must have a `particle_ref`.")
        if gas_density is None:
            raise ValueError("`gas_density` is required.")
        if sigma_z is None:
            raise ValueError("`sigma_z` is required.")
        if sigma_delta is None:
            raise ValueError("`sigma_delta` is required.")
        if n_scattering_events is None:
            raise ValueError("`n_scattering_events` is required.")

        particle_ref = line.particle_ref
        pdg_id = int(np.atleast_1d(particle_ref.pdg_id)[0])
        if pdg_id == 0:
            raise ValueError(
                "The reference particle has no PDG id set. The beam-gas "
                "models need to know whether the beam is made of electrons "
                "or positrons.")
        if pdg_id not in (PDG_ID_ELECTRON, PDG_ID_POSITRON):
            raise ValueError(
                "The beam-gas models implemented in Xcoll are only valid for "
                "electron and positron beams, but the reference particle is "
                f"a {pdg.get_name_from_pdg_id(pdg_id)} (PDG id {pdg_id}).")

        # q0 sets the sign of the McKinley-Feshbach interference term, so it
        # must agree with the PDG id rather than be trusted blindly
        q0_from_pdg = -1.0 if pdg_id == PDG_ID_ELECTRON else 1.0
        if not np.isclose(float(particle_ref.q0), q0_from_pdg):
            raise ValueError(
                f"The reference particle is a "
                f"{pdg.get_name_from_pdg_id(pdg_id)} (PDG id {pdg_id}), which "
                f"must have q0={q0_from_pdg:+.0f}, but it has "
                f"q0={float(particle_ref.q0):+g}.")

        self.line = line
        self.particle_ref = particle_ref
        self.pdg_id = pdg_id
        self.twiss = twiss
        self.p0c = float(particle_ref.p0c[0])
        self.q0 = q0_from_pdg
        self.beta0 = float(particle_ref.beta0[0])

        self.process = _resolve_process(process)
        self.elements = self._resolve_elements(elements)
        self.gas_density = self._validate_gas_density(gas_density)
        self.gas_species = [cc for cc in self.gas_density._col_names
                            if cc not in ('name', 's')]
        self.atomic_numbers = {kk: atomic_number_from_symbol(kk)
                               for kk in self.gas_species}

        self.sigma_z = float(sigma_z)
        self.sigma_delta = float(sigma_delta)
        self.bunch_intensity = float(bunch_intensity)
        if self.bunch_intensity < 0:
            raise ValueError("`bunch_intensity` must be non-negative.")
        self.n_scattering_events = _resolve_n_scattering_events(
            n_scattering_events=n_scattering_events, default=None)
        self.brems_energy_cut = float(brems_energy_cut)
        self.coulomb_theta = (float(coulomb_theta[0]),
                              float(coulomb_theta[1]))

        self.seed = seed
        self.rng = np.random.default_rng(seed)

        # Emittance validation
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
            beta0 = particle_ref.beta0[0]
            gamma0 = particle_ref.gamma0[0]
            self.gemitt_x = float(nemitt_x/(beta0*gamma0))
            self.gemitt_y = float(nemitt_y/(beta0*gamma0))
        else:
            self.gemitt_x = float(gemitt_x)
            self.gemitt_y = float(gemitt_y)

        self.kwargs = kwargs

        # Cross-section models, shared by all the scattering elements
        self.calculators = self._build_calculators()
        self.xsecs = {kk: cc.compute_xsec()
                      for kk, cc in self.calculators.items()}

    def __repr__(self):
        return (f"<BeamGasStudy process='{self.process}', "
                f"{len(self.elements)} elements, "
                f"gas={'+'.join(self.gas_species)}>")

    # ######################################################## #
    # Validation helpers
    # ######################################################## #
    def _resolve_elements(self, elements):
        """
        Resolve and validate the beam-gas scattering elements of the study.

        Parameters
        ----------
        elements : str, sequence of str, or None
            Requested elements, or ``None`` for all the
            :class:`xcoll.BeamGasScattering` elements in the line.

        Returns
        -------
        elements : list of str
            Names of the selected elements, in the order they appear in the
            line.
        """
        line = self.line
        tab = line.get_table()

        if elements is None:
            elements = [nn for nn in tab.name[:-1]
                        if isinstance(line[nn], BeamGasScattering)]
            if len(elements) == 0:
                raise ValueError(
                    "The line does not contain any BeamGasScattering. Please "
                    "add them before initialising the BeamGasStudy.")
            return elements

        if isinstance(elements, str):
            elements = [elements]
        else:
            elements = list(elements)

        if len(elements) == 0:
            raise ValueError(
                "No BeamGasScattering elements selected for this study.")

        for nn in elements:
            if nn not in line.element_names:
                raise ValueError(f"Element '{nn}' is not present in the line.")
            if not isinstance(line[nn], BeamGasScattering):
                raise TypeError(
                    f"Element '{nn}' is not a BeamGasScattering "
                    f"(got {type(line[nn]).__name__}).")

        return elements

    def _validate_gas_density(self, gas_density):
        """
        Validate the gas-density table.

        Parameters
        ----------
        gas_density : xtrack.Table
            Table with a column ``s`` and one column per gas species.

        Returns
        -------
        gas_density : xtrack.Table
            The validated table.
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

        return gas_density

    def _build_calculators(self):
        """
        Build the cross-section model of each gas species.

        Parameters
        ----------
        None

        Returns
        -------
        calculators : dict
            Mapping ``{element_symbol: calculator}``.
        """
        if self.process == 'brems':
            return {kk: BremsstrahlungCalculator(
                        Z, self.p0c, energy_cut=self.brems_energy_cut)
                    for kk, Z in self.atomic_numbers.items()}

        return {kk: CoulombScatteringCalculator(
                    Z, self.p0c, q0=self.q0, theta_lim=self.coulomb_theta)
                for kk, Z in self.atomic_numbers.items()}

    # ######################################################## #
    # Rates
    # ######################################################## #
    def _integrated_atomic_densities(self, s_start, s_end):
        """
        Integrate the atomic density of each species over a lattice section.

        The tabulated density profile is treated as piecewise linear in ``s``,
        so the trapezoidal rule evaluated on the union of the tabulated
        positions inside the section and the two section edges is exact.

        Parameters
        ----------
        s_start, s_end : float
            Limits of the lattice section [m].

        Returns
        -------
        integrals : dict
            Mapping ``{element_symbol: integral}`` with the integrated atomic
            density [atoms/m^2].
        """
        if s_end <= s_start:
            return {kk: 0.0 for kk in self.gas_species}

        s_table = np.asarray(self.gas_density['s'], dtype=float)
        s_inner = s_table[(s_table > s_start) & (s_table < s_end)]
        s_nodes = np.concatenate(([s_start], s_inner, [s_end]))

        integrals = {}
        for kk in self.gas_species:
            n_at = np.interp(s_nodes, s_table,
                             np.asarray(self.gas_density[kk], dtype=float))
            integrals[kk] = float(np.trapezoid(n_at, s_nodes))

        return integrals

    def _local_rate(self, s):
        """
        Compute the local beam-gas interaction rate per unit length.

        Parameters
        ----------
        s : float
            Longitudinal position [m].

        Returns
        -------
        local_rate : float
            Interaction rate per metre of lattice [1/(s m)], for the whole
            bunch and restricted to the generated part of the cross section.
        """
        s_table = np.asarray(self.gas_density['s'], dtype=float)
        inverse_mfp = sum(
            np.interp(s, s_table,
                      np.asarray(self.gas_density[kk], dtype=float))
            * self.xsecs[kk]
            for kk in self.gas_species)

        return float(self.bunch_intensity*self._f_rev*inverse_mfp)

    @property
    def _f_rev(self):
        """
        Revolution frequency of the reference particle.

        Parameters
        ----------
        None

        Returns
        -------
        f_rev : float
            Revolution frequency [Hz].
        """
        return self.beta0*C_LIGHT/self.line.get_length()

    # ######################################################## #
    # Initialisation
    # ######################################################## #
    def initialise_beamgas(self, element=None, verbose=True):
        """
        Compute and configure the beam-gas interaction rates in the lattice.

        For each :class:`xcoll.BeamGasScattering` element this method:

        1. Integrates the atomic density of every gas species over the lattice
           section represented by the element (from the previous scattering
           element, or the start of the line, up to the element itself).
        2. Multiplies it by the total cross section of each species to obtain
           the interaction rate of the section, for the whole bunch and per
           second.
        3. Stores that rate together with the local optics, closed orbit, beam
           parameters and gas composition on the element via
           :meth:`xcoll.BeamGasScattering._configure`, so that
           :meth:`xcoll.BeamGasScattering.scatter` can weight the Monte Carlo
           macro-particles correctly.

        Parameters
        ----------
        element : str or None, optional
            If ``None`` (default), all the scattering elements of the study
            are initialised. If a string, only the named element is
            (re-)initialised.
        verbose : bool, optional
            If ``True`` (default), print one line per initialised element.

        Returns
        -------
        None
        """
        line = self.line
        tab = line.get_table()

        if self.twiss is None:
            twiss_method = self.kwargs.get('method', '6d')
            twiss_kwargs = {kk: vv for kk, vv in self.kwargs.items()
                            if kk != 'method'}
            self.twiss = line.twiss(method=twiss_method, **twiss_kwargs)

        if element is None:
            elements = self.elements
        else:
            if not isinstance(element, str):
                raise TypeError(
                    f"`element` must be a string "
                    f"(got {type(element).__name__}).")
            if element not in self.elements:
                raise ValueError(
                    f"`element='{element}'` is not among the elements of this "
                    f"BeamGasStudy.")
            elements = [element]

        s_all = np.array([float(tab['s', nn]) for nn in self.elements])
        line_length = float(line.get_length())
        if s_all[-1] < line_length - 1e-6:
            warn(f"The last BeamGasScattering is at s={s_all[-1]:.6g} m, but "
                 f"the line is {line_length:.6g} m long: the last "
                 f"{line_length - s_all[-1]:.6g} m of the line are not "
                 f"represented by any scattering element. Place a "
                 f"BeamGasScattering at the end of the line to cover the "
                 f"whole circumference.", stacklevel=2)

        s_table = np.asarray(self.gas_density['s'], dtype=float)
        density_table = {kk: np.asarray(self.gas_density[kk], dtype=float)
                         for kk in self.gas_species}
        f_rev = self._f_rev
        twiss = self.twiss

        for nn in elements:
            ii = self.elements.index(nn)
            s = s_all[ii]
            s_start = s_all[ii - 1] if ii > 0 else 0.0

            integrated_densities = self._integrated_atomic_densities(
                s_start, s)
            interaction_rate = self.bunch_intensity*f_rev*sum(
                integrated_densities[kk]*self.xsecs[kk]
                for kk in self.gas_species)

            atomic_densities = {
                kk: float(np.interp(s, s_table, density_table[kk]))
                for kk in self.gas_species}

            if verbose:
                print(f'Initialising BeamGasScattering for {nn}')

            line[nn]._configure(
                s=float(s),
                ds=float(s - s_start),
                particle_ref=self.particle_ref,
                element_index=line.element_names.index(nn),
                alfx=twiss['alfx', nn], betx=twiss['betx', nn],
                alfy=twiss['alfy', nn], bety=twiss['bety', nn],
                dx=twiss['dx', nn], dpx=twiss['dpx', nn],
                dy=twiss['dy', nn], dpy=twiss['dpy', nn],
                x_co=twiss['x', nn], px_co=twiss['px', nn],
                y_co=twiss['y', nn], py_co=twiss['py', nn],
                zeta_co=twiss['zeta', nn], delta_co=twiss['delta', nn],
                gemitt_x=self.gemitt_x, gemitt_y=self.gemitt_y,
                sigma_z=self.sigma_z, sigma_delta=self.sigma_delta,
                n_scattering_events=self.n_scattering_events,
                interaction_rate=float(interaction_rate),
                process=self.process,
                atomic_densities=atomic_densities,
                _calculators=self.calculators,
                _xsecs=self.xsecs,
            )

    # ######################################################## #
    # Running
    # ######################################################## #
    def generate_particles(self):
        """
        Generate beam-gas-scattered particles at the configured elements.

        Parameters
        ----------
        None

        Returns
        -------
        particles_by_element : dict
            Mapping ``{element_name: xtrack.Particles}``.
        """
        return {nn: self.line[nn].scatter(rng=self.rng)
                for nn in self.elements}

    def run(self, *, track=False, n_turns=None, generate_particles=None,
            keep_particles=False, with_progress=False):
        """
        Run the configured beam-gas rate or loss study.

        The scattering rate reported in the result is the beam-gas interaction
        rate configured on the scattering elements. If ``track`` is true,
        particles are generated and tracked from each scattering element back
        to itself, and the beam-gas loss rate is obtained as the total weight
        of the particles that are lost.

        Parameters
        ----------
        track : bool, optional
            If ``True``, track the generated particles and compute the loss
            rate from the particles lost during tracking. If ``False``, only
            the interaction rate is reported.
        n_turns : int or None, optional
            Number of turns to track. Required when ``track`` is ``True``.
        generate_particles : bool or None, optional
            If ``True``, generate weighted scattered particles even when
            tracking is disabled. If ``None``, particles are generated only
            when ``track`` is ``True`` or ``keep_particles`` is ``True``.
        keep_particles : bool, optional
            If ``True``, store the generated particle samples and the
            interaction log in the returned :class:`BeamGasResult`.
        with_progress : bool, optional
            Forwarded to :meth:`xtrack.Line.track` when tracking is enabled.

        Returns
        -------
        result : BeamGasResult
            Study result containing the local rates, the interaction rate, the
            optional tracking-derived loss rate and lifetime, and optionally
            the generated and lost particle samples.
        """
        if keep_particles and generate_particles is False:
            raise ValueError(
                "`keep_particles=True` is incompatible with "
                "`generate_particles=False`.")

        if track:
            if generate_particles is False:
                raise ValueError(
                    "`generate_particles=False` is incompatible with "
                    "`track=True`.")
            if n_turns is None:
                raise ValueError("`n_turns` is required when `track=True`.")
            generate_particles = True
        elif generate_particles is None:
            generate_particles = keep_particles

        particles_by_element = {}
        merged_particles = None
        lost_particles = None
        interaction_log = None

        rate_scattering = float(sum(self.line[nn].interaction_rate
                                    for nn in self.elements))

        if generate_particles:
            for nn in self.elements:
                particles = self.line[nn].scatter(rng=self.rng)
                if track:
                    self.line.track(particles,
                                    ele_start=nn, ele_stop=nn,
                                    num_turns=n_turns,
                                    with_progress=with_progress)
                particles_by_element[nn] = particles

            merged_particles = xt.Particles.merge(
                list(particles_by_element.values()))

        rate_tracking = None
        lifetime_tracking = None
        rate_tracking_error = None
        lifetime_tracking_error = None
        cutoff_scan = None
        rate_above_theta_max = None
        if track:
            lost = merged_particles.filter(_lost_mask(merged_particles))
            rate_tracking = float(np.sum(lost.weight))
            lifetime_tracking = (
                np.inf if rate_tracking == 0
                else float(self.bunch_intensity/rate_tracking))
            lost_particles = lost

            events = self._event_losses(particles_by_element)
            cutoff_scan = self._cutoff_scan(*events)
            rate_tracking_error = float(cutoff_scan.rate_tracking_error[0])
            lifetime_tracking_error = (
                np.nan if rate_tracking == 0
                else lifetime_tracking*rate_tracking_error/rate_tracking)
            if self.process == 'coulomb':
                rate_above_theta_max = self._rate_above_theta_max()
            self._warn_on_truncation(events[0], events[1]*events[2],
                                     cutoff_scan, rate_above_theta_max)

        lifetime_scattering = (
            np.inf if rate_scattering == 0
            else float(self.bunch_intensity/rate_scattering))

        local_rates = self.local_rates(
            particles_by_element=(
                particles_by_element if generate_particles else None),
            include_tracking=track)

        if keep_particles:
            interaction_log = self.interaction_log()
        else:
            particles_by_element = None
            merged_particles = None
            lost_particles = None

        return BeamGasResult(
            element_names=list(self.elements),
            gas_density=self.gas_density,
            local_rates=local_rates,
            rate_scattering=rate_scattering,
            lifetime_scattering=lifetime_scattering,
            rate_tracking=rate_tracking,
            lifetime_tracking=lifetime_tracking,
            tracked=track,
            rate_tracking_error=rate_tracking_error,
            lifetime_tracking_error=lifetime_tracking_error,
            cutoff_scan=cutoff_scan,
            rate_above_theta_max=rate_above_theta_max,
            particles_by_element=particles_by_element,
            particles=merged_particles,
            lost_particles=lost_particles,
            interaction_log=interaction_log,
        )

    # ######################################################## #
    # Diagnostics
    # ######################################################## #
    # Resolution of BeamGasResult.cutoff_scan
    _CUTOFF_SCAN_BINS_PER_DECADE = 10
    # Warn when more than this fraction of the loss rate comes from events
    # within a factor two of the lower cut
    _LOWER_CUT_WARN_FRACTION = 1e-3
    # Warn when the Coulomb interaction rate above theta_max exceeds this
    # fraction of the tracked loss rate
    _UPPER_CUT_WARN_FRACTION = 1e-2

    @property
    def _generation_window(self):
        """
        Range of the variable that bounds the generated events.

        Parameters
        ----------
        None

        Returns
        -------
        window : tuple of float
            ``(theta_min, theta_max)`` [rad] for Coulomb scattering, or
            ``(brems_energy_cut, ekin)`` [eV] for bremsstrahlung.
        """
        if self.process == 'coulomb':
            return self.coulomb_theta
        ekin = next(iter(self.calculators.values())).ekin
        return (self.brems_energy_cut, float(ekin))

    def _event_losses(self, particles_by_element):
        """
        Join the interaction log of every element with the tracked states.

        Xtrack reorders the particles while tracking, so the generated events
        are matched to the particles by ``particle_id``.

        Parameters
        ----------
        particles_by_element : dict
            Mapping ``{element_name: xtrack.Particles}`` of tracked particles.

        Returns
        -------
        x : ndarray
            Generation variable of each event: scattering angle [rad] for
            Coulomb scattering, photon energy [eV] for bremsstrahlung.
        weight : ndarray
            Weight of each event [1/s].
        lost : ndarray of bool
            Whether the primary particle of each event is lost.
        i_element : ndarray of int
            Index in :attr:`elements` of the element that generated each
            event.
        """
        column = 'theta' if self.process == 'coulomb' else 'photon_energy'
        x, weight, lost, i_element = [], [], [], []
        for ii, nn in enumerate(self.elements):
            log = self.line[nn].scatter_log
            particles = particles_by_element[nn]
            allocated = _allocated_mask(particles)
            ids = particles.particle_id[allocated]
            order = np.argsort(ids)
            idx = order[np.searchsorted(ids, log['particle_id'], sorter=order)]
            x.append(np.asarray(log[column], dtype=float))
            weight.append(np.asarray(log['weight'], dtype=float))
            lost.append(particles.state[allocated][idx] <= 0)
            i_element.append(np.full(idx.size, ii))

        return tuple(np.concatenate(vv) for vv in (x, weight, lost, i_element))

    def _cutoff_scan(self, x, weight, lost, i_element):
        """
        Tabulate the loss rate as a function of the lower generation cut.

        The events above a cut are an unbiased importance sample of a study
        generated with its lower cut set there, so the loss rate of every
        such study follows from the single tracked sample. The loss rate of
        each element is a sum of independent event contributions, so its
        variance is the number of events times their sample variance.

        Parameters
        ----------
        x, weight, lost, i_element : ndarray
            Per-event arrays returned by :meth:`_event_losses`.

        Returns
        -------
        table : xtrack.Table
            See :attr:`BeamGasResult.cutoff_scan`.
        """
        lo, hi = self._generation_window
        n_bins = max(int(np.ceil(
            self._CUTOFF_SCAN_BINS_PER_DECADE*np.log10(hi/lo))), 1)
        edges = np.geomspace(lo, hi, n_bins + 1)
        i_bin = np.clip(np.searchsorted(edges, x, side='right') - 1,
                        0, n_bins - 1)
        contribution = weight*lost
        n_elements = len(self.elements)

        def per_bin(values, element_wise=False):
            if not element_wise:
                return np.bincount(i_bin, weights=values, minlength=n_bins)
            flat = np.bincount(i_element*n_bins + i_bin, weights=values,
                               minlength=n_elements*n_bins)
            return flat.reshape(n_elements, n_bins)

        def above(binned):
            # Sum over all the bins from the current one upwards
            return np.flip(np.cumsum(np.flip(binned, axis=-1), axis=-1),
                           axis=-1)

        num_events = per_bin(np.ones_like(x))
        weight_in_bin = per_bin(weight)
        with np.errstate(invalid='ignore', divide='ignore'):
            loss_probability = per_bin(contribution)/weight_in_bin

        s1 = above(per_bin(contribution, element_wise=True))
        s2 = above(per_bin(contribution**2, element_wise=True))
        n = np.bincount(i_element, minlength=n_elements)[:, None]
        with np.errstate(invalid='ignore', divide='ignore'):
            variance = np.where(n > 1, n*(s2 - s1**2/n)/(n - 1), 0.0)
        rate = s1.sum(axis=0)
        with np.errstate(divide='ignore'):
            lifetime = np.where(rate > 0, self.bunch_intensity/rate, np.inf)

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

    def _rate_above_theta_max(self):
        """
        Compute the Coulomb interaction rate above ``theta_max``.

        Parameters
        ----------
        None

        Returns
        -------
        rate : float
            Interaction rate of the whole study for scattering angles between
            ``theta_max`` and ``pi`` [1/s].
        """
        theta_max = self.coulomb_theta[1]
        if theta_max >= np.pi:
            return 0.0

        xsecs_above = {
            kk: CoulombScatteringCalculator(
                    Z, self.p0c, q0=self.q0,
                    theta_lim=(theta_max, np.pi)).xsec
            for kk, Z in self.atomic_numbers.items()}

        rate = 0.0
        for nn in self.elements:
            elem = self.line[nn]
            integrated_densities = self._integrated_atomic_densities(
                float(elem.s) - float(elem.ds), float(elem.s))
            rate += sum(integrated_densities[kk]*xsecs_above[kk]
                        for kk in self.gas_species)

        return float(self.bunch_intensity*self._f_rev*rate)

    def _warn_on_truncation(self, x, contribution, cutoff_scan,
                            rate_above_theta_max):
        """
        Warn when the generation window visibly truncates the loss rate.

        Parameters
        ----------
        x : ndarray
            Generation variable of each event, see :meth:`_event_losses`.
        contribution : ndarray
            Loss-rate contribution of each event, i.e. its weight if lost and
            zero otherwise [1/s].
        cutoff_scan : xtrack.Table
            Table returned by :meth:`_cutoff_scan`.
        rate_above_theta_max : float or None
            Coulomb interaction rate above ``theta_max`` [1/s].

        Returns
        -------
        None
        """
        rate = float(np.sum(contribution))
        if rate <= 0:
            return

        lo = self._generation_window[0]
        if self.process == 'coulomb':
            name = '`coulomb_theta[0]`'
            what = 'scattering angles'
            hint = (" The log-uniform angular sampling makes a lower cut "
                    "cheap: the statistics only degrade as "
                    "log(theta_max/theta_min).")
        else:
            name = '`brems_energy_cut`'
            what = 'photon energies'
            hint = ""

        near_lower = float(np.sum(contribution[x < 2*lo]))/rate
        if near_lower > self._LOWER_CUT_WARN_FRACTION:
            warn(f"{100*near_lower:.2g}% of the tracked loss rate comes from "
                 f"events generated within a factor 2 of the lower cut "
                 f"{name}={lo:.3g}. The cut lies inside the range of {what} "
                 f"that cause losses, so the losses from below it are missing "
                 f"and the lifetime is overestimated. Lower {name}, and check "
                 f"that `BeamGasResult.cutoff_scan.lifetime_tracking` is flat "
                 f"in its first rows.{hint}", stacklevel=3)

        if rate_above_theta_max is not None:
            above_upper = rate_above_theta_max/rate
            if above_upper > self._UPPER_CUT_WARN_FRACTION:
                warn(f"The Coulomb interaction rate above `coulomb_theta[1]`="
                     f"{self.coulomb_theta[1]:.3g} rad, which is not "
                     f"generated, amounts to {100*above_upper:.2g}% of the "
                     f"tracked loss rate. Up to that much loss rate is "
                     f"missing, and all of it if these events are always lost "
                     f"(loss probability in the last `cutoff_scan` row: "
                     f"{cutoff_scan.loss_probability[-1]:.2f}). Raise "
                     f"`coulomb_theta[1]`, or add `rate_above_theta_max` to "
                     f"the loss rate.", stacklevel=3)

    def local_rates(self, *, particles_by_element=None,
                    include_tracking=False):
        """
        Return an ``xt.Table`` with per-scattering-element diagnostics.

        Parameters
        ----------
        particles_by_element : dict or None, optional
            Mapping from element name to the particles generated at that
            element. If provided, Monte Carlo particle counts and weight sums
            are included in the table.
        include_tracking : bool, optional
            If ``True``, include loss-count and lost-weight columns computed
            from the tracked particle states. This option is meaningful only
            when ``particles_by_element`` is provided.

        Returns
        -------
        table : xtrack.Table
            Per-element diagnostics table. It always contains ``name``, ``s``,
            ``ds``, ``mean_free_path``, ``local_rate`` and
            ``interaction_rate``. If particles are provided, it also contains
            ``num_particles`` and ``sum_weight``. If ``include_tracking`` is
            ``True``, it additionally contains ``num_lost_particles`` and
            ``sum_lost_weight``.
        """
        data = {'name': [], 's': [], 'ds': [], 'mean_free_path': [],
                'local_rate': [], 'interaction_rate': []}
        include_particles = particles_by_element is not None
        if include_particles:
            data.update({'num_particles': [], 'sum_weight': []})
        if include_tracking:
            data.update({'num_lost_particles': [], 'sum_lost_weight': []})

        for nn in self.elements:
            elem = self.line[nn]
            data['name'].append(nn)
            data['s'].append(float(elem.s))
            data['ds'].append(float(elem.ds))
            data['mean_free_path'].append(float(elem.mean_free_path))
            data['local_rate'].append(self._local_rate(float(elem.s)))
            data['interaction_rate'].append(float(elem.interaction_rate))

            particles = (None if particles_by_element is None
                         else particles_by_element.get(nn))

            if include_particles:
                if particles is None:
                    data['num_particles'].append(0)
                    data['sum_weight'].append(np.nan)
                else:
                    allocated = _allocated_mask(particles)
                    data['num_particles'].append(int(np.sum(allocated)))
                    data['sum_weight'].append(
                        float(np.sum(particles.weight[allocated])))

            if include_tracking:
                if particles is None:
                    data['num_lost_particles'].append(0)
                    data['sum_lost_weight'].append(np.nan)
                else:
                    lost = _lost_mask(particles)
                    data['num_lost_particles'].append(int(np.sum(lost)))
                    data['sum_lost_weight'].append(
                        float(np.sum(particles.weight[lost])))

        return xt.Table({kk: np.array(vv) for kk, vv in data.items()})

    def interaction_log(self):
        """
        Return an ``xt.Table`` with the record of the generated interactions.

        Parameters
        ----------
        None

        Returns
        -------
        table : xtrack.Table or None
            Table with columns ``name``, ``s``, ``particle_id``, ``gas``,
            ``theta``, ``photon_energy`` and ``weight``, with one row per
            generated interaction. ``None`` if no particles have been
            generated yet.

        Notes
        -----
        ``particle_id`` identifies the particle within the sample generated at
        the element named in the same row. Particle ids travel with the
        particles during tracking, so the log can always be joined to
        ``BeamGasResult.particles_by_element[name]``. It cannot be joined to
        the merged ``BeamGasResult.particles``, because
        :meth:`xtrack.Particles.merge` renumbers the particles.
        """
        logs = [(nn, self.line[nn].scatter_log) for nn in self.elements
                if self.line[nn].scatter_log is not None]
        if len(logs) == 0:
            return None

        data = {'name': [], 's': [], 'particle_id': [], 'gas': [],
                'theta': [], 'photon_energy': [], 'weight': []}
        for nn, log in logs:
            n = len(log['particle_id'])
            data['name'].append(np.full(n, nn))
            data['s'].append(np.full(n, float(self.line[nn].s)))
            for cc in ('particle_id', 'gas', 'theta', 'photon_energy',
                       'weight'):
                data[cc].append(np.asarray(log[cc]))

        return xt.Table({kk: np.concatenate(vv) for kk, vv in data.items()})
