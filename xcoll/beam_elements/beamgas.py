# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2026.                 #
# ######################################### #

import numpy as np

import xobjects as xo
import xtrack as xt
from xtrack.particles import pdg

from ..general import _pkg_root


BEAMGAS_PROCESSES = ('brems', 'coulomb')

# Number of gas species that a marker can store
MAX_GAS_SPECIES = 16


def _resolve_process(process):
    """
    Validate and normalise the name of a beam-gas process.

    Parameters
    ----------
    process : str
        Name of the process, ``'brems'`` or ``'coulomb'``.

    Returns
    -------
    process : str
        Normalised process name.
    """
    if isinstance(process, str):
        process = process.lower()
        if process in BEAMGAS_PROCESSES:
            return process
        if process in ('bremsstrahlung', 'bremss'):
            return 'brems'

    raise ValueError(
        "`process` must be either 'brems' (bremsstrahlung) or 'coulomb' "
        f"(Coulomb scattering). Got {process!r} instead.")


def _resolve_n_scattering_events(*, n_scattering_events, default):
    """
    Validate the requested number of Monte Carlo scattering events.

    Parameters
    ----------
    n_scattering_events : int or None
        Requested number of events.
    default : int or None
        Value used when ``n_scattering_events`` is ``None``.

    Returns
    -------
    n_scattering_events : int
        Validated number of events.
    """
    if n_scattering_events is None:
        n_scattering_events = default

    if n_scattering_events is None:
        raise ValueError("`n_scattering_events` is required.")

    if n_scattering_events < 0:
        raise ValueError("`n_scattering_events` must be non-negative.")

    return int(n_scattering_events)


class BeamGasMarker(xt.BeamElement):
    """
    Passive marker of a beam-gas scattering location.

    Each marker represents the section of lattice that precedes it, up to the
    previous marker (or the start of the line for the first one). It holds the
    configuration needed to generate beam-gas interactions there: the local
    optics and closed orbit, the beam parameters, the gas composition, the
    process and its generation window, and the interaction rate of the
    section. All of it is stored in xofields, so a configured line can be
    saved to JSON and loaded back, e.g. by independent jobs.

    The marker does nothing during tracking. The markers are installed,
    configured and used through :class:`xcoll.BeamGasGenerator`, available as
    ``line.xcoll.beamgas``.

    Parameters
    ----------
    ds : float, optional
        Length of the lattice section represented by this marker [m].
    p0c : float, optional
        Reference momentum at configuration [eV].
    bunch_intensity : float, optional
        Number of particles in the bunch.
    alfx, betx : float, optional
        Horizontal Twiss parameters at the marker.
    alfy, bety : float, optional
        Vertical Twiss parameters at the marker.
    dx, dpx : float, optional
        Horizontal dispersion and its derivative at the marker.
    dy, dpy : float, optional
        Vertical dispersion and its derivative at the marker.
    x_co, px_co : float, optional
        Horizontal closed-orbit position and normalised momentum.
    y_co, py_co : float, optional
        Vertical closed-orbit position and normalised momentum.
    zeta_co, delta_co : float, optional
        Longitudinal closed-orbit coordinate and relative momentum deviation.
    gemitt_x, gemitt_y : float, optional
        Horizontal and vertical geometric emittances [m rad].
    sigma_z : float, optional
        RMS bunch length [m].
    sigma_delta : float, optional
        RMS relative momentum spread.
    n_scattering_events : int, optional
        Number of beam-gas interactions generated at this marker.
    interaction_rate : float, optional
        Beam-gas interaction rate of the represented section [1/s], restricted
        to the generated part of the cross section. It is the total weight of
        the particles generated at this marker.
    rate_above_theta_max : float, optional
        Coulomb-scattering rate of the represented section above
        ``coulomb_theta_max`` [1/s], which is not generated.
    mean_free_path : float, optional
        Local beam-gas mean free path [m], restricted to the generated part of
        the cross section.
    process : {'brems', 'coulomb'}, optional
        Beam-gas process generated at this marker.
    brems_energy_cut : float, optional
        Lowest photon energy generated for bremsstrahlung [eV].
    coulomb_theta_min, coulomb_theta_max : float, optional
        Range of scattering angles generated for Coulomb scattering [rad].
    atomic_densities : dict, optional
        Mapping ``{element_symbol: atomic_density}`` of the local gas
        composition [atoms/m^3], with at most ``MAX_GAS_SPECIES`` species.
    seed : int, optional
        Seed of the random streams of the marker.

    Attributes
    ----------
    process, atomic_densities, gas_species
        See the properties.
    scatter_log : xtrack.Table or None
        Record of the interactions of the last generation at this marker; see
        :meth:`xcoll.BeamGasGenerator.generate`. It is kept in memory only.
        ``None`` before the first generation.

    Notes
    -----
    **Physics summary**

    For each of the ``n_scattering_events`` generated macro-particles:

    1. Initial coordinates are drawn from the local matched Gaussian
       distribution, built from the Twiss parameters, the geometric
       emittances, ``sigma_z`` and ``sigma_delta``, and shifted onto the local
       closed orbit. Betatron and dispersive contributions are included;
       transverse coupling is not.
    2. A gas species is drawn with a probability proportional to its
       contribution ``n_gas * sigma_gas`` to the local interaction rate.
    3. One interaction of the configured ``process`` is applied, using the
       cross-section models of :mod:`xcoll.beamgas.cross_sections`.
    4. The macro-particle is given a weight
       ``interaction_rate/n_scattering_events`` times the importance weight of
       the sampled event, so that the total weight of the sample is the
       section interaction rate [1/s].

    Every generated macro-particle interacts exactly once. This *forced
    interaction* scheme is the correct estimator in the single-scattering
    regime of a dilute residual gas, where the real probability of
    interacting in one turn is many orders of magnitude below one.
    """

    _xofields = {
        'p0c': xo.Float64,
        'bunch_intensity': xo.Float64,
        'ds': xo.Float64,
        'interaction_rate': xo.Float64,
        'rate_above_theta_max': xo.Float64,
        'mean_free_path': xo.Float64,
        'n_scattering_events': xo.Int64,
        'gemitt_x': xo.Float64,
        'gemitt_y': xo.Float64,
        'alfx': xo.Float64,
        'betx': xo.Float64,
        'alfy': xo.Float64,
        'bety': xo.Float64,
        'dx': xo.Float64,
        'dpx': xo.Float64,
        'dy': xo.Float64,
        'dpy': xo.Float64,
        'x_co': xo.Float64,
        'px_co': xo.Float64,
        'y_co': xo.Float64,
        'py_co': xo.Float64,
        'zeta_co': xo.Float64,
        'delta_co': xo.Float64,
        'sigma_z': xo.Float64,
        'sigma_delta': xo.Float64,
        '_process': xo.Int64,
        'brems_energy_cut': xo.Float64,
        'coulomb_theta_min': xo.Float64,
        'coulomb_theta_max': xo.Float64,
        'n_gas_species': xo.Int64,
        'gas_atomic_numbers': xo.Int64[MAX_GAS_SPECIES],
        'gas_densities': xo.Float64[MAX_GAS_SPECIES],
        'seed': xo.Int64,
    }

    isthick = False
    behaves_like_drift = False
    allow_track = True
    skip_in_loss_location_refinement = True

    _extra_c_sources = [
        _pkg_root.joinpath('beam_elements', 'elements_src', 'beamgas.h')
    ]

    def __init__(self, ds=0.0,
                 p0c=0.0,
                 bunch_intensity=0.0,
                 alfx=0.0, betx=0.0, alfy=0.0, bety=0.0,
                 dx=0.0, dpx=0.0, dy=0.0, dpy=0.0,
                 x_co=0.0, px_co=0.0, y_co=0.0, py_co=0.0,
                 zeta_co=0.0, delta_co=0.0,
                 gemitt_x=0.0, gemitt_y=0.0,
                 sigma_z=0.0, sigma_delta=0.0,
                 n_scattering_events=None,
                 interaction_rate=0.0,
                 rate_above_theta_max=0.0,
                 mean_free_path=0.0,
                 process=None,
                 brems_energy_cut=0.0,
                 coulomb_theta_min=0.0,
                 coulomb_theta_max=0.0,
                 atomic_densities=None,
                 seed=0,
                 **kwargs):
        """
        Create a beam-gas scattering marker.

        Parameters
        ----------
        See the class docstring.

        Returns
        -------
        None
        """
        self.scatter_log = None

        if '_xobject' in kwargs.keys():
            self.xoinitialize(**kwargs)
            return

        super().__init__(**kwargs)

        self.ds = ds
        self.p0c = p0c
        self.bunch_intensity = bunch_intensity
        self.alfx = alfx
        self.betx = betx
        self.alfy = alfy
        self.bety = bety
        self.dx = dx
        self.dpx = dpx
        self.dy = dy
        self.dpy = dpy
        self.x_co = x_co
        self.px_co = px_co
        self.y_co = y_co
        self.py_co = py_co
        self.zeta_co = zeta_co
        self.delta_co = delta_co
        self.gemitt_x = gemitt_x
        self.gemitt_y = gemitt_y
        self.sigma_z = sigma_z
        self.sigma_delta = sigma_delta
        self.n_scattering_events = _resolve_n_scattering_events(
            n_scattering_events=n_scattering_events, default=0)
        self.interaction_rate = interaction_rate
        self.rate_above_theta_max = rate_above_theta_max
        self.mean_free_path = mean_free_path
        # Not overwritten when None, to keep the values loaded from a dict
        if process is not None:
            self.process = process
        self.brems_energy_cut = brems_energy_cut
        self.coulomb_theta_min = coulomb_theta_min
        self.coulomb_theta_max = coulomb_theta_max
        if atomic_densities is not None:
            self.atomic_densities = atomic_densities
        self.seed = seed

    @property
    def process(self):
        """
        Beam-gas process generated by this element.

        Parameters
        ----------
        None

        Returns
        -------
        process : {'brems', 'coulomb'}
            Name of the configured process.
        """
        return BEAMGAS_PROCESSES[int(self._process)]

    @process.setter
    def process(self, value):
        """
        Set the beam-gas process generated by this element.

        Parameters
        ----------
        value : {'brems', 'coulomb'}
            Name of the process.

        Returns
        -------
        None
        """
        self._process = BEAMGAS_PROCESSES.index(_resolve_process(value))

    @property
    def atomic_densities(self):
        """
        Local gas composition.

        Parameters
        ----------
        None

        Returns
        -------
        atomic_densities : dict
            Mapping ``{element_symbol: atomic_density}`` [atoms/m^3].
        """
        n_species = int(self.n_gas_species)
        return {pdg.get_element_name_from_Z(int(zz)): float(dd)
                for zz, dd in zip(self.gas_atomic_numbers[:n_species],
                                  self.gas_densities[:n_species])}

    @atomic_densities.setter
    def atomic_densities(self, value):
        """
        Set the local gas composition.

        Parameters
        ----------
        value : dict
            Mapping ``{element_symbol: atomic_density}`` [atoms/m^3], with at
            most ``MAX_GAS_SPECIES`` species.

        Returns
        -------
        None
        """
        value = dict(value)
        if len(value) > MAX_GAS_SPECIES:
            raise ValueError(
                f"At most {MAX_GAS_SPECIES} gas species are supported "
                f"(got {len(value)}).")
        atomic_numbers = np.zeros(MAX_GAS_SPECIES, dtype=np.int64)
        densities = np.zeros(MAX_GAS_SPECIES)
        for ii, (kk, vv) in enumerate(value.items()):
            atomic_numbers[ii] = int(pdg.get_Z_from_element_name(kk))
            densities[ii] = float(vv)
        self.n_gas_species = len(value)
        self.gas_atomic_numbers[:] = atomic_numbers
        self.gas_densities[:] = densities

    @property
    def gas_species(self):
        """
        Gas species present at this element.

        Parameters
        ----------
        None

        Returns
        -------
        species : list of str
            Chemical symbols of the local gas species, in a fixed order.
        """
        return list(self.atomic_densities.keys())

    def _configure(self, **kwargs):
        """
        Set the configuration supplied by
        :meth:`xcoll.BeamGasGenerator.configure`.

        Parameters
        ----------
        **kwargs
            Any of the constructor parameters.

        Returns
        -------
        None
        """
        config_allowed = {
            'ds', 'p0c', 'bunch_intensity',
            'alfx', 'betx', 'alfy', 'bety',
            'dx', 'dpx', 'dy', 'dpy',
            'x_co', 'px_co', 'y_co', 'py_co',
            'zeta_co', 'delta_co',
            'gemitt_x', 'gemitt_y',
            'sigma_z', 'sigma_delta',
            'n_scattering_events', 'interaction_rate',
            'rate_above_theta_max', 'mean_free_path',
            'process', 'brems_energy_cut', 'coulomb_theta_min',
            'coulomb_theta_max', 'atomic_densities', 'seed',
        }

        unknown = set(kwargs) - config_allowed
        if unknown:
            bad = ", ".join(sorted(unknown))
            raise KeyError(f"Unsupported configure() keys: {bad}")

        n_scattering_events = kwargs.pop('n_scattering_events', None)
        if n_scattering_events is not None:
            self.n_scattering_events = _resolve_n_scattering_events(
                n_scattering_events=n_scattering_events,
                default=self.n_scattering_events)

        for kk, vv in kwargs.items():
            setattr(self, kk, vv)

    def _generate_primaries(self, n, rng):
        """
        Draw macro-particles from the local matched Gaussian distribution.

        Betatron and dispersive contributions are generated relative to the
        closed orbit, which is then added so that the returned coordinates are
        the physical ones. The interaction is applied afterwards, on the
        physical momentum direction of each macro-particle.

        Parameters
        ----------
        n : int
            Number of macro-particles to generate.
        rng : numpy.random.Generator
            Random number generator.

        Returns
        -------
        x, px, y, py, zeta, delta : ndarray
            Physical coordinates of the generated macro-particles.
        """
        u = rng.standard_normal((6, n))

        delta = self.sigma_delta*u[5]
        zeta = self.sigma_z*u[4]

        x = np.sqrt(self.gemitt_x*self.betx)*u[0] + self.dx*delta
        px = (np.sqrt(self.gemitt_x/self.betx)*(u[1] - self.alfx*u[0])
              + self.dpx*delta)
        y = np.sqrt(self.gemitt_y*self.bety)*u[2] + self.dy*delta
        py = (np.sqrt(self.gemitt_y/self.bety)*(u[3] - self.alfy*u[2])
              + self.dpy*delta)

        return (x + self.x_co, px + self.px_co,
                y + self.y_co, py + self.py_co,
                zeta + self.zeta_co, delta + self.delta_co)

    def track(self, particles):
        """
        Track particles through the element without applying a kick.

        Beam-gas interactions are generated by
        :meth:`xcoll.BeamGasGenerator.generate`; during tracking this element
        is a passive marker.

        Parameters
        ----------
        particles : xtrack.Particles
            Particles tracked through the passive marker element.

        Returns
        -------
        None
        """
        super().track(particles)
