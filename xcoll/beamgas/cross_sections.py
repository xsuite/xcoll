# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2026.                 #
# ######################################### #
"""
Atomic data and single-scattering models for beam-residual-gas interactions.

This module provides the differential and total cross sections, and the
associated Monte Carlo samplers, used by :class:`xcoll.BeamGasStudy` to
simulate the interaction of a relativistic electron or positron beam with the
residual gas of an accelerator vacuum chamber.

Two processes are available:

* **Bremsstrahlung** (``'brems'``): emission of a hard photon in the field of a
  gas nucleus. The energy loss is the dominant effect; the associated angular
  deflection is small but is included.
* **Coulomb scattering** (``'coulomb'``): elastic scattering off the screened
  nuclear Coulomb field. The energy loss is negligible; the angular deflection
  is the dominant effect.

Acknowledgements
----------------
The models implemented here are transcriptions of the corresponding Geant4
models, restricted to the relativistic electron/positron-on-dilute-gas regime
relevant for storage rings. The provenance of each piece is given in the
individual class and method docstrings. In particular:

* The bremsstrahlung differential cross section, the Tsai screening functions,
  the Coulomb correction, the low-Z Dirac-Fock radiation logarithms and the
  Gauss-Legendre integration of the total cross section follow
  ``G4eBremsstrahlungRelModel`` and ``G4SeltzerBergerModel``.
* The bremsstrahlung photon angular distribution follows ``G4ModifiedTsai``.
* The screened Rutherford / Mott elastic cross section, the Moliere screening
  parameter, the McKinley-Feshbach Mott-to-Rutherford ratio and the
  exponential nuclear form factor follow ``G4ScreeningMottCrossSection`` and
  ``G4WentzelVIModel``.

References
----------
.. [1] S. Agostinelli et al. (Geant4 Collaboration), "Geant4 - a simulation
   toolkit", Nucl. Instrum. Meth. A **506**, 250 (2003).
   https://doi.org/10.1016/S0168-9002(03)01368-8
.. [2] Geant4 Collaboration, "Geant4 Physics Reference Manual", Rev. 11.x.
   https://geant4.web.cern.ch/docs/
.. [3] Y.-S. Tsai, "Pair production and bremsstrahlung of charged leptons",
   Rev. Mod. Phys. **46**, 815 (1974).
   https://doi.org/10.1103/RevModPhys.46.815
.. [4] S. M. Seltzer and M. J. Berger, "Bremsstrahlung spectra from electron
   interactions with screened atomic nuclei and orbital electrons",
   Nucl. Instrum. Meth. B **12**, 95 (1985).
   https://doi.org/10.1016/0168-583X(85)90707-4
.. [5] G. Moliere, "Theorie der Streuung schneller geladener Teilchen I.
   Einzelstreuung am abgeschirmten Coulomb-Feld", Z. Naturforsch. A **2**,
   133 (1947). https://doi.org/10.1515/zna-1947-0302
.. [6] W. A. McKinley and H. Feshbach, "The Coulomb scattering of relativistic
   electrons by nuclei", Phys. Rev. **74**, 1759 (1948).
   https://doi.org/10.1103/PhysRev.74.1759
.. [7] M. J. Boschini et al., "Nuclear and non-ionizing energy-loss for
   Coulomb scattered particles from low energy up to relativistic regime in
   space radiation environment", in Proc. ICATPP 2010, p. 9 (2011).
   https://doi.org/10.1142/9789814329033_0002
"""

from dataclasses import dataclass

import numpy as np
from scipy import constants as sc

import xtrack as xt
from xtrack.particles import pdg


# ############################################################ #
# Constants
# ############################################################ #
ELECTRON_MASS_EV = xt.ELECTRON_MASS_EV

ALPHA = sc.alpha
# Reduced Planck constant times the speed of light, in eV m (~1.9733e-7 eV m)
HBAR_C_EV_M = sc.hbar * sc.c / sc.e
CLASSICAL_ELECTRON_RADIUS = sc.value('classical electron radius')
BOHR_RADIUS = sc.value('Bohr radius')
# Atomic mass constant (unified atomic mass unit) in eV
ATOMIC_MASS_CONSTANT_EV = sc.m_u * sc.c**2 / sc.e

# Thomas-Fermi constant, a_TF = C_TF * a_0 * Z^(-1/3)
C_TF = 0.5 * (3.0*np.pi/4.0)**(2.0/3.0)

# Standard atomic weights (IUPAC 2021), for Z = 1 to 54. This range covers all
# species found in the residual gas of an accelerator vacuum chamber
# (H2, He, CH4, H2O, CO, N2, O2, Ne, Ar, CO2, Kr, Xe, ...). Keeping the table
# local avoids adding a third-party periodic-table dependency to Xcoll.
_STANDARD_ATOMIC_WEIGHT = {
     1: 1.008,    2: 4.0026,   3: 6.94,     4: 9.0122,   5: 10.81,
     6: 12.011,   7: 14.007,   8: 15.999,   9: 18.998,  10: 20.180,
    11: 22.990,  12: 24.305,  13: 26.982,  14: 28.085,  15: 30.974,
    16: 32.06,   17: 35.45,   18: 39.95,   19: 39.098,  20: 40.078,
    21: 44.956,  22: 47.867,  23: 50.942,  24: 51.996,  25: 54.938,
    26: 55.845,  27: 58.933,  28: 58.693,  29: 63.546,  30: 65.38,
    31: 69.723,  32: 72.630,  33: 74.922,  34: 78.971,  35: 79.904,
    36: 83.798,  37: 85.468,  38: 87.62,   39: 88.906,  40: 91.224,
    41: 92.906,  42: 95.95,   43: 97.0,    44: 101.07,  45: 102.91,
    46: 106.42,  47: 107.87,  48: 112.41,  49: 114.82,  50: 118.71,
    51: 121.76,  52: 127.60,  53: 126.90,  54: 131.29,
}
_MAX_Z = max(_STANDARD_ATOMIC_WEIGHT)

# Maximum number of rejection-sampling iterations before giving up. Reaching
# this limit means the majorant of the sampled distribution is wrong.
_MAX_REJECTION_ITERATIONS = 1000


def atomic_number_from_symbol(symbol):
    """
    Return the atomic number of a chemical element from its symbol.

    Parameters
    ----------
    symbol : str
        Short (``'Ar'``) or long (``'Argon'``) element name.

    Returns
    -------
    Z : int
        Atomic number of the element.
    """
    return int(pdg.get_Z_from_element_name(symbol))


# ############################################################ #
# Scattering sample
# ############################################################ #
@dataclass
class ScatteringSample:
    """
    Outcome of a batch of single-scattering events.

    All momentum components are normalised to the reference momentum ``p0c``,
    consistently with the Xsuite convention (``px = Px/p0c``,
    ``delta = (p - p0)/p0``).

    Parameters
    ----------
    px, py : ndarray
        Normalised transverse momenta after the interaction.
    delta : ndarray
        Relative momentum deviation after the interaction.
    weight : ndarray
        Importance-sampling weight of each event, normalised to the total
        cross section returned by ``compute_xsec`` so that its expectation
        value is exactly one. For an unbiased sampler every weight is one.
    theta : ndarray
        Polar angle of the interaction [rad]. For Coulomb scattering this is
        the scattering angle of the beam particle; for bremsstrahlung it is
        the emission angle of the photon.
    photon_energy : ndarray
        Energy of the radiated photon [eV]. ``NaN`` for Coulomb scattering.

    Attributes
    ----------
    px, py, delta, weight, theta, photon_energy : ndarray
        See above.
    """
    px: np.ndarray
    py: np.ndarray
    delta: np.ndarray
    weight: np.ndarray
    theta: np.ndarray
    photon_energy: np.ndarray


# ############################################################ #
# Helpers
# ############################################################ #
def _energy_from_momentum(pc):
    """
    Return the total energy of an electron/positron from its momentum.

    Parameters
    ----------
    pc : float or ndarray
        Momentum times the speed of light [eV].

    Returns
    -------
    energy : float or ndarray
        Total energy [eV].
    """
    return np.sqrt(pc**2 + ELECTRON_MASS_EV**2)


def _momentum_frame(px, py, delta):
    """
    Build a right-handed orthonormal frame attached to the particle momentum.

    The returned basis ``(u_hat, v_hat, p_hat)`` has ``p_hat`` along the
    particle momentum, so that a direction given by a polar angle ``theta``
    around ``p_hat`` and an azimuth ``phi`` reads
    ``sin(theta) cos(phi) u_hat + sin(theta) sin(phi) v_hat
    + cos(theta) p_hat``.

    The azimuthal orientation of ``u_hat`` in the plane transverse to
    ``p_hat`` is arbitrary; this is irrelevant because ``phi`` is always
    sampled uniformly in ``[0, 2 pi)``.

    Parameters
    ----------
    px, py : ndarray
        Normalised transverse momenta.
    delta : ndarray
        Relative momentum deviation.

    Returns
    -------
    u_hat, v_hat, p_hat : ndarray
        Basis vectors, each with shape ``(n, 3)``.
    p_norm : ndarray
        Momentum modulus normalised to ``p0c``, i.e. ``1 + delta``.
    """
    pz = np.sqrt(np.maximum((1.0 + delta)**2 - px**2 - py**2, 0.0))
    p_vec = np.stack((px, py, pz), axis=1)
    p_norm = np.linalg.norm(p_vec, axis=1)
    p_hat = p_vec / p_norm[:, None]

    # u_hat = z_hat x p_hat, normalised; degenerate when p_hat is along z_hat
    u_hat = np.stack((-p_hat[:, 1], p_hat[:, 0], np.zeros_like(px)), axis=1)
    u_norm = np.linalg.norm(u_hat, axis=1)
    degenerate = u_norm < 1e-12
    u_hat[degenerate] = (1.0, 0.0, 0.0)
    u_norm[degenerate] = 1.0
    u_hat /= u_norm[:, None]

    v_hat = np.cross(p_hat, u_hat)

    return u_hat, v_hat, p_hat, p_norm


def _rejection_sample(n, propose, rng, name):
    """
    Draw ``n`` accepted samples with vectorised rejection sampling.

    The proposal callable is invoked with a batch size and must return the
    array of accepted values for that batch. Batches are over-drawn based on
    the running acceptance rate, so that typically one or two calls suffice.

    Parameters
    ----------
    n : int
        Number of accepted samples required.
    propose : callable
        ``propose(batch_size, rng) -> ndarray`` returning accepted values.
    rng : numpy.random.Generator
        Random number generator.
    name : str
        Name of the sampled quantity, used in the error message.

    Returns
    -------
    samples : ndarray
        Array of ``n`` accepted samples.
    """
    if n == 0:
        return np.zeros(0)

    out = np.empty(n)
    n_filled = 0
    n_drawn = 0
    for _ in range(_MAX_REJECTION_ITERATIONS):
        n_missing = n - n_filled
        # Over-draw using the measured acceptance rate (guessed on first pass)
        efficiency = max(n_filled/n_drawn, 1e-3) if n_drawn > 0 else 0.1
        batch_size = int(min(n_missing/efficiency*1.1 + 16, 1e7))
        accepted = propose(batch_size, rng)
        n_drawn += batch_size
        n_take = min(accepted.size, n_missing)
        out[n_filled:n_filled + n_take] = accepted[:n_take]
        n_filled += n_take
        if n_filled == n:
            return out

    raise RuntimeError(
        f"Rejection sampling of `{name}` did not converge after "
        f"{_MAX_REJECTION_ITERATIONS} iterations.")


# ############################################################ #
# Atomic data
# ############################################################ #
class ElementData:
    """
    Atomic data of a gas species needed by the beam-gas cross sections.

    The radiation logarithms, the Coulomb correction and the derived ``f_Z``
    factors follow the Geant4 relativistic bremsstrahlung model
    (``G4eBremsstrahlungRelModel``), which in turn follows Tsai's complete
    screening results with the Dirac-Fock corrections of Seltzer and Berger
    for low-Z elements.

    Parameters
    ----------
    Z : int
        Atomic number of the element. Must be between 1 and 54.

    Attributes
    ----------
    Z : int
        Atomic number.
    symbol : str
        Chemical symbol.
    atomic_weight : float
        Standard atomic weight (IUPAC 2021), in unified atomic mass units.
    A : int
        Mass number, taken as the rounded standard atomic weight.
    mass : float
        Atomic mass [eV].
    nuclear_radius : float
        Effective nuclear radius ``1.27 fm * A^0.27`` [m], as used in the
        exponential nuclear form factor of ``G4ScreeningMottCrossSection``.
    f_el, f_inel : float
        Elastic and inelastic radiation logarithms.
    f_c : float
        Coulomb correction factor of Davies, Bethe and Maximon.
    f_Z_factor_1, f_Z_factor_2 : float
        Combinations of the above entering the low-Z bremsstrahlung
        differential cross section; their sum is also the majorant used when
        sampling the photon energy.
    f_Z : float
        ``log(Z)/3 + f_c``, entering the high-Z differential cross section.
    gamma_factor, epsilon_factor : float
        Prefactors of the Tsai screening variables ``gamma`` and ``epsilon``
        [eV], to be multiplied by ``y/(etot - photon_energy)`` [1/eV].
    """

    def __init__(self, Z):
        Z = int(Z)
        if Z < 1 or Z > _MAX_Z:
            raise ValueError(
                f"Atomic number Z={Z} is out of the supported range "
                f"[1, {_MAX_Z}].")

        self.Z = Z
        self.symbol = pdg.get_element_name_from_Z(Z)
        self.atomic_weight = _STANDARD_ATOMIC_WEIGHT[Z]
        self.A = int(round(self.atomic_weight))
        self.mass = self.atomic_weight * ATOMIC_MASS_CONSTANT_EV
        self.nuclear_radius = 1.27e-15 * self.A**0.27

        self.f_el, self.f_inel = self._compute_radiation_logarithms()
        self.f_c = self._compute_coulomb_factor()
        (self.f_Z_factor_1, self.f_Z_factor_2,
         self.f_Z) = self._compute_Z_factors()
        self.gamma_factor, self.epsilon_factor = \
            self._compute_gamma_epsilon_factors()

    def __repr__(self):
        return (f"<ElementData {self.symbol} (Z={self.Z}, A={self.A})>")

    def _compute_coulomb_factor(self):
        """
        Compute the Coulomb correction factor of Davies, Bethe and Maximon.

        Parameters
        ----------
        None

        Returns
        -------
        f_c : float
            Coulomb correction factor.
        """
        K1, K2, K3, K4 = 0.0083, 0.20206, 0.0020, 0.0369
        a_Z = ALPHA * self.Z

        return ((K1*a_Z**4 + K2 + 1/(1 + a_Z**2)) * a_Z**2
                - (K3*a_Z**4 + K4) * a_Z**4)

    def _compute_radiation_logarithms(self):
        """
        Compute the elastic and inelastic radiation logarithms.

        For ``Z < 5`` the Thomas-Fermi model is inaccurate, and the tabulated
        values obtained with the Dirac-Fock atomic model are used instead.

        Parameters
        ----------
        None

        Returns
        -------
        f_el : float
            Elastic radiation logarithm.
        f_inel : float
            Inelastic radiation logarithm.
        """
        F_EL_LOWZ = [0.0, 5.3104, 4.7935, 4.7402, 4.7112, 4.6694, 4.6134,
                     4.5520]
        F_INEL_LOWZ = [0.0, 5.9173, 5.6125, 5.5377, 5.4728, 5.4174, 5.3688,
                       5.3236]

        if self.Z < 5:
            return F_EL_LOWZ[self.Z], F_INEL_LOWZ[self.Z]

        return (np.log(184.15) - np.log(self.Z)/3,
                np.log(1194) - 2*np.log(self.Z)/3)

    def _compute_Z_factors(self):
        """
        Compute the Z-dependent factors of the bremsstrahlung cross section.

        Parameters
        ----------
        None

        Returns
        -------
        f_Z_factor_1, f_Z_factor_2, f_Z : float
            Z-dependent factors, see the class attributes.
        """
        f_Z_factor_1 = (self.f_el - self.f_c) + self.f_inel/self.Z
        f_Z_factor_2 = (1 + 1/self.Z)/12
        f_Z = np.log(self.Z)/3 + self.f_c

        return f_Z_factor_1, f_Z_factor_2, f_Z

    def _compute_gamma_epsilon_factors(self):
        """
        Compute the prefactors of the Tsai screening variables.

        Parameters
        ----------
        None

        Returns
        -------
        gamma_factor, epsilon_factor : float
            Prefactors of the screening variables [eV].
        """
        # The Tsai screening variables are gamma = 100 m_e k /(E (E-k) Z^(1/3))
        # and epsilon = gamma / Z^(1/3), which are dimensionless only if m_e
        # and the energies are expressed in the same unit. Geant4 works in
        # MeV; here everything is in eV, so the electron mass enters in eV to
        # pair with `dum1 = y/(etot - photon_energy)` in 1/eV.
        gamma_factor = 100*ELECTRON_MASS_EV / np.cbrt(self.Z)
        epsilon_factor = 100*ELECTRON_MASS_EV / np.cbrt(self.Z)**2

        return gamma_factor, epsilon_factor


_ELEMENT_DATA_CACHE = {}


def element_data(Z):
    """
    Return the (cached) :class:`ElementData` of an element.

    Parameters
    ----------
    Z : int
        Atomic number.

    Returns
    -------
    data : ElementData
        Atomic data of the element.
    """
    Z = int(Z)
    if Z not in _ELEMENT_DATA_CACHE:
        _ELEMENT_DATA_CACHE[Z] = ElementData(Z)
    return _ELEMENT_DATA_CACHE[Z]


# ############################################################ #
# Bremsstrahlung
# ############################################################ #
class BremsstrahlungCalculator:
    """
    Bremsstrahlung of relativistic electrons and positrons on a gas nucleus.

    The differential cross section, its Gauss-Legendre integration and the
    photon-energy sampling follow the Geant4 relativistic bremsstrahlung model
    ``G4eBremsstrahlungRelModel``, in the complete-screening approximation with
    the Tsai screening functions. The photon emission angle is sampled with
    the modified Tsai distribution of ``G4ModifiedTsai``.

    The recoil of the emitted photon is applied to the beam particle by exact
    momentum conservation, so the particle acquires both an energy loss and a
    (small) angular deflection.

    Only photon energies above ``energy_cut`` are generated; ``compute_xsec``
    returns the cross section restricted to the same range. Choosing
    ``energy_cut`` close to the momentum acceptance of the machine therefore
    concentrates the Monte Carlo statistics on the events that can actually
    lead to a loss.

    Parameters
    ----------
    Z : int
        Atomic number of the gas species.
    p0c : float
        Reference momentum of the beam times the speed of light [eV].
    energy_cut : float, optional
        Lower limit of the generated photon energy [eV]. Default 10 keV.

    Attributes
    ----------
    Z : int
        Atomic number of the gas species.
    p0c : float
        Reference momentum times the speed of light [eV].
    energy_cut : float
        Lower limit of the generated photon energy [eV].
    element_data : ElementData
        Atomic data of the gas species.
    ekin : float
        Kinetic energy of the reference beam particle [eV].

    Notes
    -----
    The density (Ter-Mikaelian) correction is set to zero: for the gas
    densities of an accelerator vacuum chamber the plasma energy is many
    orders of magnitude below the photon-energy range of interest, so the
    correction is entirely negligible.

    References
    ----------
    .. [1] Y.-S. Tsai, Rev. Mod. Phys. **46**, 815 (1974).
    .. [2] S. M. Seltzer and M. J. Berger, Nucl. Instrum. Meth. B **12**,
       95 (1985).
    .. [3] Geant4 Collaboration, "Geant4 Physics Reference Manual", sections on
       ``G4eBremsstrahlungRelModel`` and ``G4ModifiedTsai``.
    """

    # Abscissas and weights of an 8-point Gauss-Legendre quadrature on [0, 1]
    _GL_X = np.array([1.98550718e-02, 1.01666761e-01, 2.37233795e-01,
                      4.08282679e-01, 5.91717321e-01, 7.62766205e-01,
                      8.98333239e-01, 9.80144928e-01])
    _GL_W = np.array([5.06142681e-02, 1.11190517e-01, 1.56853323e-01,
                      1.81341892e-01, 1.81341892e-01, 1.56853323e-01,
                      1.11190517e-01, 5.06142681e-02])

    def __init__(self, Z, p0c, energy_cut=10e3):
        self.Z = int(Z)
        self.p0c = float(p0c)
        self.energy_cut = float(energy_cut)
        self.element_data = element_data(Z)
        self.ekin = _energy_from_momentum(self.p0c) - ELECTRON_MASS_EV

        if self.energy_cut <= 0.0:
            raise ValueError("`brems_energy_cut` must be positive.")
        if self.energy_cut >= self.ekin:
            raise ValueError(
                f"`brems_energy_cut` ({self.energy_cut:.3e} eV) must be below "
                f"the beam kinetic energy ({self.ekin:.3e} eV).")

    def __repr__(self):
        return (f"<BremsstrahlungCalculator {self.element_data.symbol} "
                f"(Z={self.Z}), p0c={self.p0c:.3e} eV, "
                f"energy_cut={self.energy_cut:.3e} eV>")

    def _compute_screening_functions(self, gamma, epsilon):
        """
        Evaluate the Tsai screening functions.

        Parameters
        ----------
        gamma, epsilon : float or ndarray
            Screening variables for the elastic and inelastic contributions.

        Returns
        -------
        phi1, phi1m2, psi1, psi1m2 : float or ndarray
            Screening functions ``phi1``, ``phi1 - phi2``, ``psi1`` and
            ``psi1 - psi2``.
        """
        phi1 = (16.863 - 2*np.log(1 + 0.311877*gamma**2)
                + 2.4*np.exp(-0.9*gamma) + 1.6*np.exp(-1.5*gamma))
        phi1m2 = 2/(3 + 19.5*gamma + 18*gamma**2)
        psi1 = (24.34 - 2*np.log(1 + 13.111641*epsilon**2)
                + 2.8*np.exp(-8*epsilon) + 1.2*np.exp(-29.2*epsilon))
        psi1m2 = 2/(3 + 120*epsilon + 1200*epsilon**2)

        return phi1, phi1m2, psi1, psi1m2

    def _compute_dxsec(self, photon_energy):
        """
        Evaluate the scaled bremsstrahlung differential cross section.

        The returned quantity is the differential cross section per unit
        ``log(k)`` divided by the constant prefactor
        ``16 alpha r_e^2 Z^2 / 3``; it is the quantity integrated by
        :meth:`compute_xsec` and sampled by :meth:`sample_photon_energies`.

        Parameters
        ----------
        photon_energy : float or ndarray
            Energy of the emitted photon [eV].

        Returns
        -------
        dxsec : float or ndarray
            Scaled differential cross section.
        """
        etot = self.ekin + ELECTRON_MASS_EV
        y = photon_energy/etot
        dum0 = (1 - y) + 0.75*y**2
        dum1 = y/(etot - photon_energy)
        gamma = dum1 * self.element_data.gamma_factor
        epsilon = dum1 * self.element_data.epsilon_factor

        if self.Z < 5:
            return (dum0*self.element_data.f_Z_factor_1
                    + (1 - y)*self.element_data.f_Z_factor_2)

        phi1, phi1m2, psi1, psi1m2 = self._compute_screening_functions(
            gamma, epsilon)

        return (dum0*((0.25*phi1 - self.element_data.f_Z)
                      + (0.25*psi1 - 2*np.log(self.Z)/3)/self.Z)
                + 0.125*(1 - y)*(phi1m2 + psi1m2/self.Z))

    def compute_xsec(self):
        """
        Compute the total bremsstrahlung cross section above the energy cut.

        The differential cross section is integrated over
        ``log(k)`` from ``energy_cut`` to the beam kinetic energy with the
        adaptive sub-interval Gauss-Legendre scheme of
        ``G4eBremsstrahlungRelModel::ComputeCrossSectionPerAtom``.

        Parameters
        ----------
        None

        Returns
        -------
        xsec : float
            Cross section per gas atom [m^2].
        """
        etot = _energy_from_momentum(self.p0c)
        alpha_min = np.log(self.energy_cut/etot)
        alpha_max = np.log(self.ekin/self.energy_cut)
        n_sub = max(int(0.45*alpha_max), 0) + 4
        delta = alpha_max/n_sub

        # Sub-interval left edges, and the quadrature nodes inside each of them
        alpha_i = alpha_min + delta*np.arange(n_sub)
        photon_energy = np.exp(
            alpha_i[:, None] + self._GL_X[None, :]*delta) * etot
        xsec = np.sum(self._GL_W[None, :] * self._compute_dxsec(photon_energy))
        xsec *= delta

        return 16*ALPHA*CLASSICAL_ELECTRON_RADIUS**2*self.Z**2/3 * xsec

    def sample_photon_energies(self, n, rng):
        """
        Sample photon energies from the bremsstrahlung spectrum.

        The sampling follows ``G4eBremsstrahlungRelModel::SampleSecondaries``:
        a trial energy is drawn uniformly in ``log(k^2)`` between the energy
        cut and the beam kinetic energy, and accepted against the majorant
        ``f_Z_factor_1 + f_Z_factor_2`` of the differential cross section.

        Parameters
        ----------
        n : int
            Number of photon energies to sample.
        rng : numpy.random.Generator
            Random number generator.

        Returns
        -------
        photon_energy : ndarray
            Sampled photon energies [eV].
        """
        func_max = (self.element_data.f_Z_factor_1
                    + self.element_data.f_Z_factor_2)
        # Density (Ter-Mikaelian) correction, negligible for a dilute gas
        f_density_corr = 0.0
        x_min = np.log(self.energy_cut**2 + f_density_corr)
        x_range = np.log(self.ekin**2 + f_density_corr) - x_min

        def propose(batch_size, rng):
            rndm = rng.random((2, batch_size))
            photon_energy = np.sqrt(np.maximum(
                np.exp(x_min + rndm[0]*x_range) - f_density_corr, 0.0))
            accept = self._compute_dxsec(photon_energy) >= func_max*rndm[1]
            return photon_energy[accept]

        return _rejection_sample(n, propose, rng, 'photon energy')

    def sample_photon_cos_theta(self, n, rng):
        """
        Sample the cosine of the photon emission angle.

        The modified Tsai angular distribution of ``G4ModifiedTsai`` is used:
        it depends only on the Lorentz factor of the emitting particle and not
        on the photon energy, which is an excellent approximation in the
        ultra-relativistic regime.

        Parameters
        ----------
        n : int
            Number of angles to sample.
        rng : numpy.random.Generator
            Random number generator.

        Returns
        -------
        cos_theta : ndarray
            Cosine of the photon emission angle with respect to the direction
            of the emitting particle.
        """
        u_max = 2*(1 + self.ekin/ELECTRON_MASS_EV)
        a1 = 1.6
        a2 = a1/3.0
        border = 0.25

        def propose(batch_size, rng):
            rndm = rng.random((3, batch_size))
            uu = -np.log(rndm[0]*rndm[1])
            u = np.where(rndm[2] < border, uu*a1, uu*a2)
            return u[u <= u_max]

        u = _rejection_sample(n, propose, rng, 'photon emission angle')

        return 1.0 - 2.0*u**2/u_max**2

    def sample_deflections(self, px, py, delta, rng):
        """
        Sample bremsstrahlung interactions for a batch of beam particles.

        A photon energy and emission direction are drawn for each particle,
        and the outgoing beam-particle momentum is obtained by exact momentum
        conservation.

        Parameters
        ----------
        px, py : ndarray
            Normalised transverse momenta before the interaction.
        delta : ndarray
            Relative momentum deviation before the interaction.
        rng : numpy.random.Generator
            Random number generator.

        Returns
        -------
        sample : ScatteringSample
            Outgoing momenta, unit weights, photon emission angles and photon
            energies.
        """
        n = px.size
        photon_energy = self.sample_photon_energies(n, rng)
        cos_theta = self.sample_photon_cos_theta(n, rng)
        sin_theta = np.sqrt(np.maximum(1.0 - cos_theta**2, 0.0))
        phi = rng.uniform(0.0, 2.0*np.pi, n)

        u_hat, v_hat, p_hat, p_norm = _momentum_frame(px, py, delta)
        photon_dir = ((sin_theta*np.cos(phi))[:, None]*u_hat
                      + (sin_theta*np.sin(phi))[:, None]*v_hat
                      + cos_theta[:, None]*p_hat)

        # Momentum conservation, in units of p0c
        p_out = (p_norm[:, None]*p_hat
                 - (photon_energy/self.p0c)[:, None]*photon_dir)

        return ScatteringSample(
            px=p_out[:, 0],
            py=p_out[:, 1],
            delta=np.linalg.norm(p_out, axis=1) - 1.0,
            weight=np.ones(n),
            theta=np.arccos(np.clip(cos_theta, -1.0, 1.0)),
            photon_energy=photon_energy,
        )


# ############################################################ #
# Coulomb scattering
# ############################################################ #
class CoulombScatteringCalculator:
    """
    Elastic Coulomb scattering of electrons and positrons on a gas nucleus.

    The differential cross section is the Rutherford cross section corrected
    by the Moliere screening factor, the McKinley-Feshbach Mott-to-Rutherford
    ratio and the exponential nuclear form factor, as implemented in
    ``G4ScreeningMottCrossSection``. In terms of ``z = 1 - cos(theta)`` it
    reads

    .. math::
        \\frac{d\\sigma}{dz} = \\frac{2 \\pi Z^2 r_e^2}{\\beta^4 \\gamma^2}
            \\, \\frac{R_{MF}(z)\\, F_N^2(z)}{(2 A_s + z)^2}

    with ``A_s`` the Moliere screening parameter. The energy transferred to
    the nucleus is negligible for a relativistic lepton on a gas atom, so
    ``delta`` is left unchanged by the interaction.

    Events are generated only within ``theta_lim``, and ``compute_xsec``
    returns the cross section restricted to the same range. Scattering
    outside ``theta_lim`` is not represented at all, so any loss it would
    cause is missing from the result.

    Within ``theta_lim`` the polar angle is drawn from a log-uniform proposal
    in ``z``, which over-samples the large angles that actually cause losses;
    the resulting bias is corrected exactly by the importance weight returned
    in :attr:`ScatteringSample.weight`. With this proposal the fraction of
    events above a given angle only decreases as ``1/log(theta_max/theta_min)``
    when ``theta_min`` is lowered, so a wide ``theta_lim`` costs little
    statistics, whereas a ``theta_lim`` that cuts into the loss-relevant range
    biases the loss rate low.

    Parameters
    ----------
    Z : int
        Atomic number of the gas species.
    p0c : float
        Reference momentum of the beam times the speed of light [eV].
    q0 : float, optional
        Charge of the beam particle in units of the elementary charge:
        ``-1`` for electrons (default), ``+1`` for positrons. It enters
        through the sign of the McKinley-Feshbach term.
    theta_lim : tuple of float, optional
        Lower and upper limit of the generated polar scattering angle [rad].

    Attributes
    ----------
    Z : int
        Atomic number of the gas species.
    p0c : float
        Reference momentum times the speed of light [eV].
    q0 : float
        Charge of the beam particle in units of the elementary charge.
    theta_lim : tuple of float
        Generated polar scattering-angle range [rad].
    element_data : ElementData
        Atomic data of the gas species.
    ekin : float
        Kinetic energy of the reference beam particle [eV].
    gamma, beta : float
        Relativistic factors of the reference beam particle.
    a_TF : float
        Thomas-Fermi screening radius of the gas atom [m].
    tmax : float
        Maximum kinetic energy transferable to the recoiling nucleus [eV].

    References
    ----------
    .. [1] G. Moliere, Z. Naturforsch. A **2**, 133 (1947).
    .. [2] W. A. McKinley and H. Feshbach, Phys. Rev. **74**, 1759 (1948).
    .. [3] M. J. Boschini et al., Proc. ICATPP 2010, p. 9 (2011).
    """

    def __init__(self, Z, p0c, q0=-1.0, theta_lim=(1e-7, 50e-3)):
        self.Z = int(Z)
        self.p0c = float(p0c)
        self.q0 = float(q0)
        self.theta_lim = (float(theta_lim[0]), float(theta_lim[1]))
        self.element_data = element_data(Z)

        if not 0.0 < self.theta_lim[0] < self.theta_lim[1] <= np.pi:
            raise ValueError(
                "`coulomb_theta` must be a pair (theta_min, theta_max) with "
                "0 < theta_min < theta_max <= pi.")

        etot = _energy_from_momentum(self.p0c)
        self.ekin = etot - ELECTRON_MASS_EV
        self.gamma = etot/ELECTRON_MASS_EV
        self.beta = np.sqrt(1.0 - 1.0/self.gamma**2)

        # Thomas-Fermi screening radius
        self.a_TF = C_TF * BOHR_RADIUS * self.Z**(-1/3)

        # Maximum kinetic energy transferred to the recoiling nucleus
        mass = self.element_data.mass
        self.tmax = (2*mass*self.ekin*(self.ekin + 2*ELECTRON_MASS_EV)
                     / (ELECTRON_MASS_EV**2 + mass**2 + 2*mass*etot))

        self._screening_As = self._compute_screening_As()
        # z = 1 - cos(theta) = 2 sin^2(theta/2); the latter form keeps full
        # relative precision at small angles, where 1 - cos(theta) cancels
        # catastrophically (it is exactly zero below theta ~ 1e-8)
        self._z_lim = (2.0*np.sin(0.5*self.theta_lim[0])**2,
                       2.0*np.sin(0.5*self.theta_lim[1])**2)

        # Integrated once here, and reused as the normalisation of the
        # importance weights, so that the sampler and the total cross section
        # are by construction the same differential cross section.
        self.xsec = self.compute_xsec()

    def __repr__(self):
        return (f"<CoulombScatteringCalculator {self.element_data.symbol} "
                f"(Z={self.Z}), p0c={self.p0c:.3e} eV, "
                f"theta_lim={self.theta_lim}>")

    def _compute_screening_As(self):
        """
        Compute the Moliere screening parameter.

        Parameters
        ----------
        None

        Returns
        -------
        As : float
            Moliere screening parameter, in units of ``1 - cos(theta)``.
        """
        return ((HBAR_C_EV_M/(2.0*self.p0c*self.a_TF))**2
                * (1.13 + 3.76*(ALPHA*self.Z/self.beta)**2))

    def _mott_ratio(self, z):
        """
        Evaluate the McKinley-Feshbach Mott-to-Rutherford ratio.

        Parameters
        ----------
        z : ndarray
            Scattering variable ``z = 1 - cos(theta) = 2 sin^2(theta/2)``.
            Working directly in ``z`` avoids the ``arccos`` round-trip, which
            loses relative precision at the small angles that dominate the
            cross section.

        Returns
        -------
        ratio : ndarray
            Ratio of the Mott to the Rutherford differential cross section.
            The sign of the interference term follows the charge of the beam
            particle: it increases the cross section for electrons and
            decreases it for positrons.
        """
        s2 = 0.5*z
        s = np.sqrt(s2)

        return (1.0 - self.beta**2*s2
                - self.q0*self.Z*ALPHA*self.beta*np.pi*s*(1.0 - s))

    def _nuclear_form_factor2(self, z):
        """
        Evaluate the squared exponential nuclear form factor.

        Follows ``G4ScreeningMottCrossSection::FormFactor2ExpHof``.

        Parameters
        ----------
        z : ndarray
            Scattering variable ``z = 1 - cos(theta) = 2 sin^2(theta/2)``.

        Returns
        -------
        form_factor2 : ndarray
            Squared nuclear form factor. It is within a few ppm of one for the
            small angles relevant to beam-gas losses in a storage ring, and
            only becomes significant at large angles, where it can suppress
            the cross section by tens of percent for heavy gas species.
        """
        t = self.tmax * 0.5*z
        q2 = t*(t + 2.0*self.element_data.mass) / HBAR_C_EV_M**2
        form_factor = 1.0/(1.0 + self.element_data.nuclear_radius**2*q2/12.0)

        return form_factor**2

    def _dxsec_dz(self, z):
        """
        Evaluate the differential cross section in ``z = 1 - cos(theta)``.

        This is the single definition of the Coulomb differential cross
        section: it is both integrated by :meth:`compute_xsec` and used to
        weight the events drawn by :meth:`sample_theta`, so the total cross
        section and the generated sample cannot drift apart.

        Parameters
        ----------
        z : ndarray
            Scattering variable ``z = 1 - cos(theta)``.

        Returns
        -------
        dxsec_dz : ndarray
            Differential cross section ``dsigma/dz`` [m^2].
        """
        return (2.0*np.pi*self.Z**2*CLASSICAL_ELECTRON_RADIUS**2
                / (self.beta**4*self.gamma**2)
                * self._mott_ratio(z) * self._nuclear_form_factor2(z)
                / (2*self._screening_As + z)**2)

    def compute_xsec(self):
        """
        Compute the total cross section within ``theta_lim``.

        Integrates :meth:`_dxsec_dz` -- the same differential cross section
        the sampler uses, i.e. screened Rutherford times the McKinley-Feshbach
        Mott ratio and the squared nuclear form factor -- over
        ``z = 1 - cos(theta)``. The integration is carried out in ``log(z)``,
        where the integrand is a smooth bump even when ``theta_lim`` spans
        many decades, with a fixed 128-node Gauss-Legendre rule (accurate to
        ~1e-14 relative over the full range of interest).

        Parameters
        ----------
        None

        Returns
        -------
        xsec : float
            Cross section per gas atom, restricted to ``theta_lim`` [m^2].
        """
        t1, t2 = np.log(self._z_lim[0]), np.log(self._z_lim[1])
        nodes, weights = np.polynomial.legendre.leggauss(128)
        t = 0.5*(t2 - t1)*nodes + 0.5*(t2 + t1)
        z = np.exp(t)

        # dz = z dt for the log substitution
        return 0.5*(t2 - t1) * float(np.sum(weights * self._dxsec_dz(z) * z))

    def sample_theta(self, n, rng):
        """
        Sample polar scattering angles with importance sampling.

        The proposal density is log-uniform in ``z = 1 - cos(theta)``,
        ``g(z) = 1/(z log(z2/z1))``, which strongly over-samples the large
        angles responsible for the losses. The returned weights restore the
        correct differential cross section, including the Mott and nuclear
        form-factor corrections.

        Parameters
        ----------
        n : int
            Number of angles to sample.
        rng : numpy.random.Generator
            Random number generator.

        Returns
        -------
        theta : ndarray
            Sampled polar scattering angles [rad].
        weight : ndarray
            Importance weights, normalised to the total cross section
            returned by :meth:`compute_xsec`, so that their expectation value
            is exactly one.
        """
        z1, z2 = self._z_lim
        log_ratio = np.log(z2/z1)

        u = rng.random(n)
        z = z1*np.exp(u*log_ratio)
        theta = 2.0*np.arcsin(np.minimum(np.sqrt(0.5*z), 1.0))

        # w = (dsigma/dz) / (g(z) * sigma), with g(z) = 1/(z log(z2/z1))
        weight = self._dxsec_dz(z) * z * log_ratio / self.xsec

        return theta, weight

    def sample_deflections(self, px, py, delta, rng):
        """
        Sample Coulomb scattering events for a batch of beam particles.

        Parameters
        ----------
        px, py : ndarray
            Normalised transverse momenta before the interaction.
        delta : ndarray
            Relative momentum deviation before the interaction.
        rng : numpy.random.Generator
            Random number generator.

        Returns
        -------
        sample : ScatteringSample
            Outgoing momenta, importance weights and scattering angles. The
            relative momentum deviation is returned unchanged, since the
            nuclear recoil energy is negligible.
        """
        n = px.size
        theta, weight = self.sample_theta(n, rng)
        phi = rng.uniform(0.0, 2.0*np.pi, n)

        sin_theta = np.sin(theta)
        cos_theta = np.cos(theta)

        u_hat, v_hat, p_hat, p_norm = _momentum_frame(px, py, delta)
        scattered_dir = ((sin_theta*np.cos(phi))[:, None]*u_hat
                         + (sin_theta*np.sin(phi))[:, None]*v_hat
                         + cos_theta[:, None]*p_hat)
        p_out = p_norm[:, None]*scattered_dir

        return ScatteringSample(
            px=p_out[:, 0],
            py=p_out[:, 1],
            delta=delta.copy(),
            weight=weight,
            theta=theta,
            photon_energy=np.full(n, np.nan),
        )
