# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2026.                 #
# ######################################### #
import warnings

import numpy as np
import pytest
from scipy import constants as sc
from scipy.integrate import quad

import xtrack as xt
import xcoll as xc
from xcoll.beamgas.cross_sections import (BremsstrahlungCalculator,
                                          CoulombScatteringCalculator,
                                          ElementData, HBAR_C_EV_M, ALPHA,
                                          CLASSICAL_ELECTRON_RADIUS,
                                          ELECTRON_MASS_EV)


#############################################################
# Shared parameters
#############################################################
NEMITT_X = 1e-5
NEMITT_Y = 1e-7
SIGMA_Z = 4e-3
SIGMA_DELTA = 1e-3
BUNCH_INTENSITY = 4e9
P0C = 1e9

# N2 at 1e-7 mbar and room temperature, as atomic density
ATOMIC_DENSITY = 2 * 1e-7 * 1e2 / (sc.Boltzmann * 293.15)


#############################################################
# Module-level fixture: toy ring with beam-gas scattering centres
#############################################################
@pytest.fixture(scope='module')
def toy_ring():
    """
    Build a FODO-like toy ring with BeamGasScattering elements and apertures.

    Yields a dict with keys ``line`` and ``gas_density``.
    """
    env = xt.Environment()
    line = env.new_line(components=[
        env.new('mqf.1', xt.Quadrupole, length=0.3, k1=0.1),
        env.new('d1.1', xt.Drift, length=1),
        env.new('mb1.1', xt.Bend, length=3, angle=np.pi/2),
        env.new('d2.1', xt.Drift, length=1),
        env.new('mqd.1', xt.Quadrupole, length=0.3, k1=-0.7),
        env.new('d3.1', xt.Drift, length=1),
        env.new('mb2.1', xt.Bend, length=3, angle=np.pi/2),
        env.new('d4.1', xt.Drift, length=1),
        env.new('mqf.2', xt.Quadrupole, length=0.3, k1=0.1),
        env.new('d1.2', xt.Drift, length=1),
        env.new('mb1.2', xt.Bend, length=3, angle=np.pi/2),
        env.new('d2.2', xt.Drift, length=1),
        env.new('mqd.2', xt.Quadrupole, length=0.3, k1=-0.7),
        env.new('d3.2', xt.Drift, length=1),
        env.new('mb2.2', xt.Bend, length=3, angle=np.pi/2),
        env.new('d4.2', xt.Drift, length=1),
    ])
    line.set_particle_ref('electron', p0c=P0C)
    line.configure_bend_model(core='full', edge=None)

    circumference = line.get_length()
    placements = []
    for ii, ss in enumerate(np.linspace(0, circumference, 9)[1:]):
        env.elements[f'BeamGasScattering.{ii}'] = xc.BeamGasScattering()
        placements.append(env.place(f'BeamGasScattering.{ii}', at=ss))
    line.insert(placements)

    tab = line.get_table()
    needs_aperture = tab.rows.match_not(
        element_type='Drift.*|Marker|BeamGasScattering|').name
    env.new('aper', xt.LimitRect, min_x=-0.04, max_x=0.04,
            min_y=-0.04, max_y=0.04)
    placements = []
    for nn in needs_aperture:
        env.new(f'{nn}_aper_entry', 'aper')
        env.new(f'{nn}_aper_exit', 'aper')
        placements.append(env.place(f'{nn}_aper_entry', at=f'{nn}@start'))
        placements.append(env.place(f'{nn}_aper_exit', at=f'{nn}@end'))
    line.insert(placements)
    line.build_tracker()

    tab = line.get_table()
    tt = tab.rows[tab.element_type == 'BeamGasScattering']
    gas_density = xt.Table({
        'name': tt.name,
        's': tt.s,
        'N': np.ones(len(tt.name))*ATOMIC_DENSITY,
    })

    yield {'line': line, 'gas_density': gas_density}


def make_study(toy_ring, **kwargs):
    """Build and initialise a BeamGasStudy with the toy-ring defaults."""
    kwargs.setdefault('process', 'brems')
    kwargs.setdefault('brems_energy_cut', 1e6)
    kwargs.setdefault('n_scattering_events', 200)
    kwargs.setdefault('seed', 1997)
    study = xc.BeamGasStudy(
        line=toy_ring['line'],
        gas_density=kwargs.pop('gas_density', toy_ring['gas_density']),
        nemitt_x=NEMITT_X, nemitt_y=NEMITT_Y,
        sigma_z=SIGMA_Z, sigma_delta=SIGMA_DELTA,
        bunch_intensity=BUNCH_INTENSITY,
        method='4d', **kwargs)
    study.initialise_beamgas(verbose=False)
    return study


#############################################################
# Atomic data
#############################################################
class TestElementData:

    def test_known_elements(self):
        for Z, symbol, A in [(1, 'H', 1), (6, 'C', 12), (7, 'N', 14),
                             (8, 'O', 16), (18, 'Ar', 40)]:
            data = ElementData(Z)
            assert data.symbol == symbol
            assert data.A == A

    def test_raises_out_of_range(self):
        with pytest.raises(ValueError):
            ElementData(0)
        with pytest.raises(ValueError):
            ElementData(55)

    def test_low_z_uses_dirac_fock_logarithms(self):
        # Below Z=5 the tabulated Dirac-Fock values must be used, and they
        # differ from the Thomas-Fermi expression
        for Z in (1, 2, 3, 4):
            data = ElementData(Z)
            assert not np.isclose(data.f_el,
                                  np.log(184.15) - np.log(Z)/3, rtol=1e-3)
        for Z in (5, 7, 18):
            data = ElementData(Z)
            assert np.isclose(data.f_el, np.log(184.15) - np.log(Z)/3)


#############################################################
# Cross sections
#############################################################
class TestCoulombCrossSection:

    @staticmethod
    def _quad_xsec(calc):
        """Integrate the calculator's own dsigma/dz independently, with quad.

        The integrand is sharply peaked at small z, so it is integrated in
        log(z), where it is a smooth bump.
        """
        z1, z2 = calc._z_lim
        value, _ = quad(lambda t: calc._dxsec_dz(np.exp(t))*np.exp(t),
                        np.log(z1), np.log(z2),
                        epsabs=0, epsrel=1e-12, limit=400)
        return value

    @pytest.mark.parametrize('Z', [1, 6, 7, 18, 54])
    @pytest.mark.parametrize('theta_lim', [(1e-5, 50e-3), (20e-3, 200e-3),
                                           (0.1, np.pi)])
    def test_xsec_matches_numerical_integral(self, Z, theta_lim):
        # compute_xsec must integrate exactly the dsigma/dz that the sampler
        # weights with -- including the Mott ratio and the nuclear form
        # factor, which suppress the cross section by tens of percent at
        # large angles for heavy species
        calc = CoulombScatteringCalculator(Z, P0C, q0=-1.0,
                                           theta_lim=theta_lim)
        assert np.isclose(calc.compute_xsec(), self._quad_xsec(calc),
                          rtol=1e-10)

    @pytest.mark.parametrize('Z,theta_lim', [(7, (1e-3, 50e-3)),
                                             (54, (20e-3, 200e-3))])
    def test_importance_weights_are_normalised(self, Z, theta_lim):
        # Weights are normalised to compute_xsec, so their mean is one
        calc = CoulombScatteringCalculator(Z, P0C, q0=-1.0,
                                           theta_lim=theta_lim)
        rng = np.random.default_rng(0)
        _, weight = calc.sample_theta(400_000, rng)
        assert np.isclose(weight.mean(), 1.0,
                          atol=5*weight.std()/np.sqrt(weight.size))

    def test_sampled_distribution_matches_cross_section(self):
        calc = CoulombScatteringCalculator(7, P0C, q0=-1.0,
                                           theta_lim=(1e-3, 50e-3))
        rng = np.random.default_rng(0)
        theta, weight = calc.sample_theta(400_000, rng)

        # The weighted sample must reproduce the fraction of the cross section
        # contained in any sub-range of theta
        z = 1.0 - np.cos(theta)
        z1, z2 = calc._z_lim
        total = self._quad_xsec(calc)
        for z_cut in (z1*10, np.sqrt(z1*z2), z2/10):
            observed = weight[z < z_cut].sum()/weight.size
            partial, _ = quad(
                lambda t: calc._dxsec_dz(np.exp(t))*np.exp(t),
                np.log(z1), np.log(z_cut), epsabs=0, epsrel=1e-12, limit=400)
            assert np.isclose(observed, partial/total, rtol=0.02)

    def test_mott_correction_sign_flips_with_charge(self):
        electron = CoulombScatteringCalculator(7, P0C, q0=-1.0,
                                               theta_lim=(8e-3, 20e-3))
        positron = CoulombScatteringCalculator(7, P0C, q0=+1.0,
                                               theta_lim=(8e-3, 20e-3))
        z = np.array([1.0 - np.cos(1e-2)])
        # McKinley-Feshbach enhances the cross section for electrons
        assert electron._mott_ratio(z)[0] > 1.0
        assert positron._mott_ratio(z)[0] < 1.0
        # ... and the total cross section now carries that charge dependence,
        # since compute_xsec integrates the Mott-corrected dsigma/dz
        assert electron.compute_xsec() > positron.compute_xsec()

    def test_deflection_conserves_momentum_modulus(self):
        calc = CoulombScatteringCalculator(7, P0C, q0=-1.0,
                                           theta_lim=(8e-3, 20e-3))
        rng = np.random.default_rng(1)
        n = 5000
        px = rng.normal(0, 1e-5, n)
        py = rng.normal(0, 1e-5, n)
        delta = rng.normal(0, 1e-3, n)
        sample = calc.sample_deflections(px, py, delta, rng)

        # Elastic scattering: delta is unchanged
        assert np.array_equal(sample.delta, delta)
        # The angle between the incoming and outgoing momentum is theta
        pz_in = np.sqrt((1 + delta)**2 - px**2 - py**2)
        pz_out = np.sqrt((1 + sample.delta)**2 - sample.px**2 - sample.py**2)
        cos_angle = ((px*sample.px + py*sample.py + pz_in*pz_out)
                     / ((1 + delta)*(1 + sample.delta)))
        assert np.allclose(np.arccos(np.clip(cos_angle, -1, 1)), sample.theta,
                           atol=1e-9)

    def test_small_angles_keep_full_precision(self):
        # 1 - cos(theta) is exactly zero below theta ~ 1e-8, so the window
        # and the sampled angles must be computed without it
        calc = CoulombScatteringCalculator(7, P0C, q0=-1.0,
                                           theta_lim=(1e-10, 1e-3))
        assert np.isclose(calc._z_lim[0], 0.5e-20, rtol=1e-12)
        assert np.isfinite(calc.xsec) and calc.xsec > 0
        theta, weight = calc.sample_theta(100_000, np.random.default_rng(0))
        assert theta.min() >= 1e-10*(1 - 1e-12)
        assert theta.max() <= 1e-3*(1 + 1e-12)
        assert np.all(np.isfinite(weight))
        # Below the screening angle the cross section is flat in z, so
        # lowering theta_min further only adds a negligible amount
        wider = CoulombScatteringCalculator(7, P0C, q0=-1.0,
                                            theta_lim=(1e-12, 1e-3))
        assert np.isclose(wider.xsec, calc.xsec, rtol=1e-6)

    def test_raises_on_invalid_theta_range(self):
        with pytest.raises(ValueError):
            CoulombScatteringCalculator(7, P0C, theta_lim=(2e-2, 1e-3))
        with pytest.raises(ValueError):
            CoulombScatteringCalculator(7, P0C, theta_lim=(0.0, 1e-3))


class TestBremsstrahlungCrossSection:

    @pytest.mark.parametrize('Z', [1, 6, 7, 18])
    def test_xsec_matches_numerical_integral(self, Z):
        calc = BremsstrahlungCalculator(Z, P0C, energy_cut=1e6)
        etot = np.sqrt(P0C**2 + xt.ELECTRON_MASS_EV**2)
        prefactor = 16*ALPHA*CLASSICAL_ELECTRON_RADIUS**2*Z**2/3

        # The Gauss-Legendre scheme integrates over alpha = log(k/etot)
        numerical, _ = quad(
            lambda a: calc._compute_dxsec(np.exp(a)*etot),
            np.log(calc.energy_cut/etot), np.log(calc.ekin/etot),
            epsabs=0, epsrel=1e-10)
        assert np.isclose(calc.compute_xsec(), prefactor*numerical, rtol=1e-6)

    def test_sampled_spectrum_matches_cross_section(self):
        calc = BremsstrahlungCalculator(7, P0C, energy_cut=1e6)
        rng = np.random.default_rng(0)
        k = calc.sample_photon_energies(200_000, rng)

        etot = np.sqrt(P0C**2 + xt.ELECTRON_MASS_EV**2)
        norm, _ = quad(lambda a: calc._compute_dxsec(np.exp(a)*etot),
                       np.log(calc.energy_cut/etot), np.log(calc.ekin/etot),
                       epsabs=0, epsrel=1e-10)
        for k_cut in (1e7, 1e8, 5e8):
            observed = np.mean(k < k_cut)
            expected, _ = quad(
                lambda a: calc._compute_dxsec(np.exp(a)*etot),
                np.log(calc.energy_cut/etot), np.log(k_cut/etot),
                epsabs=0, epsrel=1e-10)
            assert np.isclose(observed, expected/norm, rtol=0.02)

    @pytest.mark.parametrize('Z', [1, 7, 18, 54])
    def test_screening_variables_are_dimensionless(self, Z):
        # The Tsai screening variables are
        #   gamma   = 100 m_e k / (E (E-k) Z^(1/3))
        #   epsilon = gamma / Z^(1/3)
        # which are dimensionless only if m_e and the energies share a unit.
        # Everything here is in eV, so the prefactors must carry m_e in eV;
        # a stray eV->MeV conversion would make them 1e6 too small and pin
        # the screening functions to the complete-screening limit.
        calc = BremsstrahlungCalculator(Z, P0C, energy_cut=1e6)
        etot = calc.ekin + ELECTRON_MASS_EV
        k = 0.99*etot
        dum1 = (k/etot)/(etot - k)
        expected_gamma = 100*ELECTRON_MASS_EV*k/(etot*(etot - k)*np.cbrt(Z))
        assert np.isclose(dum1*calc.element_data.gamma_factor, expected_gamma,
                          rtol=1e-12)
        assert np.isclose(dum1*calc.element_data.epsilon_factor,
                          expected_gamma/np.cbrt(Z), rtol=1e-12)
        # ... and at the tip of the spectrum screening must actually bite
        assert expected_gamma > 1.0

    def test_screening_suppresses_the_hard_tip(self):
        # Screening only matters where gamma is O(1), i.e. y -> 1. Check the
        # differential cross section is pushed below its complete-screening
        # (gamma = 0) value there, and is untouched at small y.
        calc = BremsstrahlungCalculator(7, P0C, energy_cut=1e6)
        etot = calc.ekin + ELECTRON_MASS_EV
        complete_screening = (calc.element_data.f_Z_factor_1,
                              calc.element_data.f_Z_factor_2)

        def unscreened(y):
            dum0 = (1 - y) + 0.75*y**2
            return dum0*complete_screening[0] + (1 - y)*complete_screening[1]

        for y, expected_ratio in [(1e-3, 1.0), (0.99, 0.75)]:
            ratio = calc._compute_dxsec(y*etot)/unscreened(y)
            if expected_ratio == 1.0:
                assert np.isclose(ratio, 1.0, rtol=1e-3)
            else:
                assert ratio < expected_ratio

    def test_sampled_energies_within_bounds(self):
        calc = BremsstrahlungCalculator(7, P0C, energy_cut=1e6)
        rng = np.random.default_rng(2)
        k = calc.sample_photon_energies(50_000, rng)
        assert k.min() >= calc.energy_cut
        assert k.max() <= calc.ekin

    def test_deflection_conserves_energy(self):
        calc = BremsstrahlungCalculator(7, P0C, energy_cut=1e6)
        rng = np.random.default_rng(3)
        n = 5000
        px = np.zeros(n)
        py = np.zeros(n)
        delta = np.zeros(n)
        sample = calc.sample_deflections(px, py, delta, rng)

        # The particle always loses momentum, and the loss is bounded by the
        # photon energy (equality only for exactly collinear emission)
        assert np.all(sample.delta < 0)
        momentum_loss = -sample.delta*P0C
        assert np.all(momentum_loss <= sample.photon_energy*(1 + 1e-9))

        # Exact vector momentum conservation: starting along z with p = p0c,
        # |p_out| = sqrt(p0c^2 - 2 k p0c cos(theta) + k^2). Asserting this
        # rather than momentum_loss == photon_energy matters at the tip of
        # the spectrum, where p0c - k -> 0 and the transverse recoil is no
        # longer a negligible correction.
        kk = sample.photon_energy/P0C
        expected = np.sqrt(1.0 - 2.0*kk*np.cos(sample.theta) + kk**2)
        assert np.allclose(1.0 + sample.delta, expected, rtol=1e-12)

        # Away from the tip the loss is the photon energy to good accuracy
        soft = sample.photon_energy < 0.9*(calc.ekin + ELECTRON_MASS_EV)
        assert np.allclose(momentum_loss[soft], sample.photon_energy[soft],
                           rtol=1e-3)
        # Bremsstrahlung is sampled without biasing
        assert np.all(sample.weight == 1.0)

    def test_raises_on_invalid_energy_cut(self):
        with pytest.raises(ValueError):
            BremsstrahlungCalculator(7, P0C, energy_cut=-1.0)
        with pytest.raises(ValueError):
            BremsstrahlungCalculator(7, P0C, energy_cut=1e12)


#############################################################
# Study construction and validation
#############################################################
class TestBeamGasStudyValidation:

    def test_construction_with_normalised_emittances(self, toy_ring):
        study = make_study(toy_ring)
        beta0 = toy_ring['line'].particle_ref.beta0[0]
        gamma0 = toy_ring['line'].particle_ref.gamma0[0]
        assert np.isclose(study.gemitt_x, NEMITT_X/(beta0*gamma0))
        assert np.isclose(study.gemitt_y, NEMITT_Y/(beta0*gamma0))

    def test_construction_with_geometric_emittances(self, toy_ring):
        study = xc.BeamGasStudy(
            line=toy_ring['line'], gas_density=toy_ring['gas_density'],
            process='brems', gemitt_x=1e-9, gemitt_y=1e-11,
            sigma_z=SIGMA_Z, sigma_delta=SIGMA_DELTA,
            n_scattering_events=10, method='4d')
        assert np.isclose(study.gemitt_x, 1e-9)

    @pytest.mark.parametrize('missing_key', ['line', 'gas_density', 'sigma_z',
                                             'sigma_delta',
                                             'n_scattering_events'])
    def test_raises_on_missing_required_kwarg(self, toy_ring, missing_key):
        kwargs = dict(line=toy_ring['line'],
                      gas_density=toy_ring['gas_density'],
                      process='brems', nemitt_x=NEMITT_X, nemitt_y=NEMITT_Y,
                      sigma_z=SIGMA_Z, sigma_delta=SIGMA_DELTA,
                      n_scattering_events=10)
        kwargs[missing_key] = None
        with pytest.raises(ValueError):
            xc.BeamGasStudy(**kwargs)

    def test_raises_on_invalid_process(self, toy_ring):
        with pytest.raises(ValueError):
            make_study(toy_ring, process='compton')

    def test_raises_on_both_nemitt_and_gemitt(self, toy_ring):
        with pytest.raises(ValueError):
            xc.BeamGasStudy(
                line=toy_ring['line'], gas_density=toy_ring['gas_density'],
                process='brems', nemitt_x=NEMITT_X, nemitt_y=NEMITT_Y,
                gemitt_x=1e-9, gemitt_y=1e-11, sigma_z=SIGMA_Z,
                sigma_delta=SIGMA_DELTA, n_scattering_events=10)

    def test_raises_on_neither_nemitt_nor_gemitt(self, toy_ring):
        with pytest.raises(ValueError):
            xc.BeamGasStudy(
                line=toy_ring['line'], gas_density=toy_ring['gas_density'],
                process='brems', sigma_z=SIGMA_Z, sigma_delta=SIGMA_DELTA,
                n_scattering_events=10)

    def test_raises_on_wrong_gas_density_type(self, toy_ring):
        with pytest.raises(TypeError):
            make_study(toy_ring, gas_density={'s': [0.0], 'N': [1e15]})

    def test_raises_on_gas_density_without_species(self, toy_ring):
        with pytest.raises(ValueError):
            make_study(toy_ring, gas_density=xt.Table(
                {'name': np.array(['a', 'b']), 's': np.array([0.0, 1.0])}))

    @pytest.mark.parametrize('bad_value', [np.nan, np.inf, -1.0])
    def test_raises_on_bad_density_values(self, toy_ring, bad_value):
        gd = toy_ring['gas_density']
        with pytest.raises(ValueError):
            make_study(toy_ring, gas_density=xt.Table(
                {'name': gd.name, 's': gd.s,
                 'N': np.full(len(gd.s), bad_value)}))

    def test_raises_on_unknown_gas_species(self, toy_ring):
        gd = toy_ring['gas_density']
        with pytest.raises(ValueError):
            make_study(toy_ring, gas_density=xt.Table(
                {'name': gd.name, 's': gd.s,
                 'Zz': np.ones(len(gd.s))*ATOMIC_DENSITY}))

    def _study_with_particle_ref(self, toy_ring, particle_ref):
        line = toy_ring['line'].copy()
        line.particle_ref = particle_ref
        return xc.BeamGasStudy(
            line=line, gas_density=toy_ring['gas_density'],
            process='brems', nemitt_x=NEMITT_X, nemitt_y=NEMITT_Y,
            sigma_z=SIGMA_Z, sigma_delta=SIGMA_DELTA,
            n_scattering_events=10)

    def test_raises_on_non_lepton_beam(self, toy_ring):
        with pytest.raises(ValueError, match='electron and positron'):
            self._study_with_particle_ref(
                toy_ring, xt.Particles(pdg_id=2212, p0c=1e12))

    def test_raises_when_pdg_id_not_set(self, toy_ring):
        # An electron built without a PDG id: right mass, but the study
        # cannot tell an electron from a positron
        with pytest.raises(ValueError, match='no PDG id'):
            self._study_with_particle_ref(
                toy_ring, xt.Particles(mass0=xt.ELECTRON_MASS_EV, q0=-1.0,
                                       p0c=P0C))

    def test_raises_on_charge_inconsistent_with_pdg_id(self, toy_ring):
        with pytest.raises(ValueError, match='q0'):
            self._study_with_particle_ref(
                toy_ring, xt.Particles(pdg_id=11, mass0=xt.ELECTRON_MASS_EV,
                                       q0=+1.0, p0c=P0C))

    @pytest.mark.parametrize('name,pdg_id,q0', [('electron', 11, -1.0),
                                                ('positron', -11, 1.0)])
    def test_accepts_electrons_and_positrons(self, toy_ring, name, pdg_id, q0):
        line = toy_ring['line'].copy()
        line.set_particle_ref(name, p0c=P0C)
        study = xc.BeamGasStudy(
            line=line, gas_density=toy_ring['gas_density'],
            process='coulomb', nemitt_x=NEMITT_X, nemitt_y=NEMITT_Y,
            sigma_z=SIGMA_Z, sigma_delta=SIGMA_DELTA,
            n_scattering_events=10)
        assert study.pdg_id == pdg_id
        # q0 is taken from the PDG id, and drives the Mott term
        assert study.q0 == q0
        assert study.calculators['N'].q0 == q0

    def test_raises_when_no_beamgas_elements(self, toy_ring):
        line = xt.Line(elements=[xt.Drift(length=1.0)])
        line.set_particle_ref('electron', p0c=P0C)
        with pytest.raises(ValueError, match='does not contain'):
            xc.BeamGasStudy(
                line=line, gas_density=toy_ring['gas_density'],
                process='brems', nemitt_x=NEMITT_X, nemitt_y=NEMITT_Y,
                sigma_z=SIGMA_Z, sigma_delta=SIGMA_DELTA,
                n_scattering_events=10)


#############################################################
# Initialisation of the scattering elements
#############################################################
class TestInitialiseBeamGas:

    @pytest.fixture(autouse=True)
    def _study(self, toy_ring):
        self.toy_ring = toy_ring
        self.study = make_study(toy_ring)

    def test_sections_cover_the_whole_line(self):
        table = self.study.local_rates()
        assert np.isclose(table.ds.sum(), self.toy_ring['line'].get_length())

    def test_all_elements_are_configured(self):
        for nn in self.study.elements:
            elem = self.toy_ring['line'][nn]
            assert elem.p0c > 0
            assert elem.betx > 0 and elem.bety > 0
            assert elem.interaction_rate > 0
            assert elem.n_scattering_events == self.study.n_scattering_events
            assert dict(elem.atomic_densities) == {'N': pytest.approx(
                ATOMIC_DENSITY)}

    def test_interaction_rate_matches_analytic_formula(self):
        # rate = N_bunch * f_rev * integral(n ds) * sigma, and for a flat
        # profile the integral over the whole ring is n * circumference
        line = self.toy_ring['line']
        xsec = self.study.xsecs['N']
        beta0 = line.particle_ref.beta0[0]
        expected = BUNCH_INTENSITY*beta0*sc.c*ATOMIC_DENSITY*xsec
        total = sum(line[nn].interaction_rate for nn in self.study.elements)
        assert np.isclose(total, expected, rtol=1e-9)

    def test_mean_free_path(self):
        line = self.toy_ring['line']
        expected = 1.0/(ATOMIC_DENSITY*self.study.xsecs['N'])
        for nn in self.study.elements:
            assert np.isclose(line[nn].mean_free_path, expected)

    def test_partial_initialise_single_element(self):
        nn = self.study.elements[0]
        self.study.initialise_beamgas(element=nn, verbose=False)
        assert self.toy_ring['line'][nn].interaction_rate > 0

    def test_partial_initialise_raises_on_unknown_element(self):
        with pytest.raises(ValueError):
            self.study.initialise_beamgas(element='mqf.1', verbose=False)
        with pytest.raises(TypeError):
            self.study.initialise_beamgas(element=0, verbose=False)

    def test_warns_when_line_end_is_not_covered(self, toy_ring):
        study = xc.BeamGasStudy(
            line=toy_ring['line'], gas_density=toy_ring['gas_density'],
            process='brems', elements=['BeamGasScattering.0'],
            nemitt_x=NEMITT_X, nemitt_y=NEMITT_Y,
            sigma_z=SIGMA_Z, sigma_delta=SIGMA_DELTA,
            n_scattering_events=10, method='4d')
        with pytest.warns(UserWarning, match='not represented'):
            study.initialise_beamgas(verbose=False)


#############################################################
# Particle generation
#############################################################
class TestGenerateParticles:

    @pytest.fixture(autouse=True)
    def _study(self, toy_ring):
        self.toy_ring = toy_ring
        self.study = make_study(toy_ring, n_scattering_events=500)

    def test_weights_sum_to_the_interaction_rate(self):
        line = self.toy_ring['line']
        for nn, particles in self.study.generate_particles().items():
            allocated = particles.particle_id >= 0
            assert np.isclose(particles.weight[allocated].sum(),
                              line[nn].interaction_rate, rtol=1e-12)

    def test_every_particle_interacts_once(self):
        particles_by_element = self.study.generate_particles()
        log = self.study.interaction_log()
        n_generated = sum(int(np.sum(pp.particle_id >= 0))
                          for pp in particles_by_element.values())
        assert len(log.particle_id) == n_generated

    def test_generated_distribution_matches_the_optics(self, toy_ring):
        # With Coulomb scattering delta is untouched, so the generated
        # longitudinal distribution must be the matched one
        study = make_study(toy_ring, process='coulomb',
                           coulomb_theta=(1e-5, 1e-3),
                           n_scattering_events=20_000)
        nn = study.elements[0]
        particles = toy_ring['line'][nn].scatter(rng=study.rng)
        allocated = particles.particle_id >= 0
        elem = toy_ring['line'][nn]

        assert np.isclose(particles.delta[allocated].std(), SIGMA_DELTA,
                          rtol=0.05)
        assert np.isclose(particles.zeta[allocated].std(), SIGMA_Z, rtol=0.05)
        sigma_x = np.sqrt(study.gemitt_x*elem.betx
                          + (elem.dx*SIGMA_DELTA)**2)
        assert np.isclose(particles.x[allocated].std(), sigma_x, rtol=0.05)
        assert np.isclose(particles.x[allocated].mean(), elem.x_co,
                          atol=0.1*sigma_x)

    def test_particles_start_at_the_element(self):
        line = self.toy_ring['line']
        for nn, particles in self.study.generate_particles().items():
            allocated = particles.particle_id >= 0
            assert np.allclose(particles.s[allocated], line[nn].s)
            assert np.all(particles.at_element[allocated]
                          == line.element_names.index(nn))

    def test_raises_if_not_initialised(self, toy_ring):
        elem = xc.BeamGasScattering()
        with pytest.raises(ValueError, match='not been initialised'):
            elem.scatter()


#############################################################
# Running the study
#############################################################
class TestRun:

    @pytest.fixture(autouse=True)
    def _study(self, toy_ring):
        self.toy_ring = toy_ring
        self.study = make_study(toy_ring, process='coulomb',
                                coulomb_theta=(8e-3, 20e-3),
                                n_scattering_events=200)

    def test_run_without_tracking(self):
        result = self.study.run()
        assert result.tracked is False
        assert result.rate_tracking is None
        assert result.particles is None
        assert result.rate_scattering > 0
        assert np.isclose(result.lifetime_scattering,
                          BUNCH_INTENSITY/result.rate_scattering)

    def test_run_with_tracking(self):
        result = self.study.run(track=True, n_turns=20, keep_particles=True)
        assert result.tracked is True
        assert 0 < result.rate_tracking <= result.rate_scattering
        assert np.isclose(result.lifetime_tracking,
                          BUNCH_INTENSITY/result.rate_tracking)
        assert result.particles is not None
        assert result.lost_particles is not None
        assert np.all(result.lost_particles.state <= 0)
        assert np.all(result.lost_particles.particle_id >= 0)
        assert np.isclose(result.lost_particles.weight.sum(),
                          result.rate_tracking)

    def test_lifetime_is_independent_of_bunch_intensity(self, toy_ring):
        kwargs = dict(process='coulomb', coulomb_theta=(8e-3, 20e-3),
                      n_scattering_events=200, seed=7)
        a = make_study(toy_ring, **kwargs).run(track=True, n_turns=20)
        study_b = xc.BeamGasStudy(
            line=toy_ring['line'], gas_density=toy_ring['gas_density'],
            nemitt_x=NEMITT_X, nemitt_y=NEMITT_Y, sigma_z=SIGMA_Z,
            sigma_delta=SIGMA_DELTA, bunch_intensity=10*BUNCH_INTENSITY,
            method='4d', **kwargs)
        study_b.initialise_beamgas(verbose=False)
        b = study_b.run(track=True, n_turns=20)
        assert np.isclose(a.rate_tracking*10, b.rate_tracking)
        assert np.isclose(a.lifetime_tracking, b.lifetime_tracking)

    def test_local_rates_columns(self):
        result = self.study.run(track=True, n_turns=5, keep_particles=True)
        table = result.local_rates
        for col in ('name', 's', 'ds', 'mean_free_path', 'local_rate',
                    'interaction_rate', 'num_particles', 'sum_weight',
                    'num_lost_particles', 'sum_lost_weight'):
            assert col in table._col_names
        # The event weights are normalised to the configured cross section,
        # so this now agrees up to the Monte Carlo noise of the importance
        # sampling alone
        assert np.isclose(table.sum_weight.sum(), result.rate_scattering,
                          rtol=0.01)
        assert np.isclose(table.sum_lost_weight.sum(), result.rate_tracking)

    def test_interaction_log(self):
        result = self.study.run(track=True, n_turns=5, keep_particles=True)
        log = result.interaction_log
        assert set(np.unique(log.gas)) == {'N'}
        assert np.all(log.theta >= 8e-3) and np.all(log.theta <= 20e-3)
        assert np.all(np.isnan(log.photon_energy))  # Coulomb: no photon
        assert np.isclose(log.weight.sum(), result.rate_scattering, rtol=0.01)

    def test_reproducible_with_seed(self, toy_ring):
        kwargs = dict(process='coulomb', coulomb_theta=(8e-3, 20e-3),
                      n_scattering_events=200)
        a = make_study(toy_ring, seed=11, **kwargs).run(track=True, n_turns=10)
        b = make_study(toy_ring, seed=11, **kwargs).run(track=True, n_turns=10)
        c = make_study(toy_ring, seed=12, **kwargs).run(track=True, n_turns=10)
        assert a.rate_tracking == b.rate_tracking
        assert a.rate_tracking != c.rate_tracking

    def test_statistical_error(self):
        result = self.study.run(track=True, n_turns=20)
        assert result.rate_tracking_error > 0
        assert result.rate_tracking_error < result.rate_tracking
        assert np.isclose(
            result.lifetime_tracking_error/result.lifetime_tracking,
            result.rate_tracking_error/result.rate_tracking)
        assert self.study.run().rate_tracking_error is None

    def test_cutoff_scan(self):
        result = self.study.run(track=True, n_turns=20, keep_particles=True)
        scan = result.cutoff_scan
        assert np.isclose(scan.cut[0], 8e-3)
        assert np.all(np.diff(scan.cut) > 0)
        # The first row is the study itself: no collimators, so every lost
        # particle is a generated primary
        assert np.isclose(scan.rate_tracking[0], result.rate_tracking)
        assert np.isclose(scan.rate_tracking_error[0],
                          result.rate_tracking_error)
        assert scan.num_events.sum() == len(result.interaction_log.theta)
        assert scan.num_lost.sum() == len(result.lost_particles.state)
        # Raising the cut can only remove losses
        assert np.all(np.diff(scan.rate_tracking) <= 0)
        assert np.all(np.diff(scan.lifetime_tracking) >= 0)
        assert np.all((scan.loss_probability[scan.num_events > 0] >= 0)
                      & (scan.loss_probability[scan.num_events > 0] <= 1))
        # The rate above a cut is the lost weight generated above it
        lost_ids = set()
        log = result.interaction_log
        for nn, pp in result.particles_by_element.items():
            ids = pp.particle_id[(pp.particle_id >= 0) & (pp.state <= 0)]
            lost_ids |= {(nn, ii) for ii in ids}
        lost = np.array([(nn, ii) in lost_ids
                         for nn, ii in zip(log.name, log.particle_id)])
        for ii in (0, len(scan.cut)//2, len(scan.cut) - 1):
            above = log.theta >= scan.cut[ii]
            assert np.isclose(scan.rate_tracking[ii],
                              log.weight[above & lost].sum())

    def test_warns_when_the_window_truncates_the_losses(self):
        # The toy ring loses particles well below 8 mrad and virtually all
        # of them above 20 mrad, so both cuts of this window bias the result
        with pytest.warns(UserWarning) as record:
            result = self.study.run(track=True, n_turns=20)
        messages = [str(ww.message) for ww in record]
        assert any('lower cut `coulomb_theta[0]`' in mm for mm in messages)
        assert any('above `coulomb_theta[1]`' in mm for mm in messages)
        assert result.rate_above_theta_max > 0.01*result.rate_tracking

    def test_no_truncation_warning_for_a_wide_window(self, toy_ring):
        study = make_study(toy_ring, process='coulomb',
                           coulomb_theta=(1e-7, 1.0),
                           n_scattering_events=200)
        with warnings.catch_warnings():
            warnings.filterwarnings('error', message='.*coulomb_theta.*')
            result = study.run(track=True, n_turns=20)
        assert result.rate_above_theta_max < 1e-3*result.rate_tracking

    def test_brems_cutoff_scan(self, toy_ring):
        study = make_study(toy_ring, n_scattering_events=200)
        result = study.run(track=True, n_turns=5)
        assert np.isclose(result.cutoff_scan.cut[0], 1e6)
        assert np.isclose(result.cutoff_scan.rate_tracking[0],
                          result.rate_tracking)
        assert result.rate_above_theta_max is None

    def test_raises_on_inconsistent_run_arguments(self):
        with pytest.raises(ValueError):
            self.study.run(track=True)  # missing n_turns
        with pytest.raises(ValueError):
            self.study.run(track=True, n_turns=1, generate_particles=False)
        with pytest.raises(ValueError):
            self.study.run(keep_particles=True, generate_particles=False)


#############################################################
# line.xcoll.beamgas_configure facade
#############################################################
class TestLineFacade:

    def _kwargs(self, toy_ring):
        return dict(gas_density=toy_ring['gas_density'], process='coulomb',
                    coulomb_theta=(8e-3, 20e-3),
                    nemitt_x=NEMITT_X, nemitt_y=NEMITT_Y,
                    sigma_z=SIGMA_Z, sigma_delta=SIGMA_DELTA,
                    n_scattering_events=200, seed=1997, method='4d')

    def test_returns_an_initialised_study(self, toy_ring):
        study = toy_ring['line'].xcoll.beamgas_configure(
            verbose=False, **self._kwargs(toy_ring))
        assert isinstance(study, xc.BeamGasStudy)
        # initialise_beamgas() has already run, so the elements are configured
        assert len(study.elements) == len(
            [nn for nn in toy_ring['line'].element_names
             if isinstance(toy_ring['line'][nn], xc.BeamGasScattering)])
        for nn in study.elements:
            assert toy_ring['line'][nn].interaction_rate > 0
        assert study.twiss is not None

    def test_equivalent_to_direct_construction(self, toy_ring):
        kwargs = self._kwargs(toy_ring)
        via_facade = toy_ring['line'].xcoll.beamgas_configure(
            verbose=False, **kwargs)
        result_facade = via_facade.run(track=True, n_turns=20)

        direct = xc.BeamGasStudy(line=toy_ring['line'], **kwargs)
        direct.initialise_beamgas(verbose=False)
        result_direct = direct.run(track=True, n_turns=20)

        assert result_facade.rate_scattering == result_direct.rate_scattering
        assert result_facade.rate_tracking == result_direct.rate_tracking

    def test_forwards_element_selection(self, toy_ring):
        names = ['BeamGasScattering.0', 'BeamGasScattering.1']
        with pytest.warns(UserWarning, match='not represented'):
            study = toy_ring['line'].xcoll.beamgas_configure(
                elements=names, verbose=False, **self._kwargs(toy_ring))
        assert study.elements == names

    def test_propagates_validation_errors(self, toy_ring):
        kwargs = self._kwargs(toy_ring)
        kwargs['process'] = 'compton'
        with pytest.raises(ValueError, match='brems'):
            toy_ring['line'].xcoll.beamgas_configure(
                verbose=False, **kwargs)


#############################################################
# The element is passive during tracking
#############################################################
class TestPassiveTracking:

    def test_tracking_is_not_affected(self, toy_ring):
        line = toy_ring['line']
        particles = line.build_particles(x=[1e-4, 2e-4], px=[1e-6, -1e-6])
        reference = particles.copy()
        line.track(particles, num_turns=5)

        line_no_beamgas = line.copy()
        line_no_beamgas.discard_tracker()
        for nn in [nn for nn in line_no_beamgas.element_names
                   if isinstance(line_no_beamgas[nn], xc.BeamGasScattering)]:
            line_no_beamgas.element_dict[nn] = xt.Marker()
        line_no_beamgas.build_tracker()
        line_no_beamgas.track(reference, num_turns=5)

        assert np.allclose(particles.x, reference.x, rtol=0, atol=1e-15)
        assert np.allclose(particles.px, reference.px, rtol=0, atol=1e-15)

    def test_element_is_not_collective(self, toy_ring):
        assert not toy_ring['line']['BeamGasScattering.0'].iscollective
