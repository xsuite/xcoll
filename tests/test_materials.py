# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

import pytest
import numpy as np

import xtrack as xt
import xcoll as xc
from xcoll.materials import Material, CrystalMaterial, db as mdb


# TODO: make test_crystal_material_creation and test_adapt(), and expand test_db


def test_elemental_material_creation():
    with pytest.raises(ValueError, match="Invalid material definition! Use either `Z` and `A` for "
                       "elemental materials, or `components` and `n_atoms`, `mass_fractions`, "
                       "`volume_fractions`, `molar_fractions`, or `atomic_fractions` for compound materials."):
        Material()
    with pytest.raises(ValueError, match="A must be provided for an elemental Material"):
        Material(Z=4)
    with pytest.raises(ValueError, match="density must be provided for Material"):
        Material(Z=4, A=3)

    # Create a simple material
    mat = Material(A=12.01, Z=6, density=2.265, name='TestMAT')
    assert isinstance(mat, Material)
    assert mat.name == 'TestMAT'
    assert mat.short_name is None
    assert mat.fluka_name is None
    assert mat.geant4_name is None
    assert mat.temperature is None
    assert mat.pressure is None
    assert mat.state is None
    assert mat.is_elemental
    assert not mat.is_compound
    assert not mat.is_mixture
    assert not mat.full_everest_supported
    assert np.isclose(mat.A, 12.01)
    assert np.isclose(mat.Z, 6)
    assert np.isclose(mat.density, 2.265)
    assert mat.components is None
    assert mat.n_atoms is None
    assert mat.mass_fractions is None
    assert mat.volume_fractions is None
    assert mat.molar_fractions is None
    assert mat.atomic_fractions is None
    assert mat.composition is None
    assert np.isclose(mat.molar_mass, 12.01)
    assert mat.average_molar_mass is None
    assert str(mat) == 'Elemental Material(TestMAT, Z=6, A=12.010, density=2.2650 g/cm^3)'
    assert np.isclose(mat._ZA_mean, 0.4995836802664446)
    assert np.isclose(mat._Z2_eff, 36.0)
    assert np.isclose(mat.radiation_length, 0.18849567315806856)
    assert np.isclose(mat.excitation_energy, 81.0)
    assert mat._radiation_length_set_manually is False
    assert mat._excitation_energy_set_manually is False
    assert mat.nuclear_radius is None
    assert np.isclose(mat._nuclear_radius, -1.)
    assert mat.nuclear_elastic_slope is None
    assert np.isclose(mat._nuclear_elastic_slope, -1.)
    assert mat.cross_section is None
    assert np.all(mat._cross_section == np.array([-1.]*6))
    assert mat.hcut is None
    assert np.isclose(mat._hcut, -1.)
    assert np.isclose(mat.electron_density, 3.0085632439633638e+23)
    assert np.isclose(mat.plasma_energy, 30.652924826509626)
    assert np.isclose(mat.atoms_per_volume, 1.1357326245961699e+29)
    assert np.isclose(mat._atoms_per_volume, mat.atoms_per_volume)
    assert np.isclose(mat.num_effective_nucleons, 3.704356404188558)
    assert np.isclose(mat._num_nucleons_eff, mat.num_effective_nucleons)
    assert xt.line._dicts_equal(mat.to_dict(), {'__class__': 'Material', 'name': 'TestMAT',
                                                'A': 12.01, 'Z': 6, 'density': 2.265})
    assert xt.line._dicts_equal(mat._xobject._to_dict(), {
        '_density': np.float64(2.265),
        '_ZA_mean': np.float64(0.4995836802664446),
        '_Z2_eff': np.float64(36.0),
        '_atoms_per_volume': np.float64(1.1357326245961699e+29),
        '_num_nucleons_eff': np.float64(3.704356404188558),
        '_radiation_length': np.float64(0.18849567315806856),
        '_excitation_energy': np.float64(81.0),
        '_nuclear_radius': np.float64(-1.0),
        '_nuclear_elastic_slope': np.float64(-1.0),
        '_cross_section': [np.float64(-1.0),
        np.float64(-1.0),
        np.float64(-1.0),
        np.float64(-1.0),
        np.float64(-1.0),
        np.float64(-1.0)],
        '_hcut': np.float64(-1.0)
    })

    # Remove Z
    mat.Z = None
    assert mat.Z is None
    assert np.isclose(mat.A, 12.01)
    assert np.isclose(mat._ZA_mean, -1)
    assert np.isclose(mat._Z2_eff, -1)
    assert np.isclose(mat.density, 2.265)
    assert mat.radiation_length is None
    assert np.isclose(mat._radiation_length, -1)
    assert mat.excitation_energy is None
    assert np.isclose(mat._excitation_energy, -1)
    assert mat.electron_density is None
    assert mat.plasma_energy is None
    assert np.isclose(mat.atoms_per_volume, 1.1357326245961699e+29)
    assert np.isclose(mat._atoms_per_volume, mat.atoms_per_volume)
    assert np.isclose(mat.num_effective_nucleons, 3.704356404188558)
    assert np.isclose(mat._num_nucleons_eff, mat.num_effective_nucleons)

    # Restore Z
    mat.Z = 6
    assert np.isclose(mat._Z, 6)
    assert np.isclose(mat.A, 12.01)
    assert np.isclose(mat._ZA_mean, 0.4995836802664446)
    assert np.isclose(mat._Z2_eff, 36)
    assert np.isclose(mat.density, 2.265)
    assert np.isclose(mat.radiation_length, 0.18849567315806856)
    assert np.isclose(mat.excitation_energy, 81.0)
    assert np.isclose(mat.electron_density, 3.0085632439633638e+23)
    assert np.isclose(mat.plasma_energy, 30.652924826509626)
    assert np.isclose(mat.atoms_per_volume, 1.1357326245961699e+29)
    assert np.isclose(mat._atoms_per_volume, mat.atoms_per_volume)
    assert np.isclose(mat.num_effective_nucleons, 3.704356404188558)
    assert np.isclose(mat._num_nucleons_eff, mat.num_effective_nucleons)

    # Change Z
    mat.Z = 5
    assert np.isclose(mat._Z, 5)
    assert np.isclose(mat.A, 12.01)
    assert np.isclose(mat._ZA_mean, 0.41631973355537055)
    assert np.isclose(mat._Z2_eff, 25)
    assert np.isclose(mat.density, 2.265)
    assert np.isclose(mat.radiation_length, 0.25840988913517837)
    assert np.isclose(mat.excitation_energy, 76.0)
    assert np.isclose(mat.electron_density, 2.5071360366361366e+23)
    assert np.isclose(mat.plasma_energy, 27.982163968315756)
    assert np.isclose(mat.atoms_per_volume, 1.1357326245961699e+29)
    assert np.isclose(mat._atoms_per_volume, mat.atoms_per_volume)
    assert np.isclose(mat.num_effective_nucleons, 3.704356404188558)
    assert np.isclose(mat._num_nucleons_eff, mat.num_effective_nucleons)

    # Remove A
    mat.A = None
    assert mat.A is None
    assert np.isclose(mat.Z, 5)
    assert np.isclose(mat._ZA_mean, -1)
    assert np.isclose(mat._Z2_eff, 25)
    assert np.isclose(mat.density, 2.265)
    assert mat.radiation_length is None
    assert np.isclose(mat._radiation_length, -1)
    assert np.isclose(mat.excitation_energy, 76.0)
    assert mat.electron_density is None
    assert mat.plasma_energy is None
    assert mat.atoms_per_volume is None
    assert np.isclose(mat._atoms_per_volume, -1)
    assert mat.num_effective_nucleons is None
    assert np.isclose(mat._num_nucleons_eff, -1)

    # Restore A
    mat.A = 12.01
    assert np.isclose(mat._A, 12.01)
    assert np.isclose(mat.Z, 5)
    assert np.isclose(mat._ZA_mean, 0.41631973355537055)
    assert np.isclose(mat._Z2_eff, 25)
    assert np.isclose(mat.density, 2.265)
    assert np.isclose(mat.radiation_length, 0.25840988913517837)
    assert np.isclose(mat.excitation_energy, 76.0)
    assert np.isclose(mat.radiation_length, 0.25840988913517837)
    assert np.isclose(mat.excitation_energy, 76.0)
    assert np.isclose(mat.electron_density, 2.5071360366361366e+23)
    assert np.isclose(mat.plasma_energy, 27.982163968315756)
    assert np.isclose(mat.atoms_per_volume, 1.1357326245961699e+29)
    assert np.isclose(mat._atoms_per_volume, mat.atoms_per_volume)
    assert np.isclose(mat.num_effective_nucleons, 3.704356404188558)
    assert np.isclose(mat._num_nucleons_eff, mat.num_effective_nucleons)

    # Change A
    mat.A = 14.3
    assert np.isclose(mat._A, 14.3)
    assert np.isclose(mat.Z, 5)
    assert np.isclose(mat._ZA_mean, 0.34965034965034963)
    assert np.isclose(mat._Z2_eff, 25)
    assert np.isclose(mat.density, 2.265)
    assert np.isclose(mat.radiation_length, 0.3076820495114947)
    assert np.isclose(mat.excitation_energy, 76.0)
    assert np.isclose(mat.electron_density, 2.1056436223776224e+23)
    assert np.isclose(mat.plasma_energy, 25.643941771779268)
    assert np.isclose(mat.atoms_per_volume, 9.538565609370629e+28)
    assert np.isclose(mat._atoms_per_volume, mat.atoms_per_volume)
    assert np.isclose(mat.num_effective_nucleons, 3.926242391361843)
    assert np.isclose(mat._num_nucleons_eff, mat.num_effective_nucleons)

    # Remove density
    mat.density = None
    assert np.isclose(mat._density, -1)
    assert np.isclose(mat.A, 14.3)
    assert np.isclose(mat.Z, 5)
    assert np.isclose(mat._ZA_mean, 0.34965034965034963)
    assert np.isclose(mat._Z2_eff, 25)
    assert mat.radiation_length is None
    assert np.isclose(mat._radiation_length, -1)
    assert np.isclose(mat.excitation_energy, 76.0)
    assert np.isclose(mat.electron_density, 2.1056436223776224e+23)
    assert mat.plasma_energy is None
    assert mat.atoms_per_volume is None
    assert np.isclose(mat._atoms_per_volume, -1)
    assert np.isclose(mat.num_effective_nucleons, 3.926242391361843)
    assert np.isclose(mat._num_nucleons_eff, mat.num_effective_nucleons)

    # Restore density
    mat.density = 2.265
    assert np.isclose(mat._density, 2.265)
    assert np.isclose(mat.A, 14.3)
    assert np.isclose(mat.Z, 5)
    assert np.isclose(mat._ZA_mean, 0.34965034965034963)
    assert np.isclose(mat._Z2_eff, 25)
    assert np.isclose(mat.radiation_length, 0.3076820495114947)
    assert np.isclose(mat.excitation_energy, 76.0)
    assert np.isclose(mat.electron_density, 2.1056436223776224e+23)
    assert np.isclose(mat.plasma_energy, 25.643941771779268)
    assert np.isclose(mat.atoms_per_volume, 9.538565609370629e+28)
    assert np.isclose(mat._atoms_per_volume, mat.atoms_per_volume)
    assert np.isclose(mat.num_effective_nucleons, 3.926242391361843)
    assert np.isclose(mat._num_nucleons_eff, mat.num_effective_nucleons)

    # Change density
    mat.density = 4.345
    assert np.isclose(mat._density, 4.345)
    assert np.isclose(mat.A, 14.3)
    assert np.isclose(mat.Z, 5)
    assert np.isclose(mat._ZA_mean, 0.34965034965034963)
    assert np.isclose(mat._Z2_eff, 25)
    assert np.isclose(mat.radiation_length, 0.16039121798470324)
    assert np.isclose(mat.excitation_energy, 76.0)
    assert np.isclose(mat.electron_density, 2.1056436223776224e+23)
    assert np.isclose(mat.plasma_energy, 35.51776008183468)
    assert np.isclose(mat.atoms_per_volume, 1.8298043078461536e+29)
    assert np.isclose(mat._atoms_per_volume, mat.atoms_per_volume)
    assert np.isclose(mat.num_effective_nucleons, 3.926242391361843)
    assert np.isclose(mat._num_nucleons_eff, mat.num_effective_nucleons)

    # Set radiation length manually
    mat.radiation_length = 0.5
    assert np.isclose(mat._A, 14.3)
    assert np.isclose(mat.Z, 5)
    assert np.isclose(mat._ZA_mean, 0.34965034965034963)
    assert np.isclose(mat._Z2_eff, 25)
    assert np.isclose(mat.density, 4.345)
    assert np.isclose(mat.radiation_length, 0.5)
    assert np.isclose(mat.excitation_energy, 76.0)
    assert mat._radiation_length_set_manually is True
    assert mat._excitation_energy_set_manually is False

    # Remove Z again - radiation length should stay
    mat.Z = None
    assert mat.Z is None
    assert np.isclose(mat.A, 14.3)
    assert np.isclose(mat._ZA_mean, -1)
    assert np.isclose(mat._Z2_eff, -1)
    assert np.isclose(mat.density, 4.345)
    assert np.isclose(mat.radiation_length, 0.5)
    assert mat.excitation_energy is None
    assert np.isclose(mat._excitation_energy, -1)

    # Restore Z again - radiation length should stay
    mat.Z = 5
    assert np.isclose(mat._Z, 5)
    assert np.isclose(mat.A, 14.3)
    assert np.isclose(mat._ZA_mean, 0.34965034965034963)
    assert np.isclose(mat._Z2_eff, 25)
    assert np.isclose(mat.density, 4.345)
    assert np.isclose(mat.radiation_length, 0.5)
    assert np.isclose(mat.excitation_energy, 76.0)
    assert mat._radiation_length_set_manually is True
    assert mat._excitation_energy_set_manually is False

    # Set excitation energy manually
    mat.excitation_energy = 897
    assert np.isclose(mat.Z, 5)
    assert np.isclose(mat.A, 14.3)
    assert np.isclose(mat._ZA_mean, 0.34965034965034963)
    assert np.isclose(mat._Z2_eff, 25)
    assert np.isclose(mat.density, 4.345)
    assert np.isclose(mat.radiation_length, 0.5)
    assert np.isclose(mat.excitation_energy, 897.0)
    assert mat._radiation_length_set_manually is True
    assert mat._excitation_energy_set_manually is True

    # Remove Z again - radiation length and excitation energy should stay
    mat.Z = None
    assert mat.Z is None
    assert np.isclose(mat.A, 14.3)
    assert np.isclose(mat._ZA_mean, -1)
    assert np.isclose(mat._Z2_eff, -1)
    assert np.isclose(mat.density, 4.345)
    assert np.isclose(mat.radiation_length, 0.5)
    assert np.isclose(mat.excitation_energy, 897.0)
    assert mat._radiation_length_set_manually is True
    assert mat._excitation_energy_set_manually is True

    # Restore Z again - radiation length and excitation energy should stay
    mat.Z = 5
    assert np.isclose(mat._Z, 5)
    assert np.isclose(mat.A, 14.3)
    assert np.isclose(mat._ZA_mean, 0.34965034965034963)
    assert np.isclose(mat._Z2_eff, 25)
    assert np.isclose(mat.density, 4.345)
    assert np.isclose(mat.radiation_length, 0.5)
    assert np.isclose(mat.excitation_energy, 897.0)
    assert mat._radiation_length_set_manually is True
    assert mat._excitation_energy_set_manually is True

    # Remove radiation length - should revert to calculated value
    mat.radiation_length = None
    assert np.isclose(mat.Z, 5)
    assert np.isclose(mat.A, 14.3)
    assert np.isclose(mat._ZA_mean, 0.34965034965034963)
    assert np.isclose(mat._Z2_eff, 25)
    assert np.isclose(mat.density, 4.345)
    assert np.isclose(mat.radiation_length, 0.16039121798470324)
    assert np.isclose(mat.excitation_energy, 897.0)
    assert mat._radiation_length_set_manually is False
    assert mat._excitation_energy_set_manually is True

    # Remove excitation energy - should revert to calculated value
    mat.excitation_energy = None
    assert np.isclose(mat.Z, 5)
    assert np.isclose(mat.A, 14.3)
    assert np.isclose(mat._ZA_mean, 0.34965034965034963)
    assert np.isclose(mat._Z2_eff, 25)
    assert np.isclose(mat.density, 4.345)
    assert np.isclose(mat.radiation_length, 0.16039121798470324)
    assert np.isclose(mat.excitation_energy, 76.0)
    assert mat._radiation_length_set_manually is False
    assert mat._excitation_energy_set_manually is False


def test_compound_material_creation():
    with pytest.raises(ValueError, match="One of `n_atoms`, `mass_fractions`, `volume_fractions`, "
                       "or `molar_fractions` must be provided"):
        Material(components=['C', 'H', 'O'])
    with pytest.raises(ValueError, match="Variable `components` must be provided"):
        Material(n_atoms=[2, 6, 1])

    # Create a compound material
    mat = Material(components=['C', 'H', 'C', 'H', 'O', 'H'], n_atoms=[1, 3, 1, 2, 1, 1], density=0.78945,
                   name='Ethanol', state='liquid', temperature=293.15)
    assert isinstance(mat, Material)
    assert mat.name == 'Ethanol'
    assert mat.short_name is None
    assert mat.fluka_name is None
    assert mat.geant4_name is None
    assert np.isclose(mat.temperature, 293.15)
    assert mat.pressure is None
    assert mat.state == 'liquid'
    assert not mat.is_elemental
    assert mat.is_compound
    assert not mat.is_mixture
    assert not mat.full_everest_supported
    assert mat.A is None
    assert mat.Z is None
    assert np.isclose(mat.density, 0.78945)
    assert all(mat.components == [mdb['C'], mdb['H'], mdb['O']])
    assert all(mat.n_atoms == [2, 6, 1])
    assert mat._mass_fractions is None
    assert np.allclose(mat.mass_fractions, [0.52143524, 0.13128134, 0.34728342])
    assert np.allclose(mat.volume_fractions, [1.35111388e-04, 8.57234650e-01, 1.42630238e-01])
    assert np.allclose(mat.molar_fractions, [0.22222222, 0.66666667, 0.11111111])
    assert np.allclose(mat.atomic_fractions, [0.22222222, 0.66666667, 0.11111111])
    assert np.isclose(mat.molar_mass, 46.069)
    assert mat.average_molar_mass is None
    assert mat.composition == [['C', 2.0], ['H', 6.0], ['O', 1.0]]
    assert str(mat) == 'Compound Material(Ethanol, [C: 2, H: 6, O: 1], density=0.7894 g/cm^3)'
    assert np.isclose(mat._ZA_mean, 0.5643708350517702)
    assert np.isclose(mat._Z2_eff, 15.777777777777777)
    assert np.isclose(mat.radiation_length, 0.5183349281307253)
    assert np.isclose(mat.excitation_energy, 61.02615629345832)
    assert mat._radiation_length_set_manually is False
    assert mat._excitation_energy_set_manually is False
    assert mat.nuclear_radius is None
    assert np.isclose(mat._nuclear_radius, -1.)
    assert mat.nuclear_elastic_slope is None
    assert np.isclose(mat._nuclear_elastic_slope, -1.)
    assert mat.cross_section is None
    assert np.all(mat._cross_section == np.array([-1.]*6))
    assert mat.hcut is None
    assert np.isclose(mat._hcut, -1.)
    assert np.isclose(mat.electron_density, 3.398720609520502e+23)
    assert np.isclose(mat.plasma_energy, 19.23438440668927)
    assert np.isclose(mat.atoms_per_volume, 1.0319692250715231e+28)
    assert np.isclose(mat._atoms_per_volume, mat.atoms_per_volume)
    assert np.isclose(mat.num_effective_nucleons, 2.3573500951095088)
    assert np.isclose(mat._num_nucleons_eff, mat.num_effective_nucleons)
    assert xt.line._dicts_equal(mat.to_dict(), {'__class__': 'Material',
        'components': [mdb['C'], mdb['H'], mdb['O']],
        'temperature': 293.15,
        'name': 'Ethanol',
        'n_atoms': [2, 6, 1],
        'state': 'liquid',
        'density': 0.78945})
    assert xt.line._dicts_equal(mat._xobject._to_dict(), {
        '_density': np.float64(0.78945),
        '_ZA_mean': np.float64(0.5643708350517702),
        '_Z2_eff': np.float64(15.777777777777777),
        '_atoms_per_volume': np.float64(1.0319692250715231e+28),
        '_num_nucleons_eff': np.float64(2.357350095109508),
        '_radiation_length': np.float64(0.5183349281307253),
        '_excitation_energy': np.float64(61.02615629345832),
        '_nuclear_radius': np.float64(-1.0),
        '_nuclear_elastic_slope': np.float64(-1.0),
        '_cross_section': [np.float64(-1.0),
        np.float64(-1.0),
        np.float64(-1.0),
        np.float64(-1.0),
        np.float64(-1.0),
        np.float64(-1.0)],
        '_hcut': np.float64(-1.0)
    })

    # Remove density
    mat.density = None
    assert np.isclose(mat._density, -1)
    assert all(mat.components == [mdb['C'], mdb['H'], mdb['O']])
    assert all(mat.n_atoms == [2, 6, 1])
    assert np.allclose(mat.mass_fractions, [0.52143524, 0.13128134, 0.34728342])
    assert np.isclose(mat._ZA_mean, 0.5643708350517702)
    assert np.isclose(mat._Z2_eff, 15.777777777777777)
    assert mat.radiation_length is None
    assert np.isclose(mat._radiation_length, -1)
    assert np.isclose(mat.excitation_energy, 61.02615629345832)
    assert mat.plasma_energy is None
    assert mat.atoms_per_volume is None
    assert np.isclose(mat._atoms_per_volume, -1)
    assert np.isclose(mat.num_effective_nucleons, 2.3573500951095088)
    assert np.isclose(mat._num_nucleons_eff, mat.num_effective_nucleons)

    # Restore density
    mat.density = 0.78945
    assert np.isclose(mat.density, 0.78945)
    assert all(mat.components == [mdb['C'], mdb['H'], mdb['O']])
    assert all(mat.n_atoms == [2, 6, 1])
    assert np.allclose(mat.mass_fractions, [0.52143524, 0.13128134, 0.34728342])
    assert np.isclose(mat._ZA_mean, 0.5643708350517702)
    assert np.isclose(mat._Z2_eff, 15.777777777777777)
    assert np.isclose(mat.radiation_length, 0.5183349281307253)
    assert np.isclose(mat.excitation_energy, 61.02615629345832)
    assert np.isclose(mat.electron_density, 3.398720609520502e+23)
    assert np.isclose(mat.plasma_energy, 19.23438440668927)
    assert np.isclose(mat.atoms_per_volume, 1.0319692250715231e+28)
    assert np.isclose(mat._atoms_per_volume, mat.atoms_per_volume)
    assert np.isclose(mat.num_effective_nucleons, 2.3573500951095088)
    assert np.isclose(mat._num_nucleons_eff, mat.num_effective_nucleons)

    # Change density
    mat.density = 4.345
    assert np.isclose(mat._density, 4.345)
    assert all(mat.components == [mdb['C'], mdb['H'], mdb['O']])
    assert all(mat.n_atoms == [2, 6, 1])
    assert np.allclose(mat.mass_fractions, [0.52143524, 0.13128134, 0.34728342])
    assert np.isclose(mat._ZA_mean, 0.5643708350517702)
    assert np.isclose(mat._Z2_eff, 15.777777777777777)
    assert np.isclose(mat.radiation_length, 0.09417710218936735)
    assert np.isclose(mat.excitation_energy, 61.02615629345832)
    assert np.isclose(mat.electron_density, 3.398720609520502e+23)
    assert np.isclose(mat.plasma_energy, 45.12434045076062)
    assert np.isclose(mat.atoms_per_volume, 5.679785018602531e+28)
    assert np.isclose(mat._atoms_per_volume, mat.atoms_per_volume)
    assert np.isclose(mat.num_effective_nucleons, 2.357350095109508)
    assert np.isclose(mat._num_nucleons_eff, mat.num_effective_nucleons)

    # Set radiation length manually
    mat.radiation_length = 0.5
    assert np.isclose(mat._density, 4.345)
    assert all(mat.components == [mdb['C'], mdb['H'], mdb['O']])
    assert all(mat.n_atoms == [2, 6, 1])
    assert np.allclose(mat.mass_fractions, [0.52143524, 0.13128134, 0.34728342])
    assert np.isclose(mat._ZA_mean, 0.5643708350517702)
    assert np.isclose(mat._Z2_eff, 15.777777777777777)
    assert np.isclose(mat.radiation_length, 0.5)
    assert np.isclose(mat.excitation_energy, 61.02615629345832)
    assert np.isclose(mat.electron_density, 3.398720609520502e+23)
    assert np.isclose(mat.plasma_energy, 45.12434045076062)
    assert np.isclose(mat.atoms_per_volume, 5.679785018602531e+28)
    assert np.isclose(mat._atoms_per_volume, mat.atoms_per_volume)
    assert np.isclose(mat.num_effective_nucleons, 2.357350095109508)
    assert np.isclose(mat._num_nucleons_eff, mat.num_effective_nucleons)
    assert mat._radiation_length_set_manually is True
    assert mat._excitation_energy_set_manually is False

    # Set excitation energy manually
    mat.excitation_energy = 897
    assert np.isclose(mat._density, 4.345)
    assert all(mat.components == [mdb['C'], mdb['H'], mdb['O']])
    assert all(mat.n_atoms == [2, 6, 1])
    assert np.allclose(mat.mass_fractions, [0.52143524, 0.13128134, 0.34728342])
    assert np.isclose(mat._ZA_mean, 0.5643708350517702)
    assert np.isclose(mat._Z2_eff, 15.777777777777777)
    assert np.isclose(mat.radiation_length, 0.5)
    assert np.isclose(mat.excitation_energy, 897)
    assert np.isclose(mat.electron_density, 3.398720609520502e+23)
    assert np.isclose(mat.plasma_energy, 45.12434045076062)
    assert np.isclose(mat.atoms_per_volume, 5.679785018602531e+28)
    assert np.isclose(mat._atoms_per_volume, mat.atoms_per_volume)
    assert np.isclose(mat.num_effective_nucleons, 2.357350095109508)
    assert np.isclose(mat._num_nucleons_eff, mat.num_effective_nucleons)
    assert mat._radiation_length_set_manually is True
    assert mat._excitation_energy_set_manually is True

    # Remove radiation length - should revert to calculated value
    mat.radiation_length = None
    assert np.isclose(mat._density, 4.345)
    assert all(mat.components == [mdb['C'], mdb['H'], mdb['O']])
    assert all(mat.n_atoms == [2, 6, 1])
    assert np.allclose(mat.mass_fractions, [0.52143524, 0.13128134, 0.34728342])
    assert np.isclose(mat._ZA_mean, 0.5643708350517702)
    assert np.isclose(mat._Z2_eff, 15.777777777777777)
    assert np.isclose(mat.radiation_length, 0.09417710218936735)
    assert np.isclose(mat.excitation_energy, 897)
    assert np.isclose(mat.electron_density, 3.398720609520502e+23)
    assert np.isclose(mat.plasma_energy, 45.12434045076062)
    assert np.isclose(mat.atoms_per_volume, 5.679785018602531e+28)
    assert np.isclose(mat._atoms_per_volume, mat.atoms_per_volume)
    assert np.isclose(mat.num_effective_nucleons, 2.357350095109508)
    assert np.isclose(mat._num_nucleons_eff, mat.num_effective_nucleons)
    assert mat._radiation_length_set_manually is False
    assert mat._excitation_energy_set_manually is True

    # Remove excitation energy - should revert to calculated value
    mat.excitation_energy = None
    assert np.isclose(mat._density, 4.345)
    assert all(mat.components == [mdb['C'], mdb['H'], mdb['O']])
    assert all(mat.n_atoms == [2, 6, 1])
    assert np.allclose(mat.mass_fractions, [0.52143524, 0.13128134, 0.34728342])
    assert np.isclose(mat._ZA_mean, 0.5643708350517702)
    assert np.isclose(mat._Z2_eff, 15.777777777777777)
    assert np.isclose(mat.radiation_length, 0.09417710218936735)
    assert np.isclose(mat.excitation_energy, 61.02615629345832)
    assert np.isclose(mat.electron_density, 3.398720609520502e+23)
    assert np.isclose(mat.plasma_energy, 45.12434045076062)
    assert np.isclose(mat.atoms_per_volume, 5.679785018602531e+28)
    assert np.isclose(mat._atoms_per_volume, mat.atoms_per_volume)
    assert np.isclose(mat.num_effective_nucleons, 2.357350095109508)
    assert np.isclose(mat._num_nucleons_eff, mat.num_effective_nucleons)
    assert mat._radiation_length_set_manually is False
    assert mat._excitation_energy_set_manually is False


def test_mixture_material_creation():
    with pytest.raises(ValueError, match="Variable `components` must be provided"):
        Material(mass_fractions=[0.01, 0.001, 0.529107, 0.016, 0.002, 0.033872, 0.337021, 0.013, 0.044, 0.014])

    # Create a compound material
    mat = Material(components=['H', 'C', 'O', 'Na', 'Mg', 'Al', 'Si', 'K', 'Ca', 'Fe'],
                   mass_fractions=[0.01, 0.001, 0.529107, 0.016, 0.002, 0.033872, 0.337021, 0.013, 0.044, 0.014],
                   name='Concrete', density=2.35, state='solid')
    assert isinstance(mat, Material)
    assert mat.name == 'Concrete'
    assert mat.short_name is None
    assert mat.fluka_name is None
    assert mat.geant4_name is None
    assert mat.temperature is None
    assert mat.pressure is None
    assert mat.state == 'solid'
    assert not mat.is_elemental
    assert not mat.is_compound
    assert mat.is_mixture
    assert not mat.full_everest_supported
    assert mat.A is None
    assert mat.Z is None
    assert np.isclose(mat.density, 2.35)
    assert all(mat.components == [mdb[nn] for nn in ['H', 'C', 'O', 'Na', 'Mg', 'Al', 'Si', 'K', 'Ca', 'Fe']])
    assert mat.n_atoms is None
    assert np.allclose(mat.mass_fractions, [0.01, 0.001, 0.529107, 0.016, 0.002, 0.033872, 0.337021, 0.013, 0.044, 0.014])
    assert np.allclose(mat.volume_fractions, [2.30951391e-01, 9.16464063e-07, 7.68590345e-01, 3.43105967e-05,
                                              2.38871243e-06, 2.60411838e-05, 3.00302710e-04, 3.13054343e-05,
                                              5.93083172e-05, 3.69076396e-06])
    assert np.allclose(mat.molar_fractions, [0.16874747, 0.00141618, 0.56253359, 0.01183801, 0.00139969,
                                             0.02135328, 0.20411754, 0.0056557 , 0.0186743 , 0.00426424])
    assert np.allclose(mat.atomic_fractions, mat.molar_fractions)
    assert mat.molar_mass is None
    assert np.isclose(mat.average_molar_mass, 17.009744521602588)
    assert mat.composition == [['H', 0.01], ['C', 0.001], ['O', 0.529107], ['Na', 0.016], ['Mg', 0.002],
                               ['Al', 0.033872], ['Si', 0.337021], ['K', 0.013], ['Ca', 0.044], ['Fe', 0.014]]
    assert str(mat) == 'Mixture Material(Concrete, [H: 0.01, C: 0.001, O: 0.52911, Na: 0.016, Mg: 0.002, Al: 0.03387, Si: 0.33702, K: 0.013, Ca: 0.044, Fe: 0.014], density=2.3500 g/cm^3)'
    assert np.isclose(mat._ZA_mean, 0.5027459683349229)
    assert np.isclose(mat._Z2_eff, 93.86563197810295)
    assert np.isclose(mat.radiation_length, 0.11306709197802262)
    assert np.isclose(mat.excitation_energy, 121.74338267324313)
    assert mat._radiation_length_set_manually is False
    assert mat._excitation_energy_set_manually is False
    assert mat.nuclear_radius is None
    assert np.isclose(mat._nuclear_radius, -1.)
    assert mat.nuclear_elastic_slope is None
    assert np.isclose(mat._nuclear_elastic_slope, -1.)
    assert mat.cross_section is None
    assert np.all(mat._cross_section == np.array([-1.]*6))
    assert mat.hcut is None
    assert np.isclose(mat._hcut, -1.)
    assert np.isclose(mat.electron_density, 3.0276069878354084e+23)
    assert np.isclose(mat.plasma_energy, 31.321454741172133)
    assert np.isclose(mat.atoms_per_volume, 8.319954934083074e+28)
    assert np.isclose(mat._atoms_per_volume, mat.atoms_per_volume)
    assert np.isclose(mat.num_effective_nucleons, 3.9008141200039312)
    assert np.isclose(mat._num_nucleons_eff, mat.num_effective_nucleons)
    assert xt.line._dicts_equal(mat.to_dict(), {'__class__': 'Material',
        'components': [mdb[nn] for nn in ['H', 'C', 'O', 'Na', 'Mg', 'Al', 'Si', 'K', 'Ca', 'Fe']],
        'name': 'Concrete',
        'mass_fractions': [0.01, 0.001, 0.529107, 0.016, 0.002, 0.033872,
                            0.337021, 0.013, 0.044, 0.014],
        'state': 'solid',
        'density': 2.35})
    assert xt.line._dicts_equal(mat._xobject._to_dict(), {
        '_density': np.float64(2.35),
        '_ZA_mean': np.float64(0.5027459683349229),
        '_Z2_eff': np.float64(93.86563197810295),
        '_atoms_per_volume': np.float64(8.319954934083074e+28),
        '_num_nucleons_eff': np.float64(3.9008141200039312),
        '_radiation_length': np.float64(0.11306709197802262),
        '_excitation_energy': np.float64(121.74338267324313),
        '_nuclear_radius': np.float64(-1.0),
        '_nuclear_elastic_slope': np.float64(-1.0),
        '_cross_section': [np.float64(-1.0),
        np.float64(-1.0),
        np.float64(-1.0),
        np.float64(-1.0),
        np.float64(-1.0),
        np.float64(-1.0)],
        '_hcut': np.float64(-1.0)
    })


def test_different_fractions():
    Ethanol_v0 = Material(components=['C', 'H', 'O'], n_atoms=[2, 6, 1], density=0.78945, name='Ethanol_v0')
    mf = Ethanol_v0.mass_fractions
    vf = Ethanol_v0.volume_fractions
    molf = Ethanol_v0.molar_fractions
    Ethanol_v1 = Material(components=['C', 'H', 'O'], mass_fractions=mf, density=0.78945, name='Ethanol_v1')
    assert np.allclose(Ethanol_v0.mass_fractions, Ethanol_v1.mass_fractions)
    assert np.allclose(Ethanol_v0.volume_fractions, Ethanol_v1.volume_fractions)
    assert np.allclose(Ethanol_v0.molar_fractions, Ethanol_v1.molar_fractions)
    Ethanol_v2 = Material(components=['C', 'H', 'O'], volume_fractions=vf, density=0.78945, name='Ethanol_v2')
    assert np.allclose(Ethanol_v0.mass_fractions, Ethanol_v2.mass_fractions)
    assert np.allclose(Ethanol_v0.volume_fractions, Ethanol_v2.volume_fractions)
    assert np.allclose(Ethanol_v0.molar_fractions, Ethanol_v2.molar_fractions)
    Ethanol_v3 = Material(components=['C', 'H', 'O'], molar_fractions=molf, density=0.78945, name='Ethanol_v3')
    assert np.allclose(Ethanol_v0.mass_fractions, Ethanol_v3.mass_fractions)
    assert np.allclose(Ethanol_v0.volume_fractions, Ethanol_v3.volume_fractions)
    assert np.allclose(Ethanol_v0.molar_fractions, Ethanol_v3.molar_fractions)
    StrongBooze = Material(components=['Ethanol_v0', 'Water'], volume_fractions=[0.76, 0.24], density=1.2)
    assert StrongBooze.composition == [['C', 0.3724505492185721],
                                       ['H', 0.12574562181545387],
                                       ['O', 0.5018038289659741]]
    assert np.allclose(StrongBooze.mass_fractions, [0.37245055, 0.12574562, 0.50180383])
    assert np.allclose(StrongBooze.volume_fractions, [9.39448227e-05, 7.99286076e-01, 2.00619979e-01])
    assert np.allclose(StrongBooze.molar_fractions, [0.16571654, 0.66666667, 0.1676168 ])


def test_crystal_material_creation():
    pass

def test_adapt():
    pass

def test_db():
    # Create a material that is not in the database
    this_mat = Material(A=183.84, Z=74, density=19.25)
    assert this_mat.name is None
    assert this_mat.fluka_name is None
    assert this_mat.geant4_name is None
    assert this_mat not in mdb
    assert this_mat not in mdb.values()
    assert this_mat not in mdb.fluka.values()
    assert this_mat not in mdb.geant4.values()
    assert xc.materials.CarbonFibreCarbon in mdb.values()
    assert xc.materials.CarbonFibreCarbon in mdb.fluka.values()
    assert xc.materials.CarbonFibreCarbon not in mdb.geant4.values()

    this_mat_fluka  = Material(name='Yawa_FLUKA', A=1.84, Z=74, density=19.25, fluka_name='Yawa')
    this_mat_geant4 = Material(name='Yawa_GEANT4', A=111.84, Z=74, density=19.25, geant4_name='Yawa')
    assert this_mat_fluka in mdb.values()
    assert this_mat_fluka in mdb.fluka.values()
    assert this_mat_fluka not in mdb.geant4.values()
    assert this_mat_geant4 in mdb.values()
    assert this_mat_geant4 not in mdb.fluka.values()
    assert this_mat_geant4 in mdb.geant4.values()
