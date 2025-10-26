# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

import pytest
import numpy as np

import xcoll as xc
from xcoll.materials import Material, CrystalMaterial, db as mdb


def _check_to_dict(mat, expected_keys):
    dct = mat.to_dict()
    assert set(dct.keys()) == set(expected_keys)
    assert dct['__class__'] == 'Material'
    for kk in dct:
        if kk == '__class__':
            continue
        elif kk.endswith('_set_manually'):
            assert dct[kk] == getattr(mat, kk)
        else:
            assert np.allclose(dct[kk], getattr(mat, kk))
    mat2 = Material.from_dict(dct)
    assert isinstance(mat2, Material)
    assert mat2 == mat


def test_material_creation():
    with pytest.raises(ValueError, match="Z must be provided for Material"):
        Material()
    with pytest.raises(ValueError, match="A must be provided for Material"):
        Material(Z=4)
    with pytest.raises(ValueError, match="density must be provided for Material"):
        Material(Z=4, A=3)

    # Create a simple material
    mat = Material(A=12.01, Z=6, density=2.265)
    assert isinstance(mat, Material)
    assert mat.name is None
    assert mat.fluka_name is None
    assert mat.geant4_name is None
    assert mat.temperature is None
    assert mat.pressure is None
    assert mat.state is None
    assert mat.A == 12.01
    assert mat.Z == 6
    assert mat.density == 2.265
    assert np.isclose(mat._ZA_mean, 0.4995836802664446)
    assert np.isclose(mat._Z_eff, 6.0)
    assert np.isclose(mat.radiation_length, 0.18849567315806856)
    assert np.isclose(mat.excitation_energy, 81.0)
    assert mat._radiation_length_set_manually is False
    assert mat._excitation_energy_set_manually is False
    assert mat.nuclear_radius is None
    assert mat.nuclear_elastic_slope is None
    assert np.all(mat.cross_section == np.array([0.0]*6))
    assert mat.hcut == 0.02
    assert np.isclose(mat.electron_density, 3.0085632439633638e+23)
    assert np.isclose(mat.plasma_energy, 30.652924826509626)
    assert np.isclose(mat.atoms_per_volume, 1.1357326245961699e+29)
    _check_to_dict(mat, {'A', 'Z', '__class__', '_excitation_energy_set_manually',
                         '_radiation_length_set_manually', 'cross_section',
                         'density', 'excitation_energy', 'hcut',
                         'radiation_length'})

    # Remove Z
    mat.Z = None
    assert mat._Z == -1
    assert mat.A == 12.01
    assert mat._ZA_mean == -1
    assert mat._Z_eff == -1
    assert mat.density == 2.265
    assert mat.radiation_length is None
    assert mat._radiation_length == -1
    assert mat.excitation_energy is None
    assert mat._excitation_energy == -1
    with pytest.raises(ValueError, match="Z must be provided for Material"):
        _check_to_dict(mat, {'A', '__class__', '_excitation_energy_set_manually',
                            '_radiation_length_set_manually', 'cross_section',
                            'density', 'hcut'})

    # Restore Z
    mat.Z = 7
    assert mat._Z == 7
    assert mat.A == 12.01
    assert np.isclose(mat._ZA_mean, 0.5828476269775188)
    assert mat._Z_eff == 7
    assert mat.density == 2.265
    assert np.isclose(mat.radiation_length, 0.1438080404832462)
    assert np.isclose(mat.excitation_energy, 82.0)
    _check_to_dict(mat, {'A', 'Z', '__class__', '_excitation_energy_set_manually',
                         '_radiation_length_set_manually', 'cross_section',
                         'density', 'excitation_energy', 'hcut',
                         'radiation_length'})

    # Change Z   TODO
    mat.Z = 7
    assert mat._Z == 7
    assert mat.A == 12.01
    assert np.isclose(mat._ZA_mean, 0.5828476269775188)
    assert mat._Z_eff == 7
    assert mat.density == 2.265
    assert np.isclose(mat.radiation_length, 0.1438080404832462)
    assert np.isclose(mat.excitation_energy, 82.0)
    _check_to_dict(mat, {'A', 'Z', '__class__', '_excitation_energy_set_manually',
                         '_radiation_length_set_manually', 'cross_section',
                         'density', 'excitation_energy', 'hcut',
                         'radiation_length'})

    # Remove A
    mat.A = None
    assert mat._A == -1
    assert mat.Z == 7
    assert mat._ZA_mean == -1
    assert mat._Z_eff == 7
    assert mat.density == 2.265
    assert mat.radiation_length is None
    assert mat._radiation_length == -1
    assert np.isclose(mat.excitation_energy, 82.0)
    with pytest.raises(ValueError, match="A must be provided for Material"):
        _check_to_dict(mat, {'Z', '__class__', '_excitation_energy_set_manually',
                            '_radiation_length_set_manually', 'cross_section',
                            'density', 'excitation_energy', 'hcut'})

    # Restore A
    mat.A = 14.3
    assert mat._A == 14.3
    assert mat.Z == 7
    assert np.isclose(mat._ZA_mean, 0.4895104895104895)
    assert mat._Z_eff == 7
    assert mat.density == 2.265
    assert np.isclose(mat.radiation_length, 0.17122855777772028)
    assert np.isclose(mat.excitation_energy, 82.0)
    dct = mat.to_dict()
    _check_to_dict(mat, {'A', 'Z', '__class__', '_excitation_energy_set_manually',
                         '_radiation_length_set_manually', 'cross_section',
                         'density', 'excitation_energy', 'hcut',
                         'radiation_length'})

    # Change A   TODO
    mat.A = 14.3
    assert mat._A == 14.3
    assert mat.Z == 7
    assert np.isclose(mat._ZA_mean, 0.4895104895104895)
    assert mat._Z_eff == 7
    assert mat.density == 2.265
    assert np.isclose(mat.radiation_length, 0.17122855777772028)
    assert np.isclose(mat.excitation_energy, 82.0)
    dct = mat.to_dict()
    _check_to_dict(mat, {'A', 'Z', '__class__', '_excitation_energy_set_manually',
                         '_radiation_length_set_manually', 'cross_section',
                         'density', 'excitation_energy', 'hcut',
                         'radiation_length'})

    # Remove density
    mat.density = None
    assert mat._density == -1
    assert mat.A == 14.3
    assert mat.Z == 7
    assert np.isclose(mat._ZA_mean, 0.4895104895104895)
    assert mat._Z_eff == 7
    assert mat.radiation_length is None
    assert mat._radiation_length == -1
    assert np.isclose(mat.excitation_energy, 82.0)
    with pytest.raises(ValueError, match="density must be provided for Material"):
        _check_to_dict(mat, {'A', 'Z', '__class__', '_excitation_energy_set_manually',
                            '_radiation_length_set_manually', 'cross_section',
                            'excitation_energy', 'hcut'})

    # Restore density
    mat.density = 4.345
    assert mat._density == 4.345
    assert mat.A == 14.3
    assert mat.Z == 7
    assert np.isclose(mat._ZA_mean, 0.4895104895104895)
    assert mat._Z_eff == 7
    assert np.isclose(mat.radiation_length, 0.08925953587262059)
    assert np.isclose(mat.excitation_energy, 82.0)
    _check_to_dict(mat, {'A', 'Z', '__class__', '_excitation_energy_set_manually',
                         '_radiation_length_set_manually', 'cross_section',
                         'density', 'excitation_energy', 'hcut',
                         'radiation_length'})

    # Change density   TODO
    mat.density = 4.345
    assert mat._density == 4.345
    assert mat.A == 14.3
    assert mat.Z == 7
    assert np.isclose(mat._ZA_mean, 0.4895104895104895)
    assert mat._Z_eff == 7
    assert np.isclose(mat.radiation_length, 0.08925953587262059)
    assert np.isclose(mat.excitation_energy, 82.0)
    _check_to_dict(mat, {'A', 'Z', '__class__', '_excitation_energy_set_manually',
                         '_radiation_length_set_manually', 'cross_section',
                         'density', 'excitation_energy', 'hcut',
                         'radiation_length'})

    # Set radiation length manually
    mat.radiation_length = 0.5
    assert mat._A == 14.3
    assert mat.Z == 7
    assert np.isclose(mat._ZA_mean, 0.4895104895104895)
    assert mat._Z_eff == 7
    assert mat.density == 4.345
    assert np.isclose(mat.radiation_length, 0.5)
    assert np.isclose(mat.excitation_energy, 82.0)
    assert mat._radiation_length_set_manually is True
    assert mat._excitation_energy_set_manually is False
    _check_to_dict(mat, {'A', 'Z', '__class__', '_excitation_energy_set_manually',
                         '_radiation_length_set_manually', 'cross_section',
                         'density', 'excitation_energy', 'hcut',
                         'radiation_length'})

    # Remove Z again - radiation length should stay
    mat.Z = None
    assert mat._Z == -1
    assert mat.A == 14.3
    assert mat._ZA_mean == -1
    assert mat._Z_eff == -1
    assert mat.density == 4.345
    assert np.isclose(mat.radiation_length, 0.5)
    assert mat.excitation_energy is None
    assert mat._excitation_energy == -1
    with pytest.raises(ValueError, match="Z must be provided for Material"):
        _check_to_dict(mat, {'A', '__class__', '_excitation_energy_set_manually',
                            '_radiation_length_set_manually', 'cross_section',
                            'density', 'hcut', 'radiation_length'})

    # Restore Z again - radiation length should stay
    mat.Z = 7
    assert mat._Z == 7
    assert mat.A == 14.3
    assert np.isclose(mat._ZA_mean, 0.4895104895104895)
    assert mat._Z_eff == 7
    assert mat.density == 4.345
    assert np.isclose(mat.radiation_length, 0.5)
    assert np.isclose(mat.excitation_energy, 82.0)
    assert mat._radiation_length_set_manually is True
    assert mat._excitation_energy_set_manually is False
    _check_to_dict(mat, {'A', 'Z', '__class__', '_excitation_energy_set_manually',
                         '_radiation_length_set_manually', 'cross_section',
                         'density', 'excitation_energy', 'hcut',
                         'radiation_length'})

    # Set excitation energy manually
    mat.excitation_energy = 897
    assert mat._Z == 7
    assert mat.A == 14.3
    assert np.isclose(mat._ZA_mean, 0.4895104895104895)
    assert mat._Z_eff == 7
    assert mat.density == 4.345
    assert np.isclose(mat.radiation_length, 0.5)
    assert np.isclose(mat.excitation_energy, 897.0)
    assert mat._radiation_length_set_manually is True
    assert mat._excitation_energy_set_manually is True
    _check_to_dict(mat, {'A', 'Z', '__class__', '_excitation_energy_set_manually',
                         '_radiation_length_set_manually', 'cross_section',
                         'density', 'excitation_energy', 'hcut',
                         'radiation_length'})

    # Remove Z again - radiation length and excitation energy should stay
    mat.Z = None
    assert mat._Z == -1
    assert mat.A == 14.3
    assert mat._ZA_mean == -1
    assert mat._Z_eff == -1
    assert mat.density == 4.345
    assert np.isclose(mat.radiation_length, 0.5)
    assert np.isclose(mat.excitation_energy, 897.0)
    assert mat._radiation_length_set_manually is True
    assert mat._excitation_energy_set_manually is True
    with pytest.raises(ValueError, match="Z must be provided for Material"):
        _check_to_dict(mat, {'A', '__class__', '_excitation_energy_set_manually',
                            '_radiation_length_set_manually', 'cross_section',
                            'density', 'excitation_energy', 'hcut',
                            'radiation_length'})

    # Restore Z again - radiation length and excitation energy should stay
    mat.Z = 7
    assert mat._Z == 7
    assert mat.A == 14.3
    assert np.isclose(mat._ZA_mean, 0.4895104895104895)
    assert mat._Z_eff == 7
    assert mat.density == 4.345
    assert np.isclose(mat.radiation_length, 0.5)
    assert np.isclose(mat.excitation_energy, 897.0)
    assert mat._radiation_length_set_manually is True
    assert mat._excitation_energy_set_manually is True
    _check_to_dict(mat, {'A', 'Z', '__class__', '_excitation_energy_set_manually',
                         '_radiation_length_set_manually', 'cross_section',
                         'density', 'excitation_energy', 'hcut',
                         'radiation_length'})

    # Remove radiation length - should revert to calculated value
    mat.radiation_length = None
    assert mat.Z == 7
    assert mat.A == 14.3
    assert np.isclose(mat._ZA_mean, 0.4895104895104895)
    assert mat._Z_eff == 7
    assert mat.density == 4.345
    assert np.isclose(mat.radiation_length, 0.08925953587262059)
    assert np.isclose(mat.excitation_energy, 897.0)
    assert mat._radiation_length_set_manually is False
    assert mat._excitation_energy_set_manually is True
    _check_to_dict(mat, {'A', 'Z', '__class__', '_excitation_energy_set_manually',
                         '_radiation_length_set_manually', 'cross_section',
                         'density', 'excitation_energy', 'hcut',
                         'radiation_length'})

    # Remove excitation energy - should revert to calculated value
    mat.excitation_energy = None
    assert mat.Z == 7
    assert mat.A == 14.3
    assert np.isclose(mat._ZA_mean, 0.4895104895104895)
    assert mat._Z_eff == 7
    assert mat.density == 4.345
    assert np.isclose(mat.radiation_length, 0.08925953587262059)
    assert np.isclose(mat.excitation_energy, 82.0)
    assert mat._radiation_length_set_manually is False
    assert mat._excitation_energy_set_manually is False
    _check_to_dict(mat, {'A', 'Z', '__class__', '_excitation_energy_set_manually',
                         '_radiation_length_set_manually', 'cross_section',
                         'density', 'excitation_energy', 'hcut',
                         'radiation_length'})

    # Test copy
    mat_cp = mat.copy()
    assert mat_cp == mat
    assert mat_cp is not mat
    # TODO check dicts!


def test_db():
    # Create a material that is not in the database
    this_mat = Material(A=183.84, Z=74, density=19.25)
    assert this_mat.name is None
    assert this_mat.fluka_name is None
    assert this_mat.geant4_name is None
    assert this_mat not in mdb.values()
    assert this_mat not in mdb.fluka.values()
    assert this_mat not in mdb.geant4.values()
    assert xc.materials.Carbon in mdb.values()
    assert xc.materials.Carbon in mdb.fluka.values()     # Only if has fluka_name
    assert xc.materials.Carbon in mdb.geant4.values()    # Only if has geant4_name

    this_mat_fluka  = Material(A=1.84, Z=74, density=19.25, fluka_name='Yawa')
    this_mat_geant4 = Material(A=111.84, Z=74, density=19.25, geant4_name='Yawa')
    assert this_mat_fluka not in mdb.values() # Because conflicting A
    assert this_mat_fluka in mdb.fluka.values()
    assert this_mat_fluka not in mdb.geant4.values()
    assert this_mat_geant4 not in mdb.values() # Because conflicting A
    assert this_mat_geant4 not in mdb.fluka.values()
    assert this_mat_geant4 in mdb.geant4.values()

    # assert isinstance(mdb['tungsten'], Material)
    # assert mdb['tungsten'].name == 'Tungsten'
    # assert mdb['tungsten'].density == 19.25
    # assert mdb['tungsten'].n_elements == 1
    # assert mdb['tungsten'].elements == ['W']
    # assert mdb['tungsten'].Z == [74]
    # assert mdb['tungsten'].A == [183.84]
    # assert mdb['tungsten'].I == [727.0e-6]
    # assert mdb['tungsten'].state == 'solid'
    # assert mdb['tungsten'].radiation_length == 0.35
    # assert mdb['tungsten'].interaction_length == 9.6
    # assert mdb['tungsten'].cross_section[0] == 0.0
    # assert mdb['tungsten'].cross_section[1] == 0.0
    # assert mdb['tungsten'].cross_section[2] == 0.0
    # assert mdb['tungsten'].cross_section[3] == 0.0
    # assert mdb['tungsten'].cross_section[4] == 0.0
    # assert mdb['tungsten'].cross_section[5] == 0.0
    # assert mdb['tungsten'].cross_section[6] == 0.0
    # assert mdb['tungsten'].cross_section[7] == 0.0
    # assert mdb['tungsten'].cross_section[8] == 0.0
    # assert mdb['tungsten'].cross_section[9] == 0.0
    # assert mdb['tungsten'].fluka_name == 'W'
    # assert mdb['tungsten'].geant4_name == 'G4_W'
    # assert not mdb['tungsten'].is_compound
    # assert not mdb['tungsten'].is_crystal

    # assert isinstance(mdb['silicon_crystal'], CrystalMaterial)
    # assert mdb['silicon_crystal'].name == 'SiliconCrystal'
    # assert mdb['silicon_crystal'].density == 2.33
    # assert mdb['silicon_crystal'].n_elements == 1
    # assert




# test_compoundmaterial.py
import math
import numpy as np
import pytest

# ---- Helpers ---------------------------------------------------------------

TOL = 1e-10

def approx(a, rel=1e-12, abs=1e-12):
    return pytest.approx(a, rel=rel, abs=abs)

def inv_weighted_mean(values, weights):
    # 1 / sum(w / v)
    return 1.0 / sum(w / v for v, w in zip(values, weights))

def log_mix(values, weights):
    # exp(sum(w * ln v)), weights assumed normalised
    s = sum(weights)
    w = [wi / s for wi in weights]
    return math.exp(sum(wi * math.log(vi) for wi, vi in zip(w, values)))

# ---- Test fixtures ---------------------------------------------------------

@pytest.fixture
def materials_db(monkeypatch):
    """
    Provide a minimal database with a few elemental materials.
    Fields used by CompoundMaterial:
      name, Z, A, density, radiation_length (MASS X0 in g/cm^2), excitation_energy (I in eV)
    """
    from xcoll.materials import Material

    # Mass radiation lengths (g/cm^2) are PDG-ish; densities are typical room-T values.
    # These are sufficient for consistency checks (not for metrology).
    elements = {
        "H":  Material(name="H",  Z=1,  A=1.00794, density=0.0708,  radiation_length=63.04e3, excitation_energy=19.2),
        "C":  Material(name="C",  Z=6,  A=12.011,  density=2.265,   radiation_length=42.70,   excitation_energy=81.0),
        "O":  Material(name="O",  Z=8,  A=15.999,  density=1.141,   radiation_length=34.24,   excitation_energy=95.0),
        "Na": Material(name="Na", Z=11, A=22.989,  density=0.971,   radiation_length=63.01,   excitation_energy=149.0),
        "Cl": Material(name="Cl", Z=17, A=35.453,  density=3.214,   radiation_length=23.97,   excitation_energy=180.0),
        "Ni": Material(name="Ni", Z=28, A=58.6934, density=8.90,    radiation_length=13.07,   excitation_energy=135.0),
        "Cu": Material(name="Cu", Z=29, A=63.546,  density=8.96,    radiation_length=12.86,   excitation_energy=322.0),
        "W":  Material(name="W",  Z=74, A=183.84,  density=19.25,   radiation_length=6.76,    excitation_energy=727.0),
        "Zn": Material(name="Zn", Z=30, A=65.38,   density=7.14,    radiation_length=12.86,   excitation_energy=330.0),
    }

    # Monkeypatch the central DB dict used by CompoundMaterial
    # (the class does 'from xcoll.materials.database import db as mdb')
    import types
    fake_db_module = types.SimpleNamespace(db=elements)

    monkeypatch.setitem(sys.modules, 'xcoll.materials.database', fake_db_module)

    return elements

# ---- Tests -----------------------------------------------------------------

def test_by_counts_water(materials_db):
    from xcoll.materials import CompoundMaterial

    H2O = CompoundMaterial(components=["H", "O"], n_atoms=np.array([2, 1]), name="H2O", density=0.997)
    # Mass fractions from stoichiometry
    M = 2 * materials_db["H"].A + materials_db["O"].A
    wH = 2 * materials_db["H"].A / M
    wO = materials_db["O"].A / M

    assert H2O.mass_fractions[0] == approx(wH)
    assert H2O.mass_fractions[1] == approx(wO)
    # Effective molar mass must equal the true molar mass for a true compound
    assert H2O.molar_mass == approx(M)
    assert H2O.effective_molar_mass is None  # you only define it for mass-defined compounds
    # Atomic fractions are the normalised counts
    af = H2O.atomic_fractions
    assert af[0] == approx(2/3)
    assert af[1] == approx(1/3)
    # X0 mass mix: 1/X0 = sum w_i/X0_i
    X0_expected = inv_weighted_mean(
        [materials_db["H"].radiation_length, materials_db["O"].radiation_length],
        [wH, wO]
    )
    assert H2O.radiation_length == approx(X0_expected)

def test_by_mass_w_nicu(materials_db):
    from xcoll.materials import CompoundMaterial

    comp = CompoundMaterial(
        components=["W", "Ni", "Cu"],
        mass_fractions=np.array([0.95, 0.035, 0.015]),
        name="Inermet-like"
    )
    w = comp.mass_fractions
    # Effective molar mass per atom 1/sum(w/A)
    Aeff = 1.0 / sum(wi / ai for wi, ai in zip(w, [materials_db[s].A for s in ["W","Ni","Cu"]]))
    assert comp.effective_molar_mass == approx(Aeff)

    # Atomic fractions from mass fractions: x_i ∝ w_i/A_i
    x = np.array([w[0]/materials_db["W"].A,
                  w[1]/materials_db["Ni"].A,
                  w[2]/materials_db["Cu"].A])
    x /= x.sum()
    af = comp.atomic_fractions
    assert af[0] == approx(x[0])
    assert af[1] == approx(x[1])
    assert af[2] == approx(x[2])

    # X0 mass mixing check
    X0_expected = inv_weighted_mean(
        [materials_db["W"].radiation_length,
         materials_db["Ni"].radiation_length,
         materials_db["Cu"].radiation_length],
        w
    )
    assert comp.radiation_length == approx(X0_expected)

def test_nested_counts_and_mass_paths_equivalence(materials_db):
    from xcoll.materials import CompoundMaterial

    # Define true compounds
    NaCl = CompoundMaterial(components=["Na","Cl"], n_atoms=np.array([1,1]), name="NaCl")
    CO2  = CompoundMaterial(components=["C","O"],  n_atoms=np.array([1,2]), name="CO2")
    H2O  = CompoundMaterial(components=["H","O"],  n_atoms=np.array([2,1]), name="H2O")

    # Mixture by counts of formula units: 3 NaCl + 2 CO2 + 4 H2O
    mix_counts = CompoundMaterial(
        components=[NaCl, CO2, H2O],
        n_atoms=np.array([3,2,4]),
        name="mix-by-counts"
    )

    # Recompute expected elemental mass fractions manually:
    def elem_w_of(comp):
        return dict(zip([e.name for e in comp.components], comp.mass_fractions))

    Ms = [NaCl.molar_mass, CO2.molar_mass, H2O.molar_mass]
    W = np.array([3*Ms[0], 2*Ms[1], 4*Ms[2]])
    W /= W.sum()

    # element mass fractions from the three compounds
    e_all = {"Na":0.0,"Cl":0.0,"C":0.0,"O":0.0,"H":0.0}
    for comp, Wi in zip([NaCl, CO2, H2O], W):
        efrac = elem_w_of(comp)
        for sym, wi in efrac.items():
            e_all[sym] += Wi * wi

    # Compare to class result
    got = dict(zip([e.name for e in mix_counts.components], mix_counts.mass_fractions))
    # Merge duplicate element names
    from collections import defaultdict
    agg = defaultdict(float)
    for k, v in got.items():
        agg[k] += v
    for sym in e_all:
        assert agg[sym] == approx(e_all[sym], rel=1e-10, abs=1e-12)

def test_log_mix_excitation_energy(materials_db, monkeypatch):
    from xcoll.materials import CompoundMaterial

    # Two-element mix with known I-values; check PDG logarithmic rule with weights ∝ w*Z/A
    comp = CompoundMaterial(components=["C","O"], mass_fractions=np.array([0.7, 0.3]), name="C/O mix")
    w = comp.mass_fractions
    Z = np.array([materials_db["C"].Z, materials_db["O"].Z])
    A = np.array([materials_db["C"].A, materials_db["O"].A])
    I = np.array([materials_db["C"].excitation_energy, materials_db["O"].excitation_energy])

    weights = w * Z / A
    weights /= weights.sum()
    I_expected = log_mix(I.tolist(), weights.tolist())
    assert comp.excitation_energy == approx(I_expected)

def test_density_fallback_inverse_weighted(materials_db):
    from xcoll.materials import CompoundMaterial

    comp = CompoundMaterial(components=["C","O"], mass_fractions=np.array([0.6, 0.4]))
    # Fallback density uses inverse weighted mean of component densities
    rho_expected = inv_weighted_mean([materials_db["C"].density, materials_db["O"].density], [0.6, 0.4])
    assert comp.density == approx(rho_expected)

def test_error_both_counts_and_mass_raises(materials_db):
    from xcoll.materials import CompoundMaterial
    with pytest.raises(ValueError):
        CompoundMaterial(components=["C","O"], n_atoms=np.array([1,1]), mass_fractions=np.array([0.5,0.5]))

def test_error_nonpositive_shares(materials_db):
    from xcoll.materials import CompoundMaterial
    with pytest.raises(ValueError):
        CompoundMaterial(components=["C","O"], mass_fractions=np.array([0.6, 0.0]))

def test_error_nested_counts_without_child_counts(materials_db):
    from xcoll.materials import CompoundMaterial
    # Child defined by mass (no n_atoms)
    child = CompoundMaterial(components=["C","O"], mass_fractions=np.array([0.6, 0.4]), name="child-mass")
    with pytest.raises(ValueError):
        CompoundMaterial(components=[child, "H"], n_atoms=np.array([2, 1]))

def test_atomic_fractions_sum_to_one(materials_db):
    from xcoll.materials import CompoundMaterial
    comp = CompoundMaterial(components=["W","Cu","Ni"], mass_fractions=np.array([0.95,0.015,0.035]))
    s = comp.atomic_fractions.sum()
    assert s == approx(1.0, abs=1e-12)

def test_effective_A_matches_definition(materials_db):
    from xcoll.materials import CompoundMaterial
    comp = CompoundMaterial(components=["W","Cu","Ni"], mass_fractions=np.array([0.95,0.015,0.035]))
    Aeff = 1.0 / sum(comp.mass_fractions[i] / materials_db[s].A
                     for i, s in enumerate(["W","Cu","Ni"]))
    assert comp.effective_molar_mass == approx(Aeff)
