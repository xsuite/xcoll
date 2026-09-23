# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2026.                 #
# ######################################### #

import numpy as np
import pytest

import xtrack as xt
import xtrack.particles.pdg as pdg
from xcoll.scattering_routines.engine import BaseEngine
from xcoll.scattering_routines.physics_settings import PhysicsSettingsHelper


RETURN_GROUPS = {
    "leptons": ["electrons", "muons", "tauons", "neutrinos"],
    "baryons": ["protons", "neutrons", "other_baryons"],
    "mesons": ["pions", "kaons", "other_mesons"],
}

INCLUDE_PROCESSES = [
    "showers",
    "single_coulomb",
    "multiple_coulomb",
    "ionisation_fluctuations",
    "pair_production",
    "bremsstrahlung",
    "elastic",
    "inelastic",
]

CUTS = [
    "hadron_lower_momentum_cut",
    "photon_lower_momentum_cut",
    "electron_lower_momentum_cut",
    "relative_energy_cut",
]

EXPECTED_ALL_FLAGS = [
    "return_all",
    "return_all_charged",
    "return_none",
    "return_photons",
    "return_leptons",
    "return_baryons",
    "return_mesons",
    "return_ions",
    "return_neutral",
    "return_electrons",
    "return_muons",
    "return_tauons",
    "return_neutrinos",
    "return_protons",
    "return_neutrons",
    "return_other_baryons",
    "return_pions",
    "return_kaons",
    "return_other_mesons",
    "hadron_lower_momentum_cut",
    "photon_lower_momentum_cut",
    "electron_lower_momentum_cut",
    "relative_energy_cut",
    "include_showers",
    "include_single_coulomb",
    "include_multiple_coulomb",
    "include_ionisation_fluctuations",
    "include_pair_production",
    "include_bremsstrahlung",
    "include_elastic",
    "include_inelastic",
]


class DummyEngine(BaseEngine):
    name = "dummy"
    _element_classes = []

    def __init__(self, particle_ref, veto=()):
        super().__init__()
        self.particle_ref = particle_ref
        self._physics_settings_veto_list = list(veto)
        self.stop_calls = 0
        self.messages = []

    def stop(self, *args, **kwargs):
        self.stop_calls += 1

    def _print(self, *args, **kwargs):
        self.messages.append(" ".join(str(arg) for arg in args))


def make_settings(particle="proton", p0c=7e12, veto=()):
    particle_ref = xt.Particles(particle, p0c=p0c)
    engine = DummyEngine(particle_ref, veto=veto)
    settings = PhysicsSettingsHelper(engine)
    return engine, settings

def snapshot(settings):
    return {
        flag: getattr(settings, flag)
        for flag in settings.all_flags
    }


def test_all_flags():
    # Assert we have all flags and no duplicates (to ensure test is
    # updated when new flags are added)
    _, settings = make_settings()
    assert settings.all_flags == EXPECTED_ALL_FLAGS
    assert len(settings.all_flags) == len(set(settings.all_flags))


@pytest.mark.parametrize(
    "particle,p0c,is_lepton,is_proton,is_ion,A",
    [
        ("proton",   7e12,    False, True,  False, None),
        ("electron", 7e12,    True,  False, False, None),
        ("Pu-239",   94*7e12, False, False, True,  239),
    ],
)
def test_defaults(particle, p0c, is_lepton, is_proton, is_ion, A):
    _, settings = make_settings(particle, p0c)

    assert settings.return_neutral is False
    assert settings.return_photons is False
    assert settings.return_pions is False
    assert settings.return_kaons is False
    assert settings.return_other_mesons is False

    assert settings.return_electrons is is_lepton
    assert settings.return_muons is is_lepton
    assert settings.return_tauons is is_lepton
    assert settings.return_neutrinos is False

    assert settings.return_protons is (is_proton or is_ion)
    assert settings.return_neutrons is False
    assert settings.return_other_baryons is False

    assert settings.return_ions is is_ion

    # Cuts
    expected_hadron_cut = p0c / 10
    if is_ion:
        expected_hadron_cut /= A

    assert np.isclose(
        settings.hadron_lower_momentum_cut,
        expected_hadron_cut,
    )

    assert np.isclose(
        settings.photon_lower_momentum_cut,
        p0c * 1e-3,
    )

    expected_electron_cut = p0c / 10 if is_lepton else p0c * 1e-3
    assert np.isclose(
        settings.electron_lower_momentum_cut,
        expected_electron_cut,
    )

    assert settings.relative_energy_cut == 0.1

    # Physics processes
    assert settings.include_showers is is_lepton

    for process in INCLUDE_PROCESSES:
        if process != "showers":
            assert getattr(settings, f"include_{process}") is True


@pytest.mark.parametrize("group,children", RETURN_GROUPS.items())
def test_return_groups(group, children):
    _, settings = make_settings()
    settings.return_none = True
    attr = f"return_{group}"
    setattr(settings, attr, True)
    assert getattr(settings, attr) is True
    assert all(
        getattr(settings, f"return_{child}")
        for child in children
    )

    # Mixed group -> None
    setattr(settings, f"return_{children[0]}", False)
    assert getattr(settings, attr) is None

    setattr(settings, attr, False)
    assert getattr(settings, attr) is False
    assert not any(
        getattr(settings, f"return_{child}")
        for child in children
    )


def test_return_none():
    _, settings = make_settings()

    settings.return_none = True

    assert settings.return_none
    assert not settings.return_all
    assert not settings.return_all_charged

    for flag in [
        "photons", "electrons", "muons", "tauons", "neutrinos",
        "protons", "neutrons", "other_baryons",
        "pions", "kaons", "other_mesons", "ions",
    ]:
        assert getattr(settings, f"return_{flag}") is False


def test_return_all():
    _, settings = make_settings()

    settings.return_all = True

    assert settings.return_all
    assert settings.return_neutral

    for flag in [
        "photons", "electrons", "muons", "tauons", "neutrinos",
        "protons", "neutrons", "other_baryons",
        "pions", "kaons", "other_mesons", "ions",
    ]:
        assert getattr(settings, f"return_{flag}") is True


def test_return_all_charged():
    _, settings = make_settings()

    settings.return_all_charged = True

    assert settings.return_all_charged
    assert not settings.return_neutral
    assert not settings.return_photons
    assert not settings.return_neutrinos
    assert not settings.return_neutrons

    for flag in [
        "electrons", "muons", "tauons",
        "protons", "other_baryons",
        "pions", "kaons", "other_mesons", "ions",
    ]:
        assert getattr(settings, f"return_{flag}") is True


@pytest.mark.parametrize("process", INCLUDE_PROCESSES)
def test_include_process(process):
    _, settings = make_settings("proton")
    attr = f"include_{process}"
    default = False if process == "showers" else True
    assert getattr(settings, attr) is default
    setattr(settings, attr, not default)
    assert getattr(settings, attr) is not default
    setattr(settings, attr, default)
    assert getattr(settings, attr) is default
    # None means "use default"
    setattr(settings, attr, None)
    assert getattr(settings, attr) is default


def test_showers_default_depends_on_reference():
    engine, settings = make_settings("proton")
    assert settings.include_showers is False
    engine.particle_ref = xt.Particles("electron", p0c=7e12)
    assert settings.include_showers is True


def test_showers_default_depends_on_return_all():
    _, settings = make_settings("proton")
    assert settings.include_showers is False
    settings.return_all = True
    assert settings.include_showers is True


def test_showers_follows_return_all_when_default():
    _, settings = make_settings("proton")
    assert settings.include_showers is False
    settings.return_all = True
    assert settings.include_showers is True
    settings.return_none = True
    assert settings.include_showers is False


def test_showers_explicit_value_does_not_follow_return_all():
    _, settings = make_settings("proton")
    settings.include_showers = False
    settings.return_all = True
    assert settings.return_all is True
    assert settings.include_showers is False


@pytest.mark.parametrize("attr,value",
    [
        ("hadron_lower_momentum_cut",   2e9),
        ("photon_lower_momentum_cut",   2e6),
        ("electron_lower_momentum_cut", 20e6),
        ("relative_energy_cut",         1e-3),
    ],
)
def test_explicit_cut(attr, value):
    _, settings = make_settings()
    setattr(settings, attr, value)
    assert getattr(settings, attr) == value


@pytest.mark.parametrize(
    "attr,value",
    [
        ("hadron_lower_momentum_cut",   -1),
        ("photon_lower_momentum_cut",   -1),
        ("electron_lower_momentum_cut", -1),
        ("relative_energy_cut",         -1),
    ],
)
def test_invalid_cut(attr, value):
    engine, settings = make_settings()
    with pytest.raises(ValueError):
        setattr(settings, attr, value)
    assert engine.stop_calls == 1


@pytest.mark.parametrize(
    "attr,value",
    [
        ("hadron_lower_momentum_cut",   1e8),
        ("photon_lower_momentum_cut",   1e2),
        ("electron_lower_momentum_cut", 1e5),
        ("relative_energy_cut",         1e-7),
    ],
)
def test_low_cut_warning(attr, value):
    engine, settings = make_settings()
    setattr(settings, attr, value)
    assert len(engine.messages) == 1
    assert "very low" in engine.messages[0]


def test_update_recomputes_defaults():
    engine, settings = make_settings("proton", 7e12)
    engine.particle_ref = xt.Particles("electron", p0c=5e12)
    _, reference = make_settings("electron", 5e12)
    assert snapshot(settings) == snapshot(reference)


def test_update_preserves_explicit_settings():
    engine, settings = make_settings("proton", 7e12)
    settings.return_electrons = False
    settings.hadron_lower_momentum_cut = 123e9
    settings.include_elastic = False
    engine.particle_ref = xt.Particles("electron", p0c=5e12)
    assert settings.return_electrons is False
    assert settings.hadron_lower_momentum_cut == 123e9
    assert settings.include_elastic is False


def test_reset():
    _, settings = make_settings("proton", 7e12)
    _, reference = make_settings("proton", 7e12)
    settings.return_all = True
    settings.hadron_lower_momentum_cut = 2e9
    settings.photon_lower_momentum_cut = 2e6
    settings.electron_lower_momentum_cut = 20e6
    settings.relative_energy_cut = 1e-3
    settings.include_single_coulomb = False
    settings.include_multiple_coulomb = False
    settings.include_elastic = False
    settings.include_inelastic = False
    settings.reset()
    assert snapshot(settings) == snapshot(reference)


PARTICLES = [
    ("photon",            22,            0),
    ("electron",          11,           -1),
    ("positron",         -11,            1),
    ("nu_e",              12,            0),
    ("anti_nu_e",        -12,            0),
    ("muon-",             13,           -1),
    ("muon+",            -13,            1),
    ("nu_mu",             14,            0),
    ("anti_nu_mu",       -14,            0),
    ("tau-",              15,           -1),
    ("tau+",             -15,            1),
    ("nu_tau",            16,            0),
    ("anti_nu_tau",      -16,            0),
    ("pi+",              211,            1),
    ("pi-",             -211,           -1),
    ("pi0",              111,            0),
    ("K+",               321,            1),
    ("K-",              -321,           -1),
    ("K_L",              130,            0),
    ("K_S",              310,            0),
    ("K0",               311,            0),
    # Generic mesons
    ("rho+",             213,            1),
    ("rho-",            -213,           -1),
    ("rho0",             113,            0),
    ("proton",          2212,            1),
    ("antiproton",     -2212,           -1),
    ("neutron",         2112,            0),
    ("antineutron",    -2112,            0),
    # Generic baryons
    ("sigma+",          3222,            1),
    ("antisigma",      -3222,           -1),
    ("lambda",          3122,            0),
    ("antilambda",     -3122,            0),
    ("C12",       1000060120,            6),
    ("anti-C12", -1000060120,           -6),
    ("unknown",          999,            1),
]

NAMES = np.array([pp[0] for pp in PARTICLES])
PDG_IDS = np.array([pp[1] for pp in PARTICLES], dtype=np.int64)
CHARGES = np.array([pp[2] for pp in PARTICLES], dtype=float)

def selected(settings):
    mask = settings.mask_particle_return_types(PDG_IDS, CHARGES)
    return set(NAMES[mask])


@pytest.mark.parametrize(
    "flag,expected",
    [
        ("photons", {"photon"}),
        ("electrons", {"electron", "positron"}),
        ("muons", {"muon-", "muon+"}),
        ("tauons", {"tau-", "tau+"}),
        # Neutrinos additionally require the corresponding lepton family.
        ("neutrinos", set()),
        ("protons", {"proton", "antiproton"}),
        ("neutrons", {"neutron", "antineutron"}),
        ("pions", {"pi+", "pi-"}),
        ("kaons", {"K+", "K-"}),
        ("other_mesons", {"rho+", "rho-"}),
        ("other_baryons", {"sigma+", "antisigma"}),
        ("ions", {"C12", "anti-C12"}),
    ],
)
def test_mask_single_return_type(flag, expected):
    _, settings = make_settings()
    settings.return_none = True
    setattr(settings, f"return_{flag}", True)
    assert selected(settings) == expected


def test_mask_leptons():
    _, settings = make_settings()
    settings.return_none = True
    settings.return_leptons = True
    assert selected(settings) == {
        "electron", "positron", "nu_e", "anti_nu_e",
        "muon-", "muon+", "nu_mu", "anti_nu_mu",
        "tau-", "tau+", "nu_tau", "anti_nu_tau",
    }


def test_mask_mesons_with_neutral():
    _, settings = make_settings()
    settings.return_none = True
    settings.return_neutral = True
    settings.return_mesons = True
    assert selected(settings) == {
        "pi+", "pi-", "pi0",
        "K+", "K-", "K_L", "K_S", "K0",
        "rho+", "rho-", "rho0",
    }


def test_mask_baryons_with_neutral():
    _, settings = make_settings()
    settings.return_none = True
    settings.return_neutral = True
    settings.return_baryons = True
    assert selected(settings) == {
        "proton", "antiproton",
        "neutron", "antineutron",
        "sigma+", "antisigma",
        "lambda", "antilambda",
    }


def test_all_ions_are_classified_as_ions():
    _, settings = make_settings()
    A = []
    Z = []
    for aa in range(2, 240):
        for zz in range(1, min(aa, 95) + 1):
            A.append(aa)
            Z.append(zz)
    A = np.array(A)
    Z = np.array(Z)
    ions = pdg.get_pdg_id_ion(A, Z)
    # Also exercise anti-ions.
    pdg_ids = np.concatenate([ions, -ions])
    charges = np.concatenate([Z, -Z]).astype(float)
    settings.return_none = True
    settings.return_ions = True
    mask = settings.mask_particle_return_types(pdg_ids, charges)
    assert np.all(mask)


def test_ions_do_not_alias_mesons_or_baryons():
    _, settings = make_settings()
    A = []
    Z = []
    for aa in range(2, 240):
        for zz in range(1, min(aa, 95) + 1):
            A.append(aa)
            Z.append(zz)
    A = np.array(A)
    Z = np.array(Z)
    pdg_ids = pdg.get_pdg_id_ion(A, Z)
    charges = Z.astype(float)
    settings.return_none = True
    settings.return_other_mesons = True
    mask = settings.mask_particle_return_types(pdg_ids, charges)
    assert not np.any(mask)
    settings.return_none = True
    settings.return_other_baryons = True
    mask = settings.mask_particle_return_types(pdg_ids, charges)
    assert not np.any(mask)


def test_vetoed_setting():
    engine, settings = make_settings(
        veto=["include_elastic"],
    )
    with pytest.raises(
        AttributeError,
        match="does not support physics setting 'include_elastic'",
    ):
        _ = settings.include_elastic
    with pytest.raises(
        AttributeError,
        match="does not support physics setting 'include_elastic'",
    ):
        settings.include_elastic = False


@pytest.mark.parametrize(
    "attr",
    [
        "return_dinosaurs",
        "banana_cut",
        "include_gravity",
    ],
)
def test_unknown_setting(attr):
    engine, settings = make_settings()
    with pytest.raises(AttributeError):
        getattr(settings, attr)
    assert engine.stop_calls == 1


