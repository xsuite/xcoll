# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2026.                 #
# ######################################### #

import numpy as np
import pytest

import xtrack as xt
from xtrack.particles.pdg import get_properties_from_pdg_id

from xcoll.scattering_routines.physics_settings import PhysicsSettingsHelper
from xcoll.scattering_routines.fluka.includes import _scoring_include_file
from xcoll.scattering_routines.fluka.reference_names import fluka_names_meta


FLUKA_NAME_TO_PDG = {
    value["value"]: value["pdg_id"]
    for value in fluka_names_meta.values()
}


# Every concrete FLUKA particle name appearing in the USRBDX block.
#
# PHOTON also controls OPTIPHOT and RAY.
# HEAVYION is treated separately because it does not represent one PDG ID.
SCORING_PARTICLES = [
    "PHOTON",
    "ELECTRON",
    "POSITRON",
    "MUON+",
    "MUON-",
    "TAU+",
    "TAU-",
    "NEUTRIE",
    "ANEUTRIE",
    "NEUTRIM",
    "ANEUTRIM",
    "NEUTRIT",
    "ANEUTRIT",
    "PION+",
    "PION-",
    "PIZERO",
    "KAON+",
    "KAON-",
    "KAONZERO",
    "AKAONZER",
    "KAONLONG",
    "KAONSHRT",
    "PROTON",
    "APROTON",
    "NEUTRON",
    "ANEUTRON",
    "DEUTERON",
    "TRITON",
    "3-HELIUM",
    "4-HELIUM",
    "LAMBDA",
    "ALAMBDA",
    "LAMBDAC+",
    "ALAMBDC-",
    "SIGMA-",
    "SIGMA+",
    "SIGMAZER",
    "ASIGMAZE",
    "ASIGMA-",
    "ASIGMA+",
    "XSI-",
    "AXSI+",
    "XSIZERO",
    "AXSIZERO",
    "XSIC+",
    "AXSIC-",
    "XSIC0",
    "AXSIC0",
    "XSIPC+",
    "AXSIPC-",
    "XSIPC0",
    "AXSIPC0",
    "OMEGA-",
    "AOMEGA+",
    "OMEGAC0",
    "AOMEGAC0",
    "D+",
    "D-",
    "D0",
    "D0BAR",
    "DS+",
    "DS-",
]


class DummyEngine:
    name = "dummy"
    _physics_settings_veto_list = []

    def __init__(self):
        self.particle_ref = xt.Particles(
            "proton",
            p0c=7e12,
        )

    def stop(self, *args, **kwargs):
        pass

    def _print(self, *args, **kwargs):
        pass


@pytest.fixture
def settings():
    return PhysicsSettingsHelper(DummyEngine())


def _generate_scoring(tmp_path, monkeypatch, settings, *, touches=False, crystals=False,):
    monkeypatch.chdir(tmp_path)
    path = _scoring_include_file(
        verbose=False,
        return_list=settings,
        get_touches=touches,
        use_crystals=crystals,
    )
    return path.read_text()


def _active_usrbdx_particles(text):
    particles = [
        line.split()[2]
        for line in text.splitlines()
        if line.lstrip().startswith("USRBDX")
    ]
    assert len(particles) == len(set(particles)), (
        "Duplicate USRBDX scorers found:\n"
        + "\n".join(particles)
    )
    return set(particles)

def _scorer_is_enabled(text, particle):
    return particle in _active_usrbdx_particles(text)

def _global_scorer_is_enabled(text, particle):
    return _scorer_is_enabled(text, particle)


def test_scoring_matches_physics_settings(tmp_path, monkeypatch, settings):
    """Every concrete FLUKA scorer must agree with PhysicsSettingsHelper."""
    settings.return_none = True
    settings.return_neutral = True
    settings.return_leptons = True
    settings.return_baryons = True
    settings.return_mesons = True
    settings.return_ions = True
    # Exercise exact overrides in both directions as well.
    settings.dont_return_pdg_id(211)
    settings.dont_return_pdg_id(421)
    settings.return_pdg_id(-411)
    text = _generate_scoring(tmp_path, monkeypatch, settings)
    for fluka_name in SCORING_PARTICLES:
        pdg_id = FLUKA_NAME_TO_PDG[fluka_name]
        assert _scorer_is_enabled(
            text,
            fluka_name,
        ) == settings.pdg_id_is_returned(pdg_id), (
            f"Mismatch for FLUKA particle {fluka_name} "
            f"(PDG ID {pdg_id})"
        )


def test_return_none_disables_all_particle_scoring(tmp_path, monkeypatch, settings):
    settings.return_none = True
    text = _generate_scoring(tmp_path, monkeypatch, settings)
    assert not _global_scorer_is_enabled(text, "ALL-PART")
    assert not _global_scorer_is_enabled(text, "ALL-CHAR")
    for particle in SCORING_PARTICLES:
        assert not _scorer_is_enabled(text, particle)
    assert not _scorer_is_enabled(text, "HEAVYION")


def test_return_all_uses_only_all_part(tmp_path, monkeypatch, settings):
    settings.return_all = True
    text = _generate_scoring(tmp_path, monkeypatch, settings)
    assert _global_scorer_is_enabled(text, "ALL-PART")
    assert not _global_scorer_is_enabled(text, "ALL-CHAR")
    # Individual cards are redundant when ALL-PART is active.
    for particle in SCORING_PARTICLES:
        assert not _scorer_is_enabled(text, particle)
    assert not _scorer_is_enabled(text, "HEAVYION")


def test_return_all_charged(tmp_path, monkeypatch, settings):
    settings.return_all_charged = True
    text = _generate_scoring(tmp_path, monkeypatch, settings)
    assert not _global_scorer_is_enabled(text, "ALL-PART")
    assert _global_scorer_is_enabled(text, "ALL-CHAR")
    # Charged particles are already covered by ALL-CHAR.
    # Neutral particles are not requested.
    for particle in SCORING_PARTICLES:
        assert not _scorer_is_enabled(text, particle)
    assert not _scorer_is_enabled(text, "HEAVYION")


@pytest.mark.parametrize(
    "flag,particles",
    [
        ("return_pions", [
            "PION+",
            "PION-",
            "PIZERO",
        ]),
        ("return_kaons", [
            "KAON+",
            "KAON-",
            "KAONZERO",
            "AKAONZER",
            "KAONLONG",
            "KAONSHRT",
        ]),
        ("return_other_mesons", [
            "D+",
            "D-",
            "DS+",
            "DS-",
            "D0",
            "D0BAR",
        ]),
        ("return_protons", [
            "PROTON",
            "APROTON",
        ]),
        ("return_neutrons", [
            "NEUTRON",
            "ANEUTRON",
        ]),
        ("return_other_baryons", [
            "LAMBDA",
            "ALAMBDA",
            "LAMBDAC+",
            "ALAMBDC-",
            "SIGMA-",
            "SIGMA+",
            "SIGMAZER",
            "ASIGMAZE",
            "ASIGMA-",
            "ASIGMA+",
            "XSI-",
            "AXSI+",
            "XSIZERO",
            "AXSIZERO",
            "XSIC+",
            "AXSIC-",
            "XSIC0",
            "AXSIC0",
            "XSIPC+",
            "AXSIPC-",
            "XSIPC0",
            "AXSIPC0",
            "OMEGA-",
            "AOMEGA+",
            "OMEGAC0",
            "AOMEGAC0",
        ]),
    ],
)
def test_family_scoring_matches_physics_settings(
    tmp_path,
    monkeypatch,
    settings,
    flag,
    particles,
):
    settings.return_none = True
    setattr(settings, flag, True)
    text = _generate_scoring(tmp_path, monkeypatch, settings)
    for particle in particles:
        pdg_id = FLUKA_NAME_TO_PDG[particle]
        assert _scorer_is_enabled(
            text,
            particle,
        ) == settings.pdg_id_is_returned(pdg_id), (
            f"Mismatch for {particle} / PDG {pdg_id}"
        )


@pytest.mark.parametrize(
    "pdg_id,fluka_name",
    [
        (411, "D+"),
        (-411, "D-"),
        (421, "D0"),
        (-421, "D0BAR"),
        (111, "PIZERO"),
        (3122, "LAMBDA"),
        (2112, "NEUTRON"),
        (22, "PHOTON"),
    ],
)
def test_explicit_pdg_return(tmp_path, monkeypatch, settings, pdg_id, fluka_name):
    settings.return_none = True
    settings.return_pdg_id(pdg_id)
    text = _generate_scoring(tmp_path, monkeypatch, settings)
    assert _scorer_is_enabled(text, fluka_name)
    # No global scorer should have been enabled as a side effect.
    assert not _global_scorer_is_enabled(text, "ALL-PART")
    assert not _global_scorer_is_enabled(text, "ALL-CHAR")


def test_explicit_neutral_pdg_with_all_charged(tmp_path, monkeypatch, settings):
    settings.return_all_charged = True
    settings.return_pdg_id(421)  # D0
    text = _generate_scoring(tmp_path, monkeypatch, settings)
    assert _global_scorer_is_enabled(text, "ALL-CHAR")
    # Charged particles remain covered by ALL-CHAR.
    assert not _scorer_is_enabled(text, "D+")
    # Explicit neutral exception needs its own scorer.
    assert _scorer_is_enabled(text, "D0")


def test_explicit_kill_disables_individual_card(tmp_path, monkeypatch, settings):
    settings.return_none = True
    settings.return_neutral = True
    settings.return_mesons = True
    settings.dont_return_pdg_id(421)  # D0
    text = _generate_scoring(
        tmp_path,
        monkeypatch,
        settings,
    )
    assert _scorer_is_enabled(text, "D+")
    assert not _scorer_is_enabled(text, "D0")


def test_ions_enable_heavyion(tmp_path, monkeypatch, settings):
    settings.return_none = True
    settings.return_ions = True
    text = _generate_scoring(
        tmp_path,
        monkeypatch,
        settings,
    )
    assert _scorer_is_enabled(text, "DEUTERON")
    assert _scorer_is_enabled(text, "TRITON")
    assert _scorer_is_enabled(text, "3-HELIUM")
    assert _scorer_is_enabled(text, "4-HELIUM")
    assert _scorer_is_enabled(text, "HEAVYION")


def test_explicit_heavy_ion_enables_heavyion(tmp_path, monkeypatch, settings):
    settings.return_none = True
    settings.return_pdg_id(1000060120)  # C12
    text = _generate_scoring(
        tmp_path,
        monkeypatch,
        settings,
    )
    assert _scorer_is_enabled(text, "HEAVYION")
    # No light-ion cards are needed just because one heavy ion was requested.
    assert not _scorer_is_enabled(text, "DEUTERON")
    assert not _scorer_is_enabled(text, "TRITON")
    assert not _scorer_is_enabled(text, "3-HELIUM")
    assert not _scorer_is_enabled(text, "4-HELIUM")


@pytest.mark.parametrize(
    "pdg_id,fluka_name",
    [
        (1000010020, "DEUTERON"),
        (1000010030, "TRITON"),
        (1000020030, "3-HELIUM"),
        (1000020040, "4-HELIUM"),
    ],
)
def test_explicit_light_ion_only_enables_its_card(
    tmp_path,
    monkeypatch,
    settings,
    pdg_id,
    fluka_name
):
    settings.return_none = True
    settings.return_pdg_id(pdg_id)
    text = _generate_scoring(tmp_path, monkeypatch, settings)
    assert _scorer_is_enabled(text, fluka_name)
    if fluka_name != "DEUTERON":
        assert not _scorer_is_enabled(text, "DEUTERON")
    if fluka_name != "TRITON":
        assert not _scorer_is_enabled(text, "TRITON")
    if fluka_name != "3-HELIUM":
        assert not _scorer_is_enabled(text, "3-HELIUM")
    if fluka_name != "4-HELIUM":
        assert not _scorer_is_enabled(text, "4-HELIUM")
    assert not _scorer_is_enabled(text, "HEAVYION")


def test_touches_scoring(tmp_path, monkeypatch, settings):
    settings.return_none = True
    text = _generate_scoring(
        tmp_path,
        monkeypatch,
        settings,
        touches=True,
    )
    line = next(
        line
        for line in text.splitlines()
        if "100.0" in line
        and "USERDUMP" in line
    )
    assert not line.lstrip().startswith("*USERDUMP")


def test_crystal_scoring(tmp_path, monkeypatch, settings):
    settings.return_none = True
    text = _generate_scoring(tmp_path, monkeypatch, settings, crystals=True)
    line = next(
        line
        for line in text.splitlines()
        if "CRYSTAL" in line
        and "USRICALL" in line
    )
    assert not line.lstrip().startswith("*USRICALL")


def test_photon_scoring_includes_fluka_aliases(
    tmp_path,
    monkeypatch,
    settings,
):
    settings.return_none = True
    settings.return_photons = True
    text = _generate_scoring(tmp_path, monkeypatch, settings)
    assert _active_usrbdx_particles(text) == {"PHOTON", "OPTIPHOT", "RAY"}


def test_return_all_with_explicit_kill_still_uses_all_part(
    tmp_path,
    monkeypatch,
    settings,
):
    settings.return_all = True
    settings.dont_return_pdg_id(411)
    text = _generate_scoring(tmp_path, monkeypatch, settings)
    assert _active_usrbdx_particles(text) == {"ALL-PART"}
    # Filtering happens after FLUKA returns the particle.
    assert not settings.pdg_id_is_returned(411)


def test_return_all_charged_with_explicit_kill_still_uses_all_char(
    tmp_path,
    monkeypatch,
    settings,
):
    settings.return_all_charged = True
    settings.dont_return_pdg_id(411)
    text = _generate_scoring(tmp_path, monkeypatch, settings)
    assert _active_usrbdx_particles(text) == {"ALL-CHAR"}
    assert not settings.pdg_id_is_returned(411)


@pytest.mark.parametrize(
    "pdg_id,score",
    [
        (2114, "ALL-NEUT"),   # Delta0 is known to Xtrack
        (2224, "ALL-CHAR"),   # Delta++ is known to Xtrack
        (521,  "ALL-PART"),   # B+ is unknown to Xtrack (so we cannot deduce the charge)
    ],
)
def test_unrepresented_explicit_pdg_falls_back_to_all_part(
    tmp_path,
    monkeypatch,
    settings,
    pdg_id,
    score
):
    assert pdg_id not in FLUKA_NAME_TO_PDG.values()
    settings.return_none = True
    settings.return_pdg_id(pdg_id)
    with pytest.warns(RuntimeWarning, match=str(pdg_id)):
        text = _generate_scoring(tmp_path, monkeypatch, settings)
    assert _active_usrbdx_particles(text) == {score}
    assert not settings.return_all
    assert not settings.return_all_charged
    try:
        assert settings.pdg_id_is_returned(pdg_id)
    except ValueError:
        assert settings.pdg_id_is_returned(pdg_id, 1)
