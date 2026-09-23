# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2026.                 #
# ######################################### #

import pytest
import numpy as np
from pathlib import Path
from types import SimpleNamespace

import xtrack as xt
import xcoll as xc
from xcoll.scattering_routines.fluka.environment import format_fluka_float
from xcoll.scattering_routines.fluka.includes import (
    _physics_include_file,
    _scoring_include_file,
)
from xcoll.scattering_routines.fluka.fluka_input import get_collimators_from_input_file


@pytest.mark.serial
@pytest.mark.fluka
@pytest.mark.parametrize("el_type", ['collimator', 'crystal'])
def test_fluka_input_single(el_type, register_cleanup):
    print(f"\nTesting FLUKA input generation for single {el_type}... in {Path.cwd()}")
    if xc.fluka.engine.is_running():
        xc.fluka.engine.stop(clean=True)
    if el_type == 'collimator':
        coll = xc.FlukaCollimator(length=0.456, angle=32, jaw=[0.01, -0.02], tilt=[10e-6, -8.7e-6], material='Yttrium')
    else:
        coll = xc.FlukaCrystal(length=0.002, side='-', angle=90, jaw=-0.01, tilt=43e-6, material='Yttrium', bending_radius=65)
    with pytest.raises(ValueError, match="Need to provide either a line with a reference particle, or `particle_ref`."):
        input_file = xc.fluka.engine.generate_input_file(elements=coll, names='TestColl', clean=False)
    path_tmp = Path.cwd() / f'fluka_run_temp_test_single_{el_type}'
    register_cleanup(path_tmp)
    particle_ref = xt.Particles('proton', p0c=7e12)
    input_file = xc.fluka.engine.generate_input_file(elements=coll, names='TestColl', clean=False,
                        particle_ref=particle_ref, cwd=path_tmp, filename=path_tmp / 'fluka_input_test.inp')
    assert hasattr(input_file, '__iter__') and not isinstance(input_file, str)
    assert input_file == [path_tmp / 'fluka_input_test.inp', path_tmp / 'insertion.txt']
    assert input_file[0].exists()
    assert input_file[1].exists()
    assert path_tmp.exists()
    for file in ['fluka_input_orig.inp', 'include_custom_biasing.inp', 'include_define.inp',
                 'include_settings_physics.inp', 'fluka_input.log', 'include_custom_assignmat.inp',
                 'include_custom_scoring.inp', 'include_settings_beam.inp', 'linebuilder.log',
                 'prototypes.lbp']:
        assert (path_tmp / file).exists()

    new_coll = get_collimators_from_input_file(input_file[0])
    name = 'TestColl'.lower()
    assert name in new_coll
    if el_type == 'collimator':
        assert np.isclose(new_coll[name]['length'],  0.558)
        assert np.isclose(new_coll[name]['angle'],  32)
        assert np.isclose(new_coll[name]['jaw'][0], 0.01)
        assert np.isclose(new_coll[name]['jaw'][1], -0.02)
        assert np.isclose(new_coll[name]['tilt'][0], 10e-6)
        assert np.isclose(new_coll[name]['tilt'][1], -8.7e-6)
    else:
        assert np.isclose(new_coll[name]['length'],  0.1044)
        assert np.isclose(new_coll[name]['angle'],  90)
        assert new_coll[name]['jaw'][0] is None
        assert np.isclose(new_coll[name]['jaw'][1], -0.01)
        assert new_coll[name]['tilt'][0] is None
        assert np.isclose(new_coll[name]['tilt'][1], 43e-6)

    # Check assembly
    prototypes = coll.assembly.prototypes
    assert len(prototypes) == 2
    jaw = [pp for pp in coll.assembly.prototypes if pp.fedb_tag.endswith('B')]
    assert len(jaw) == 1
    jaw = jaw[0]
    xc.FlukaPrototype.inspect_prototypes_file(coll, path_tmp / 'prototypes.lbp')
    found = False
    with input_file[0].open('r') as fp:
        for line in fp:
            if f"RPP {jaw.fedb_tag}   0.0 20.0 -10.0 10.0 -22.8 22.8" in line \
            and el_type == 'collimator':
                found = True
            if f"RPP {jaw.fedb_tag}   0.0 2.2 -2.75 2.75 -0.24 0.24" in line \
            and el_type == 'crystal':
                found = True
    assert found
    found = False
    with input_file[1].open('r') as fp:
        for line in fp:
            if "     1          INROT_1              INROT_1             0.051000" in line \
            and el_type == 'collimator':
                found = True
            if "     1          INROT_1              INROT_1             0.051200" in line \
            and el_type == 'crystal':
                found = True
        print(line)
    assert found

    # Check material assignment
    assert coll.material == xc.materials.Yttrium
    assert coll.assembly.material == xc.materials.Yttrium
    assert coll.material.fluka_name is not None
    assert jaw.material == xc.materials.Yttrium
    found_1 = False
    found_2 = False
    found_3 = False
    with input_file[0].open('r') as fp:
        for line in fp:
            if f"ASSIGNMA    {coll.material.fluka_name}  {jaw.fedb_tag}" in line:
                found_1 = True
            if f"MATERIAL        39.0               4.469                              {coll.material.fluka_name}" in line:
                found_2 = True
            if f"MAT-PROP                           379.0  {coll.material.fluka_name}" in line:
                found_3 = True
        print(line)
    assert found_1
    assert found_2
    assert found_3

    # Check crystal definition
    if el_type == 'crystal':
        found_1 = False
        found_2 = False
        found_3 = False
        with input_file[0].open('r') as fp:
            for line in fp:
                if f"CRYSTAL     {coll.assembly.prototypes[1].fedb_tag}  0.030769       0.2       0.0       0.0     300.0 110" in line:
                    found_1 = True
                if "CRYSTAL          0.0      -1.0       0.0       0.0       0.0       1.0 &" in line:
                    found_2 = True
                if "USRICALL        50.0                                                  CRYSTAL" in line:
                    found_3 = True
        assert found_1
        assert found_2
        assert found_3


@pytest.mark.serial
@pytest.mark.fluka
@pytest.mark.parametrize("ignore_crystals", [True, False], ids=['no_crystals', 'with_crystals'])
def test_fluka_input_line(ignore_crystals, register_cleanup):
    print(f"\nTesting FLUKA input generation for line (ignore_crystals={ignore_crystals})... in {Path.cwd()}")
    if xc.fluka.engine.is_running():
        xc.fluka.engine.stop(clean=True)
    beam = 1
    path = Path(__file__).parent
    env = xt.load(path / 'data' / f'sequence_lhc_run3_b{beam}.json')
    line = env.lines[f'lhcb{beam}']
    colldb = xc.CollimatorDatabase.from_yaml(path / 'data' / f'colldb_lhc_run3_ir7.yaml', beam=beam,
                                             ignore_crystals=ignore_crystals)
    colldb.install_fluka_collimators(line=line, verbose=True)
    tt_colls = line.get_table().rows.match(
        element_type='|'.join(cc.__name__ for cc in xc.collimator_classes)
    )
    colls = [line[name] for name in tt_colls.name]
    line.build_tracker()
    line.xcoll.collimators.assign_optics()
    if not ignore_crystals:
        line.xcoll.collimators.align_to_beam_divergence()
    path_tmp = Path.cwd() / f'fluka_run_temp_test_line_{ignore_crystals}'
    register_cleanup(path_tmp)
    particle_ref = xt.Particles('proton', p0c=7e12)
    input_file = xc.fluka.engine.generate_input_file(line=line, clean=False, cwd=path_tmp,
                        particle_ref=particle_ref, filename=path_tmp / 'fluka_input_test.inp')
    assert input_file == [path_tmp / 'fluka_input_test.inp', path_tmp / 'insertion.txt']
    assert input_file[0].exists()
    assert input_file[1].exists()
    assert path_tmp.exists()
    for file in ['fluka_input_orig.inp', 'include_custom_biasing.inp', 'include_define.inp',
                 'include_settings_physics.inp', 'fluka_input.log', 'include_custom_assignmat.inp',
                 'include_custom_scoring.inp', 'include_settings_beam.inp', 'linebuilder.log',
                 'prototypes.lbp']:
        assert (path_tmp / file).exists()

    new_coll_dct = get_collimators_from_input_file(input_file[0])
    if ignore_crystals:
        this_coll_dct = {
            'tcp.d6l7.b1':   {'length': 1.482, 'angle': 90.0,  'tilt': [0.0, 0.0], 'jaw': [0.0009238348691078535, -0.0009219591452698239]},
            'tcp.c6l7.b1':   {'length': 1.482, 'angle': 0.0,   'tilt': [2.5e-06, -2.5e-06], 'jaw': [0.0013138622735122674, -0.0013103666451246276]},
            'tcp.b6l7.b1':   {'length': 1.482, 'angle': 127.5, 'tilt': [0.0, 0.0], 'jaw': [0.0010997388044970968, -0.001100319238767078]},
            'tcsg.a6l7.b1':  {'length': 1.482, 'angle': 141.1, 'tilt': [0.0, 0.0], 'jaw': [0.0014749366037034584, -0.0014740571819316095]},
            'tcsg.b5l7.b1':  {'length': 1.482, 'angle': 143.5, 'tilt': [0.0, 0.0], 'jaw': [0.001814181900021694, -0.0018086857446322213]},
            'tcsg.a5l7.b1':  {'length': 1.102, 'angle': 40.7,  'tilt': [0.0, 0.0], 'jaw': [0.0018463405948736522, -0.0018520014167018317]},
            'tcsg.d4l7.b1':  {'length': 1.482, 'angle': 90.0,  'tilt': [0.0, 0.0], 'jaw': [0.0011923168368670467, -0.0011930737268484037]},
            'tcspm.b4l7.b1': {'length': 1.482, 'angle': 0.0,   'tilt': [0.0, 0.0], 'jaw': [0.0016543438662348642, -0.0016608919799403488]},
            'tcsg.a4l7.b1':  {'length': 1.482, 'angle': 134.6, 'tilt': [0.0, 0.0], 'jaw': [0.0016554613686823316, -0.0016540730713363594]},
            'tcsg.a4r7.b1':  {'length': 1.482, 'angle': 46.3,  'tilt': [0.0, 0.0], 'jaw': [0.0016562125801096172, -0.0016637353535826627]}
        }
    else:
        this_coll_dct = {
            'tcp.d6l7.b1':   {'length': 1.482,  'angle': 90.0,  'tilt': [0.0, 0.0], 'jaw': [0.0009238348692113263, -0.0009219591453732967]},
            'tcp.c6l7.b1':   {'length': 1.482,  'angle': 0.0,   'tilt': [2.5e-06, -2.5e-06], 'jaw': [0.0013138622732093985, -0.0013103666448213147]},
            'tcp.b6l7.b1':   {'length': 1.482,  'angle': 127.5, 'tilt': [0.0, 0.0], 'jaw': [0.0010997388044393652, -0.0011003192387097904]},
            'tcsg.a6l7.b1':  {'length': 1.482,  'angle': 141.1, 'tilt': [0.0, 0.0], 'jaw': [0.001474936603826471, -0.0014740571820550663]},
            'tcpcv.a6l7.b1': {'length': 0.1068, 'angle': 90.0,  'tilt': [1.66577e-05, None], 'jaw': [0.0017621469072926072, None]},
            'tcsg.b5l7.b1':  {'length': 1.482,  'angle': 143.5, 'tilt': [0.0, 0.0], 'jaw': [0.001814181900504419, -0.0018086857451153904]},
            'tcsg.a5l7.b1':  {'length': 1.102,  'angle': 40.7,  'tilt': [0.0, 0.0], 'jaw': [0.0018463405953683676, -0.001852001417195659]},
            'tcsg.d4l7.b1':  {'length': 1.482,  'angle': 90.0,  'tilt': [0.0, 0.0], 'jaw': [0.0011923168370637782, -0.0011930737270451353]},
            'tcpch.a4l7.b1': {'length': 0.1068, 'angle': 0.0,   'tilt': [1.21108e-05, None], 'jaw': [0.0020205086648474287, None]},
            'tcspm.b4l7.b1': {'length': 1.482,  'angle': 0.0,   'tilt': [0.0, 0.0], 'jaw': [0.0016543438666660748, -0.0016608919803706712]},
            'tcsg.a4l7.b1':  {'length': 1.482,  'angle': 134.6, 'tilt': [0.0, 0.0], 'jaw': [0.001655461368718747, -0.0016540730713732188]},
            'tcsg.a4r7.b1':  {'length': 1.482,  'angle': 46.3,  'tilt': [0.0, 0.0], 'jaw': [0.001656212580093186, -0.0016637353535662314]}
        }
    for name, params in this_coll_dct.items():
        assert name in new_coll_dct
        assert np.isclose(new_coll_dct[name]['length'],  params['length'])
        assert np.isclose(new_coll_dct[name]['angle'],   params['angle'])
        assert np.isclose(new_coll_dct[name]['jaw'][0],  params['jaw'][0])
        if params['jaw'][1] is None:
            assert new_coll_dct[name]['jaw'][1] is None
        else:
            assert np.isclose(new_coll_dct[name]['jaw'][1],  params['jaw'][1])
        assert np.isclose(new_coll_dct[name]['tilt'][0], params['tilt'][0])
        if params['tilt'][1] is None:
            assert new_coll_dct[name]['tilt'][1] is None
        else:
            assert np.isclose(new_coll_dct[name]['tilt'][1], params['tilt'][1])

    # Check assembly
    xc.FlukaPrototype.inspect_prototypes_file(colls, path_tmp / 'prototypes.lbp')
    for coll in colls:
        _ = coll.assembly.prototypes # Check that prototypes were created
    with input_file[1].open('r') as fp:
        insertion_txt = fp.read()
    if ignore_crystals:
        assert """\
     1          INROT_1              INROT_1             -0.041500  
     2          INROT_2              INROT_2             0.441000   
     3          INROT_3              INROT_3             0.441000   
     4          INROT_4              INROT_4             0.241000   
     5          INROT_5              INROT_5             0.241000   
     6          INROT_6              INROT_6             0.051000   
     7          INROT_7              INROT_7             0.241000   
     8          INROT_8              INROT_8             0.241000   
     9          INROT_9              INROT_9             0.241000   
     10         INROT_10             INROT_10            0.241000""" in insertion_txt
    else:
        assert """\
     1          INROT_1              INROT_1             -0.041500  
     2          INROT_2              INROT_2             0.441000   
     3          INROT_3              INROT_3             0.441000   
     4          INROT_4              INROT_4             0.241000   
     5          INROT_5              INROT_5             0.241000   
     6          INROT_6              INROT_6             0.051000   
     7          INROT_7              INROT_7             0.241000   
     8          INROT_8              INROT_8             0.241000   
     9          INROT_9              INROT_9             0.241000   
     10         INROT_10             INROT_10            0.241000   
     11         INROT_11             INROT_11            0.051400   
     12         INROT_12             INROT_12            0.051400""" in insertion_txt



def _active_cards(text, card):
    """Return uncommented FLUKA cards of the requested type."""
    result = []
    for line in text.splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith("*"):
            continue
        if stripped.split()[0] == card:
            result.append(line)
    return result

def _make_physics_input(tmp_path, monkeypatch, **kwargs):
    monkeypatch.chdir(tmp_path)
    params = {
        "verbose": False,
        "particle_ref": xt.Particles("proton", p0c=7e12),
        "hadron_lower_momentum_cut": 2e9,
        "photon_lower_momentum_cut": 3e6,
        "electron_lower_momentum_cut": 4e6,
        "include_showers": True,
        "include_single_coulomb": True,
        "include_multiple_coulomb": True,
        "include_elastic": True,
        "include_inelastic": True,
        "include_pair_production": True,
        "include_bremsstrahlung": True,
        "include_ionisation_fluctuations": True,
    }
    params.update(kwargs)
    filename = _physics_include_file(**params)
    return filename.read_text()


@pytest.mark.parametrize("include_showers", [True, False])
def test_fluka_input_showers(tmp_path, monkeypatch, include_showers):
    text = _make_physics_input(
        tmp_path, monkeypatch,
        include_showers=include_showers,
    )
    emfcut = _active_cards(text, "EMFCUT")
    emf = _active_cards(text, "EMF")
    deltaray = _active_cards(text, "DELTARAY")
    if include_showers:
        assert len(emfcut) == 2
        assert len(emf) == 0
        assert len(deltaray) == 0
    else:
        assert len(emfcut) == 0
        assert len(emf) == 1
        assert len(deltaray) == 1


@pytest.mark.parametrize(
    "single,multiple,expected",
    [
        (
            True, True,
            [
                "MULSOPT                                        "
                "1.0       1.0       1.0GLOBAL",
            ],
        ),
        (
            True, False,
            [
                "MULSOPT          0.0       0.0       0.0       "
                "1.0       1.0 99999999.GLOBAL",
            ],
        ),
        (
            False, True,
            [
                "MULSOPT                                       "
                "-1.0      -1.0      -1.0GLOBAL",
            ],
        ),
        (
            False, False,
            [
                "MULSOPT                                       "
                "-1.0      -1.0      -1.0GLOBAL",
                "MULSOPT                    3.0       3.0  "
                "BLCKHOLE  @LASTMAT",
            ],
        ),
    ],
)
def test_fluka_input_coulomb(
    tmp_path, monkeypatch, single, multiple, expected
):
    text = _make_physics_input(
        tmp_path, monkeypatch,
        include_single_coulomb=single,
        include_multiple_coulomb=multiple,
    )
    assert _active_cards(text, "MULSOPT") == expected


@pytest.mark.parametrize(
    "pair,brem,what1",
    [
        (True,  True,   3.0),
        (True,  False,  1.0),
        (False, True,   2.0),
        (False, False, -3.0),
    ],
)
def test_fluka_input_pair_bremsstrahlung(
    tmp_path, monkeypatch, pair, brem, what1
):
    text = _make_physics_input(
        tmp_path, monkeypatch,
        include_pair_production=pair,
        include_bremsstrahlung=brem,
    )
    cards = _active_cards(text, "PAIRBREM")
    assert len(cards) == 1
    assert float(cards[0].split()[1]) == what1


@pytest.mark.parametrize("fluctuations", [True, False])
def test_fluka_input_ionisation_fluctuations(
    tmp_path, monkeypatch, fluctuations
):
    text = _make_physics_input(
        tmp_path, monkeypatch,
        include_ionisation_fluctuations=fluctuations,
    )
    cards = _active_cards(text, "IONFLUCT")
    if fluctuations:
        assert cards == []
    else:
        assert len(cards) == 1
        assert "-1.0" in cards[0]


@pytest.mark.parametrize(
    "elastic,inelastic,comment,n_threshold_values",
    [
        (True,  True,  None,                                         0),
        (False, True,  "Deactivate elastic hadronic interactions",  1),
        (True,  False, "Deactivate inelastic hadronic interactions",1),
        (False, False, "Deactivate hadronic interactions",          2),
    ],
)
def test_fluka_input_hadronic_interactions(
    tmp_path, monkeypatch,
    elastic, inelastic, comment, n_threshold_values,
):
    particle_ref = xt.Particles("proton", p0c=7e12)
    text = _make_physics_input(
        tmp_path, monkeypatch,
        particle_ref=particle_ref,
        include_elastic=elastic,
        include_inelastic=inelastic,
    )
    cards = _active_cards(text, "THRESHOLd")
    if elastic and inelastic:
        assert cards == []
        return
    assert len(cards) == 1
    assert comment in text
    threshold = format_fluka_float(
        5 * particle_ref.energy0[0] / 1e9
    ).strip()
    assert cards[0].count(threshold) == n_threshold_values


def test_fluka_input_momentum_cuts(tmp_path, monkeypatch):
    hadron_cut = 2e9
    photon_cut = 3e6
    electron_cut = 4e6
    text = _make_physics_input(
        tmp_path, monkeypatch,
        hadron_lower_momentum_cut=hadron_cut,
        photon_lower_momentum_cut=photon_cut,
        electron_lower_momentum_cut=electron_cut,
        include_showers=True,
    )
    # Hadron thresholds: generic, D, T, He3, He4
    part_thr = _active_cards(text, "PART-THR")
    assert len(part_thr) == 5
    values = [float(line.split()[1]) for line in part_thr]
    assert values == [2.0, 4.0, 6.0, 6.0, 8.0]
    # EM production cuts.
    emfcut = _active_cards(text, "EMFCUT")
    assert len(emfcut) == 2
    photon_gev = format_fluka_float(photon_cut / 1e9).strip()
    electron_energy = np.sqrt(electron_cut**2 + (511e3)**2)
    electron_gev = format_fluka_float(electron_energy / 1e9).strip()
    for line in emfcut:
        assert electron_gev in line
        assert photon_gev in line


@pytest.mark.parametrize(
    "particle,p0c,active",
    [
        ("proton", 7e12, False),
        ("Pu-239", 94*7e12, True),
    ],
)
def test_fluka_input_em_dissociation(
    tmp_path, monkeypatch, particle, p0c, active
):
    text = _make_physics_input(
        tmp_path, monkeypatch,
        particle_ref=xt.Particles(particle, p0c=p0c),
    )
    cards = [
        line for line in _active_cards(text, "PHYSICS")
        if "EM-DISSO" in line
    ]
    assert bool(cards) is active


_RETURN_DEFAULTS = {
    "return_all": False,
    "return_all_charged": False,
    "return_neutral": False,
    "return_photons": False,
    "return_electrons": False,
    "return_muons": False,
    "return_tauons": False,
    "return_neutrinos": False,
    "return_protons": False,
    "return_neutrons": False,
    "return_other_baryons": False,
    "return_pions": False,
    "return_kaons": False,
    "return_other_mesons": False,
    "return_ions": False,
}

def _return_settings(**kwargs):
    values = _RETURN_DEFAULTS.copy()
    values.update(kwargs)
    return SimpleNamespace(**values)

def _make_scoring_input(tmp_path, monkeypatch, **kwargs):
    monkeypatch.chdir(tmp_path)
    filename = _scoring_include_file(
        verbose=False,
        return_list=_return_settings(**kwargs),
    )
    return filename.read_text()

def _active_usrbdx_particles(text):
    return {
        line.split()[2]
        for line in _active_cards(text, "USRBDX")
    }


@pytest.mark.parametrize(
    "flag,particles",
    [
        (
            "return_photons",
            {"PHOTON", "OPTIPHOT", "RAY"},
        ),
        (
            "return_electrons",
            {"ELECTRON", "POSITRON"},
        ),
        (
            "return_muons",
            {"MUON+", "MUON-"},
        ),
        (
            "return_tauons",
            {"TAU+", "TAU-"},
        ),
        (
            "return_neutrinos",
            {
                "NEUTRIE", "ANEUTRIE",
                "NEUTRIM", "ANEUTRIM",
                "NEUTRIT", "ANEUTRIT",
            },
        ),
        (
            "return_protons",
            {"PROTON", "APROTON"},
        ),
        (
            "return_neutrons",
            {"NEUTRON", "ANEUTRON"},
        ),
        (
            "return_pions",
            {"PION+", "PION-"},
        ),
        (
            "return_kaons",
            {"KAON+", "KAON-"},
        ),
        (
            "return_ions",
            {"DEUTERON", "TRITON", "3-HELIUM", "4-HELIUM", "HEAVYION"},
        ),
        (
            "return_other_mesons",
            {"D+", "D-", "DS+", "DS-"},
        ),
        (
            "return_other_baryons",
            {
                "LAMBDAC+", "ALAMBDC-",
                "SIGMA-", "SIGMA+",
                "ASIGMA-", "ASIGMA+",
                "XSI-", "AXSI+",
                "XSIC+", "AXSIC-",
                "XSIPC+", "AXSIPC-",
                "OMEGA-", "AOMEGA+",
            },
        ),
    ],
)
def test_fluka_input_return_types(
    tmp_path, monkeypatch, flag, particles
):
    text = _make_scoring_input(
        tmp_path, monkeypatch,
        **{flag: True},
    )
    assert _active_usrbdx_particles(text) == particles


@pytest.mark.parametrize(
    "flag,particles",
    [
        (
            "return_pions",
            {"PION+", "PION-", "PIZERO"},
        ),
        (
            "return_kaons",
            {
                "KAON+", "KAON-",
                "KAONZERO", "AKAONZER", "KAONLONG", "KAONSHRT",
            },
        ),
        (
            "return_other_mesons",
            {"D+", "D-", "DS+", "DS-", "D0", "D0BAR"},
        ),
        (
            "return_other_baryons",
            {
                "LAMBDAC+", "ALAMBDC-",
                "SIGMA-", "SIGMA+",
                "ASIGMA-", "ASIGMA+",
                "XSI-", "AXSI+",
                "XSIC+", "AXSIC-",
                "XSIPC+", "AXSIPC-",
                "OMEGA-", "AOMEGA+",

                "LAMBDA", "ALAMBDA",
                "SIGMAZER", "ASIGMAZE",
                "XSIZERO", "AXSIZERO",
                "XSIC0", "AXSIC0",
                "XSIPC0", "AXSIPC0",
                "OMEGAC0", "AOMEGAC0",
            },
        ),
    ],
)
def test_fluka_input_return_types_neutral(
    tmp_path, monkeypatch, flag, particles
):
    text = _make_scoring_input(
        tmp_path, monkeypatch,
        return_neutral=True,
        **{flag: True},
    )
    assert _active_usrbdx_particles(text) == particles


def test_fluka_input_return_all(tmp_path, monkeypatch):
    text = _make_scoring_input(
        tmp_path, monkeypatch,
        return_all=True,
    )
    assert _active_usrbdx_particles(text) == {"ALL-PART"}


def test_fluka_input_return_all_charged(tmp_path, monkeypatch):
    text = _make_scoring_input(
        tmp_path, monkeypatch,
        return_all_charged=True,
    )
    assert _active_usrbdx_particles(text) == {"ALL-CHAR"}

