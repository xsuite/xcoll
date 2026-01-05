# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2026.                 #
# ######################################### #

import numpy as np
from pathlib import Path
import xtrack as xt
import xcoll as xc
import pytest
import shutil

from xcoll.scattering_routines.fluka.fluka_input import get_collimators_from_input_file


@pytest.mark.fluka
@pytest.mark.parametrize("el_type", ['collimator', 'crystal'])
def test_fluka_input_single(el_type):
    if el_type == 'collimator':
        coll = xc.FlukaCollimator(length=0.456, angle=32, jaw=[0.01, -0.02], tilt=[10e-6, -8.7e-6], material='Yttrium')
    else:
        coll = xc.FlukaCrystal(length=0.002, side='-', angle=90, jaw=-0.01, tilt=43e-6, material='Yttrium', bending_radius=65)
    with pytest.raises(ValueError, match="Need to provide either a line with a reference particle, or `particle_ref`."):
        input_file = xc.fluka.engine.generate_input_file(elements=coll, names='TestColl', clean=False)
    path_tmp = Path.cwd() / f'temp_fluka_test_single_{el_type}'
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
        assert np.isclose(new_coll[name]['length'],  1.014)
        assert np.isclose(new_coll[name]['angle'],  32)
        assert np.isclose(new_coll[name]['jaw'][0], 0.01)
        assert np.isclose(new_coll[name]['jaw'][1], -0.02)
        assert np.isclose(new_coll[name]['tilt'][0], 10e-6)
        assert np.isclose(new_coll[name]['tilt'][1], -8.7e-6)
    else:
        assert np.isclose(new_coll[name]['length'],  0.106)
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
    xc.FlukaPrototype.inspect_prototypes_file(path_tmp / 'prototypes.lbp')
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
            if "     1          INROT_1              INROT_1             0.279000" in line \
            and el_type == 'collimator':
                found = True
            if "     1          INROT_1              INROT_1             0.052000" in line \
            and el_type == 'crystal':
                found = True
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
    assert found_1 and found_2 and found_3

    # Check crystal definition
    if el_type == 'crystal':
        found_1 = False
        found_2 = False
        found_3 = False
        found_4 = False
        with input_file[0].open('r') as fp:
            for line in fp:
                if f"CRYSTAL     {jaw.fedb_tag}  0.030769       0.2       0.0       0.0     300.0 110" in line:
                    found_1 = True
                if "CRYSTAL          0.0      -1.0       0.0       0.0       0.0       1.0 &" in line:
                    found_2 = True
                if "CRYSTAL   -1049.9999   -3000.0    1000.0                              &&" in line:
                    found_3 = True
                if "USRICALL        50.0                                                  CRYSTAL" in line:
                    found_4 = True
        assert found_1 and found_2 and found_3 and found_4

    # Clean up
    shutil.rmtree(path_tmp)


@pytest.mark.fluka
@pytest.mark.parametrize("ignore_crystals", [True, False])
def test_fluka_input_line(ignore_crystals):
    beam = 1
    path = Path(__file__).parent
    env = xt.load(path / 'data' / f'sequence_lhc_run3_b{beam}.json')
    line = env.lines[f'lhcb{beam}']
    colldb = xc.CollimatorDatabase.from_yaml(path / 'data' / f'colldb_lhc_run3_ir7.yaml', beam=beam,
                                             ignore_crystals=ignore_crystals)
    colldb.install_fluka_collimators(line=line, verbose=True)
    colls, _ = line.get_elements_of_type(xc.collimator_classes)
    line.build_tracker()
    line.collimators.assign_optics()
    line.collimators.align_to_beam_divergence()
    path_tmp = Path.cwd() / f'temp_fluka_test_line_{ignore_crystals}'
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
            'tcsg.a5l7.b1':  {'length': 1.482, 'angle': 40.7,  'tilt': [0.0, 0.0], 'jaw': [0.0018463405948736522, -0.0018520014167018317]},
            'tcsg.d4l7.b1':  {'length': 1.482, 'angle': 90.0,  'tilt': [0.0, 0.0], 'jaw': [0.0011923168368670467, -0.0011930737268484037]},
            'tcspm.b4l7.b1': {'length': 1.482, 'angle': 0.0,   'tilt': [0.0, 0.0], 'jaw': [0.0016543438662348642, -0.0016608919799403488]},
            'tcsg.a4l7.b1':  {'length': 1.482, 'angle': 134.6, 'tilt': [0.0, 0.0], 'jaw': [0.0016554613686823316, -0.0016540730713363594]},
            'tcsg.a4r7.b1':  {'length': 1.482, 'angle': 46.3,  'tilt': [0.0, 0.0], 'jaw': [0.0016562125801096172, -0.0016637353535826627]}
        }
    else:
        this_coll_dct = {
            'tcp.d6l7.b1':   {'length': 1.482, 'angle': 90.0,  'tilt': [0.0, 0.0], 'jaw': [0.0009238348692113263, -0.0009219591453732967]},
            'tcp.c6l7.b1':   {'length': 1.482, 'angle': 0.0,   'tilt': [2.5e-06, -2.5e-06], 'jaw': [0.0013138622732093985, -0.0013103666448213147]},
            'tcp.b6l7.b1':   {'length': 1.482, 'angle': 127.5, 'tilt': [0.0, 0.0], 'jaw': [0.0010997388044393652, -0.0011003192387097904]},
            'tcsg.a6l7.b1':  {'length': 1.482, 'angle': 141.1, 'tilt': [0.0, 0.0], 'jaw': [0.001474936603826471, -0.0014740571820550663]},
            'tcpcv.a6l7.b1': {'length': 0.11,  'angle': 90.0,  'tilt': [1.66577e-05, None], 'jaw': [0.0017621469072926072, None]},
            'tcsg.b5l7.b1':  {'length': 1.482, 'angle': 143.5, 'tilt': [0.0, 0.0], 'jaw': [0.001814181900504419, -0.0018086857451153904]},
            'tcsg.a5l7.b1':  {'length': 1.482, 'angle': 40.7,  'tilt': [0.0, 0.0], 'jaw': [0.0018463405953683676, -0.001852001417195659]},
            'tcsg.d4l7.b1':  {'length': 1.482, 'angle': 90.0,  'tilt': [0.0, 0.0], 'jaw': [0.0011923168370637782, -0.0011930737270451353]},
            'tcpch.a4l7.b1': {'length': 0.11,  'angle': 0.0,   'tilt': [1.21108e-05, None], 'jaw': [0.0020205086648474287, None]},
            'tcspm.b4l7.b1': {'length': 1.482, 'angle': 0.0,   'tilt': [0.0, 0.0], 'jaw': [0.0016543438666660748, -0.0016608919803706712]},
            'tcsg.a4l7.b1':  {'length': 1.482, 'angle': 134.6, 'tilt': [0.0, 0.0], 'jaw': [0.001655461368718747, -0.0016540730713732188]},
            'tcsg.a4r7.b1':  {'length': 1.482, 'angle': 46.3,  'tilt': [0.0, 0.0], 'jaw': [0.001656212580093186, -0.0016637353535662314]}
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
    xc.FlukaPrototype.inspect_prototypes_file(path_tmp / 'prototypes.lbp')
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
     6          INROT_6              INROT_6             0.241000   
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
     5          INROT_5              INROT_5             0.053000   
     6          INROT_6              INROT_6             0.241000   
     7          INROT_7              INROT_7             0.241000   
     8          INROT_8              INROT_8             0.241000   
     9          INROT_9              INROT_9             0.053000   
     10         INROT_10             INROT_10            0.241000   
     11         INROT_11             INROT_11            0.241000   
     12         INROT_12             INROT_12            0.241000""" in insertion_txt

    # Clean up
    shutil.rmtree(path_tmp)
