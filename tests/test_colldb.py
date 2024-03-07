# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

import json
from pathlib import Path
import xcoll as xc


path =  Path(__file__).parent / 'data'


def test_loading():
    colldb_1 = xc.CollimatorDatabase.from_yaml(path / 'colldb_lhc_run3.yaml', beam=1, ignore_crystals=False)
    colldb_2 = xc.CollimatorDatabase.from_yaml(path / 'colldb_lhc_run3_b1.yaml', ignore_crystals=False)
    colldb_3 = xc.CollimatorDatabase.from_json(path / 'colldb_lhc_run3_b1.json', ignore_crystals=False)
    with open(path / 'colldb_lhc_run3_b1.json', 'r') as fid:
        colldb_4 = xc.CollimatorDatabase.from_dict(json.load(fid), ignore_crystals=False)
    df1 = colldb_1._colldb.reindex(sorted(colldb_1._colldb.columns), axis=1)
    df2 = colldb_2._colldb.reindex(sorted(colldb_2._colldb.columns), axis=1)
    df3 = colldb_3._colldb.reindex(sorted(colldb_3._colldb.columns), axis=1)
    df4 = colldb_4._colldb.reindex(sorted(colldb_4._colldb.columns), axis=1)
    assert df1.equals(df2)
    assert df2.equals(df3)
    assert df3.equals(df4)


def test_loading_no_families():
    colldb_1 = xc.CollimatorDatabase.from_yaml(path / 'colldb_lhc_run3.yaml', beam=1, ignore_crystals=False)
    colldb_2 = xc.CollimatorDatabase.from_yaml(path / 'colldb_lhc_run3_b1_no_families.yaml', nemitt_x=3.5e-6, nemitt_y=3.5e-6, ignore_crystals=False)
    df1 = colldb_1._colldb.reindex(sorted(colldb_1._colldb.columns), axis=1)
    df2 = colldb_2._colldb.reindex(sorted(colldb_2._colldb.columns), axis=1)
    df1_no_fam = df1.drop(['family', 'overwritten_keys'], axis=1)
    df2_no_fam = df2.drop(['family', 'overwritten_keys'], axis=1)
    assert df1_no_fam.equals(df2_no_fam)


def test_loading_no_merge():
    colldb_1 = xc.CollimatorDatabase.from_yaml(path / 'colldb_lhc_run3.yaml', beam=1, ignore_crystals=False)
    colldb_2 = xc.CollimatorDatabase.from_yaml(path / 'colldb_lhc_run3_b1_no_merge.yaml', ignore_crystals=False)
    df1 = colldb_1._colldb.reindex(sorted(colldb_1._colldb.columns), axis=1)
    df2 = colldb_2._colldb.reindex(sorted(colldb_2._colldb.columns), axis=1)
    assert df1.equals(df2)

def test_loading_crystals():
    colldb_1 = xc.CollimatorDatabase.from_yaml(path / 'colldb_lhc_run3.yaml', beam=1)
    colldb_2 = xc.CollimatorDatabase.from_yaml(path / 'colldb_lhc_run3.yaml', beam=1, ignore_crystals=False)
    df1 = colldb_1._colldb.reindex(sorted(colldb_1._colldb.columns), axis=1)
    df2 = colldb_2._colldb.reindex(sorted(colldb_2._colldb.columns), axis=1)
    df2_no_cry = df2.drop(['tcpch.a4l7.b1', 'tcpcv.a6l7.b1'], axis=0)
    assert df2_no_cry.equals(df1)


def test_loading_SixTrack():
    import numpy as np
    colldb_1 = xc.CollimatorDatabase.from_yaml(path / 'colldb_lhc_run3.yaml', beam=1)
    colldb_2 = xc.CollimatorDatabase.from_SixTrack(path / 'colldb_lhc_run3_b1.dat', nemitt_x=3.5e-6, nemitt_y=3.5e-6)
    df1 = colldb_1._colldb.reindex(sorted(colldb_1._colldb.columns), axis=1)
    df2 = colldb_2._colldb.reindex(sorted(colldb_2._colldb.columns), axis=1)
    # In SixTrack loader, families are not (yet) loaded
    df1 = df1.drop(['family', 'overwritten_keys'], axis=1)
    df2 = df2.drop(['family', 'overwritten_keys'], axis=1)
    # In SixTrack loader, non-active collimators default to active
    df1_only_active = df1.drop(['tcsg.b4l7.b1', 'tcsg.e5r7.b1', 'tcsg.6r7.b1'], axis=0)
    df1_only_active = df1_only_active.drop('parking', axis=1)
    df2_only_active = df2.drop(['tcsg.b4l7.b1', 'tcsg.e5r7.b1', 'tcsg.6r7.b1'], axis=0)
    df2_only_active = df2_only_active.drop('parking', axis=1)
    assert df1_only_active.equals(df2_only_active)


def test_loading_SixTrack_crystals():
    import numpy as np
    colldb_1 = xc.CollimatorDatabase.from_yaml(path / 'colldb_lhc_run3.yaml', beam=1, ignore_crystals=False)
    colldb_2 = xc.CollimatorDatabase.from_SixTrack(path / 'colldb_lhc_run3_b1.dat', nemitt_x=3.5e-6, nemitt_y=3.5e-6, ignore_crystals=False)
    df1 = colldb_1._colldb.reindex(sorted(colldb_1._colldb.columns), axis=1)
    df2 = colldb_2._colldb.reindex(sorted(colldb_2._colldb.columns), axis=1)
    # In SixTrack loader, families are not (yet) loaded
    df1 = df1.drop(['family', 'overwritten_keys'], axis=1)
    df2 = df2.drop(['family', 'overwritten_keys'], axis=1)
    # In SixTrack loader, non-active collimators default to active
    df1_only_active = df1.drop(['tcsg.b4l7.b1', 'tcsg.e5r7.b1', 'tcsg.6r7.b1'], axis=0)
    df1_only_active = df1_only_active.drop('parking', axis=1)
    df2_only_active = df2.drop(['tcsg.b4l7.b1', 'tcsg.e5r7.b1', 'tcsg.6r7.b1'], axis=0)
    df2_only_active = df2_only_active.drop('parking', axis=1)
    df2_only_active.sort_index(axis=0, inplace=True)
    df1_only_active.sort_index(axis=0, inplace=True)
    assert df1_only_active.equals(df2_only_active)


def test_dumping():
    colldb_1 = xc.CollimatorDatabase.from_yaml(path / 'colldb_lhc_run3.yaml', beam=1)
    colldb_2 = xc.CollimatorDatabase.from_yaml(path / 'colldb_lhc_run3.yaml', beam=2)
    colldb_1.write_to_yaml('out1')
    colldb_2.write_to_yaml('out2')
    colldb_3 = xc.CollimatorDatabase.from_yaml('out1.yaml', beam=1)
    colldb_4 = xc.CollimatorDatabase.from_yaml('out2.yaml', beam=2)
    df1 = colldb_1._colldb.reindex(sorted(colldb_1._colldb.columns), axis=1)
    df2 = colldb_2._colldb.reindex(sorted(colldb_2._colldb.columns), axis=1)
    df3 = colldb_3._colldb.reindex(sorted(colldb_3._colldb.columns), axis=1)
    df4 = colldb_4._colldb.reindex(sorted(colldb_4._colldb.columns), axis=1)
    assert df1.equals(df3)
    assert df2.equals(df4)
    (Path.cwd() / 'out1.yaml').unlink()
    (Path.cwd() / 'out2.yaml').unlink()

    