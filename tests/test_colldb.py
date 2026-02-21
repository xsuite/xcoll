# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

import json
from pathlib import Path

import xcoll as xc


path = Path(__file__).parent / 'data'


def test_loading():
    colldb_1 = xc.CollimatorDatabase.from_yaml(path / 'colldb_lhc_run3.yaml', beam=1, ignore_crystals=False)
    colldb_2 = xc.CollimatorDatabase.from_yaml(path / 'colldb_lhc_run3_b1.yaml', ignore_crystals=False)
    colldb_3 = xc.CollimatorDatabase.from_json(path / 'colldb_lhc_run3_b1.json', ignore_crystals=False)
    with open(path / 'colldb_lhc_run3_b1.json', 'r') as fid:
        colldb_4 = xc.CollimatorDatabase.from_dict(json.load(fid), ignore_crystals=False)
    df1 = colldb_1.to_pandas().reindex(sorted(colldb_1.to_pandas().columns), axis=1)
    df2 = colldb_2.to_pandas().reindex(sorted(colldb_2.to_pandas().columns), axis=1)
    df3 = colldb_3.to_pandas().reindex(sorted(colldb_3.to_pandas().columns), axis=1)
    df4 = colldb_4.to_pandas().reindex(sorted(colldb_4.to_pandas().columns), axis=1)
    assert df1.equals(df2)
    assert df2.equals(df3)
    assert df3.equals(df4)


def test_loading_no_families():
    colldb_1 = xc.CollimatorDatabase.from_yaml(path / 'colldb_lhc_run3.yaml', beam=1, ignore_crystals=False)
    colldb_2 = xc.CollimatorDatabase.from_yaml(path / 'colldb_lhc_run3_b1_no_families.yaml', nemitt_x=3.5e-6, nemitt_y=3.5e-6, ignore_crystals=False)
    df1 = colldb_1.to_pandas().reindex(sorted(colldb_1.to_pandas().columns), axis=1)
    df2 = colldb_2.to_pandas().reindex(sorted(colldb_2.to_pandas().columns), axis=1)
    df1_no_fam = df1.drop(['family', 'overwritten_keys'], axis=1)
    df2_no_fam = df2.drop(['family', 'overwritten_keys'], axis=1)
    assert df1_no_fam.equals(df2_no_fam)


def test_loading_no_merge():
    colldb_1 = xc.CollimatorDatabase.from_yaml(path / 'colldb_lhc_run3.yaml', beam=1, ignore_crystals=False)
    colldb_2 = xc.CollimatorDatabase.from_yaml(path / 'colldb_lhc_run3_b1_no_merge.yaml', ignore_crystals=False)
    df1 = colldb_1.to_pandas().reindex(sorted(colldb_1.to_pandas().columns), axis=1)
    df2 = colldb_2.to_pandas().reindex(sorted(colldb_2.to_pandas().columns), axis=1)
    assert df1.equals(df2)


def test_loading_crystals():
    colldb_1 = xc.CollimatorDatabase.from_yaml(path / 'colldb_lhc_run3.yaml', beam=1)
    colldb_2 = xc.CollimatorDatabase.from_yaml(path / 'colldb_lhc_run3.yaml', beam=1, ignore_crystals=False)
    df1 = colldb_1.to_pandas().reindex(sorted(colldb_1.to_pandas().columns), axis=1)
    df2 = colldb_2.to_pandas().reindex(sorted(colldb_2.to_pandas().columns), axis=1)
    df2_no_cry = df2.drop(['tcpch.a4l7.b1', 'tcpcv.a6l7.b1'], axis=0)
    assert df2_no_cry.equals(df1)


def test_loading_SixTrack():
    colldb_1 = xc.CollimatorDatabase.from_yaml(path / 'colldb_lhc_run3.yaml', beam=1)
    colldb_2 = xc.CollimatorDatabase.from_SixTrack(path / 'colldb_lhc_run3_b1.dat', nemitt_x=3.5e-6, nemitt_y=3.5e-6)
    df1 = colldb_1.to_pandas().reindex(sorted(colldb_1.to_pandas().columns), axis=1)
    df2 = colldb_2.to_pandas().reindex(sorted(colldb_2.to_pandas().columns), axis=1)
    # In SixTrack loader, families, and tilts are not loaded
    df1 = df1.drop(['overwritten_keys', 'tilt'], axis=1)
    df2 = df2.drop(['overwritten_keys', 'tilt'], axis=1)
    # In SixTrack loader, if a collimator has an individual gap setting, no family nor stage is assigned
    df2.loc['tcla.a6r7.b1', 'family'] = 'tcla7'
    df2.loc['tcla.a6r7.b1', 'stage'] = 'tertiary'
    # In SixTrack loader, non-active collimators default to active (but miss family keys etc)
    df2.loc[['tcsg.b4l7.b1', 'tcsg.e5r7.b1', 'tcsg.6r7.b1'], 'active'] = False
    df2.loc[['tcsg.b4l7.b1', 'tcsg.e5r7.b1', 'tcsg.6r7.b1'], 'family'] = 'tcsg7'
    df2.loc[['tcsg.b4l7.b1', 'tcsg.e5r7.b1', 'tcsg.6r7.b1'], 'gap'] = 6.5
    df2.loc[['tcsg.b4l7.b1', 'tcsg.e5r7.b1', 'tcsg.6r7.b1'], 'stage'] = 'secondary'
    # In SixTrack loader, information for parking is missing
    df2.loc[['tctph.4l1.b1', 'tctpv.4l1.b1', 'tctph.4l5.b1', 'tctpv.4l5.b1'], 'parking'] = 0.02
    df1.sort_index(axis=0, inplace=True)
    df2.sort_index(axis=0, inplace=True)
    assert df1.equals(df2)


def test_loading_SixTrack_crystals():
    colldb_1 = xc.CollimatorDatabase.from_yaml(path / 'colldb_lhc_run3.yaml', beam=1, ignore_crystals=False)
    colldb_2 = xc.CollimatorDatabase.from_SixTrack(path / 'colldb_lhc_run3_b1.dat', nemitt_x=3.5e-6, nemitt_y=3.5e-6, ignore_crystals=False)
    df1 = colldb_1.to_pandas().reindex(sorted(colldb_1.to_pandas().columns), axis=1)
    df2 = colldb_2.to_pandas().reindex(sorted(colldb_2.to_pandas().columns), axis=1)
    # In SixTrack loader, families, and tilts are not loaded
    df1 = df1.drop(['overwritten_keys', 'tilt'], axis=1)
    df2 = df2.drop(['overwritten_keys', 'tilt'], axis=1)
    # In SixTrack loader, if a collimator has an individual gap setting, no family nor stage is assigned
    df2.loc['tcla.a6r7.b1', 'family'] = 'tcla7'
    df2.loc['tcla.a6r7.b1', 'stage'] = 'tertiary'
    # In SixTrack loader, non-active collimators default to active (but miss family keys etc)
    df2.loc[['tcsg.b4l7.b1', 'tcsg.e5r7.b1', 'tcsg.6r7.b1'], 'active'] = False
    df2.loc[['tcsg.b4l7.b1', 'tcsg.e5r7.b1', 'tcsg.6r7.b1'], 'family'] = 'tcsg7'
    df2.loc[['tcsg.b4l7.b1', 'tcsg.e5r7.b1', 'tcsg.6r7.b1'], 'gap'] = 6.5
    df2.loc[['tcsg.b4l7.b1', 'tcsg.e5r7.b1', 'tcsg.6r7.b1'], 'stage'] = 'secondary'
    # In SixTrack loader, information for parking is missing
    df2.loc[['tctph.4l1.b1', 'tctpv.4l1.b1', 'tctph.4l5.b1', 'tctpv.4l5.b1'], 'parking'] = 0.02
    df1.sort_index(axis=0, inplace=True)
    df2.sort_index(axis=0, inplace=True)
    assert df1.equals(df2)


def test_dumping():
    colldb_1 = xc.CollimatorDatabase.from_yaml(path / 'colldb_lhc_run3.yaml', beam=1)
    colldb_2 = xc.CollimatorDatabase.from_yaml(path / 'colldb_lhc_run3.yaml', beam=2)
    colldb_1.to_yaml('out1')
    colldb_2.to_yaml('out2')
    colldb_3 = xc.CollimatorDatabase.from_yaml('out1.yaml', beam=1)
    colldb_4 = xc.CollimatorDatabase.from_yaml('out2.yaml', beam=2)
    df1 = colldb_1.to_pandas().reindex(sorted(colldb_1.to_pandas().columns), axis=1)
    df2 = colldb_2.to_pandas().reindex(sorted(colldb_2.to_pandas().columns), axis=1)
    df3 = colldb_3.to_pandas().reindex(sorted(colldb_3.to_pandas().columns), axis=1)
    df4 = colldb_4.to_pandas().reindex(sorted(colldb_4.to_pandas().columns), axis=1)
    assert df1.equals(df3)
    assert df2.equals(df4)
    (Path.cwd() / 'out1.yaml').unlink()
    (Path.cwd() / 'out2.yaml').unlink()
