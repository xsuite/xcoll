# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

import json
from pathlib import Path
import xcoll as xc
import pandas as pd


path =  Path(__file__).parent / 'data'


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
    import numpy as np
    colldb_1 = xc.CollimatorDatabase.from_yaml(path / 'colldb_lhc_run3.yaml', beam=1)
    colldb_2 = xc.CollimatorDatabase.from_SixTrack(path / 'colldb_lhc_run3_b1.dat', nemitt_x=3.5e-6, nemitt_y=3.5e-6)
    df1 = colldb_1.to_pandas().reindex(sorted(colldb_1.to_pandas().columns), axis=1)
    df2 = colldb_2.to_pandas().reindex(sorted(colldb_2.to_pandas().columns), axis=1)
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
    df1 = colldb_1.to_pandas().reindex(sorted(colldb_1.to_pandas().columns), axis=1)
    df2 = colldb_2.to_pandas().reindex(sorted(colldb_2.to_pandas().columns), axis=1)
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


def test_dumping_from_Sixtrack():
    colldb_yaml = xc.CollimatorDatabase.from_yaml(path / 'colldb_lhc_run3.yaml', beam=1,ignore_crystals=False)
    colldb_dat = xc.CollimatorDatabase.from_SixTrack(path / 'colldb_lhc_run3_b1.dat', nemitt_x=3.5e-6, nemitt_y=3.5e-6, ignore_crystals=False)
    colldb_dat.to_yaml('ne.to_pandas().lhc_run3')

    colldb_yaml_new = xc.CollimatorDatabase.from_yaml('ne.to_pandas().lhc_run3.yaml', beam=1,ignore_crystals=False)
    df_yaml = colldb_yaml.to_pandas().reindex(sorted(colldb_yaml.to_pandas().columns), axis=1)
    df_dat = colldb_dat.to_pandas().reindex(sorted(colldb_dat.to_pandas().columns), axis=1)
    df_yaml_new = colldb_yaml_new.to_pandas().reindex(sorted(colldb_yaml_new.to_pandas().columns), axis=1)

    df_yaml = df_yaml.drop(['family', 'overwritten_keys', 'parking'], axis=1)
    df_yaml_new = df_yaml_new.drop(['family', 'overwritten_keys', 'parking'], axis=1)
    df_dat = df_dat.drop(['family', 'overwritten_keys', 'parking'], axis=1)

    df_yaml_new = df_yaml_new.drop(['tcsg.b4l7.b1', 'tcsg.e5r7.b1', 'tcsg.6r7.b1'], axis=0)
    df_yaml = df_yaml.drop(['tcsg.b4l7.b1', 'tcsg.e5r7.b1', 'tcsg.6r7.b1'], axis=0)
    df_dat = df_dat.drop(['tcsg.b4l7.b1', 'tcsg.e5r7.b1', 'tcsg.6r7.b1'], axis=0)

    df_yaml_new.sort_index(axis=0, inplace=True)
    df_dat.sort_index(axis=0, inplace=True)
    df_yaml.sort_index(axis=0, inplace=True)

    assert df_dat.equals(df_yaml)
    assert df_yaml_new.equals(df_yaml)
    assert df_yaml_new.equals(df_dat)
    (Path.cwd() / 'ne.to_pandas().lhc_run3.yaml').unlink()