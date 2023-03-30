import json
from pathlib import Path
import xcoll as xc


path = Path.cwd() / 'data'


def test_loading():
    colldb_1 = xc.CollDB.from_yaml(path / 'colldb_lhc_run3.yaml', beam=1)
    colldb_2 = xc.CollDB.from_yaml(path / 'colldb_lhc_run3_b1.yaml')
    colldb_3 = xc.CollDB.from_yaml(path / 'colldb_lhc_run3_b1_no_merge.yaml')
    colldb_4 = xc.CollDB.from_json(path / 'colldb_lhc_run3_b1.json')
    with open(path / 'colldb_lhc_run3_b1.json', 'r') as fid:
        colldb_5 = xc.CollDB.from_dict(json.load(fid))
    colldb_6 = xc.CollDB.from_yaml(path / 'colldb_lhc_run3_b1_no_families.yaml', nemitt_x=3.5e-6, nemitt_y=3.5e-6)
    colldb_7 = xc.CollDB.from_SixTrack(path / 'colldb_lhc_run3_b1.dat', nemitt_x=3.5e-6, nemitt_y=3.5e-6)
    df1 = colldb_1._colldb.reindex(sorted(colldb_1._colldb.columns), axis=1)
    df2 = colldb_2._colldb.reindex(sorted(colldb_2._colldb.columns), axis=1)
    df3 = colldb_3._colldb.reindex(sorted(colldb_3._colldb.columns), axis=1)
    df4 = colldb_4._colldb.reindex(sorted(colldb_4._colldb.columns), axis=1)
    df5 = colldb_5._colldb.reindex(sorted(colldb_5._colldb.columns), axis=1)
    df6 = colldb_6._colldb.reindex(sorted(colldb_6._colldb.columns), axis=1)
    df7 = colldb_7._colldb.reindex(sorted(colldb_7._colldb.columns), axis=1)
    assert df1.equals(df2)
    assert df2.equals(df3)
    assert df3.equals(df4)
    assert df4.equals(df5)
    df5_no_fam = df5.drop(['family', 'overwritten_keys'], axis=1)
    df6_no_fam = df6.drop(['family', 'overwritten_keys'], axis=1)
    assert df5_no_fam.equals(df6_no_fam)
    # Crystals not supported in SixTrack loader, and non-active collimators default to active
    df6_only_active_no_cry = df6.drop(['tcsg.b4l7.b1', 'tcsg.e5r7.b1', 'tcsg.6r7.b1', 'tcpch.a4l7.b1', 'tcpcv.a6l7.b1'], axis=0)
    df6_only_active_no_cry = df6_only_active_no_cry.drop('parking', axis=1)
    df7_only_active        = df7.drop(['tcsg.b4l7.b1', 'tcsg.e5r7.b1', 'tcsg.6r7.b1'], axis=0)
    df7_only_active        = df7_only_active.drop('parking', axis=1)
    assert df6_only_active_no_cry.equals(df7_only_active)

