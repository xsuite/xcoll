import json
import numpy as np
from pathlib import Path
import xtrack as xt
import xcoll as xc



def _run_lossmap(beam, plane, npart, interpolation):

    with open(Path(Path.cwd(),'data',f'lhc_run3_b{beam}.json'), 'r') as fid:
        loaded_dct = json.load(fid)
    line = xt.Line.from_dict(loaded_dct)

    line_is_reversed = True if beam=='2' else False
    coll_manager = xc.CollimatorManager(
        line=line, line_is_reversed=line_is_reversed,
        colldb=xc.load_SixTrack_colldb(Path(Path.cwd(),'data',f'lhc_run3_b{beam}.dat'), emit=3.5e-6)
        )

    coll_manager.install_everest_collimators()
    tracker = coll_manager.build_tracker()
    coll_manager.set_openings()

    tcp  = f"tcp.{'c' if plane=='H' else 'd'}6{'l' if beam=='1' else 'r'}7.b{beam}"
    part = coll_manager.generate_pencil_on_collimator(tcp, num_particles=npart, zeta=0.0, delta=0.0)

    coll_manager.enable_scattering()
    tracker.track(part, num_turns=2)
    coll_manager.disable_scattering()

    coll_manager.create_lossmap(part, interpolation=interpolation)

    summ = coll_manager.summary
    assert list(summ.columns) == ['collname', 'nabs', 'length', 's']
    assert len(summ) == 10
    assert summ.loc[summ.collname==tcp,'nabs'].values[0] > 0.1*npart

    lm = coll_manager.lossmap
    assert list(lm.keys()) == ['collimator', 'aperture', 'machine_length', 'interpolation', 'reversed']
    assert list(lm['collimator'].keys()) == ['s', 'name', 'length']
    assert len(lm['collimator']['s']) == summ.nabs.sum()
    assert len(lm['collimator']['name']) == summ.nabs.sum()
    assert len(lm['collimator']['length']) == summ.nabs.sum()
    assert np.all([ s == summ.loc[summ.collname==name,'s'].values[0]
              for name,s in list(zip(lm['collimator']['name'],lm['collimator']['s'])) ])
    assert np.all([ l == summ.loc[summ.collname==name,'length'].values[0]
              for name,l in list(zip(lm['collimator']['name'],lm['collimator']['length'])) ])
    assert np.all([nn[:3] in ['tcp', 'tcs'] for nn in lm['collimator']['name']])
    assert np.all([s < lm['machine_length'] for s in lm['collimator']['s']])
    assert list(lm['aperture'].keys()) == ['s', 'name']
    assert np.all([s < lm['machine_length'] for s in lm['aperture']['s']])
    assert lm['interpolation'] == interpolation
    assert lm['reversed'] == line_is_reversed


def test_lossmap_B1H():
    _run_lossmap('1', 'H', 25000, 0.2)

def test_lossmap_B2V():
    _run_lossmap('2', 'V', 25000, 0.3)

#     # Save to json
#     # These files can be loaded, combined (for more statistics), and plotted with the 'lossmaps' package
#     with open(Path(path_out,f'lossmap_B{beam+plane}.json'), 'w') as fid:
#         json.dump(coll_manager.lossmap, fid, indent=True)

#     # Save a summary of the collimator losses to a text file
#     with open(Path(path_out,f'coll_summary_B{beam+plane}.out'), 'w') as fid:
#         fid.write(coll_manager.summary.__repr__())
#     print(summary)

#     exit()
