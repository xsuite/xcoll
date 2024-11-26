# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

import json
import numpy as np
from pathlib import Path
import xtrack as xt
import xcoll as xc
import xpart as xp
import pytest
from xpart.test_helpers import flaky_assertions, retry
from xobjects.test_helpers import for_all_test_contexts

#path = Path(__file__).parent / '..' / 'data'
path = xc._pkg_root.parent / 'examples'

@pytest.mark.parametrize("beam, plane, npart, interpolation, ignore_crystals", [
                            [1, 'H', 10, 0.2, True],
                            [1, 'V', 100, 0.3, True]
                            # [2, 'H', 1000, 0.15, True],
                            # [2, 'V', 10000, 0.1, True]
                        ], ids=["B1H", "B1V"])#, "B2H", "B2V"])
def test_run_lossmap(beam, plane, npart, interpolation, ignore_crystals):
    
    ### XXX Crushes at the second attempt since the particle reference is not updated properly.
    import xcoll as xc
    # import importlib

    # importlib.reload(xc)
    # importlib.reload(xp)

    #xc.FlukaEngine = xc.FlukaEngine()

    # Load from json
    line = xt.Line.from_json(path / 'machines' / f'lhc_run3_b{beam}.json')

    colldb = xc.CollimatorDatabase.from_yaml(path / 'colldb' / f'lhc_run3_fluka.yaml', beam=beam)


    colldb.install_fluka_collimators(line=line)
    df_with_coll = line.check_aperture()
    assert not np.any(df_with_coll.has_aperture_problem)

    line.build_tracker()
    line.collimators.assign_optics()

    # Connect to FLUKA
    xc.FlukaEngine.start(line=line, _capacity=2*npart, cwd='run_fluka_temp', debug_level=1)
    particle_ref = xp.Particles.reference_from_pdg_id(pdg_id='proton', p0c=6.8e12)
    # import pdb; pdb.set_trace()
    particle_ref.mass0 = 938272310.0 # 
    line.particle_ref.mass0 = 938272310.0
    line.particle_ref = particle_ref
    print(f"Particle reference: {particle_ref}")
    xc.FlukaEngine.set_particle_ref(particle_ref=particle_ref, line=line)
    print(f"Particle reference: {particle_ref}")

    # Generate initial pencil distribution
    tcp  = f"tcp.{'c' if plane=='H' else 'd'}6{'l' if beam==1 else 'r'}7.b{beam}"
    part = line[tcp].generate_pencil(npart)

    # Track!
    line.scattering.enable()
    line.track(part, num_turns=2, time=True, with_progress=1)
    line.scattering.disable()
    print(f"Done tracking in {line.time_last_track:.1f}s.")

    xc.FlukaEngine.stop()

    line_is_reversed = True if beam == 2 else False
    with flaky_assertions():

        ThisLM = xc.LossMap(line, line_is_reversed=line_is_reversed, part=part,
                         interpolation=interpolation)

        ThisLM.to_json("lossmap.json")
        assert Path("lossmap.json").exists()
        with Path("lossmap.json").open('r') as fid:
            dct = json.load(fid)
            assert xt.line._dicts_equal(dct, ThisLM.lossmap)
        Path("lossmap.json").unlink()
        ThisLM.save_summary("coll_summary.txt")
        assert Path("coll_summary.txt").exists()
        Path("coll_summary.txt").unlink()

        # TODO: check the lossmap quantitatively: rough amount of losses at given positions
        summ = ThisLM.summary
        assert list(summ.columns) == ['collname', 'nabs', 'length', 's', 'type']

        # We want at least 5% absorption on the primary
        assert summ.loc[summ.collname==tcp,'nabs'].values[0] > 0.05*npart

        lm = ThisLM.lossmap
        summ = summ[summ.nabs > 0]
        assert list(lm.keys()) == ['collimator', 'aperture', 'machine_length', \
                                   'interpolation', 'reversed']
        assert list(lm['collimator'].keys()) == ['s', 'name', 'length', 'n']
        assert len(lm['collimator']['s']) == len(summ)
        assert len(lm['collimator']['name']) == len(summ)
        assert len(lm['collimator']['length']) == len(summ)
        assert len(lm['collimator']['n']) == len(summ)
        assert np.all(lm['collimator']['s'] == summ.s)
        assert np.all(lm['collimator']['name'] == summ.collname)
        assert np.all(lm['collimator']['length'] == summ.length)
        assert np.all(lm['collimator']['n'] == summ.nabs)
        assert np.all([nn[:3] in ['tcp', 'tcs'] for nn in lm['collimator']['name']])
        assert np.all([s < lm['machine_length'] for s in lm['collimator']['s']])
        assert list(lm['aperture'].keys()) == ['s', 'name', 'n']
        if npart > 5000:
            assert len(lm['aperture']['s']) > 0
            assert len(lm['aperture']['s']) == len(lm['aperture']['name'])
            assert len(lm['aperture']['s']) == len(lm['aperture']['n'])
            assert np.all([s < lm['machine_length'] for s in lm['aperture']['s']])
        assert lm['interpolation'] == interpolation
        # assert lm['reversed'] == line_is_reversed