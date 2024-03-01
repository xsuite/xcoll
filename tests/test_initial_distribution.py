import json
import numpy as np
from pathlib import Path
import xtrack as xt
import xcoll as xc
import pytest
import math
from xpart.test_helpers import flaky_assertions, retry
from xobjects.test_helpers import for_all_test_contexts

path = Path(__file__).parent / 'data'

@for_all_test_contexts(
    excluding=('ContextCupy', 'ContextPyopencl')  # Rutherford RNG not on GPU
)
@retry()
@pytest.mark.parametrize("beam, npart, longitudinal, longitudinal_betatron_cut", [
                        [2, 40000, None, None],
                        [1, 40000, 'bucket', None],
                        [1, 40000, 'matched_dispersion', 5]]
                        , ids=["B1", "B1_bucket", "B2_matched"])

def test_create_initial_distribution(beam, npart,longitudinal, longitudinal_betatron_cut, test_context):
    
    line = xt.Line.from_json(path / f'sequence_lhc_run3_b{beam}.json')
    coll_manager = xc.CollimatorManager.from_yaml(path / f'colldb_lhc_run3_ir7.yaml',
                                                  line=line, beam=beam,_context=test_context)
    
    coll_manager.install_everest_collimators()
    coll_manager.build_tracker()
    coll_manager.set_openings()

    tcp_conv = f"tcp.c6{'l' if beam == 1 else 'r'}7.b{beam}"
    tcp_div = f"tcp.d6{'l' if beam == 1 else 'r'}7.b{beam}"

    part_conv = xc.generate_pencil_on_collimator(line,tcp_conv, num_particles=npart,
                                                nemitt_x=3.5e-6, nemitt_y=3.5e-6, longitudinal=longitudinal,
                                                longitudinal_betatron_cut=longitudinal_betatron_cut)
    part_div = xc.generate_pencil_on_collimator(line,tcp_div, num_particles=npart,
                                                nemitt_x=3.5e-6, nemitt_y=3.5e-6, longitudinal=longitudinal,
                                                longitudinal_betatron_cut=longitudinal_betatron_cut)
    with flaky_assertions():
        # Diverging beam --------------------------------------------------------------------------
        coll_div = line[tcp_div]
        drift = xt.Drift(length=coll_div.active_length)
        drift.track(part_div)
        mask_div_L = part_div.y > 0
        mask_div_R = part_div.y < 0


        # do something else: duct tape + magical numbers solution atm
        # if longitudinal == 'matched_dispersion': 
        #     assert math.isclose(np.mean(part_div.y),0.0,abs_tol=1e-5)     # unsure about this
        #     assert math.isclose(np.mean(part_div.py),0.0,abs_tol=1e-5)
        #     assert abs(np.std(part_div.x) - np.mean(part_div.x)) < (part_div.x.max() - np.mean(part_div.x))/2
        #     assert abs(np.std(part_div.px) - np.mean(part_div.px)) < (part_div.px.max() - np.mean(part_div.px))/2
        # else:
        #     assert math.isclose(np.mean(part_div.x),0.0,abs_tol=2e-6)
        #     assert math.isclose(np.mean(part_div.px),0.0,abs_tol=2e-6)
        #     assert abs(np.std(part_div.x) - np.mean(part_div.x)) < (part_div.x.max() - np.mean(part_div.x))/3 # nah?
        #     assert abs(np.std(part_div.px) - np.mean(part_div.px)) < (part_div.px.max() - np.mean(part_div.px))/3

        # left jaw 
        pos_jawL_div = coll_div.ref_y + coll_div.jaw_L 
        pos_partL_div = part_div.y[mask_div_L].min()
        pencil_spread_divL = part_div.y[mask_div_L].max() - part_div.y[mask_div_L].min()
   
        assert math.isclose(pencil_spread_divL, 1e-6, abs_tol=1e-7)
        assert abs(pos_jawL_div - pos_partL_div) < 1e-9

        # right jaw
        pos_jawR_div = coll_div.ref_y + coll_div.jaw_R 
        pos_partR_div = part_div.y[mask_div_R].max()
        pencil_spread_divR = part_div.y[mask_div_R].max() - part_div.y[mask_div_R].min()

        assert math.isclose(pencil_spread_divR, 1e-6, abs_tol=1e-7)
        assert abs(pos_jawR_div - pos_partR_div) < 1e-8

        # Converging beam --------------------------------------------------------------------------
        coll_conv = line[tcp_conv]
        mask_conv_L = part_conv.x > 0
        mask_conv_R = part_conv.x < 0

        if longitudinal == 'matched_dispersion':
             assert math.isclose(np.mean(part_conv.x),0.0,abs_tol=1e-5)
             assert math.isclose(np.mean(part_conv.px),0.0,abs_tol=1e-5)
        else:
             assert math.isclose(np.mean(part_conv.y),0.0,abs_tol=2e-6)
             assert math.isclose(np.mean(part_conv.py),0.0,abs_tol=2e-6)     
        assert abs(np.std(part_conv.y) - np.mean(part_conv.y)) < (part_conv.y.max() - np.mean(part_conv.y))/3 
        assert abs(np.std(part_conv.py) - np.mean(part_conv.py)) < (part_conv.py.max() - np.mean(part_conv.py))/3

        # left jaw
        pos_jawL_conv = coll_conv.ref_x + coll_conv.jaw_L
        pos_partL_conv = part_conv.x[mask_conv_L].min()
        pencil_spread_convL = part_conv.x[mask_conv_L].max() - part_conv.x[mask_conv_L].min()

        assert math.isclose(pencil_spread_convL, 1e-6, abs_tol=1e-7)
        assert abs(pos_jawL_conv - pos_partL_conv) < 1e-9   
        
        # right jaw
        pos_jawR_conv = coll_conv.ref_x + coll_conv.jaw_R 
        pos_partR_conv = part_conv.x[mask_conv_R].max()
        pencil_spread_convR = part_conv.x[mask_conv_R].max() - part_conv.x[mask_conv_R].min()
        
        assert math.isclose(pencil_spread_convR, 1e-6, abs_tol=1e-7)
        assert abs(pos_jawR_conv - pos_partR_conv) < 1e-8
        