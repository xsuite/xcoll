import json
import numpy as np
from pathlib import Path
import xtrack as xt
import xcoll as xc
import pytest
import math
from scipy import stats
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

    tw = line.twiss()
    tcp_conv = f"tcp.c6{'l' if beam == 1 else 'r'}7.b{beam}"
    tcp_div = f"tcp.d6{'l' if beam == 1 else 'r'}7.b{beam}"

    # Generate particles on a collimator
    part_conv = xc.generate_pencil_on_collimator(line,tcp_conv, num_particles=npart,
                                                nemitt_x=3.5e-6, nemitt_y=3.5e-6, longitudinal=longitudinal,
                                                longitudinal_betatron_cut=longitudinal_betatron_cut)
    part_div = xc.generate_pencil_on_collimator(line,tcp_div, num_particles=npart,
                                                nemitt_x=3.5e-6, nemitt_y=3.5e-6, longitudinal=longitudinal,
                                                longitudinal_betatron_cut=longitudinal_betatron_cut)

    # Normalize coordinates
    part_norm_conv = tw.get_normalized_coordinates(part_conv, nemitt_x=3.5e-6, nemitt_y=3.5e-6)
    part_norm_div = tw.get_normalized_coordinates(part_div, nemitt_x=3.5e-6, nemitt_y=3.5e-6)

    with flaky_assertions():
        # Diverging beam --------------------------------------------------------------------------
        coll_div = line[tcp_div]
        drift = xt.Drift(length=coll_div.length)
        drift.track(part_div)
        mask_div_L = part_div.y > 0
        mask_div_R = part_div.y < 0

        # mean
        assert math.isclose(np.mean(part_norm_div.x_norm),0.0,abs_tol=2e-4)
        assert math.isclose(np.mean(part_norm_div.px_norm),0.0,abs_tol=2e-4)

        # std 
        assert math.isclose(np.std(part_norm_div.x_norm),0.01,abs_tol=1e-4)
        assert math.isclose(np.std(part_norm_div.px_norm),0.01,abs_tol=1e-4)

        # delta w/ Chi-Square test or Kolmogorov-Smirnov test
        if longitudinal == 'matched_dispersion':
            count,_ = np.histogram(part_div.delta, bins=50)
            expected_counts = np.full_like(count, np.mean(count))
            _, p = stats.chisquare(count, expected_counts)
            assert p > 0.05
        elif longitudinal == 'bucket':
            count,_ = np.histogram(part_div.delta, bins=50)
            _, p = stats.kstest((count - np.mean(count))/np.std(count), 'norm')
            assert p > 0.05


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

        # mean
        assert math.isclose(np.mean(part_norm_conv.y_norm),0.0,abs_tol=2e-4)
        assert math.isclose(np.mean(part_norm_conv.y_norm),0.0,abs_tol=2e-4)

        # std 
        assert math.isclose(np.std(part_norm_conv.y_norm),0.01,abs_tol=1e-4)
        assert math.isclose(np.std(part_norm_conv.py_norm),0.01,abs_tol=1e-4)

        # delta w/ Chi-Square test or Kolmogorov-Smirnov test
        if longitudinal == 'matched_dispersion':
            count,_ = np.histogram(part_conv.delta, bins=50)
            expected_counts = np.full_like(count, np.mean(count))
            _, p = stats.chisquare(count, expected_counts)
            assert p > 0.05
        elif longitudinal == 'bucket':
            count,_ = np.histogram(part_conv.delta, bins=50)
            _, p = stats.kstest((count - np.mean(count))/np.std(count), 'norm')
            assert p > 0.05

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
        