import numpy as np
import sys, os, contextlib
import time
import math

from pathlib import Path

import xobjects as xo
import xtrack as xt
import xcoll as xc
import xpart as xp

import argparse
import pandas as pd


path_in = Path("/eos/project-c/collimation-team/machine_configurations/LHC_run3/2025/xsuite")
path_out = Path.cwd()

beam          = 1

import matplotlib.pyplot as plt

def FlukaBeamBeamSource(line, colldb, int_type, ip, num_particles, pdg_id_b1, p0c_b1, Z_b2, A_b2,
                        sigma_p_x2 = 0, sigma_p_y2 = 0, sigma_z = 0, sigma_dpp = 0):

    _capacity = num_particles * 100
    tw = line.twiss()
    if ip == "ip1":
        xsx = line["on_x1"]
        xsv = 0
        phi_ir = line['phi_ir1']
    elif ip == "ip2":
        xsx = line["on_x2h"]
        phi_ir = line['phi_ir2']
        xsv = line["on_x2v"]
    elif ip == "ip5":
        xsx = line["on_x5"]
        xsv = 0
        phi_ir = line['phi_ir5']
    elif ip == "ip8":
        xsx = line["on_x8h"]
        xsv = line["on_x8v"]
        phi_ir = line['phi_ir8']

    theta2 = (math.sqrt(xsx**2 + xsv**2) * 1e-6)

    # if theta2 is negative, make it positive and change sing to phi_ir
    if xsx < 0 or xsv < 0:
        theta2 = -theta2
        phi_ir = -phi_ir

    geo_emittance_x = tw.beta0*tw.gamma0*colldb.nemitt_x
    geo_emittance_y = tw.beta0*tw.gamma0*colldb.nemitt_y

    sigma_p_x2 = (geo_emittance_x*(tw.rows[ip].gamx)[0])**0.5
    sigma_p_y2 = (geo_emittance_y*(tw.rows[ip].gamy)[0])**0.5
    


    bb_int = {
            'int_type': int_type,
            'theta2': theta2, # Angle between b2 and -Z
            'xs': phi_ir,
            'sigma_p_x2': sigma_p_x2,
            'sigma_p_y2': sigma_p_y2,
            'Z': Z_b2,
            'A': A_b2,
            'betx': (tw.rows[ip].betx)[0],
            'bety': (tw.rows[ip].bety)[0],
            'dx': (tw.rows[ip].dx)[0],
            'dy': (tw.rows[ip].dy)[0],
            'alfx': (tw.rows[ip].alfx)[0],
            'alfy': (tw.rows[ip].alfy)[0],
            'dpx': (tw.rows[ip].dpx)[0],
            'dpy': (tw.rows[ip].dpy)[0],
            'rms_emx': colldb.nemitt_x, #/(tw.beta0* tw.gamma0),
            'rms_emy': colldb.nemitt_y, #/(tw.beta0* tw.gamma0)
            'offset': [tw.rows[ip].x[0], tw.rows[ip].px[0], tw.rows[ip].y[0], tw.rows[ip].py[0]],
            'sigma_z': sigma_z * 1e2,
            'sigma_dpp': sigma_dpp
            }

    coll_dummy = xc.FlukaCollimator(assembly="lhc_ippipe", length=0.6, jaw=0.001)

    xc.fluka.engine.particle_ref = xt.Particles.reference_from_pdg_id(pdg_id=pdg_id_b1, p0c=p0c_b1)
    xc.fluka.engine.capacity = _capacity
    xc.fluka.engine.relative_capacity = 20 #100
    # xc.fluka.engine.seed = 5656565
    xc.fluka.engine.start(elements=coll_dummy, clean=False , verbose=True, include_showers=False, return_ions=True, bb_int=bb_int, touches=False)


    line.build_tracker()

    nemitt_x = colldb.nemitt_x
    nemitt_y = colldb.nemitt_y
    bunch_intensity = 1.8E11
    particles = xp.generate_matched_gaussian_bunch(
        num_particles=num_particles, total_intensity_particles=bunch_intensity,
        nemitt_x=nemitt_x, nemitt_y=nemitt_y, sigma_z=sigma_z,
        line=line, particle_ref=xc.fluka.engine.particle_ref)
    part_init = line.build_particles(x=particles.x, y=particles.y,
                            px=particles.px,
                            py=particles.py,
                            zeta=particles.zeta,
                            delta=particles.delta,
                            nemitt_x=colldb.nemitt_x,
                            nemitt_y=colldb.nemitt_y,
                            at_element=ip,
                            match_at_s=line.get_s_position(ip),
                            particle_ref=xc.fluka.engine.particle_ref,
                            _capacity=xc.fluka.engine.capacity,
    mode="normalized_transverse")
    part = part_init.copy()

    print(f"Tracking {num_particles} particles (FLUKA)...     ", end='')
    start = time.time()
    coll_dummy.track(part)
    print(f"Done in {round(time.time()-start, 3)}s.")
    xc.fluka.engine.stop(clean=False)

    line.discard_tracker()

    return part


start_time = time.time()

# Path to inputs
path_out = Path('./')
path_out.mkdir(exist_ok=True)

num_particles = int(100) #0

line = xt.load(path_in / f'levelling.23_b{beam}.json')

colldb = xc.CollimatorDatabase.from_yaml(path_in / ".." / "colldbs" / f'levelling.23.yaml', beam=beam)

part = FlukaBeamBeamSource(
    line=line,
    colldb=colldb,
    int_type=1.0, # 1.0 inelastic, 10.0 elastic, 100.0 Emd
    ip="ip1",
    num_particles=num_particles,
    pdg_id_b1='proton',
    p0c_b1=6.8e12,
    Z_b2=1,
    A_b2=1,
    sigma_z = 8e-2,
    sigma_dpp = 1.1e-4
)
import pdb; pdb.set_trace()

import pickle
part_pkl_init = path_out / f'part_init.pkl'
part_dpp_cut_init = part.filter(part.pdg_id != -999999999)
with open(part_pkl_init, 'wb') as fid_out:
    pickle.dump(part_dpp_cut_init.to_dict(), fid_out)

