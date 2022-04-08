import numpy as np
from pathlib import Path
import xpart as xp
import xtrack as xt
import xcoll as xc

from xcoll.beam_elements.pyk2 import pyk2_run, materials
from xcoll.general import _pkg_root

path = Path(_pkg_root / ".." / "tests" / "pyk2_data")

k2_engine = xc.K2Engine(n_alloc=200000, random_generator_seed=7569)

x_test = np.loadtxt("pyk2_data/rcx.dump")
xp_test = np.loadtxt("pyk2_data/rcxp.dump")
y_test = np.loadtxt("pyk2_data/rcy.dump")
yp_test = np.loadtxt("pyk2_data/rcyp.dump")
s_test = np.loadtxt("pyk2_data/rcs.dump")
p_test = np.loadtxt("pyk2_data/rcp.dump")

npart = len(x_test)
# part = xp.Particles(p0c=7e12, x=x_test, y=y_test)
# part.ptau = (p_test*1e9 - part.energy0)/(part.energy0*part.beta0)
# rpp_test = part.rpp.copy()
# part.px = xp_test/part.rpp
# part.py = yp_test/part.rpp
# part_test = part.copy()

part_hit = np.zeros(npart, dtype=np.int32)
part_abs = np.zeros(npart, dtype=np.int32)
part_impact = np.zeros(npart, dtype=float)
part_indiv = np.zeros(npart, dtype=float)
part_linteract = np.zeros(npart, dtype=float)
nhit_stage = np.zeros(npart, dtype=np.int32)
nabs_type = np.zeros(npart, dtype=np.int32)
linside = np.zeros(npart, dtype=np.int32)

matID = materials['MoGR']['ID']
exenergy = materials['MoGR']['exenergy']
anuc = materials['MoGR']['anuc']
zatom = materials['MoGR']['zatom']
rho = materials['MoGR']['rho']
hcut = materials['MoGR']['hcut']
bnref = materials['MoGR']['bnref']

pyk2_run(x_particles=x_test,
          xp_particles=xp_test,
          y_particles=y_test,
          yp_particles=yp_test,
          s_particles=s_test,
          p_particles=p_test,              # confusing: this is ENERGY not momentum
          part_hit=part_hit,
          part_abs=part_abs,
          part_impact=part_impact,         # impact parameter
          part_indiv=part_indiv,           # particle divergence
          part_linteract=part_linteract,   # interaction length
          nhit_stage=nhit_stage,
          nabs_type=nabs_type,
          linside=linside,
          matid=matID,
          run_exenergy=exenergy,
          run_anuc=anuc,
          run_zatom=zatom,
          run_rho=rho,
          run_hcut=hcut,
          run_bnref=bnref,
          is_crystal=False,
          c_length=0.59999999999999998,
          c_rotation=0,
          c_aperture=0.0025711021962573095,
          c_offset=0,
          c_tilt=np.array([0,0], dtype=np.float64),
          c_enom=7000000, # Reference energy
          onesided=False,
          random_generator_seed=-1 # skips rng re-initlization
          )

x_ref = np.loadtxt("pyk2_data/rcx.dump_after_REF")
xp_ref = np.loadtxt("pyk2_data/rcxp.dump_after_REF")
y_ref = np.loadtxt("pyk2_data/rcy.dump_after_REF")
yp_ref = np.loadtxt("pyk2_data/rcyp.dump_after_REF")
s_ref = np.loadtxt("pyk2_data/rcs.dump_after_REF")
p_ref = np.loadtxt("pyk2_data/rcp.dump_after_REF")

assert np.allclose(x_test, x_ref, atol=1e-9, rtol=0)
assert np.allclose(xp_test, xp_ref, atol=5e-9, rtol=0)
assert np.allclose(y_test, y_ref, atol=1e-9, rtol=0)
assert np.allclose(yp_test, yp_ref, atol=5e-9, rtol=0)
assert np.allclose(s_test, s_ref, atol=2e-4, rtol=0)
assert np.allclose(p_test, p_ref, atol=0, rtol=1e-7)

part_hit_pos_ref = np.loadtxt("pyk2_data/part_hit_pos.dump_after_REF")
part_abs_pos_ref = np.loadtxt("pyk2_data/part_abs_pos.dump_after_REF")
part_impact_ref = np.loadtxt("pyk2_data/part_impact.dump_after_REF")
part_indiv_ref = np.loadtxt("pyk2_data/part_indiv.dump_after_REF")
part_linteract_ref = np.loadtxt("pyk2_data/part_linteract.dump_after_REF")
nhit_stage_ref = np.loadtxt("pyk2_data/nhit_stage.dump_after_REF")
nabs_type_ref = np.loadtxt("pyk2_data/nabs_type.dump_after_REF")

assert np.allclose(part_hit, part_hit_pos_ref, atol=1e-9, rtol=0)
assert np.allclose(part_abs, part_abs_pos_ref, atol=1e-9, rtol=0)
assert np.allclose(part_impact, part_impact_ref, atol=1e-9, rtol=0)
assert np.allclose(part_indiv, part_indiv_ref, atol=1e-9, rtol=0)
assert np.allclose(part_linteract, part_linteract_ref, atol=2e-4, rtol=0)
assert np.allclose(nhit_stage, nhit_stage_ref, atol=1e-9, rtol=0)
assert np.allclose(nabs_type, nabs_type, atol=1e-9, rtol=0)
