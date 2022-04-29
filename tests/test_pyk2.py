import numpy as np
from pathlib import Path
import xpart as xp
import xtrack as xt
import xcoll as xc

from scipy.constants import m_p

from xcoll.beam_elements.pyk2 import k2_track
from xcoll.general import _pkg_root

path = Path(_pkg_root / ".." / "tests" / "pyk2_data")

k2_engine = xc.K2Engine(n_alloc=200000, random_generator_seed=7569)

x_test = np.loadtxt("pyk2_data/rcx.dump")
xp_test = np.loadtxt("pyk2_data/rcxp.dump")
y_test = np.loadtxt("pyk2_data/rcy.dump")
yp_test = np.loadtxt("pyk2_data/rcyp.dump")
s_test = np.loadtxt("pyk2_data/rcs.dump")
E_test = np.loadtxt("pyk2_data/rcp.dump")*1e9

e0 = 7e12
m0 = m_p
p0c = np.sqrt(e0**2 - m0**2)

one_plus_delta = np.sqrt(E_test**2 - m0**2) / p0c
px_test = xp_test * one_plus_delta
py_test = yp_test * one_plus_delta

# enom = 7000 # Reference energy
npart = len(x_test)
# part = xp.Particles(p0c=7e12, x=x_test, y=y_test)
# part.ptau = (p_test*1e9 - part.energy0)/(part.energy0*part.beta0)
# rpp_test = part.rpp.copy()
# part.px = xp_test/part.rpp
# part.py = yp_test/part.rpp
# part_test = part.copy()

particles = xp.Particles(
        x=x_test, px=px_test, y=y_test, py=py_test, s=s_test, delta=one_plus_delta-1, p0c=p0c
)

x_test, xp_test, y_test, yp_test, p_test, s_test, part_hit, part_abs = k2_track(
                            material='MoGR',
                            particles=particles,
                            closed_orbit=[0,0,0,0],
                            angle=0,
                            jaws=[0.0025711021962573095,0,0.0025711021962573095,0],
                            offset=-0.0025711021962573095/2,
                            npart=npart,
                            length=0.59999999999999998,
                            is_crystal=False,
                            onesided=False
                        )

x_ref = np.loadtxt("pyk2_data/rcx.dump_after_REF")
xp_ref = np.loadtxt("pyk2_data/rcxp.dump_after_REF")
y_ref = np.loadtxt("pyk2_data/rcy.dump_after_REF")
yp_ref = np.loadtxt("pyk2_data/rcyp.dump_after_REF")
s_ref = np.loadtxt("pyk2_data/rcs.dump_after_REF")
p_ref = np.loadtxt("pyk2_data/rcp.dump_after_REF")

print(part_abs)
print(x_test[part_abs <= 0])
print(x_ref[part_abs <= 0])

assert np.allclose(x_test, x_ref, atol=1e-9, rtol=0)
assert np.allclose(xp_test, xp_ref, atol=5e-9, rtol=0)
assert np.allclose(y_test, y_ref, atol=1e-9, rtol=0)
assert np.allclose(yp_test, yp_ref, atol=5e-9, rtol=0)
assert np.allclose(s_test, s_ref, atol=2e-4, rtol=0)
assert np.allclose(p_test, p_ref, atol=0, rtol=1e-7)

part_hit_pos_ref = np.loadtxt("pyk2_data/part_hit_pos.dump_after_REF")
part_abs_pos_ref = np.loadtxt("pyk2_data/part_abs_pos.dump_after_REF")
# part_impact_ref = np.loadtxt("pyk2_data/part_impact.dump_after_REF")
# part_indiv_ref = np.loadtxt("pyk2_data/part_indiv.dump_after_REF")
# part_linteract_ref = np.loadtxt("pyk2_data/part_linteract.dump_after_REF")
# nhit_stage_ref = np.loadtxt("pyk2_data/nhit_stage.dump_after_REF")
# nabs_type_ref = np.loadtxt("pyk2_data/nabs_type.dump_after_REF")

assert np.allclose(part_hit, part_hit_pos_ref, atol=1e-9, rtol=0)
assert np.allclose(part_abs, part_abs_pos_ref, atol=1e-9, rtol=0)
# assert np.allclose(part_impact, part_impact_ref, atol=1e-9, rtol=0)
# assert np.allclose(part_indiv, part_indiv_ref, atol=1e-9, rtol=0)
# assert np.allclose(part_linteract, part_linteract_ref, atol=2e-4, rtol=0)
# assert np.allclose(nhit_stage, nhit_stage_ref, atol=1e-9, rtol=0)
# assert np.allclose(nabs_type, nabs_type, atol=1e-9, rtol=0)
