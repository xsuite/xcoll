#!/usr/bin/env python3

import numpy as np
import xtrack as xt
import xcoll as xc
import time
import os

# ============================================================
# Load Mathematica reference (Case 0)
# ============================================================

ref_file = "bent_channelling_reference_functions/solution_case00.csv"
if not os.path.exists(ref_file):
    raise FileNotFoundError(ref_file)

ref = np.loadtxt(ref_file, delimiter=",", skiprows=1)
s_ref  = ref[:, 0]
x_ref  = ref[:, 1]
px_ref = ref[:, 2]

x0 = x_ref[0]
px0 = px_ref[0]
xF_ref = x_ref[-1]
pxF_ref = px_ref[-1]
L = s_ref[-1]

print("\n================ MATHEMATICA REFERENCE (CASE 0) =================")
print(f"sF = {L:.6e}")
print(f"x0 = {x0:.6e}, px0 = {px0:.6e}")
print(f"xF = {xF_ref:.6e}, pxF = {pxF_ref:.6e}")

def rel_err(v, vref):
    return abs(v - vref) / max(abs(vref), 1e-30)

# ============================================================
# CRYSTAL A — DEFAULTS
# ============================================================

print("\n================ CRYSTAL A — DEFAULTS for M4 =================")


crystal_a = xc.BentChannellingDev(length=L)

print("Crystal parameters:")
print(f"  length   = {crystal_a.length:.6e}")
print(f"  R        = {crystal_a.R}")
print(f"  U0       = {crystal_a.U0}")
print(f"  Umax     = {crystal_a.Umax}")
print(f"  dp       = {crystal_a.dp:.3e}")
print(f"  aTF      = {crystal_a.aTF:.3e}")
print(f"  uT       = {crystal_a.uT:.3e}")
print(f"  alpha_i  = {crystal_a.alpha_i}")
print(f"  beta_i   = {crystal_a.beta_i}")
print(f"  method   = {crystal_a.method}")
print(f"  order    = {crystal_a.order}")
print(f"  variant  = {crystal_a.variant}")
print(f"  n_steps  = {crystal_a.n_steps}")

particle_a = xt.Particles(x=x0, px=px0, p0c=150e9, delta=0.0)
bpcA = (
    particle_a.beta0[0]
    * particle_a.rvv[0]
    * (1.0 + particle_a.delta[0])
    * particle_a.p0c[0]
    * particle_a.charge_ratio[0]
)

print("\nParticle parameters:")
print(f"  beta0        = {particle_a.beta0[0]}")
print(f"  rvv          = {particle_a.rvv[0]}")
print(f"  delta        = {particle_a.delta[0]}")
print(f"  p0c          = {particle_a.p0c[0]:.6e}")
print(f"  charge_ratio = {particle_a.charge_ratio[0]}")
print(f"  --> bpc      = {bpcA:.6e}")

print("\n--- BEFORE ---")
print(f"x = {particle_a.x[0]:.6e}, px = {particle_a.px[0]:.6e}, s = {particle_a.s[0]:.6e}")

crystal_a.track(particle_a.copy())
t0A = time.perf_counter()
crystal_a.track(particle_a)
tA = time.perf_counter() - t0A

print("--- AFTER ---")
print(f"x = {particle_a.x[0]:.6e}, px = {particle_a.px[0]:.6e}, s = {particle_a.s[0]:.6e}")
print(f"Runtime = {tA*1e3:.6f} ms")

# ============================================================
# CRYSTAL B — DEFAULTS FOR M2
# ============================================================
print("\n================ CRYSTAL B — DEFAULTS FOR M2 =================")


crystal_b = xc.BentChannellingDev(length=L,
method = 2)

print("Crystal parameters:")
print(f"  length   = {crystal_b.length:.6e}")
print(f"  R        = {crystal_b.R}")
print(f"  U0       = {crystal_b.U0}")
print(f"  Umax     = {crystal_b.Umax}")
print(f"  dp       = {crystal_b.dp:.3e}")
print(f"  aTF      = {crystal_b.aTF:.3e}")
print(f"  uT       = {crystal_b.uT:.3e}")
print(f"  alpha_i  = {crystal_b.alpha_i}")
print(f"  beta_i   = {crystal_b.beta_i}")
print(f"  method   = {crystal_b.method}")
print(f"  order    = {crystal_b.order}")
print(f"  variant  = {crystal_b.variant}")
print(f"  n_steps  = {crystal_b.n_steps}")

particle_b = xt.Particles(x=x0, px=px0, p0c=150e9, delta=0.0)

bpcB = (
    particle_b.beta0[0]
    * particle_b.rvv[0]
    * (1.0 + particle_b.delta[0])
    * particle_b.p0c[0]
    * particle_b.charge_ratio[0]
)

print("\nParticle parameters:")
print(f"  beta0        = {particle_b.beta0[0]}")
print(f"  rvv          = {particle_b.rvv[0]}")
print(f"  delta        = {particle_b.delta[0]}")
print(f"  p0c          = {particle_b.p0c[0]:.6e}")
print(f"  charge_ratio = {particle_b.charge_ratio[0]}")
print(f"  --> bpc      = {bpcB:.6e}")

print("\n--- BEFORE ---")
print(f"x = {particle_b.x[0]:.6e}, px = {particle_b.px[0]:.6e}, s = {particle_b.s[0]:.6e}")

crystal_b.track(particle_b.copy())
t0B = time.perf_counter()
crystal_b.track(particle_b)
tB = time.perf_counter() - t0B

print("--- AFTER ---")
print(f"x = {particle_b.x[0]:.6e}, px = {particle_b.px[0]:.6e}, s = {particle_b.s[0]:.6e}")
print(f"Runtime = {tB*1e3:.6f} ms")

# ============================================================
# CRYSTAL C — p0c = 150 GeV, delta = 0
# ============================================================
print("\n================ CRYSTAL C — DEFAULTS FOR M3 =================")


crystal_c = xc.BentChannellingDev(length=L,
method = 3)

print("Crystal parameters:")
print(f"  length   = {crystal_c.length:.6e}")
print(f"  R        = {crystal_c.R}")
print(f"  U0       = {crystal_c.U0}")
print(f"  Umax     = {crystal_c.Umax}")
print(f"  dp       = {crystal_c.dp:.3e}")
print(f"  aTF      = {crystal_c.aTF:.3e}")
print(f"  uT       = {crystal_c.uT:.3e}")
print(f"  alpha_i  = {crystal_c.alpha_i}")
print(f"  beta_i   = {crystal_c.beta_i}")
print(f"  method   = {crystal_c.method}")
print(f"  order    = {crystal_c.order}")
print(f"  variant  = {crystal_c.variant}")
print(f"  n_steps  = {crystal_c.n_steps}")

particle_c = xt.Particles(x=x0, px=px0, p0c=150e9, delta=0.0)

bpcC = (
    particle_c.beta0[0]
    * particle_c.rvv[0]
    * (1.0 + particle_c.delta[0])
    * particle_c.p0c[0]
    * particle_c.charge_ratio[0]
)

print("\nParticle parameters:")
print(f"  beta0        = {particle_c.beta0[0]}")
print(f"  rvv          = {particle_c.rvv[0]}")
print(f"  delta        = {particle_c.delta[0]}")
print(f"  p0c          = {particle_c.p0c[0]:.6e}")
print(f"  charge_ratio = {particle_c.charge_ratio[0]}")
print(f"  --> bpc      = {bpcC:.6e}")

print("\n--- BEFORE ---")
print(f"x = {particle_c.x[0]:.6e}, px = {particle_c.px[0]:.6e}, s = {particle_c.s[0]:.6e}")

crystal_c.track(particle_c.copy())
t0C = time.perf_counter()
crystal_c.track(particle_c)
tC = time.perf_counter() - t0C

print("--- AFTER ---")
print(f"x = {particle_c.x[0]:.6e}, px = {particle_c.px[0]:.6e}, s = {particle_c.s[0]:.6e}")
print(f"Runtime = {tC*1e3:.6f} ms")

# ============================================================
# FINAL ERRORS
# ============================================================

print("\n================ FINAL  =================")
print(f"CRYSTAL A: rel_err_x={rel_err(particle_a.x[0], xF_ref):.5f}, rel_err_px={rel_err(particle_a.px[0], pxF_ref):.5f}")
print(f"CRYSTAL B: rel_err_x={rel_err(particle_b.x[0], xF_ref):.5f}, rel_err_px={rel_err(particle_b.px[0], pxF_ref):.5f}")
print(f"CRYSTAL C: rel_err_x={rel_err(particle_c.x[0], xF_ref):.5f}, rel_err_px={rel_err(particle_c.px[0], pxF_ref):.5f}")
