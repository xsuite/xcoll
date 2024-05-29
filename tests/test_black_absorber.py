# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

import numpy as np
import pytest

import xobjects as xo
import xpart as xp
import xtrack as xt
import xcoll as xc

from xpart.test_helpers import flaky_assertions, retry
from xobjects.test_helpers import for_all_test_contexts


n_part = int(2.e6)

@for_all_test_contexts(
    excluding=('ContextCupy', 'ContextPyopencl')  # BlackAbsorber not on GPU
)
def test_horizontal_parallel(test_context):
    jaw_L, jaw_R, _, _, L, coll = _make_absorber(_context=test_context)
    part, x, _, _, _ = _generate_particles(_context=test_context)
    coll.track(part)
    part.move(_context=xo.ContextCpu())
    part.sort(interleave_lost_particles=True)
    # As the angles are zero, only particles that started in front of the jaw are lost
    mask_hit = (x >= jaw_L) | (x <= jaw_R)

    lost = np.unique(part.state[mask_hit])
    alive = np.unique(part.state[~mask_hit])
    assert len(lost)==1 and lost[0]==-340
    assert len(alive)==1 and alive[0]==1
    s_lost = np.unique(part.s[mask_hit])
    s_alive = np.unique(part.s[~mask_hit])
    assert len(s_lost)==1 and s_lost[0]==0
    assert len(s_alive)==1 and s_alive[0]==L
    assert np.allclose(part.x, x, atol=1e-12, rtol=0)


@for_all_test_contexts(
    excluding=('ContextCupy', 'ContextPyopencl')  # BlackAbsorber not on GPU
)
def test_horizontal(test_context):
    jaw_L, jaw_R, _, _, L, coll = _make_absorber(_context=test_context)
    part, x, _, xp, _ = _generate_particles(four_dim=True, _context=test_context)
    coll.track(part)
    part.move(_context=xo.ContextCpu())
    part.sort(interleave_lost_particles=True)
    # Particles in front of the jaw are lost, ...
    mask_hit_front = (x >= jaw_L) | (x <= jaw_R)
    # but also those in the opening with an angle that would kick them on the jaw halfway
    mask_hit_angle_L = ~mask_hit_front & (x + xp*L >= jaw_L)
    mask_hit_angle_R = ~mask_hit_front & (x + xp*L <= jaw_R)
    mask_hit = mask_hit_front | mask_hit_angle_L | mask_hit_angle_R

    lost = np.unique(part.state[mask_hit])
    alive = np.unique(part.state[~mask_hit])
    assert len(lost)==1 and lost[0]==-340
    assert len(alive)==1 and alive[0]==1
    s_lost_front = np.unique(part.s[mask_hit_front])
    s_alive = np.unique(part.s[~mask_hit])
    assert len(s_lost_front)==1 and s_lost_front[0]==0
    assert len(s_alive)==1 and s_alive[0]==L
    assert np.allclose(part.x[mask_hit_front], x[mask_hit_front], atol=1e-12, rtol=0)
    assert np.allclose(part.x[~mask_hit],      x[~mask_hit] + xp[~mask_hit]*L, atol=1e-12, rtol=0)
    s_hit_L = (jaw_L - x[mask_hit_angle_L]) / xp[mask_hit_angle_L]
    assert np.allclose(part.s[mask_hit_angle_L], s_hit_L, atol=1e-12, rtol=0)
    assert np.allclose(part.x[mask_hit_angle_L], jaw_L, atol=1e-12, rtol=0)
    s_hit_R = (jaw_R - x[mask_hit_angle_R]) / xp[mask_hit_angle_R]
    assert np.allclose(part.s[mask_hit_angle_R], s_hit_R, atol=1e-12, rtol=0)
    assert np.allclose(part.x[mask_hit_angle_R], jaw_R, atol=1e-12, rtol=0)


@for_all_test_contexts(
    excluding=('ContextCupy', 'ContextPyopencl')  # BlackAbsorber not on GPU
)
@retry()
def test_horizontal_with_tilts(test_context):
    jaw_LU, jaw_RU, jaw_LD, jaw_RD, L, coll = _make_absorber(tilts=[0.0005, -0.00015], _context=test_context)
    part, x, _, xp, _ = _generate_particles(four_dim=True, _context=test_context)
    coll.track(part)
    part.move(_context=xo.ContextCpu())
    part.sort(interleave_lost_particles=True)
    # Particles in front of the jaw are lost, ...
    mask_hit_front = (x >= jaw_LU) | (x <= jaw_RU)
    # but also those in the opening with an angle that would kick them on the jaw halfway
    mask_hit_angle_L = ~mask_hit_front & (x + xp*L >= jaw_LD)
    mask_hit_angle_R = ~mask_hit_front & (x + xp*L <= jaw_RD)
    mask_hit = mask_hit_front | mask_hit_angle_L | mask_hit_angle_R

    lost = np.unique(part.state[mask_hit])
    alive = np.unique(part.state[~mask_hit])
    assert len(lost)==1 and lost[0]==-340
    assert len(alive)==1 and alive[0]==1
    cL = np.cos(coll.tilt_L)
    cR = np.cos(coll.tilt_R)
    s_lost_front_L = part.s[mask_hit_front & (x > (jaw_LU + jaw_RU)/2)]
    assert np.all([s <= coll.length/2*(1 - cL) for s in s_lost_front_L])
    s_lost_front_R = part.s[mask_hit_front & (x < (jaw_LU + jaw_RU)/2)]
    assert np.all([s <= coll.length/2*(1 - cR) for s in s_lost_front_R])
    s_alive = np.unique(part.s[~mask_hit])
    assert len(s_alive)==1 and s_alive[0]==L
    with flaky_assertions():
        assert np.allclose(part.x[mask_hit_front] - part.px[mask_hit_front]*part.s[mask_hit_front], x[mask_hit_front], atol=1e-12, rtol=0)
        assert np.allclose(part.x[~mask_hit], x[~mask_hit] + xp[~mask_hit]*L, atol=1e-12, rtol=0)
        s_hit_L = (jaw_LU - x[mask_hit_angle_L] - (jaw_LD-jaw_LU)/2*(1-cL)/cL) / (xp[mask_hit_angle_L] - (jaw_LD-jaw_LU)/(L*cL))
        assert np.allclose(part.s[mask_hit_angle_L], s_hit_L, atol=1e-12, rtol=0)
        assert np.allclose(part.x[mask_hit_angle_L], jaw_LU + (jaw_LD-jaw_LU)/(L*cL)*(s_hit_L - L/2*(1-cL)), atol=1e-12, rtol=0)
        s_hit_R = (jaw_RU - x[mask_hit_angle_R] - (jaw_RD-jaw_RU)/2*(1-cR)/cR) / (xp[mask_hit_angle_R] - (jaw_RD-jaw_RU)/(L*cR))
        assert np.allclose(part.s[mask_hit_angle_R], s_hit_R, atol=1e-12, rtol=0)
        assert np.allclose(part.x[mask_hit_angle_R], jaw_RU + (jaw_RD-jaw_RU)/(L*cR)*(s_hit_R - L/2*(1-cR)), atol=1e-12, rtol=0)


@for_all_test_contexts(
    excluding=('ContextCupy', 'ContextPyopencl')  # BlackAbsorber not on GPU
)
def test_vertical_parallel(test_context):
    jaw_L, jaw_R, _, _, L, coll = _make_absorber(angle=90, _context=test_context)
    part, _, y, _, _ = _generate_particles(_context=test_context)
    coll.track(part)
    part.move(_context=xo.ContextCpu())
    part.sort(interleave_lost_particles=True)
    # As the angles are zero, only particles that started in front of the jaw are lost
    mask_hit = (y >= jaw_L) | (y <= jaw_R)

    lost = np.unique(part.state[mask_hit])
    alive = np.unique(part.state[~mask_hit])
    assert len(lost)==1 and lost[0]==-340
    assert len(alive)==1 and alive[0]==1
    s_lost = np.unique(part.s[mask_hit])
    s_alive = np.unique(part.s[~mask_hit])
    assert len(s_lost)==1 and s_lost[0]==0
    assert len(s_alive)==1 and s_alive[0]==L
    assert np.allclose(part.y, y, atol=1e-12, rtol=0)


@for_all_test_contexts(
    excluding=('ContextCupy', 'ContextPyopencl')  # BlackAbsorber not on GPU
)
def test_vertical(test_context):
    jaw_L, jaw_R, _, _, L, coll = _make_absorber(angle=90, _context=test_context)
    part, _, y, _, yp = _generate_particles(four_dim=True, _context=test_context)
    coll.track(part)
    part.move(_context=xo.ContextCpu())
    part.sort(interleave_lost_particles=True)
    # Particles in front of the jaw are lost, ...
    mask_hit_front = (y >= jaw_L) | (y <= jaw_R)
    # but also those in the opening with an angle that would kick them on the jaw halfway
    mask_hit_angle_L = ~mask_hit_front & (y + yp*L >= jaw_L)
    mask_hit_angle_R = ~mask_hit_front & (y + yp*L <= jaw_R)
    mask_hit = mask_hit_front | mask_hit_angle_L | mask_hit_angle_R

    lost = np.unique(part.state[mask_hit])
    alive = np.unique(part.state[~mask_hit])
    assert len(lost)==1 and lost[0]==-340
    assert len(alive)==1 and alive[0]==1
    s_lost_front = np.unique(part.s[mask_hit_front])
    s_alive = np.unique(part.s[~mask_hit])
    assert len(s_lost_front)==1 and s_lost_front[0]==0
    assert len(s_alive)==1 and s_alive[0]==L
    assert np.allclose(part.y[mask_hit_front], y[mask_hit_front], atol=1e-12, rtol=0)
    assert np.allclose(part.y[~mask_hit],      y[~mask_hit] + yp[~mask_hit]*L, atol=1e-12, rtol=0)
    s_hit_L = (jaw_L - y[mask_hit_angle_L]) / yp[mask_hit_angle_L]
    assert np.allclose(part.s[mask_hit_angle_L], s_hit_L, atol=1e-12, rtol=0)
    assert np.allclose(part.y[mask_hit_angle_L], jaw_L, atol=1e-12, rtol=0)
    s_hit_R = (jaw_R - y[mask_hit_angle_R]) / yp[mask_hit_angle_R]
    assert np.allclose(part.s[mask_hit_angle_R], s_hit_R, atol=1e-12, rtol=0)
    assert np.allclose(part.y[mask_hit_angle_R], jaw_R, atol=1e-12, rtol=0)


@for_all_test_contexts(
    excluding=('ContextCupy', 'ContextPyopencl')  # BlackAbsorber not on GPU
)
@retry()
def test_vertical_with_tilts(test_context):
    jaw_LU, jaw_RU, jaw_LD, jaw_RD, L, coll = _make_absorber(angle=90, tilts=[0.0005, -0.00015], _context=test_context)
    part, _, y, _, yp = _generate_particles(four_dim=True, _context=test_context)
    coll.track(part)
    part.move(_context=xo.ContextCpu())
    part.sort(interleave_lost_particles=True)
    # Particles in front of the jaw are lost, ...
    mask_hit_front = (y >= jaw_LU) | (y <= jaw_RU)
    # but also those in the opening with an angle that would kick them on the jaw halfway
    mask_hit_angle_L = ~mask_hit_front & (y + yp*L >= jaw_LD)
    mask_hit_angle_R = ~mask_hit_front & (y + yp*L <= jaw_RD)
    mask_hit = mask_hit_front | mask_hit_angle_L | mask_hit_angle_R

    lost = np.unique(part.state[mask_hit])
    alive = np.unique(part.state[~mask_hit])
    assert len(lost)==1 and lost[0]==-340
    assert len(alive)==1 and alive[0]==1
    cL = np.cos(coll.tilt_L)
    cR = np.cos(coll.tilt_R)
    s_lost_front_L = part.s[mask_hit_front & (y > (jaw_LU + jaw_RU)/2)]
    assert np.all([s <= coll.length/2*(1 - cL) for s in s_lost_front_L])
    s_lost_front_R = part.s[mask_hit_front & (y < (jaw_LU + jaw_RU)/2)]
    assert np.all([s <= coll.length/2*(1 - cR) for s in s_lost_front_R])
    s_alive = np.unique(part.s[~mask_hit])
    assert len(s_alive)==1 and s_alive[0]==L
    with flaky_assertions():
        assert np.allclose(part.y[mask_hit_front] - part.py[mask_hit_front]*part.s[mask_hit_front], y[mask_hit_front], atol=1e-12, rtol=0)
        assert np.allclose(part.y[~mask_hit], y[~mask_hit] + yp[~mask_hit]*L, atol=1e-12, rtol=0)
        s_hit_L = (jaw_LU - y[mask_hit_angle_L] - (jaw_LD-jaw_LU)/2*(1-cL)/cL) / (yp[mask_hit_angle_L] - (jaw_LD-jaw_LU)/(L*cL))
        assert np.allclose(part.s[mask_hit_angle_L], s_hit_L, atol=1e-12, rtol=0)
        assert np.allclose(part.y[mask_hit_angle_L], jaw_LU + (jaw_LD-jaw_LU)/(L*cL)*(s_hit_L - L/2*(1-cL)), atol=1e-12, rtol=0)
        s_hit_R = (jaw_RU - y[mask_hit_angle_R] - (jaw_RD-jaw_RU)/2*(1-cR)/cR) / (yp[mask_hit_angle_R] - (jaw_RD-jaw_RU)/(L*cR))
        assert np.allclose(part.s[mask_hit_angle_R], s_hit_R, atol=1e-12, rtol=0)
        assert np.allclose(part.y[mask_hit_angle_R], jaw_RU + (jaw_RD-jaw_RU)/(L*cR)*(s_hit_R - L/2*(1-cR)), atol=1e-12, rtol=0)


@for_all_test_contexts(
    excluding=('ContextCupy', 'ContextPyopencl')  # BlackAbsorber not on GPU
)
def test_angled_parallel(test_context):
    angle=17.8
    jaw_L, jaw_R, _, _, L, coll = _make_absorber(angle=angle, _context=test_context)
    part, x, _, _, _ = _generate_particles(angle=angle, _context=test_context)
    coll.track(part)
    part.move(_context=xo.ContextCpu())
    part.sort(interleave_lost_particles=True)
    # As the angles are zero, only particles that started in front of the jaw are lost
    mask_hit = (x >= jaw_L) | (x <= jaw_R)

    lost = np.unique(part.state[mask_hit])
    alive = np.unique(part.state[~mask_hit])
    assert len(lost)==1 and lost[0]==-340
    assert len(alive)==1 and alive[0]==1
    s_lost = np.unique(part.s[mask_hit])
    s_alive = np.unique(part.s[~mask_hit])
    assert len(s_lost)==1 and s_lost[0]==0
    assert len(s_alive)==1 and s_alive[0]==L
    part_x_rot = part.x * np.cos(angle/180.*np.pi) + part.y * np.sin(angle/180.*np.pi)
    assert np.allclose(part_x_rot, x, atol=1e-12, rtol=0)


@for_all_test_contexts(
    excluding=('ContextCupy', 'ContextPyopencl')  # BlackAbsorber not on GPU
)
def test_angled(test_context):
    angle=17.8
    jaw_L, jaw_R, _, _, L, coll = _make_absorber(angle=angle, _context=test_context)
    part, x, _, xp, _ = _generate_particles(four_dim=True, angle=angle, _context=test_context)
    coll.track(part)
    part.move(_context=xo.ContextCpu())
    part.sort(interleave_lost_particles=True)
    # Particles in front of the jaw are lost, ...
    mask_hit_front = (x >= jaw_L) | (x <= jaw_R)
    # but also those in the opening with an angle that would kick them on the jaw halfway
    mask_hit_angle_L = ~mask_hit_front & (x + xp*L >= jaw_L)
    mask_hit_angle_R = ~mask_hit_front & (x + xp*L <= jaw_R)
    mask_hit = mask_hit_front | mask_hit_angle_L | mask_hit_angle_R

    lost = np.unique(part.state[mask_hit])
    alive = np.unique(part.state[~mask_hit])
    assert len(lost)==1 and lost[0]==-340
    assert len(alive)==1 and alive[0]==1
    s_lost_front = np.unique(part.s[mask_hit_front])
    s_alive = np.unique(part.s[~mask_hit])
    assert len(s_lost_front)==1 and s_lost_front[0]==0
    assert len(s_alive)==1 and s_alive[0]==L
    part_x_rot = part.x * np.cos(angle/180.*np.pi) + part.y * np.sin(angle/180.*np.pi)
    assert np.allclose(part_x_rot[mask_hit_front], x[mask_hit_front], atol=1e-12, rtol=0)
    assert np.allclose(part_x_rot[~mask_hit],      x[~mask_hit] + xp[~mask_hit]*L, atol=1e-12, rtol=0)
    s_hit_L = (jaw_L - x[mask_hit_angle_L]) / xp[mask_hit_angle_L]
    assert np.allclose(part.s[mask_hit_angle_L], s_hit_L, atol=1e-12, rtol=0)
    assert np.allclose(part_x_rot[mask_hit_angle_L], jaw_L, atol=1e-12, rtol=0)
    s_hit_R = (jaw_R - x[mask_hit_angle_R]) / xp[mask_hit_angle_R]
    assert np.allclose(part.s[mask_hit_angle_R], s_hit_R, atol=1e-12, rtol=0)
    assert np.allclose(part_x_rot[mask_hit_angle_R], jaw_R, atol=1e-12, rtol=0)


@for_all_test_contexts(
    excluding=('ContextCupy', 'ContextPyopencl')  # BlackAbsorber not on GPU
)
@retry()
def test_angled_with_tilts(test_context):
    angle=17.8
    jaw_LU, jaw_RU, jaw_LD, jaw_RD, L, coll = _make_absorber(angle=angle, tilts=[0.0005, -0.00015], _context=test_context)
    part, x, _, xp, _ = _generate_particles(four_dim=True,angle=angle, _context=test_context)
    coll.track(part)
    part.move(_context=xo.ContextCpu())
    part.sort(interleave_lost_particles=True)
    # Particles in front of the jaw are lost, ...
    mask_hit_front = (x >= jaw_LU) | (x <= jaw_RU)
    # but also those in the opening with an angle that would kick them on the jaw halfway
    mask_hit_angle_L = ~mask_hit_front & (x + xp*L >= jaw_LD)
    mask_hit_angle_R = ~mask_hit_front & (x + xp*L <= jaw_RD)
    mask_hit = mask_hit_front | mask_hit_angle_L | mask_hit_angle_R

    lost = np.unique(part.state[mask_hit])
    alive = np.unique(part.state[~mask_hit])
    assert len(lost)==1 and lost[0]==-340
    assert len(alive)==1 and alive[0]==1
    cL = np.cos(coll.tilt_L)
    cR = np.cos(coll.tilt_R)
    s_lost_front_L = part.s[mask_hit_front & (x > (jaw_LU + jaw_RU)/2)]
    assert np.all([s <= coll.length/2*(1 - cL) for s in s_lost_front_L])
    s_lost_front_R = part.s[mask_hit_front & (x < (jaw_LU + jaw_RU)/2)]
    assert np.all([s <= coll.length/2*(1 - cR) for s in s_lost_front_R])
    s_alive = np.unique(part.s[~mask_hit])
    assert len(s_alive)==1 and s_alive[0]==L
    with flaky_assertions():
        part_x_rot  = part.x * np.cos(np.deg2rad(angle))  + part.y * np.sin(np.deg2rad(angle))
        part_px_rot = part.px * np.cos(np.deg2rad(angle)) + part.py * np.sin(np.deg2rad(angle))
        assert np.allclose(part_x_rot[mask_hit_front] - part_px_rot[mask_hit_front]*part.s[mask_hit_front], x[mask_hit_front], atol=1e-12, rtol=0)
        assert np.allclose(part_x_rot[~mask_hit], x[~mask_hit] + xp[~mask_hit]*L, atol=1e-12, rtol=0)
        s_hit_L = (jaw_LU - x[mask_hit_angle_L] - (jaw_LD-jaw_LU)/2*(1-cL)/cL) / (xp[mask_hit_angle_L] - (jaw_LD-jaw_LU)/(L*cL))
        assert np.allclose(part.s[mask_hit_angle_L], s_hit_L, atol=1e-12, rtol=0)
        assert np.allclose(part_x_rot[mask_hit_angle_L], jaw_LU + (jaw_LD-jaw_LU)/(L*cL)*(s_hit_L - L/2*(1-cL)), atol=1e-12, rtol=0)
        s_hit_R = (jaw_RU - x[mask_hit_angle_R] - (jaw_RD-jaw_RU)/2*(1-cR)/cR) / (xp[mask_hit_angle_R] - (jaw_RD-jaw_RU)/(L*cR))
        assert np.allclose(part.s[mask_hit_angle_R], s_hit_R, atol=1e-12, rtol=0)
        assert np.allclose(part_x_rot[mask_hit_angle_R], jaw_RU + (jaw_RD-jaw_RU)/(L*cR)*(s_hit_R - L/2*(1-cR)), atol=1e-12, rtol=0)


def _make_absorber(angle=0, tilts=[0,0], _context=None):
    if _context is None:
        _context = xo.ContextCpu()
    co = [0.0075, -0.089]
    s = co[0] * np.cos(np.deg2rad(angle)) + co[1] * np.sin(np.deg2rad(angle))
    jaws = [s + 0.03,s - 0.01]
    L = 0.873
    coll = xc.BlackAbsorber(length=L, angle=angle, jaw=jaws, tilt=tilts, _context=_context)
    return coll.jaw_LU, coll.jaw_RU, coll.jaw_LD, coll.jaw_RD, L, coll

def _generate_particles(four_dim=False, angle=0, _context=None):
    if _context is None:
        _context = xo.ContextCpu()
    # Make particles
    x = np.random.uniform(-0.1, 0.1, n_part)
    y = np.random.uniform(-0.1, 0.1, n_part)
    if four_dim:
        px = np.random.uniform(-0.1, 0.1, n_part)
        py = np.random.uniform(-0.1, 0.1, n_part)
    else:
        px = 0
        py = 0
    ref = xp.Particles(mass0=xp.PROTON_MASS_EV, q0=1, p0c=7e12, _context=_context)
    part = xp.build_particles(x=x, y=y, px=px, py=py, particle_ref=ref, _context=_context)
    part_init = part.copy()
    part_init.move(_context=xo.ContextCpu())
    x_rot = part_init.x * np.cos(np.deg2rad(angle)) + part_init.y * np.sin(np.deg2rad(angle))
    y_rot = part_init.x * np.sin(np.deg2rad(angle)) + part_init.y * np.cos(np.deg2rad(angle))
    xp_rot = part_init.px * part_init.rpp * np.cos(np.deg2rad(angle)) + part_init.py * part_init.rpp * np.sin(np.deg2rad(angle))
    yp_rot = part_init.px * part_init.rpp * np.sin(np.deg2rad(angle)) + part_init.py * part_init.rpp * np.cos(np.deg2rad(angle))
    return part, x_rot, y_rot, xp_rot, yp_rot


@for_all_test_contexts(
    excluding=('ContextCupy', 'ContextPyopencl')  # BlackAbsorber not on GPU
)
@pytest.mark.parametrize("side, sign_R", [
                        ['+', 1], ['-', 1], ['+', -1], ['-', -1]]
                        , ids=["L R>0", "R R>0", "L R<0", "R R<0"])
def test_black_crystal(side, sign_R, test_context):
    ref = xp.Particles(mass0=xp.PROTON_MASS_EV, q0=1, p0c=7e12)
    x = np.random.uniform(-1, 1, n_part)
    px = np.random.uniform(-1, 1, n_part)
    y = np.random.uniform(-1, 1, n_part)
    py = np.random.uniform(-1, 1, n_part)
    part_init = xp.build_particles(x=x, px=px, y=y, py=py, particle_ref=ref)
    dri = xt.Drift(length=1) # To send the surviving particles further out

    for R in [0.87, 3.3, 17]:
        R = sign_R*R
        for tilt in [-89, -65, -42, -22, -8, 0, 11, 27, 48, 72, 89]:
            print(f"{R=}{tilt=}")
            length = 0.87
            sign = 1 if side == '+' else -1
            jaw = sign*0.13
            width = 0.2
            height = 0.5
            coll = xc.BlackCrystal(length=length, bending_radius=R, jaw=jaw, width=width, height=height, side=side, tilt=np.deg2rad(tilt))
            part = part_init.copy()
            coll.track(part)
            dri.track(part)

            Rpos = R - (sign_R-1)/2*width # full radius when positive
            Rneg = R - (sign_R+1)/2*width # shorter radius when positive (full radius when negative)

            # Four corners:
            x_TU = jaw + (1+sign)/2*width
            s_TD = Rneg*length/R
            x_TD = x_TU + Rneg*(1 - np.sqrt(1 - length**2/R**2))
            x_BU = jaw - (1-sign)/2*width
            s_BD = Rpos*length/R
            x_BD = x_BU + Rpos*(1 - np.sqrt(1 - length**2/R**2))
            Rs   = 0
            Rx   = x_BU + R if sign_R == 1 else x_TU + R

            mask_alive = part.state > 0
            assert np.allclose(part.s[mask_alive], length+1)

            # Rotate particles to tilted frame
            part_s =  part.s * np.cos(np.deg2rad(tilt)) + (part.x - x_BU) * np.sin(np.deg2rad(tilt))
            part_x = -part.s * np.sin(np.deg2rad(tilt)) + (part.x - x_BU)  * np.cos(np.deg2rad(tilt)) + x_BU

            mask_not_end = part.s < length+1
            mask_height = (part.y < height/2) & (part.y > -height/2)
            mask_front = np.isclose(part_s, 0) & mask_height
            assert np.all(part.state[mask_front] < 1)
            assert np.all(part_x[mask_front] <= x_TU)
            assert np.all(part_x[mask_front] >= x_BU)
            mask_upper_curve = np.isclose((part_x - Rx)**2 + (part_s - Rs)**2, Rneg**2) & mask_not_end & mask_height
            assert np.all(part.state[mask_upper_curve] < 1)
            mask_back = np.isclose(part_x, (x_TD - x_BD)/(s_TD - s_BD) * (part_s - s_BD) + x_BD) & mask_not_end & mask_height
            assert np.all(part.state[mask_back] < 1)
            mask_lower_curve = np.isclose((part_x - Rx)**2 + (part_s - Rs)**2, Rpos**2) & mask_not_end & mask_height
            assert np.all(part.state[mask_lower_curve] < 1)
            mask_top_face = np.isclose(part.y, height/2) & mask_not_end
            assert np.all(part.state[mask_top_face] < 1)
            mask_bottom_face = np.isclose(part.y, -height/2) & mask_not_end
            assert np.all(part.state[mask_bottom_face] < 1)
            assert np.all(mask_alive == ~mask_front & ~mask_upper_curve & ~mask_back & ~mask_lower_curve & ~mask_top_face & ~mask_bottom_face)

            # _, ax = plt.subplots(3, 2, figsize=(10, 5))
            # ax[0][0].scatter(part.s[~mask_alive], part.x[~mask_alive], s=0.1)
            # # ax[0][0].set_aspect('equal')
            # ax[0][1].scatter(part.s[mask_front], part.x[mask_front], s=0.1)
            # ax[0][1].scatter(part.s[mask_upper_curve], part.x[mask_upper_curve], s=0.1)
            # ax[0][1].scatter(part.s[mask_back], part.x[mask_back], s=0.1)
            # ax[0][1].scatter(part.s[mask_lower_curve], part.x[mask_lower_curve], s=0.1)
            # ax[1][0].scatter(part.s[mask_front], part.x[mask_front], s=0.1)
            # ax[1][1].scatter(part.s[mask_upper_curve], part.x[mask_upper_curve], s=0.1)
            # ax[2][0].scatter(part.s[mask_back], part.x[mask_back], s=0.1)
            # ax[2][1].scatter(part.s[mask_lower_curve], part.x[mask_lower_curve], s=0.1)
