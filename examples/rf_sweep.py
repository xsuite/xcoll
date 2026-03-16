# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

import time
start_time = time.time()
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt

import xobjects as xo
import xtrack as xt
import xcoll as xc
from xaux import plot_multi


# Setup
# =====

# path_in = Path(__file__).parent
path_in = xc._pkg_root.parent / 'examples'
path_out = Path.cwd() / 'plots' / 'rf_sweep'

beam = 1
env = xt.load(path_in / 'machines' / f'lhc_run3_b{beam}.json')
line = env[f'lhcb{beam}']

# Prepare cavities
tw = line.twiss()
L = line.get_length()
for cav, nn in zip(*line.get_elements_of_type(xt.Cavity)):
    frequency = round(cav.frequency * tw.T_rev0) / tw.T_rev0  # exactly correct frequency
    cav.frequency = frequency
    cav.absolute_time = True


# Reference: bucket without sweep
# ===============================

num_turns = 2000
num_particles = 21
delta0 = 2.122911e-7
zeta0 = 4.5858755e-6
line.discard_tracker()
line.build_tracker(_context=xo.ContextCpu(omp_num_threads=12))

delta_values = np.linspace(-7.5e-4, 7.5e-4, num_particles) + delta0
part_init = line.build_particles(delta=delta_values, zeta=zeta0, x_norm=0,
                                 px_norm=0, y_norm=0, py_norm=0)
part = part_init.copy()
monitor_ref = xt.ParticlesMonitor(start_at_turn=1, stop_at_turn=num_turns,
                                  num_particles=num_particles)
line.track(particles=part, num_turns=num_turns, with_progress=1,
           turn_by_turn_monitor=monitor_ref)


# Fixed 50Hz shift of the RF bucket
# =================================

df = 50
for cav, nn in zip(*line.get_elements_of_type(xt.Cavity)):
    env.elements[f"{nn}_shift"] = xc.beam_elements.SweepCavity(cavity=cav,
                                                               df=df)
    line.replace(nn, f"{nn}_shift")

line.discard_tracker()
line.build_tracker(_context=xo.ContextCpu(omp_num_threads=12))

part = part_init.copy()
monitor_shift = xt.ParticlesMonitor(start_at_turn=1, stop_at_turn=num_turns,
                                    num_particles=num_particles)
line.track(particles=part, num_turns=num_turns, with_progress=1,
           turn_by_turn_monitor=monitor_shift)

fig, axes = plt.subplots(1, 2, figsize=(12,6))
plot_multi(monitor_ref.zeta.T, monitor_ref.delta.T, ax=axes[0], cmap='berlin')
dzeta = - df/frequency * L * np.arange(1, num_turns)
dzeta_per_particle = np.array([dzeta for _ in range(num_particles)]).T # Accumulated dzeta per particle
plot_multi(monitor_shift.zeta.T + dzeta_per_particle,
           monitor_shift.delta.T, ax=axes[1], cmap='berlin')
axes[0].set_xlim(-0.7, 0.7)
axes[0].set_ylim(-9e-4, 9e-4)
axes[1].set_xlim(-0.7, 0.7)
axes[1].set_ylim(-9e-4, 9e-4)
plt.tight_layout()
plt.savefig(path_out / 'fixed_shift.png', dpi=300)
plt.close()


# Sweep -300Hz over 6000 turns
# ============================

num_turns = 6000
df_per_turn = -300/num_turns
tt = line.get_table()
previous_sweeps = tt.rows['.*_shift'].name
for nn in previous_sweeps:
    nn_new = nn.replace('_shift', '_sweep')
    env.elements[nn_new] = xc.beam_elements.SweepCavity(cavity=line[nn].cavity,
                                            df_per_turn=df_per_turn)
    line.replace(nn, nn_new)

line.discard_tracker()
line.build_tracker(_context=xo.ContextCpu(omp_num_threads=12))

part = part_init.copy()
monitor_sweep = xt.ParticlesMonitor(start_at_turn=1, stop_at_turn=num_turns, num_particles=num_particles)
line.track(particles=part, num_turns=num_turns, turn_by_turn_monitor=monitor_sweep, with_progress=1)

fig, axes = plt.subplots(1, 2, figsize=(12,6))
plot_multi(monitor_ref.zeta.T, monitor_ref.delta.T, ax=axes[0], cmap='berlin')
df = df_per_turn * np.arange(1, num_turns)
dzeta = - df/frequency * L * np.arange(1, num_turns)
dzeta_per_particle = np.array([dzeta for _ in range(num_particles)]).T # Accumulated dzeta per particle
zeta_sweep = monitor_sweep.zeta.T + dzeta_per_particle
delta_sweep = monitor_sweep.delta.T
dead = monitor_sweep.state.T < 1  # Do not plot lost particles
zeta_sweep[dead] = np.nan
delta_sweep[dead] = np.nan
plot_multi(zeta_sweep, delta_sweep, ax=axes[1], cmap='berlin')
axes[0].set_xlim(-3.3, 0.7)
axes[0].set_ylim(-9e-4, 15e-4)
axes[1].set_xlim(-3.3, 0.7)
axes[1].set_ylim(-9e-4, 15e-4)
plt.tight_layout()
plt.savefig(path_out / 'simple_sweep.png', dpi=300)
plt.close()


# Complicated example: sweep with pauses
# ======================================

num_turns = 0
num_particles = 11
df = np.array([0], dtype=np.float64)
step_turns = 5000
for _ in range(3):
    df = np.concatenate((df, df[-1] + np.array(-100/step_turns * np.arange(1, step_turns+1))))
    df = np.concatenate((df, df[-1]*np.ones(step_turns)))
    num_turns += 2*step_turns
df = df[1:]  # Remove initial zero

tt = line.get_table()
previous_sweeps = tt.rows['.*_sweep'].name
for nn in previous_sweeps:
    nn_new = nn.replace('_sweep', '_sweep_with_pauses')
    env.elements[nn_new] = xc.beam_elements.SweepCavity(cavity=line[nn].cavity, df=df)
    line.replace(nn, nn_new)

line.discard_tracker()
line.build_tracker(_context=xo.ContextCpu(omp_num_threads=12))

delta_values = np.linspace(-1e-4, 1e-4, num_particles) + delta0
part = line.build_particles(delta=delta_values, zeta=zeta0, x_norm=0, px_norm=0, y_norm=0, py_norm=0)
monitor_pause = xt.ParticlesMonitor(start_at_turn=1, stop_at_turn=num_turns, num_particles=num_particles)
line.track(particles=part, num_turns=num_turns, turn_by_turn_monitor=monitor_pause, with_progress=1)

fig, ax = plt.subplots(1, 1, figsize=(12,6))

# dzeta = -df/frequency * L * np.arange(num_turns)
dzeta = -np.cumsum(df[:-1])/frequency * L
dzeta_per_particle = np.array([dzeta for _ in range(num_particles)]).T # Accumulated dzeta per particle
zeta_pause = monitor_pause.zeta.T + dzeta_per_particle
delta_pause = monitor_pause.delta.T
dead = monitor_pause.state.T < 1  # Do not plot lost particles
zeta_pause[dead] = np.nan
delta_pause[dead] = np.nan
plot_multi(zeta_pause, delta_pause, ax=ax, cmap='berlin')
# axes[0].set_xlim(-3.3, 0.7)
# axes[0].set_ylim(-9e-4, 15e-4)
# axes[1].set_xlim(-3.3, 0.7)
# axes[1].set_ylim(-9e-4, 15e-4)
plt.tight_layout()
plt.savefig(path_out / 'sweep_with_pauses.png', dpi=300)
plt.close()

plot_multi(zeta_pause, ax=ax, cmap='berlin')
# axes[0].set_xlim(-3.3, 0.7)
# axes[0].set_ylim(-9e-4, 15e-4)
# axes[1].set_xlim(-3.3, 0.7)
# axes[1].set_ylim(-9e-4, 15e-4)
plt.tight_layout()
plt.savefig(path_out / 'sweep_with_pauses_zeta.png', dpi=300)
plt.close()
