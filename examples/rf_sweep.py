from pathlib import Path
import time
start_time = time.time()
import matplotlib.pyplot as plt

import xtrack as xt
import xcoll as xc

beam = 1
plane = 'DPpos'

sweep = 300
sweep = -abs(sweep) if plane == 'DPpos' else abs(sweep)
num_turns = int(20*abs(sweep))

path_in  = xc._pkg_root.parent / 'examples'
path_out = Path.cwd()

# Load from json
line = xt.Line.from_json(path_in / 'machines' / f'lhc_run3_b{beam}.json')

line.build_tracker()
line.optimize_for_tracking()

num_particles = 2
part = line.build_particles(delta=[-2e-4, 2e-4], x_norm=0, px_norm=0, y_norm=0, py_norm=0)
monitor = xt.ParticlesMonitor(start_at_turn=1, stop_at_turn=num_turns, num_particles=2)

rf_sweep = xc.RFSweep(line)
rf_sweep.info(sweep=sweep, num_turns=num_turns)

rf_sweep.track(sweep=sweep, num_turns=num_turns, particles=part, time=True, turn_by_turn_monitor=monitor)

plt.figure(figsize=(12,8))
plt.plot(monitor.zeta.T,monitor.delta.T)

print(f"Total calculation time {time.time()-start_time}s")
plt.show()
