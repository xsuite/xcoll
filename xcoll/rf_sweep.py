# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2023.                 #
# ######################################### #

import numpy as np
import pandas as pd

import xtrack as xt


class RFSweep:

    def __init__(self, line):
        self.line = line
        self._get_base_frequency()
        self._install_zeta_shift()

    def _get_base_frequency(self):
        self.cavities = self.line.get_elements_of_type(xt.Cavity)[1]
        freq = np.unique([round(self.line[cav].frequency, 9) for cav in cavities])
        if len(freq) > 1:
            raise NotImplementedError(f"Cannot sweep multiple cavities with different frequencies!")
        elif abs(freq[0]) < 1e-9 :
            raise ValueError(f"Cavity frequency not set!")
        self.f0 = freq[0]

    def _install_zeta_shift(self):
        if 'rf_sweep' in self.line.element_names:
            raise ValueError(f"Found existing RF sweep in line; can only use one at a time!")
        s_cav = min([self.line.get_s_position(cav) for cav in cavities])
        if self.line.tracker is not None:
            self.line.unfreeze()
            line_was_built = True
        else:
            line_was_built = False
        self.line.insert_element(element=xt.ZetaShift(dzeta=0), name='rf_sweep', at_s=s_cav)
        if line_was_built:
            self.line.build_tracker()

    @property
    def current_sweep_value(self):
        dzeta = self.line['rf_sweep'].dzeta
        return round(self.f0 * dzeta/(self.machine_length-dzeta), 6)


    def track(self, sweep=0, num_turns=0, particles=None, verbose=True, *args, **kwargs):

        # Was there a previous sweep?
        # If yes, we do not overwrite it but continue from there
        existing_sweep = self.current_sweep_value

        # Some info
        scattering_enabled = False
        if self.line.tracker is not None:
            if self.scattering_enabled:
                scattering_enabled = True
                self.disable_scattering()
            tw = self.line.twiss()
            V = np.array([self.line[cav].voltage for cav in cavities]).sum()
            beta0 = self.line.particle_ref.beta0
            q = self.line.particle_ref.q0
            h = freq * tw.T_rev0
            eta = tw.slip_factor
            E = self.line.particle_ref.energy0
            phi = np.array([self.line[cav].lag for cav in cavities])[0]*np.pi/180
            bucket_height = np.sqrt(abs(q*V*beta0**2 / (np.pi*h*eta*E) * (2*np.cos(phi) +(2*phi-np.pi)*np.sin(phi))))[0]
            delta_shift = -sweep / freq / tw.slip_factor
            bucket_shift = delta_shift / bucket_height / 2
            if verbose:
                print(f"This sweep will move the center of the bucket with \u0394\u03B4 = "
                    + f"{delta_shift} ({bucket_shift} buckets).")

        # Just set the new RF frequency, no tracking
        if num_turns == 0:
            sweep += existing_sweep
            if verbose:
                print(f"The current frequency is {freq + existing_sweep}Hz, moving to {freq + sweep}Hz."
                     + "No tracking performed.")
            self.line['rf_sweep'].dzeta = self.machine_length * sweep / (freq + sweep)

        # Sweep and track
        else:
            if self.line.tracker is None:
                raise ValueError("Need to build tracker first!")
            if particles is None:
                raise ValueError("Need particles to track!")
            rf_shift_per_turn = sweep / num_turns
            if verbose:
                print(f"The current frequency is {freq + existing_sweep}Hz, sweeping {rf_shift_per_turn}Hz "
                    + f"per turn until {freq + existing_sweep + sweep} (for {num_turns} turns).")
            if num_turns < 3*bucket_shift/tw.qs:
                print(f"Warning: This is a very fast sweep, moving ~{round(bucket_shift,2)} buckets in "
                    + f"~{round(num_turns*tw.qs,2)} synchrotron oscillations (on average). If the "
                    + f"bucket moves faster than a particle can follow, that particle will move out of "
                    + f"the bucket and remain uncaptured.")
            if scattering_enabled:
                self.enable_scattering()
            if 'time' in kwargs and ['time']:
                self.line.tracker.time_last_track = 0
            for i in range(num_turns):
                sweep = existing_sweep + i*rf_shift_per_turn
                self.line['rf_sweep'].dzeta = self.machine_length * sweep / (freq + sweep)
#                 for cav in cavities:
#                     self.line[cav].frequency = freq + sweep
                if 'time' in kwargs and ['time']:
                    prev_time = self.line.time_last_track
                self.line.track(particles, num_turns=1, *args, **kwargs)
                if 'time' in kwargs and ['time']:
                    self.line.tracker.time_last_track += prev_time
                if not np.any(particles.state == 1):
                    if verbose:
                        print(f"All particles lost at turn {i}, stopped sweep at {i*rf_shift_per_turn}Hz.")
                    break