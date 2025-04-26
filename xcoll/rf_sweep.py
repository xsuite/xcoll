# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

import numpy as np
import scipy.constants as sc
from time import perf_counter

import xtrack as xt
from xtrack.progress_indicator import progress

class RFSweep:

    def __init__(self, line):
        self.line = line
        self._get_base_frequency()
        self._install_zeta_shift()
        self.reset()

    def _get_base_frequency(self):
        self.cavities = self.line.get_elements_of_type(xt.Cavity)[1]
        freq = np.unique([round(self.line[cav].frequency, 9) for cav in self.cavities])
        if len(freq) > 1:
            raise NotImplementedError(f"Cannot sweep multiple cavities with different frequencies!")
        elif abs(freq[0]) < 1e-9 :
            raise ValueError(f"Cavity frequency not set!")
        self.f_RF = freq[0]
        phi = np.unique(np.array([self.line[cav].lag for cav in self.cavities])*np.pi/180.)
        if len(phi) > 1:
            raise NotImplementedError(f"Cannot sweep multiple cavities with different phases!")
        self.phi = phi[0]
        self.V = np.array([self.line[cav].voltage for cav in self.cavities]).sum()
        self.L = self.line.get_length()

    def _install_zeta_shift(self):
        if 'rf_sweep' in self.line.element_names:
            print(f"Found existing RF sweep in line.")
        else:
            s_cav = min([self.line.get_s_position(cav) for cav in self.cavities])
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
        return round(self.f_RF * dzeta/(self.L-dzeta), 6)


    def info(self, sweep=0, num_turns=0):
        if abs(sweep) > 1e-9:
            existing_sweep = self.current_sweep_value

            beta0 = self.line.particle_ref.beta0[0]
            E = self.line.particle_ref.energy0[0]
            q = self.line.particle_ref.q0
            h = self.f_RF * self.L / beta0 / sc.c
            tw = self.line.twiss()
            eta = tw.slip_factor

            bucket_height = np.sqrt(abs(q*self.V*beta0**2 / (np.pi*h*eta*E) * (
                                2*np.cos(self.phi) + (2*self.phi-np.pi)*np.sin(self.phi)
                            )))
            delta_shift = -sweep / self.f_RF / eta
            bucket_shift = delta_shift / bucket_height / 2
            if num_turns > 0:
                rf_shift_per_turn = sweep / num_turns
                print(f"The current frequency is {self.f_RF + existing_sweep}Hz, adding {rf_shift_per_turn}Hz "
                    + f"per turn until {self.f_RF + existing_sweep + sweep} (for {num_turns} turns).")
            print(f"This sweep will move the center of the bucket with \u0394\u03B4 = "
                + f"{delta_shift} ({bucket_shift} buckets).")
            if num_turns > 0 and num_turns < 3*bucket_shift/tw.qs:
                print(f"Warning: This is a very fast sweep, moving ~{round(bucket_shift, 2)} buckets in "
                    + f"~{round(num_turns*tw.qs, 2)} synchrotron oscillations (on average). If the "
                    + f"bucket moves faster than a particle can follow, that particle will move out of "
                    + f"the bucket and remain uncaptured.")


    def track(self, sweep=0, particles=None, num_turns=0, verbose=True, *args, **kwargs):

        # Was there a previous sweep?
        # If yes, we do not overwrite it but continue from there
        existing_sweep = self.current_sweep_value

        # Just set the new RF frequency, no tracking
        if num_turns == 0:
            sweep += existing_sweep
            if verbose:
                print(f"The current frequency is {self.f_RF + existing_sweep}Hz, moving to {self.f_RF + sweep}Hz."
                     + "No tracking performed.")
            self.line['rf_sweep'].dzeta = self.L * sweep / (self.f_RF + sweep)

        # Sweep and track
        else:
            if self.line.tracker is None:
                raise ValueError("Need to build tracker first!")
            if particles is None:
                raise ValueError("Need particles to track!")
            time = kwargs.pop('time', False)
            with_progress = kwargs.pop('with_progress', False)
            rf_shift_per_turn = sweep / num_turns

            # This is taken from xtrack.tracker.Tracker._track
            if time:
                t0 = perf_counter()
            if with_progress:
                if self.line.tracker.enable_pipeline_hold:
                    raise ValueError("Progress indicator is not supported with pipeline hold")
                if num_turns < 2:
                    raise ValueError('Tracking with progress indicator is only '
                                    'possible over more than one turn.')
                if with_progress is True:
                    batch_size = scaling = 100
                else:
                    batch_size = int(with_progress)
                    scaling = with_progress if batch_size > 1 else None
                if kwargs.get('turn_by_turn_monitor') is True:
                    ele_start = kwargs.get('ele_start') or 0
                    ele_stop = kwargs.get('ele_stop')
                    if ele_stop is None:
                        ele_stop = len(self.line)
                    if ele_start >= ele_stop:
                        # we need an additional turn and space in the monitor for
                        # the incomplete turn
                        num_turns += 1
                    _, monitor, _, _ = self.line.tracker._get_monitor(particles, True, num_turns)
                    kwargs['turn_by_turn_monitor'] = monitor

                for ii in progress(
                        range(0, num_turns, batch_size),
                        desc='Tracking',
                        unit_scale=scaling,
                ):
                    one_turn_kwargs = kwargs.copy()
                    is_first_batch = ii == 0
                    is_last_batch = ii + batch_size >= num_turns

                    if is_first_batch and is_last_batch:
                        # This is the only batch, we track as normal
                        pass
                    elif is_first_batch:
                        # Not the last batch, so track until the last element
                        one_turn_kwargs['ele_stop'] = None
                        one_turn_kwargs['num_turns'] = batch_size
                    elif is_last_batch:
                        # Not the first batch, so track from the first element
                        one_turn_kwargs['ele_start'] = None
                        remaining_turns = num_turns % batch_size
                        if remaining_turns == 0:
                            remaining_turns = batch_size
                        one_turn_kwargs['num_turns'] = remaining_turns
                        one_turn_kwargs['_reset_log'] = False
                    elif not is_first_batch and not is_last_batch:
                        # A 'middle batch', track from first to last element
                        one_turn_kwargs['num_turns'] = batch_size
                        one_turn_kwargs['ele_start'] = None
                        one_turn_kwargs['ele_stop'] = None
                        one_turn_kwargs['_reset_log'] = False
                    self._tracking_func(particles, rf_shift_per_turn, **one_turn_kwargs)
                    if not np.any(particles.state == 1):
                        break

            else:
                self._tracking_func(particles, rf_shift_per_turn, num_turns=num_turns, *args, **kwargs)

            if not np.any(particles.state == 1):
                print(f"All particles lost at turn {particles.at_turn.max()}, stopped sweep at "
                    + f"{self.current_sweep_value}Hz.")

            if time:
                t1 = perf_counter()
                self.line.tracker._context.synchronize()
                self.line.tracker.time_last_track = t1 - t0
            else:
                self.line.tracker.time_last_track = None


    def reset(self):
        self.line['rf_sweep'].dzeta = 0


    def _tracking_func(self, particles, rf_shift_per_turn, num_turns=1, *args, **kwargs):
        existing_sweep = self.current_sweep_value
        for i in range(num_turns):
            sweep = existing_sweep + i*rf_shift_per_turn
            self.line['rf_sweep'].dzeta = self.L * sweep / (self.f_RF + sweep)
#             for cav in cavities:
#                 self.line[cav].frequency = freq + sweep
            self.line.track(particles, num_turns=1, *args, **kwargs)
            if not np.any(particles.state == 1):
                break
