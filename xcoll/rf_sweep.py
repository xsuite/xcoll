# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

import numpy as np
from warnings import warn
import scipy.constants as sc

import xtrack as xt


def prepare_rf_sweep(line, *, cavities=None, sweep=None, sweep_per_turn=None, num_turns=None):
    rf_sweep = RFSweep(line=line, cavities=cavities)
    rf_sweep.prepare(sweep=sweep, sweep_per_turn=sweep_per_turn, num_turns=num_turns)


class RFSweep:
    def __init__(self, line, cavities=None):
        self.line = line
        if line.env.ref_manager is None:
            raise ValueError("Environment must have a `ref_manager` to use RFSweep!" \
                           + "Do not use `optimize_for_tracking` as it will "
                           + "disable expressions, which are needed for RFSweep.")
        line.enable_time_dependent_vars = False  # To be able to twiss in the prepare step
        self.env = line.env
        self.cavities = cavities
        self.env['rf_sweep_df'] = 0
        self._prepared = False

    @property
    def current_sweep_value(self):
        return self.env['rf_sweep_df']

    def prepare(self, *, sweep=None, num_turns=None, sweep_per_turn=None):
        # TODO: if using sweep and num_turns, the sweep should stop after num_turns,
        #       but if using sweep_per_turn, it should continue indefinitely.
        if not self._prepared:
            if sweep_per_turn is not None:
                if np.isclose(sweep_per_turn, 0):
                    raise ValueError("Variable `sweep_per_turn` must be non-zero!")
                if sweep is not None:
                    raise ValueError("Provide either `sweep` or `sweep_per_turn`, not both.")
                if num_turns is not None:
                    raise ValueError("Variable `num_turns` cannot be set when using `sweep_per_turn`.")
                self.sweep_per_turn = sweep_per_turn
            elif sweep is not None:
                if np.isclose(sweep, 0):
                    raise ValueError("Variable `sweep` must be non-zero!")
                if num_turns is None:
                    num_turns = int(abs(sweep))
                if num_turns <= 0:
                    raise ValueError("When using `sweep`, `num_turns` must be a positive integer.")
                self.sweep_per_turn = sweep / num_turns
            else:
                raise ValueError("Either `sweep` or `sweep_per_turn` must be provided.")
            self._get_cavity_data()
            self._install_zeta_shift()
            self.tw = self.line.twiss()
            self.reset() # Initialize rf_sweep_df
            self.env['rf_sweep'].dzeta = f"{self.L} * rf_sweep_df / ({self.f_RF} + rf_sweep_df)"
            for cavs in self.cavities:
                for cav in cavs['names']:
                    scale_factor = int(np.round(cavs['freq'] / self.f_RF))
                    if scale_factor == 1:
                        self.env.ref[cav].frequency += self.env.ref["rf_sweep_df"]
                    else:
                        self.env.ref[cav].frequency += scale_factor * self.env.ref["rf_sweep_df"]
            self.line.enable_time_dependent_vars = True
            print("Enabled time-dependent variables in the line.")
            self._prepared = True


    def reset(self):
        if self.sweep_per_turn is None:
            raise ValueError("RFSweep not prepared. Call `prepare` first.")
        t_turn = self.tw.T_rev0   # To start sweeping from turn 0
        current_time = self.env['t_turn_s']
        if current_time == 0:
            self.env['rf_sweep_df'] = f"(t_turn_s + {t_turn}) / {t_turn} * {self.sweep_per_turn}"
        else:
            self.env['rf_sweep_df'] = f"(t_turn_s - {current_time} + {t_turn}) / {t_turn} * {self.sweep_per_turn}"


    def info(self, sweep=None, num_turns=None):
        if sweep is not None or num_turns is not None:
            if sweep is not None:
                warn("The `sweep` argument is deprecated.", DeprecationWarning)
            if num_turns is not None:
                warn("The `num_turns` argument is deprecated.", DeprecationWarning)
            self.prepare(sweep=sweep, num_turns=num_turns)
        if self.sweep_per_turn is None:
            raise ValueError("RFSweep not prepared. Call `prepare` first.")
        existing_sweep = self.current_sweep_value

        beta0 = self.line.particle_ref.beta0[0]
        E = self.line.particle_ref.energy0[0]
        q = self.line.particle_ref.q0
        h = self.f_RF * self.L / beta0 / sc.c
        eta = self.tw.slip_factor

        bucket_height = np.sqrt(abs(q*self.V*beta0**2 / (np.pi*h*eta*E) * (
                            2*np.cos(self.phi) + (2*self.phi-np.pi)*np.sin(self.phi)
                        )))
        delta_shift = -self.sweep_per_turn / self.f_RF / eta
        bucket_turns = bucket_height * 2 / abs(delta_shift)
        print(f"The current frequency is {self.f_RF + existing_sweep}Hz, adding "
            + f"{self.sweep_per_turn}Hz per turn.")
        print(f"This sweep will move the center of the bucket with \u0394\u03B4 = "
            + f"{delta_shift:.4} per turn.")
        print(f"The bucket height is {bucket_height:.4}, so this implies the sweep "
            + f"will shift one bucket every {round(bucket_turns, 2)} turns.")
        if self.tw.qs * bucket_turns < 3:
            print(f"Warning: This is a very fast sweep, moving ~1 bucket in "
                + f"~{round(self.tw.qs * bucket_turns, 2)} synchrotron oscillations "
                + f"(on average). If the bucket moves faster than a particle "
                + f"can follow, that particle will move out of the bucket and "
                + f"remain uncaptured.")


    # For backward compatibility
    def track(self, sweep=0, particles=None, num_turns=0, verbose=True, *args, **kwargs):
        warn("Using RFSweep.track() is deprecated. Please use RFSweep.prepare() "
           + "and then xt.Line.track().", DeprecationWarning)
        self.prepare(sweep=sweep, num_turns=num_turns)
        self.line.track(particles=particles, num_turns=num_turns, *args, **kwargs)


    def _get_cavity_data(self):
        # Cavity data not yet extracted
        tt = self.line.get_table()
        tt_c = tt.rows[[nn for nn, ttt in zip(tt.name, tt.element_type)
                        if 'Cavity' in ttt and 'CrabCavity' not in ttt]]
        mask = tt_c.parent_name == None  # Unsliced cavities
        if self.cavities is None:
            self.cavities = np.unique(list(tt_c.name[mask]) + list(tt_c.parent_name[~mask]))
        else:
            if isinstance(self.cavities, str):
                self.cavities = [self.cavities]
            for cav in self.cavities:
                if cav not in self.env.elements:
                    raise ValueError(f"Cavity `{cav}` not found in environment!")
        if len(self.cavities) == 0:
            raise ValueError("No cavities found in the line!")
        freq = np.array([self.env[cav].frequency for cav in self.cavities])
        volt = np.array([self.env[cav].voltage for cav in self.cavities])
        lag  = np.array([self.env[cav].lag for cav in self.cavities])
        s_pos = []
        name_first_in_line = []
        for cav in self.cavities:
            if cav in tt_c.name:
                s_pos.append(tt_c.rows[cav].s_start[0])
                name_first_in_line.append(cav)
            elif cav in tt_c.parent_name:
                s_pos.append(tt_c.rows[tt_c.parent_name == cav].s_start.min())
                name_first_in_line.append(tt_c.rows[tt_c.parent_name == cav].name[0])
            else:
                raise ValueError(f"Cavity `{cav}` not found in the line!")
        idx  = np.argsort(freq)
        boundaries = np.where(~np.isclose(freq[idx][1:], freq[idx][:-1],
                                        rtol=1e-8, atol=1e-12)
                            )[0] + 1
        groups = []
        for gidx in np.split(idx, boundaries):
            groups.append({
                'names': self.cavities[gidx],
                'freq': np.mean(freq[gidx]),
                'voltage': volt[gidx].sum(),
                's_pos': np.array(s_pos)[gidx],
                'lag': lag[gidx],
                'name_first_in_line': np.array(name_first_in_line)[gidx],
            })
        self.cavities = sorted(groups, key=lambda g: (g['voltage'], g['freq']), reverse=True)
        if np.isclose(self.cavities[0]['voltage'], 0):
            raise ValueError("No active cavity found!")
        self.f_RF = self.cavities[0]['freq']
        self.phi = np.deg2rad(self.cavities[0]['lag'].mean())
        self.V = self.cavities[0]['voltage']
        self.L = self.line.get_length()
        if len(self.cavities) > 1:
            print("Found multiple cavities with different frequencies:")
            for g in self.cavities:
                print(f"{g['freq']}Hz  at {g['voltage']}V: {g['names']}")
            print(r"The sweep will be performed with respect to the highest "
                + f"voltage cavity at {self.f_RF}Hz. The other cavities will "
                + f"be shifted accordingly.")


    def _install_zeta_shift(self):
        if 'rf_sweep' not in self.env.elements:
            idx = np.argsort(self.cavities[0]['s_pos'])
            first_cavity = self.cavities[0]['name_first_in_line'][idx[0]]
            self.env.elements['rf_sweep'] = xt.ZetaShift(dzeta=0)
            self.line.insert('rf_sweep', at=f'{first_cavity}@start')
