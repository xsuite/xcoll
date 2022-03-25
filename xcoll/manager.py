import numpy as np
import pandas as pd

from .beam_elements import BlackAbsorber
from .colldb import CollDB
from pyk2 import K2Collimator

class CollimatorManager:
    def __init__(self, *, line, colldb: CollDB):        
        self.colldb = colldb
        self.line = line

    @property
    def collimator_names(self):
        return list(self.colldb.name)

    @property
    def s_start(self):
        return self.colldb.s_center - self.colldb.active_length/2 - self.colldb.inactive_front

    @property
    def s_start_active(self):
        return self.colldb.s_center - self.colldb.active_length/2

    @property
    def s_center(self):
        return self.colldb.s_center

    @property
    def s_end_active(self):
        return self.colldb.s_center + self.colldb.active_length/2

    @property
    def s_end(self):
        return self.colldb.s_center + self.colldb.active_length/2 + self.colldb.inactive_back


    def install_black_absorbers(self, names=None, *, verbose=False):
        line = self.line
        df = self.colldb._colldb
        if names is None:
            names = self.collimator_names
        mask = df.index.isin(names)
        for name in names:
            if name not in line.element_names:
                raise Exception(f"Collimator {name} not found in line!")
            elif name not in self.collimator_names:
                print(f"Warning: Collimator {name} not found in CollDB! Ignoring...")
        if line.tracker is not None:
            raise Exception("Tracker already built! Please install collimators before building tracker!")

        df.loc[mask,'s_center'] = line.get_s_position(names)
        df.loc[mask,'collimator_type'] = 'BlackAbsorber'

        for name in names:
            if isinstance(line[name], BlackAbsorber):
                if verbose:
                    print(f"Collimator {name} already installed. Skipping...")
            elif isinstance(line[name], K2Collimator):
                raise ValueError(f"Collimator {name} already installed as K2Collimator! Please reconstruct the line.")
            else:
                if verbose:
                    print(f"Installing {name}")
                thiscoll = df.loc[name]
                newcoll = BlackAbsorber(
                            inactive_front=thiscoll['inactive_front'],
                            inactive_back=thiscoll['inactive_back'],
                            active_length=thiscoll['active_length'],
                            angle=thiscoll['angle'],
                            jaw_R=-1, jaw_L=1,
                            jaw_D=-1, jaw_U=1
                          )
                s_install = thiscoll['s_center'] - thiscoll['active_length']/2 - thiscoll['inactive_front']
                line.insert_element(element=newcoll, name=name, at_s=s_install)


    def install_k2_collimators(self, names=None, *, k2engine, verbose=False):
        line = self.line
        df = self.colldb._colldb
        icolls = { name: icoll for icoll, name in enumerate(self.collimator_names, start=1) }
        
        if names is None:
            names = self.collimator_names
        mask = df.index.isin(names)
        for name in names:
            if name not in line.element_names:
                raise Exception(f"Collimator {name} not found in line!")
            elif name not in self.collimator_names:
                print(f"Warning: Collimator {name} not found in CollDB! Ignoring...")
        if line.tracker is not None:
            raise Exception("Tracker already built! Please install collimators before building tracker!")

        df.loc[mask,'s_center'] = line.get_s_position(names)
        df.loc[mask,'type'] = 'K2Collimator'

        for name in names:
            if isinstance(line[name], K2Collimator):
                if verbose:
                    print(f"Collimator {name} already installed. Skipping...")
            elif isinstance(line[name], BlackAbsorber):
                raise ValueError(f"Collimator {name} already installed as BlackAbsorber! Please reconstruct the line.")
            else:
                if verbose:
                    print(f"Installing {name}")
                thiscoll = df.loc[name]
                side = 2 if thiscoll['onesided'] == 'right' else 1
                opening = thiscoll['opening_R'] if thiscoll['onesided'] == 'right' else thiscoll['opening_L']
                newcoll = K2Collimator(
                        k2engine=k2engine,
                        icoll=icolls[name],
                        aperture=opening,
                        onesided=side,
                        dx=0,
                        dy=0,
                        dpx=0,
                        dpy=0,
                        inactive_front=thiscoll['inactive_front'],
                        inactive_back=thiscoll['inactive_back'],
                        active_length=thiscoll['active_length'],
                        angle=thiscoll['angle'],
                        jaw_R=-1, jaw_L=1,
                        jaw_D=-1, jaw_U=1
                        )
                s_install = thiscoll['s_center'] - thiscoll['active_length']/2 - thiscoll['inactive_front']
                line.insert_element(element=newcoll, name=name, at_s=s_install)


    def _compute_optics(self, recompute=False):
        line = self.line
        if line is None or line.tracker is None:
            raise Exception("Please build tracker before calling this method!")
        opt_funcs = ['betx', 'bety', 'x', 'px', 'y', 'py'] 

        df = self.colldb._colldb
        incomplete = np.any([ np.any([ x is None for x in df[opt] ]) for opt in opt_funcs ])
        if recompute or incomplete:
            tracker = line.tracker
            df['opening_L'] = 1
            df['opening_R'] = 1
            tw = tracker.twiss(at_s=df['s_center'])
            for opt in opt_funcs:
                df[opt] = tw[opt]
            self.colldb.gamma_rel = tracker.particle_ref._xobject.gamma0[0]


    # The variable 'gaps' is meant to specify temporary settings that will overrule the CollDB.
    # As such, its settings will be applied to the collimator elements in the line, but not
    # written to the CollDB. Hence two successive calls to set_openings will not be combined,
    # and only the last call will be applied to the line.
    # The variable 'to_parking' will send all collimators that are not listed in 'gaps' to parking.
    # Similarily, the variable 'full_open' will set all openings of the collimators that are not
    # listed in 'gaps' to 1m.
    def set_openings(self, gaps={}, *, recompute_optics=False, to_parking=False, full_open=False):
        line = self.line
        if line is None or line.tracker is None:
            raise Exception("Please build tracker before calling this method!")
        colldb = self.colldb
        if any([ x is None for x in colldb.collimator_type ]):
            raise ValueError("Some collimators have not yet been installed. "
                             + "Please install all collimators before setting the openings.")
        if to_parking and full_open:
            raise ValueError("Cannot send collimators to parking and open them fully at the same time!")

        gaps_OLD = colldb.gap
        names = self.collimator_names
        # Override gap if sending to parking
        if to_parking:
            gaps = { **{ name: None for name in names }, **gaps }
        colldb.gap = gaps

        # Get the optics (to compute the opening)
        self._compute_optics(recompute=recompute_optics)

        # Configure collimators
        for name in names:
            # Override openings if opening fully
            if full_open and name not in gaps.keys():
                colldb._colldb.loc[name,'opening_L'] = 1
                colldb._colldb.loc[name,'opening_R'] = 1

            # Apply settings to element
            if isinstance(line[name], BlackAbsorber):
                line[name].dx = colldb.x[name]
                line[name].dy = colldb.y[name]
                line[name].dpx = 0
                line[name].dpy = 0
                line[name].angle = colldb.angle[name]
                line[name].jaw_R = -colldb._colldb['opening_R'][name] + colldb.offset[name]
                line[name].jaw_L = colldb._colldb['opening_L'][name] + colldb.offset[name]
            elif isinstance(line[name], K2Collimator):
                pass
            else:
                raise ValueError(f"Missing implementation for element type of collimator {name}!")
        colldb.gap = gaps_OLD


