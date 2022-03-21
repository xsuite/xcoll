import numpy as np
import pandas as pd
from pyk2 import K2Collimator

from .beam_elements import Collimator
from .colldb import CollDB

class CollimatorManager:
    def __init__(self, *, line, coll_db_txt_file, emitx, emity):

        colldb = _load_colldb(coll_db_txt_file)
        colldb.insert(1,'s_center', None)
        colldb.insert(1,'opening', None)
        colldb['type'] = None
        
        self.colldb = colldb
        self.line = line
        self.emitx = emitx
        self.emity = emity

    @property
    def collimator_names(self):
        return list(self.colldb.index.values)


    def install_black_absorbers(self, names=None, verbose=False):
        line = self.line
        colldb = self.colldb
        if names is None:
            names = colldb.index.values
        mask = colldb.index.isin(names)
        for name in names:
            if name not in line.element_names:
                raise Exception(f"Collimator {name} not found in line!")
        if line.tracker is not None:
            raise Exception("Tracker already built! Please install collimators before building tracker!")

        colldb.loc[mask,'s_center'] = line.get_s_position(names)
        colldb.loc[mask,'type'] = 'BlackAbsorber'

        for name in names:
            if isinstance(line[name], Collimator):
                if verbose:
                    print(f"Collimator {name} already installed. Skipping...")
            elif isinstance(line[name], K2Collimator):
                raise ValueError(f"Collimator {name} already installed as K2Collimator! Please reconstruct the line.")
            else:
                if verbose:
                    print(f"Installing {name}")
                thiscoll = colldb.loc[name]
                newcoll = Collimator(
                        inactive_front=thiscoll['inactive_front'],
                        inactive_back=thiscoll['inactive_back'],
                        active_length=thiscoll['active_length'],
                        angle=thiscoll['angle'],
                        jaw_R=-1, jaw_L=1,
                        jaw_D=-1, jaw_U=1
                        )
                s_install = thiscoll['s_center'] - thiscoll['active_length']/2 - thiscoll['inactive_front']
                line.insert_element(element=newcoll, name=name, at_s=s_install)


    def install_k2_collimators(self, names=None, verbose=False):
        line = self.line
        colldb = self.colldb
        if names is None:
            names = colldb.index.values
        mask = colldb.index.isin(names)
        for name in names:
            if name not in line.element_names:
                raise Exception(f"Collimator {name} not found in line!")
        if line.tracker is not None:
            raise Exception("Tracker already built! Please install collimators before building tracker!")

        colldb.loc[mask,'s_center'] = line.get_s_position(names)
        colldb.loc[mask,'type'] = 'K2Collimator'

        for name in names:
            if isinstance(line[name], K2Collimator):
                if verbose:
                    print(f"Collimator {name} already installed. Skipping...")
            elif isinstance(line[name], Collimator):
                raise ValueError(f"Collimator {name} already installed as BlackAbsorber! Please reconstruct the line.")
            else:
                if verbose:
                    print(f"Installing {name}")
                thiscoll = colldb.loc[name]
                newcoll = K2Collimator(
                        k2_engine,
                        length,
                        rotation,
                        icoll,
                        aperture,
                        onesided,
                        dx,
                        dy,
                        dpx,
                        dpy,
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

        colldb = self.colldb
        incomplete = np.any([ np.any([ x is None for x in colldb[opt] ]) for opt in opt_funcs + ['sigmax', 'sigmay'] ])
        if recompute or incomplete:
            tracker = line.tracker
            tw = tracker.twiss(at_s=colldb['s_center'])
            for opt in opt_funcs:
                colldb[opt] = tw[opt]
            beta0_gamma0 = tracker.particle_ref._xobject.beta0[0] * tracker.particle_ref._xobject.gamma0[0]
            colldb['sigmax'] = np.sqrt(colldb['betx']*self.emitx/beta0_gamma0)
            colldb['sigmay'] = np.sqrt(colldb['bety']*self.emity/beta0_gamma0)


    def set_openings(self, gaps={}, recompute_optics=False):
        line = self.line
        if line is None or line.tracker is None:
            raise Exception("Please build tracker before calling this method!")
        colldb = self.colldb
        if any([ x is None for x in colldb.type ]):
            raise ValueError("Some collimators have not yet been installed. "
                             + "Please install all collimators before setting the openings.")
        for name, gap in gaps.items():
            if hasattr(gap, '__iter__'):
                if isinstance(gap, str):
                    raise ValueError("The gap setting has to be a number!")
                elif len(gap) > 2:
                    raise ValueError("The gap setting can maximally have two values (for the left and the right jaw)!")
                elif len(gap) == 1:
                    colldb.loc[name, 'gap'] = gap
                elif len(gap) == 2:
                    colldb.loc[name, 'gap'] = gap
                colldb.loc[name, 'gap'] = gap
            
        # Compute halfgap
        self._compute_optics(recompute=recompute_optics)
        colldb['opening'] = colldb['gap'] * np.sqrt(
            (colldb['sigmax']*np.cos(np.float_(colldb['angle'].values)*np.pi/180))**2
            + (colldb['sigmay']*np.sin(np.float_(colldb['angle'].values)*np.pi/180))**2
        )

        # Configure collimators
        for name in colldb.index:
            if isinstance(line[name], Collimator):
                line[name].dx = colldb['x'][name]
                line[name].dy = colldb['y'][name]
                line[name].jaw_R = -colldb['opening'][name] + colldb['offset'][name]
                line[name].jaw_L = colldb['opening'][name] + colldb['offset'][name]
            elif isinstance(line[name], K2Collimator):
                pass
            else:
                raise ValueError(f"Missing implementation for element type of collimator {name}!")


