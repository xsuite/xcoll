import numpy as np
import pandas as pd

from .beam_elements import BlackAbsorber, K2Collimator, K2Engine
from .colldb import CollDB

_all_collimator_types = { BlackAbsorber, K2Collimator }

class CollimatorManager:
    def __init__(self, *, line, colldb: CollDB):        
        self.colldb = colldb
        self.line = line
        self._k2engine = None   # only needed for FORTRAN K2Collimator

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
        def install_func(thiscoll, name):
            return BlackAbsorber(
                    inactive_front=thiscoll['inactive_front'],
                    inactive_back=thiscoll['inactive_back'],
                    active_length=thiscoll['active_length'],
                    angle=thiscoll['angle'],
                    jaw_L=1, jaw_R=-1, jaw_U=1, jaw_D=-1
                   )
        self._install_collimators(names, collimator_class=BlackAbsorber, install_func=install_func, verbose=verbose)


    def install_k2_collimators(self, names=None, *, colldb_filename, max_part=50000, seed=None, verbose=False):        
        # Check for the existence of a K2Engine; warn if settings are different
        # (only one instance of K2Engine should exist).
        if self._k2engine is None:
            self._k2engine = K2Engine(n_alloc=max_part, colldb_filename=colldb_filename, random_generator_seed=seed)
        else:
            if self._k2engine.n_alloc != max_part:
                print(f"Warning: K2 already initiated with a maximum allocation of {self._k2engine.n_alloc} particles.\n"
                      + f"Ignoring the requested max_part={max_part}.")
            if self._k2engine.colldb_filename != colldb_filename:
                print(f"Warning: K2 already initiated from the file {self._k2engine.colldb_filename}.\n"
                      + f"Ignoring the requested colldb_filename={colldb_filename}.")
            if self._k2engine.random_generator_seed != seed:
                print(f"Warning: K2 already initiated with seed {self._k2engine.random_generator_seed}.\n"
                      + f"Ignoring the requested seed={seed}.")
        
        # Enumerate the collimators as expected by K2
        icolls = { name: icoll for icoll, name in enumerate(self.collimator_names, start=1) }
        
        # Do the installation
        def install_func(thiscoll, name):
            return K2Collimator(
                    k2engine=self._k2engine,
                    icoll=icolls[name],
                    inactive_front=thiscoll['inactive_front'],
                    inactive_back=thiscoll['inactive_back'],
                    active_length=thiscoll['active_length'],
                    angle=thiscoll['angle'],
                    jaw=1
                   )
        self._install_collimators(names, collimator_class=K2Collimator, install_func=install_func, verbose=verbose)
        

    def _install_collimators(self, names, *, collimator_class, install_func, verbose):
        # Check that collimator marker exists in Line and CollDB,
        # and that tracker is not yet built
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
            raise Exception("Tracker already built!\nPlease install collimators before building tracker!")

        # Get collimator centers
        positions = dict(zip(names,line.get_s_position(names)))

        # Loop over collimators to install
        for name in names:
            
            # Check that collimator is not installed as different type
            # TODO: automatically replace collimator type and print warning
            for other_coll_class in _all_collimator_types - {collimator_class}:
                if isinstance(line[name], other_coll_class):
                    raise ValueError(f"Trying to install {name} as {collimator_class.__name__},"
                                     + f" but it is already installed as {other_coll_class.__name__}!\n"
                                     + "Please reconstruct the line.")

            # Check that collimator is not installed previously
            if isinstance(line[name], collimator_class):
                if verbose:
                    print(f"Collimator {name} already installed. Skipping...")
            else:
                if verbose:
                    print(f"Installing {name}")
                # Get the settings from the CollDB
                thiscoll = df.loc[name]
                # Create the collimator element
                newcoll = install_func(thiscoll, name)
                # Update the position and type in the CollDB
                df.loc[name,'s_center'] = positions[name]
                df.loc[name,'collimator_type'] = collimator_class.__name__
                # Do the installation
                s_install = df.loc[name,'s_center'] - thiscoll['active_length']/2 - thiscoll['inactive_front']
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
            df['opening_upstr_L'] = 1
            df['opening_upstr_R'] = 1
            df['opening_downstr_L'] = 1
            df['opening_downstr_R'] = 1
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
            raise ValueError("Some collimators have not yet been installed.\n"
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
                colldb._colldb.loc[name,'opening_upstr_L'] = 1
                colldb._colldb.loc[name,'opening_upstr_R'] = 1
                colldb._colldb.loc[name,'opening_downstr_L'] = 1
                colldb._colldb.loc[name,'opening_downstr_R'] = 1

            # Apply settings to element
            if isinstance(line[name], BlackAbsorber):
                line[name].dx = colldb.x[name]
                line[name].dy = colldb.y[name]
                line[name].dpx = 0
                line[name].dpy = 0
                line[name].angle = colldb.angle[name]
                line[name].jaw_R = -colldb._colldb.opening_upstr_R[name] + colldb.offset[name]
                line[name].jaw_L = colldb._colldb.opening_upstr_L[name] + colldb.offset[name]
            elif isinstance(line[name], K2Collimator):
                line[name].dx = colldb.x[name]
                line[name].dy = colldb.y[name]
                line[name].dpx = colldb.px[name]
                line[name].dpy = colldb.py[name]
                line[name].angle = colldb.angle[name]
                line[name].jaw = colldb._colldb.opening_upstr_L[name]
                if colldb.onesided[name] == 'both':
                    line[name].onesided = False
                elif colldb.onesided[name] == 'left':
                    line[name].onesided = True
                elif colldb.onesided[name] == 'right':
                    raise ValueError(f"Right-sided collimators not implemented for K2Collimator {name}!")
                line[name].offset = colldb.offset[name]
            else:
                raise ValueError(f"Missing implementation for element type of collimator {name}!")
        colldb.gap = gaps_OLD


