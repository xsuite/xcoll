import numpy as np
import pandas as pd

from .beam_elements import BlackAbsorber, K2Collimator, K2Engine
from .colldb import CollDB
from .tables import CollimatorImpacts

import xtrack as xt
import xobjects as xo

_all_collimator_types = { BlackAbsorber, K2Collimator }


class CollimatorManager:
    def __init__(self, *, line, colldb: CollDB, _context=None, _buffer=None, storage_capacity=1e6, record_impacts=False):
        if not isinstance(colldb, CollDB):
            raise ValueError("The variable 'colldb' needs to be an xcoll CollDB object!")
        else:
            self.colldb = colldb
        if not isinstance(line, xt.Line):
            raise ValueError("The variable 'line' needs to be an xtrack Line object!")
        else:
            self.line = line
        self._k2engine = None   # only needed for FORTRAN K2Collimator

        if _buffer is None:
            if _context is None:
                _context = xo.ContextCpu()
            _buffer = _context.new_buffer()
        elif _context is not None and _buffer.context != _context:
            raise ValueError("The provided buffer and context do not match! "
                             + "Make sure the buffer is generated inside the provided context, or alternatively, "
                             + "only pass one of _buffer or _context.")
        self._buffer = _buffer
        self.storage_capacity = storage_capacity
        self.record_impacts = record_impacts


    @property
    def impacts(self):
        interactions = {
            -1: 'Black', 1: 'Nuclear-Inelastic', 2: 'Nuclear-Elastic', 3: 'pp-Elastic', 4: 'Single-Diffractive', 5: 'Coulomb'
        }
        n_rows = self._impacts._row_id + 1
        df = pd.DataFrame({
                'collimator':        [self.line.element_names[elemid] for elemid in self._impacts.at_element[:n_rows]],
                's':                 self._impacts.s[:n_rows],
                'turn':              self._impacts.turn[:n_rows],
                'interaction_type':  [ interactions[int_id] for int_id in self._impacts.interaction_type[:n_rows] ],
            })
        cols = ['id', 'x', 'px', 'y', 'py', 'zeta', 'delta', 'energy']
        for particle in ['parent', 'child']:
            multicols = pd.MultiIndex.from_tuples([(particle, col) for col in cols])
            newdf = pd.DataFrame(index=df.index, columns=multicols)
            for col in cols:
                newdf[particle, col] = getattr(self._impacts,col + '_' + particle)[:n_rows]
            df = pd.concat([df, newdf], axis=1)
        return df

    @property
    def record_impacts(self):
        return self._record_impacts

    @record_impacts.setter
    def record_impacts(self, record_impacts):
        self._record_impacts = record_impacts
        if record_impacts:
            self._impacts = CollimatorImpacts(_capacity=self.storage_capacity, _buffer=self._buffer)
        else:
            self._impacts = None
        # Update the xo.Ref to the CollimatorImpacts for the installed collimators
        for name in self.collimator_names:
            if self.colldb._colldb.loc[name,'collimator_type'] is not None:
                self.line[name].impacts = self._impacts

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

    @property
    def s_match(self):
        return self.colldb.s_match

    def install_black_absorbers(self, names=None, *, verbose=False):
        def install_func(thiscoll, name):
            return BlackAbsorber(
                    _buffer=self._buffer,
                    impacts=self._impacts,
                    inactive_front=thiscoll['inactive_front'],
                    inactive_back=thiscoll['inactive_back'],
                    active_length=thiscoll['active_length'],
                    angle=thiscoll['angle'],
                    is_active=False
                   )
        self._install_collimators(names, collimator_class=BlackAbsorber, install_func=install_func, verbose=verbose)


    def install_k2_collimators(self, names=None, *, max_part=50000, seed=None, verbose=False):
        # Check for the existence of a K2Engine; warn if settings are different
        # (only one instance of K2Engine should exist).
        if self._k2engine is None:
            self._k2engine = K2Engine(n_alloc=max_part, random_generator_seed=seed)
        else:
            if self._k2engine.n_alloc != max_part:
                print(f"Warning: K2 already initiated with a maximum allocation of {self._k2engine.n_alloc} particles.\n"
                      + f"Ignoring the requested max_part={max_part}.")
            if self._k2engine.random_generator_seed != seed:
                print(f"Warning: K2 already initiated with seed {self._k2engine.random_generator_seed}.\n"
                      + f"Ignoring the requested seed={seed}.")

        # Do the installation
        def install_func(thiscoll, name):
            return K2Collimator(
                    k2engine=self._k2engine,
                    impacts=self._impacts,
                    inactive_front=thiscoll['inactive_front'],
                    inactive_back=thiscoll['inactive_back'],
                    active_length=thiscoll['active_length'],
                    angle=thiscoll['angle'],
                    is_active=False
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
                if df.loc[name,'collimator_type'] != collimator_class.__name__:
                    raise Exception(f"Something is wrong: Collimator {name} already installed in line "
                                    + f"as {collimator_class.__name__} element, but registered in CollDB "
                                    + f"as {df.loc[name,'collimator_type']}. Please reconstruct the line.")
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
                if name+'_aper' in line.element_names:
                    coll_aper = line[name+'_aper']
                    assert coll_aper.__class__.__name__.startswith('Limit')
                    if np.any([name+'_aper_tilt_' in nn for nn in line.element_names]):
                        raise NotImplementedError("Collimator apertures with tilt not implemented!")
                    if np.any([name+'_aper_offset_' in nn for nn in line.element_names]):
                        raise NotImplementedError("Collimator apertures with offset not implemented!")
                else:
                    coll_aper = None

                line.insert_element(element=newcoll, name=name, at_s=s_install)

                if coll_aper is not None:
                    line.insert_element(element=coll_aper, name=name+'_aper_front', index=name)
                    line.insert_element(element=coll_aper, name=name+'_aper_back',
                                        index=line.element_names.index(name)+1)



    def align_collimators_to(self, align):
        if any([ x is None for x in self.colldb.collimator_type ]):
            raise ValueError("Some collimators have not yet been installed.\n"
                             + "Please install all collimators before aligning the collimators.")
        self.colldb.align_to = align


    def build_tracker(self, **kwargs):
        kwargs.setdefault('_buffer', self._buffer)
        if kwargs['_buffer'] != self._buffer:
            raise ValueError("Cannot build tracker with different buffer than the CollimationManager buffer!")
        if '_context' in kwargs and kwargs['_context'] != self._buffer.context:
            raise ValueError("Cannot build tracker with different context than the CollimationManager context!")
        self.tracker = self.line.build_tracker(**kwargs)
        return self.tracker


    def _compute_optics(self, recompute=False):
        line = self.line
        if line is None or line.tracker is None:
            raise Exception("Please build tracker before computing the optics for the openings!")
        if np.any(
            [x is None for x in self.colldb._colldb.s_align_front]
            + [ x is None for x in self.colldb._colldb.s_align_back]
        ):
            raise Exception("Not all collimators are aligned! Please call 'align_collimators_to' "
                            + "on the CollimationManager before computing the optics for the openings!")

        pos = self.colldb._optics_positions_to_calculate
        if recompute or pos != {}:
            tracker = self.line.tracker
            # Calculate optics without collimators
            old_val = {}
            for name in self.collimator_names:
                old_val[name] = line[name].is_active
                line[name].is_active = False
            tw = tracker.twiss(at_s=pos)
            self.colldb._optics = pd.concat([
                                    self.colldb._optics,
                                    pd.DataFrame({
                                        opt: tw[opt] for opt in self.colldb._optics.columns
                                    },index=pos)
                                ])
            for name in self.collimator_names:
                line[name].is_active = old_val[name]
            self.colldb._optics_positions_to_calculate = {}
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
        if not self.colldb._optics_is_ready:
            raise Exception("Something is wrong: not all optics needed for the jaw openings are calculated!")

        # Configure collimators
        for name in names:
            # Override openings if opening fully
            if full_open and name not in gaps.keys():
                line[name].is_active = False
            # Apply settings to element
            elif isinstance(line[name], BlackAbsorber):
                line[name].dx = colldb.x[name]
                line[name].dy = colldb.y[name]
                line[name].angle = colldb.angle[name]
                line[name].jaw_F_L = colldb._colldb.jaw_F_L[name]
                line[name].jaw_F_R = colldb._colldb.jaw_F_R[name]
                line[name].jaw_B_L = colldb._colldb.jaw_B_L[name]
                line[name].jaw_B_R = colldb._colldb.jaw_B_R[name]
                line[name].is_active = colldb.is_active[name]
            elif isinstance(line[name], K2Collimator):
                line[name].material = colldb.material[name]
                line[name].dx = colldb.x[name]
                line[name].dy = colldb.y[name]
                line[name].dpx = colldb.px[name]   # This is a K2 curiosity; we don't want it in our future code
                line[name].dpy = colldb.py[name]   # This is a K2 curiosity; we don't want it in our future code
                line[name].angle = colldb.angle[name]
                line[name].jaw_F_L = colldb._colldb.jaw_F_L[name]
                line[name].jaw_F_R = colldb._colldb.jaw_F_R[name]
                line[name].jaw_B_L = colldb._colldb.jaw_B_L[name]
                line[name].jaw_B_R = colldb._colldb.jaw_B_R[name]
                if colldb.onesided[name] == 'both':
                    line[name].onesided = False
                elif colldb.onesided[name] == 'left':
                    line[name].onesided = True
                elif colldb.onesided[name] == 'right':
                    raise ValueError(f"Right-sided collimators not implemented for K2Collimator {name}!")
                line[name].is_active = colldb.is_active[name]
            else:
                raise ValueError(f"Missing implementation for element type of collimator {name}!")
        colldb.gap = gaps_OLD


    def track(self, *args, **kwargs):
        self.tracker.track(*args, **kwargs)


