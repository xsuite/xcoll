import numpy as np
import pandas as pd

from .beam_elements import BaseCollimator, BlackAbsorber, EverestCollimator, EverestCrystal, PyEverestCollimator, PyEverestCrystal, K2Collimator, K2Crystal, _all_collimator_types
from .scattering_routines.k2.engine import K2Engine
from .colldb import CollDB
from .tables import CollimatorImpacts
from .scattering_routines.everest.materials import SixTrack_to_xcoll

import xobjects as xo
import xpart as xp
import xtrack as xt



class CollimatorManager:
    def __init__(self, *, line, line_is_reversed=False, colldb: CollDB, capacity=1e6, record_impacts=False, \
                 _context=None, _buffer=None, io_buffer=None):

        if not isinstance(colldb, CollDB):
            raise ValueError("The variable 'colldb' needs to be an xcoll CollDB object!")
        else:
            self.colldb = colldb
        if not isinstance(line, xt.Line):
            raise ValueError("The variable 'line' needs to be an xtrack Line object!")
        else:
            self.line = line
        self._k2engine = None   # only needed for FORTRAN K2Collimator

        # Create _buffer, _context, and _io_buffer
        if _buffer is None:
            if _context is None:
                _context = xo.ContextCpu()
            _buffer = _context.new_buffer()
        elif _context is not None and _buffer.context != _context:
            raise ValueError("The provided buffer and context do not match! "
                             + "Make sure the buffer is generated inside the provided context, or alternatively, "
                             + "only pass one of _buffer or _context.")
        self._buffer = _buffer

        # TODO: currently capacity is only for io_buffer (hence for _impacts). Do we need it in the _buffer as well?
        self._capacity = int(capacity)
        if io_buffer is None:
            io_buffer = xt.new_io_buffer(_context=self._buffer.context, capacity=self.capacity)
        elif self._buffer.context != io_buffer._context:
            raise ValueError("The provided io_buffer lives on a different context than the buffer!")
        self._io_buffer = io_buffer

        # Initialise impacts table
        self._record_impacts = []
        self._impacts = None
        self.record_impacts = record_impacts

        self.tracker = None
        self._losmap = None
        self._coll_summary = None
        self._line_is_reversed = line_is_reversed


    @property
    def impacts(self):
        interactions = {
            -1: 'Black', 1: 'Nuclear-Inelastic', 2: 'Nuclear-Elastic', 3: 'pp-Elastic', 4: 'Single-Diffractive', 5: 'Coulomb'
        }
        n_rows = self._impacts._index + 1
        df = pd.DataFrame({
                'collimator':        [self.line.element_names[elemid] for elemid in self._impacts.at_element[:n_rows]],
                's':                 self._impacts.s[:n_rows],
                'turn':              self._impacts.at_turn[:n_rows],
                'interaction_type':  [ interactions[int_id] for int_id in self._impacts.interaction_type[:n_rows] ],
            })
        cols = ['id', 'x', 'px', 'y', 'py', 'zeta', 'delta', 'energy', 'mass', 'charge', 'z', 'a', 'pdgid']
        for particle in ['parent', 'child']:
            multicols = pd.MultiIndex.from_tuples([(particle, col) for col in cols])
            newdf = pd.DataFrame(index=df.index, columns=multicols)
            for col in cols:
                newdf[particle, col] = getattr(self._impacts,particle + '_' + col)[:n_rows]
            df = pd.concat([df, newdf], axis=1)
        return df

    @property
    def record_impacts(self):
        return self._record_impacts

    @record_impacts.setter
    def record_impacts(self, record_impacts):
        # TODO: how to get impacts if different collimator types in line?
        if record_impacts is True:
            record_impacts = self.collimator_names
        elif record_impacts is False or record_impacts is None:
            record_impacts = []
        record_start = set(record_impacts) - set(self._record_impacts)
        record_stop = set(self._record_impacts) - set(record_impacts)
        if record_start:
            if self._impacts is None:
                self._impacts = xt.start_internal_logging(io_buffer=self._io_buffer, capacity=self.capacity, \
                                                          elements=record_start)
            else:
                xt.start_internal_logging(io_buffer=self._io_buffer, capacity=self.capacity, \
                                          record=self._impacts, elements=record_start)
        if record_stop:
            if self.tracker is not None:
                self.tracker._check_invalidated()
            xt.stop_internal_logging(elements=record_stop)
        self._record_impacts = record_impacts

    @property
    def capacity(self):
        return self._capacity

    @capacity.setter
    def capacity(self, capacity):
        capacity = int(capacity)
        if capacity < self.capacity:
            raise NotImplementedError("Shrinking of capacity not yet implemented!")
        elif capacity == self.capacity:
            return
        else:
            self._io_buffer.grow(capacity-self.capacity)
            if self._impacts is not None:
                # TODO: increase capacity of iobuffer AND of _impacts
                raise NotImplementedError

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
                    inactive_front=thiscoll['inactive_front'],
                    inactive_back=thiscoll['inactive_back'],
                    active_length=thiscoll['active_length'],
                    angle=thiscoll['angle'],
                    is_active=False
                   )
        self._install_collimators(names, collimator_class=BlackAbsorber, install_func=install_func, verbose=verbose)


    def install_k2_collimators(self, names=None, *, max_part=50000, seed=None, verbose=False):
        # Check for the existence of a K2Engine; warn if settings are different
        self._k2engine = K2Engine(_capacity=max_part, random_generator_seed=seed)
        if self._k2engine._capacity != max_part:
            print(f"Warning: K2 already initiated with a maximum allocation of {self._k2engine._capacity} particles.\n"
                  + f"Ignoring the requested max_part={max_part}.")
        if seed is not None and self._k2engine.random_generator_seed != seed:
            print(f"Warning: K2 already initiated with seed {self._k2engine.random_generator_seed}.\n"
                  + f"Ignoring the requested seed={seed}.")

        # Do the installation
        def install_func(thiscoll, name):
            return K2Collimator(
                    inactive_front=thiscoll['inactive_front'],
                    inactive_back=thiscoll['inactive_back'],
                    active_length=thiscoll['active_length'],
                    angle=thiscoll['angle'],
                    material=SixTrack_to_xcoll[thiscoll['material']][0],
                    is_active=False
                   )
        self._install_collimators(names, collimator_class=K2Collimator, install_func=install_func, verbose=verbose)


    def install_everest_collimators(self, names=None, *, verbose=False):
        # Do the installation
        def install_func(thiscoll, name):
            return EverestCollimator(
                    inactive_front=thiscoll['inactive_front'],
                    inactive_back=thiscoll['inactive_back'],
                    active_length=thiscoll['active_length'],
                    angle=thiscoll['angle'],
                # TODO: Need to choose second element if crystal !
                # TODO: we should not use sixtrack materials here!!!
                    material=SixTrack_to_xcoll[thiscoll['material']][0],
                    is_active=False,
                    _tracking=False
                   )
        self._install_collimators(names, collimator_class=EverestCollimator, install_func=install_func, verbose=verbose)


    def install_pyeverest_collimators(self, names=None, *, verbose=False, random_seed=None):
        from .scattering_routines.pyeverest import set_random_seed
        set_random_seed(random_seed)

        # Do the installation
        def install_func(thiscoll, name):
            return PyEverestCollimator(
                    inactive_front=thiscoll['inactive_front'],
                    inactive_back=thiscoll['inactive_back'],
                    active_length=thiscoll['active_length'],
                    angle=thiscoll['angle'],
                    material=SixTrack_to_xcoll[thiscoll['material']][0],
                    is_active=False
                   )
        self._install_collimators(names, collimator_class=PyEverestCollimator, install_func=install_func, verbose=verbose)


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
                raise Exception(f"Warning: Collimator {name} not found in CollDB!...")
        if self.tracker_ready:
            raise Exception("Tracker already built!\nPlease install collimators before building tracker!")

        # Get collimator centers
        positions = dict(zip(names,line.get_s_position(names)))

        # Loop over collimators to install
        for name in names:

            # Check that collimator is not installed as different type
            # TODO: automatically replace collimator type and print warning
            if isinstance(line[name], tuple(_all_collimator_types - {collimator_class})):
                raise ValueError(f"Trying to install {name} as {collimator_class.__name__},"
                                 + f" but it is already installed as {type(line[name]).__name__}!\n"
                                 + "Please reconstruct the line.")

            # Check that collimator is not installed previously
            elif isinstance(line[name], collimator_class):
                if df.loc[name,'collimator_type'] != collimator_class.__name__:
                    raise Exception(f"Something is wrong: Collimator {name} already installed in line "
                                    + f"as {collimator_class.__name__} element, but registered in CollDB "
                                    + f"as {df.loc[name,'collimator_type']}. Please reconstruct the line.")
                if verbose:
                    print(f"Collimator {name} already installed. Skipping...")

            # TODO: only allow Marker elements, no Drifts!!
            #       How to do this with importing a line for MAD-X or SixTrack...?
            elif not isinstance(line[name], (xt.Marker, xt.Drift)):
                raise ValueError(f"Trying to install {name} as {collimator_class.__name__},"
                                 + f" but the line element to replace is not an xtrack.Marker (or xtrack.Drift)!\n"
                                 + "Please check the name, or correcft the element.")

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
    @property
    def installed(self):
        return not any([ x is None for x in self.colldb.collimator_type ])


    def align_collimators_to(self, align):
        if not self.installed:
            raise ValueError("Some collimators have not yet been installed.\n"
                             + "Please install all collimators before aligning the collimators.")
        self.colldb.align_to = align

    @property
    def aligned(self):
        return not np.any(
                    [x is None for x in self.colldb._colldb.s_align_front]
                    + [ x is None for x in self.colldb._colldb.s_align_back]
                )


    def build_tracker(self, **kwargs):
        kwargs.setdefault('_buffer', self._buffer)
        kwargs.setdefault('io_buffer', self._io_buffer)
        if kwargs['_buffer'] != self._buffer:
            raise ValueError("Cannot build tracker with different buffer than the CollimationManager buffer!")
        if kwargs['io_buffer'] != self._io_buffer:
            raise ValueError("Cannot build tracker with different io_buffer than the CollimationManager io_buffer!")
        if '_context' in kwargs and kwargs['_context'] != self._buffer.context:
            raise ValueError("Cannot build tracker with different context than the CollimationManager context!")
        self.tracker = self.line.build_tracker(**kwargs)
        return self.tracker

    @property
    def tracker_ready(self):
        return self.line is not None and self.line.tracker is not None


    def _compute_optics(self, recompute=False):
        if not self.tracker_ready:
            raise Exception("Please build tracker before computing the optics for the openings!")
        if not self.aligned:
            self.align_collimators_to('front')
        if recompute:
            self.colldb._optics_positions_to_calculate = { *set(self.colldb._colldb.s_align_front.values),\
                                                           *set(self.colldb._colldb.s_align_back.values) }
            self.colldb._optics = self.colldb._optics.iloc[0:0]

        pos = list(self.colldb._optics_positions_to_calculate)
        if pos != {}:
            tracker = self.line.tracker
#         TODO:
#         This is how we'll do it once twiss bug fixed (needed once inactive_front/back != 0):
#         (it is however much slower...)
#             tw = tracker.twiss(at_s=pos)
            tw = tracker.twiss()
            self.colldb._optics = pd.concat([
                                    self.colldb._optics,
                                    pd.DataFrame({
#                                     This is how we'll do it once twiss bug fixed:
#                                         opt: tw[opt] for opt in self.colldb._optics.columns
#                                     Current hack: getting the optics of the element within 1pm
#                                     (this will fail if the collimator is aligned to the centre)
                                        opt: [ np.array(tw[opt])[abs(tw['s']-thispos) < 1e-12][0] for thispos in pos ]
                                        for opt in self.colldb._optics.columns
                                    },index=pos)
                                ])
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
        if not self.tracker_ready:
            raise Exception("Please build tracker before setting the openings!")
        colldb = self.colldb
        if not self.installed:
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
        line = self.line
        for name in names:
            # Override openings if opening fully
            if full_open and name not in gaps.keys():
                line[name].is_active = False
            # Apply settings to element
            elif isinstance(line[name], BaseCollimator):
                line[name].dx = colldb.x[name]
                line[name].dy = colldb.y[name]
                line[name].angle = colldb.angle[name]
                line[name].jaw_F_L = colldb._colldb.jaw_F_L[name]
                line[name].jaw_F_R = colldb._colldb.jaw_F_R[name]
                line[name].jaw_B_L = colldb._colldb.jaw_B_L[name]
                line[name].jaw_B_R = colldb._colldb.jaw_B_R[name]
                line[name].is_active = colldb.is_active[name]
                if isinstance(line[name], (EverestCollimator, PyEverestCollimator, K2Collimator)):
                    line[name].material = colldb.material[name]
                    if colldb.onesided[name] == 'both':
                        line[name].onesided = False
                    elif colldb.onesided[name] == 'left':
                        line[name].onesided = True
                    elif colldb.onesided[name] == 'right':
                        raise ValueError(f"Right-sided collimators not implemented for Collimator {name}!")
                    line[name].is_active = colldb.is_active[name]
            else:
                raise ValueError(f"Missing implementation for element type of collimator {name}!")
        colldb.gap = gaps_OLD

    @property
    def openings_set(self):
        # TODO: need to delete jaw positions if some parameters (that would influence it) are changed
        return not np.any(
                    [x is None for x in self.colldb._colldb.jaw_F_L]
                    + [ x is None for x in self.colldb._colldb.jaw_F_R]
                    + [ x is None for x in self.colldb._colldb.jaw_B_L]
                    + [ x is None for x in self.colldb._colldb.jaw_B_R]
                )


    def generate_pencil_on_collimator(self, collimator, num_particles, *, side='+-', impact_parameter=1e-6, pencil_spread=1e-9,
                                     transverse_impact_parameter=0., transverse_spread_sigma=0.01, sigma_z=7.55e-2, zeta=None, delta=None):
        if not self.openings_set:
            raise ValueError("Need to set collimator openings before generating pencil distribution!")
        if not self.tracker_ready:
            raise Exception("Please build tracker before generating pencil distribution!")
        if transverse_impact_parameter != 0.:
            raise NotImplementedError

        if side == '+-':
            num_plus = int(num_particles/2)
            num_min  = int(num_particles - num_plus)
            part_plus = self.generate_pencil_on_collimator(collimator, num_plus, side='+',
                            impact_parameter=impact_parameter, pencil_spread=pencil_spread,
                            transverse_impact_parameter=transverse_impact_parameter,
                            transverse_spread_sigma=transverse_spread_sigma, sigma_z=sigma_z)
            part_min = self.generate_pencil_on_collimator(collimator, num_min, side='-',
                            impact_parameter=impact_parameter, pencil_spread=pencil_spread,
                            transverse_impact_parameter=transverse_impact_parameter,
                            transverse_spread_sigma=transverse_spread_sigma, sigma_z=sigma_z)
            part = xp.Particles.merge([part_plus, part_min])
            part.start_tracking_at_element = part_plus.start_tracking_at_element
            return part

        nemitt_x   = self.colldb.emittance[0]
        nemitt_y   = self.colldb.emittance[1]
        tracker    = self.tracker
        line       = self.line
        match_at_s = self.s_match[collimator]
        angle      = self.colldb.angle[collimator]
        sigma      = self.colldb._beam_size_front[collimator]
        dr_sigmas  = pencil_spread/sigma

        if abs(np.mod(angle-90,180)-90) < 1e-6:
            plane = 'x'
            co_pencil     = line[collimator].dx
    #         co_transverse = line[collimator].dy
        elif abs(np.mod(angle,180)-90) < 1e-6:
            plane = 'y'
            co_pencil     = line[collimator].dy
    #         co_transverse = line[collimator].dx
        else:
            raise NotImplementedError("Pencil beam on a skew collimator not yet supported!")

        if side == '+':
            absolute_cut = line[collimator].jaw_F_L + co_pencil + impact_parameter
        elif side == '-':
            absolute_cut = line[collimator].jaw_F_R + co_pencil - impact_parameter

        # Collimator plane: generate pencil distribution
        pencil, p_pencil = xp.generate_2D_pencil_with_absolute_cut(num_particles,
                        plane=plane, absolute_cut=absolute_cut, dr_sigmas=dr_sigmas,
                        side=side, tracker=tracker,
                        nemitt_x=nemitt_x, nemitt_y=nemitt_y,
                        at_element=collimator, match_at_s=match_at_s
        )

        # Other plane: generate gaussian distribution in normalized coordinates
        transverse_norm   = np.random.normal(scale=transverse_spread_sigma, size=num_particles)
        p_transverse_norm = np.random.normal(scale=transverse_spread_sigma, size=num_particles)

        # Longitudinal plane
        if zeta is None and delta is None:
            zeta, delta = xp.generate_longitudinal_coordinates(
                    num_particles=num_particles, distribution='gaussian', sigma_z=sigma_z, tracker=tracker
            )
        elif zeta is None:
            zeta = 0.0
#             delta = (self.line[collimator].jaw_F_L + self.line[collimator].dx + impact_parameter)/self.colldb.dx[collimator]
        elif delta is None:
            delta = 0.0

        if plane == 'x':
            part = xp.build_particles(
                    x=pencil, px=p_pencil, y_norm=transverse_norm, py_norm=p_transverse_norm,
                    zeta=zeta, delta=delta, nemitt_x=nemitt_x, nemitt_y=nemitt_y,
                    tracker=tracker, at_element=collimator, match_at_s=match_at_s
            )
        else:
            part = xp.build_particles(
                    x_norm=transverse_norm, px_norm=p_transverse_norm, y=pencil, py=p_pencil, 
                    zeta=zeta, delta=delta, nemitt_x=nemitt_x, nemitt_y=nemitt_y,
                    tracker=tracker, at_element=collimator, match_at_s=match_at_s
            )

        part._init_random_number_generator()

        return part



    def track(self, *args, **kwargs):
        if not self.tracker_ready:
            raise Exception("Please build tracker before tracking!")

        # Check if random generator is set
        if len(args) > 0:
            part = args[0]
        else:
            part = kwargs['particles']
        r1 = np.unique(part._rng_s1)
        r2 = np.unique(part._rng_s2)
        r3 = np.unique(part._rng_s3)
        # r4 = np.unique(part._rng_s4)  # not used
        if (len(r1)==1 and r1[0]==0) or (len(r2)==1 and r2[0]==0) or (len(r3)==1 and r3[0]==0):
            part._init_random_number_generator()

        # Prepare collimators for tracking
        for coll in self.collimator_names:
            self.line[coll]._tracking = True
        self.tracker.track(*args, **kwargs)
        # Flag collimators as no longer tracking
        for coll in self.collimator_names:
            self.line[coll]._tracking = False


    @property
    def lossmap(self):
        return self._lossmap

    def coll_summary(self, part):

        coll_s, coll_names, coll_length = self._get_collimator_losses(part)

        names = dict(zip(coll_s, coll_names))
        lengths = dict(zip(coll_s, coll_length))
        s = sorted(list(names.keys()))
        collname    =  [ names[pos] for pos in s ]
        colllengths =  [ lengths[pos] for pos in s ]
        nabs = []
        for pos in s:
            nabs.append(coll_s.count(pos))

        return pd.DataFrame({
            "collname": collname,
            "nabs":     nabs,
            "length":   colllengths,
            "s":        s
        })


    def create_lossmap(self, part, interpolation=0.1):
        # Loss location refinement
        if interpolation is not None:
            print("Performing the aperture losses refinement.")
            loss_loc_refinement = xt.LossLocationRefinement(self.tracker,
                    n_theta = 360, # Angular resolution in the polygonal approximation of the aperture
                    r_max = 0.5, # Maximum transverse aperture in m
                    dr = 50e-6, # Transverse loss refinement accuracy [m]
                    ds = interpolation, # Longitudinal loss refinement accuracy [m]
                    # save_refine_trackers=True # Diagnostics flag
                    )
            loss_loc_refinement.refine_loss_location(part)

        coll_s, coll_names, coll_length = self._get_collimator_losses(part)
        aper_s, aper_names              = self._get_aperture_losses(part)

        self._lossmap = {
            'collimator': {
                's':      coll_s,
                'name':   coll_names,
                'length': coll_length
            }
            ,
            'aperture': {
                's':    aper_s,
                'name': aper_names
            }
            ,
            'machine_length': self.line.get_length()
            ,
            'interpolation': interpolation
            ,
            'reversed': self._line_is_reversed
        }

        return self.lossmap

    def _get_collimator_losses(self, part):
        coll_names = [self.line.element_names[i] for i in part.at_element[part.state==-333]]
        # TODO: this way to get the collimator positions is a hack that needs to be cleaner with the new API
        coll_positions = dict(zip(self.collimator_names, self.s_center))
        coll_s = [coll_positions[name] for name in coll_names]
        coll_length = [self.line[i].active_length for i in part.at_element[part.state==-333]]
        machine_length = self.line.get_length()
        if self._line_is_reversed:
            coll_s = [ machine_length - s for s in coll_s ]

        return coll_s, coll_names, coll_length


    def _get_aperture_losses(self, part):
        aper_s = list(part.s[part.state==0])
        aper_names = [self.line.element_names[i] for i in part.at_element[part.state==0]]
        machine_length = self.line.get_length()
        if self._line_is_reversed:
            aper_s = [ machine_length - s for s in aper_s ]

        return aper_s, aper_names

