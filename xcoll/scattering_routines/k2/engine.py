# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

import numpy as np
import os
from pathlib import Path
import shutil

import xtrack as xt
import xpart as xp


_FILE_NAME = "temp_k2_colldb.dat"

class K2Engine:
    # The engine is a singleton
    def __new__(cls, **kwargs):
        if not hasattr(cls, 'instance'):
            cls.instance = super().__new__(cls)
            cls.instance._initialised = False
        return cls.instance

    def __init__(self, **kwargs):
        if(self._initialised):
            for kk, vv in kwargs.items():
                if not hasattr(self, kk):
                    raise ValueError(f"Invalid attribute {kk} for {self.__class__.__name__}!")
                if kk == 'capacity' or kk == '_capacity':
                    self.capacity = vv   # This ensure the setter method (which casts to int32) is called
                else:
                    setattr(self, kk, vv)
            return
        self._cwd = False
        self._old_cwd = False
        self._initialised = True
        self._warning_given = False
        self._file = None
        self._particle_ref = xp.Particles()
        self._collimator_dict = {}
        self.capacity = kwargs.get('_capacity', 50000)
        self.seed = kwargs.get('seed', None)

    def __del__(self, *args, **kwargs):
        self.stop()

    def _warn(self):
        if not self.instance._warning_given:
            print("Warning: Failed to import pyK2 (did you compile?).\n"
                + "_K2collimators will be installed but are not trackable.")
            self._warning_given = True


    @classmethod
    def start(cls, *, line=None, elements=None, names=None, seed=None, cwd=None, particle_ref=None,
              p0c=None, **kwargs):
        from ...beam_elements.k2 import _K2Collimator, _K2Crystal
        from .sixtrack_input import create_dat_file
        cls(**kwargs)
        this = cls.instance
        try:
            from .pyk2f import pyk2_init
        except ImportError:
           this._warn()
           return

        # Get reference particle
        if line is None or not hasattr(line, 'particle_ref') or line.particle_ref is None:
            if particle_ref is None and not this._has_particle_ref():
                raise ValueError("Need to provide either `line` or `particle_ref`.")
            if not isinstance(particle_ref, xp.Particles):
                if p0c is None:
                    raise ValueError("When providing `particle_ref`, it should be an "
                                        "xp.Particles object or a PDG ID. In the latter case, "
                                        "provide `p0c` as well.")
                particle_ref = xt.Particles.reference_from_pdg_id('proton', p0c=p0c)
            else:
                if p0c is not None:
                    particle_ref.p0c = p0c
            if line is not None:
                line.particle_ref = particle_ref
        else:
            if particle_ref is not None or p0c is not None:
                raise ValueError("Cannot provide both `particle_ref` or `p0c`, and a line "
                               + "with a reference particle.")
            particle_ref = line.particle_ref
        this.particle_ref = particle_ref

        # Set seed
        if seed is None:
            if this._seed is None:
                this._seed = np.int32(np.random.randint(1, int(1.e9)))
        else:
            this._seed = np.int32(abs(np.int32(seed)))
            if this._seed != int(seed):
                print(f"Warning: overflow for seed {seed}.")
            # Setting a seed here does not count as manually setting it
            this._seed_set_manually = False
        print(f"Using seed {this._seed}.")

        # Set working directory
        if cwd is None:
            from xaux.tools import ranID
            cwd = 'k2_run_' + ranID()
        cwd = Path(cwd).expanduser().resolve()
        cwd.mkdir(parents=True, exist_ok=True)
        this._old_cwd = Path.cwd()
        os.chdir(cwd)
        this._cwd = cwd

        # Get collimators
        if line is None:
            if elements is None:
                raise ValueError("Need to provide either `line` or `elements`.")
            if not hasattr(elements, '__iter__') or isinstance(elements, str):
                elements = [elements]
            if names is None:
                names = [f"coll_{i}" for i, _ in enumerate(elements)]
            else:
                if not hasattr(names, '__iter__') or isinstance(names, str):
                    names = [names]
                if len(names) != len(elements):
                    raise ValueError("Length of `elements` and `names` doesn't match.")
        else:
            if elements is not None:
                raise ValueError("Cannot provide both `line` and `elements`.")
            if names is None:
                elements, names = line.get_elements_of_type((_K2Collimator, _K2Crystal))
            else:
                if not hasattr(names, '__iter__') or isinstance(names, str):
                    names = [names]
                elements = [line[name] for name in names]
        for el in elements:
            if not isinstance(el, (_K2Collimator, _K2Crystal)):
                raise ValueError(f"Element {el} is not a _K2Collimator or a"
                                + "K2Crystal.")
        elements = [el for el in elements if el.active]
        names    = [name for el, name in zip(elements,names) if el.active]
        elements = [el for el in elements if el.jaw is not None or el.gap is not None]
        names    = [name for el, name in zip(elements,names)
                    if el.jaw is not None or el.gap is not None]

        # Check optics and emittances
        for el, name in zip(elements, names):
            if not hasattr(el, 'emittance') or el.emittance is None:
                raise ValueError(f"Collimator {name} has no emittance assigned.")
            if not hasattr(el, 'optics') or el.optics is None:
                if el.gap is not None:
                    raise ValueError(f"Collimator {name} has no optics assigned but gap is set.")
                el._optics = {}
                el._optics[el.align] = {}
                el._optics[el.align]['alfx'] = [1]
                el._optics[el.align]['alfy'] = [1]
                # Guesstimate betx and bety such that the gap is ~5 sigma
                beta_gamma_rel = particle_ref._xobject.gamma0[0] * particle_ref._xobject.beta0[0]
                if el.__class__ == _K2Crystal:
                    off2 = (el.jaw/5)**2
                    centre_x = 0
                    centre_y = 0
                    betx = 10000 if abs(el._cos_z) < 1.e-12 else off2*beta_gamma_rel/el.nemitt_x /(el._cos_z**2)
                    bety = 10000 if abs(el._sin_z) < 1.e-12 else off2*beta_gamma_rel/el.nemitt_y /(el._sin_z**2)
                else:
                    off2 = ((el.jaw_L - el.jaw_R)/10)**2
                    centre_x = (el.jaw_L + el.jaw_R)/2 * el._cos_zL
                    centre_y = (el.jaw_L + el.jaw_R)/2 * el._sin_zL
                    betx = 10000 if abs(el._cos_zL) < 1.e-12 else off2*beta_gamma_rel/el.nemitt_x /(el._cos_zL**2)
                    bety = 10000 if abs(el._sin_zL) < 1.e-12 else off2*beta_gamma_rel/el.nemitt_y /(el._sin_zL**2)
                el._optics[el.align]['betx'] = [betx]  # TODO: understand
                el._optics[el.align]['bety'] = [bety]
                el._optics[el.align]['x']  = [centre_x]
                el._optics[el.align]['y']  = [centre_y]
                el._optics[el.align]['px'] = [0]
                el._optics[el.align]['py'] = [0]
                el._optics['beta_gamma_rel'] = beta_gamma_rel
        nemitt_x = np.unique([el.nemitt_x for el in elements])
        if len(nemitt_x) > 1:
            raise ValueError("Not all collimators have the same horizontal "
                            + "emittance. This is not supported.")
        nemitt_x = nemitt_x[0]
        nemitt_y = np.unique([el.nemitt_y for el in elements])
        if len(nemitt_y) > 1:
            raise ValueError("Not all collimators have the same vertical "
                            + "emittance. This is not supported.")
        nemitt_y = nemitt_y[0]

        # Write SixTrack collimator file
        for i, el in enumerate(elements):
            el._k2_id = i + 1  # FORTRAN is 1-indexed
        this._file = create_dat_file(elements=elements, names=names, file=cwd / _FILE_NAME)
        assert this._file.exists()
        num_coll = np.int32(len(names))
        e0       = particle_ref.energy0[0]
        p0       = particle_ref.p0c[0]
        m0       = particle_ref.mass0
        beta0    = particle_ref.beta0[0]
        gamma0   = particle_ref.gamma0[0]
        alfx     = np.array([el.optics[el.align]['alfx'][0] for el in elements])
        alfy     = np.array([el.optics[el.align]['alfy'][0] for el in elements])
        betx     = np.array([el.optics[el.align]['betx'][0] for el in elements])
        bety     = np.array([el.optics[el.align]['bety'][0] for el in elements])
        x        = np.array([el.optics[el.align]['x'][0] for el in elements])
        y        = np.array([el.optics[el.align]['y'][0] for el in elements])
        px       = np.array([el.optics[el.align]['px'][0] for el in elements])
        py       = np.array([el.optics[el.align]['py'][0] for el in elements])

        pyk2_init(n_alloc=this._capacity, colldb_input_fname=this._file.name, \
                  random_generator_seed=this._seed, num_coll=num_coll, alphax=alfx, \
                  alphay=alfy, betax=betx, betay=bety, orbx=x, orby=y, orbxp=px, \
                  orbyp=py, nemitt_x=nemitt_x, nemitt_y=nemitt_y, e_ref=e0, p_ref=p0, \
                  m_ref=m0, beta_ref=beta0, gamma_ref=gamma0)


    @classmethod
    def stop(cls):
        cls()
        this = cls.instance
        if this._old_cwd is not None:
            os.chdir(this._old_cwd)
            this._old_cwd = None
        if this._seed_set_manually:
            print("Warning: seed has not been reset. Do this manually before starting "
                  "the server again to avoid repeating the same random numbers (unless "
                  "this is intentional).")
        else:
            this._seed = None
        if this._file is not None:
            if this._file.exists():
                shutil.rmtree(this._cwd)
            this._file = None
            this._cwd = None


    @classmethod
    def is_running(cls):
        # TODO: this can probably be done better; the file might be a leftover from a crash
        cls()
        this = cls.instance
        return this._file is not None and this._file.exists()


    @property
    def particle_ref(self):
        return self._particle_ref

    @particle_ref.setter
    def particle_ref(self, val):
        if val.mass0 != xt.Particles.reference_from_pdg_id('proton').mass0:
            raise ValueError("K2 only supports proton beams!")
        self._particle_ref = val

    def _has_particle_ref(self):
        initial = xp.Particles().to_dict()
        current = self.particle_ref.to_dict()
        return not xt.line._dicts_equal(initial, current)


    @property
    def capacity(self):
        return self._capacity

    @capacity.setter
    def capacity(self, val):
        self._capacity = np.int32(val)


    @property
    def seed(self):
        return self._seed

    @seed.setter
    def seed(self, val):
        if val is None:
            self._seed_set_manually = False
        else:
            self._seed_set_manually = True
            new_val = np.int32(abs(np.int32(val)))
            if new_val != int(val):
                print(f"Warning: overflow for seed {val}. Using {new_val}.")
        self._seed = val
