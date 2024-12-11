# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2023.                 #
# ######################################### #

import numpy as np

import xobjects as xo
import xtrack as xt
import xpart as xp

class Geant4Engine(xo.HybridClass):

    _xofields = {
        'particle_ref':        xp.Particles,
        'seed':                xo.Int64,
        'relative_energy_cut': xo.Float64,
        'bdsim_config_file':   xo.String
#         'random_freeze_state':    xo.Int64,  # to be implemented; number of randoms already sampled, such that this can be taken up again later
    }

    # The engine is a singleton
    def __new__(cls, *args, **kwargs):
        if not hasattr(cls, 'instance'):
            cls.instance = super().__new__(cls)
            cls.instance._initialised = False
        return cls.instance

    def __init__(self, **kwargs):
        # Apply kwargs
        if self._initialised:
            for kk, vv in kwargs.items():
                if not hasattr(self, kk):
                    raise ValueError(f"Invalid attribute {kk} for FlukaEngine!")
                setattr(self, kk, vv)
            return
        if '_xobject' not in kwargs:
            kwargs.setdefault('particle_ref', xp.Particles())
            kwargs.setdefault('seed', -1)
            kwargs.setdefault('relative_energy_cut', -1)
            kwargs.setdefault('bdsim_config_file', ''.ljust(256))  # Limit to pathnames of 256 characters
        super().__init__(**kwargs)
        if not hasattr(self, 'g4link'):
                    self.g4link = None
        self._initialised = True

    def __del__(self, *args, **kwargs):
        self.stop()


    @classmethod
    def start(cls, *, bdsim_config_file=None, line=None, elements=None, cwd=None,
              relative_energy_cut=0.15, seed=None, batch_mode=True,
              particle_ref=None, p0c=None, **kwargs):
        from ...beam_elements.geant4 import Geant4Collimator

        cls(**kwargs)
        this = cls.instance
        if this.is_running():
            print("Geant4Engine already running.", flush=True)
            return

        if bdsim_config_file is None:
            raise NotImplementedError
        # if cwd is not None:
        #     cwd = Path(cwd).expanduser().resolve()
        #     cwd.mkdir(parents=True, exist_ok=True)
        #     this._old_cwd = Path.cwd()
        #     os.chdir(cwd)
        # else:
        #     cwd = Path.cwd()
        # this._cwd = cwd

        this.bdsim_config_file = bdsim_config_file.as_posix()
        cls.set_particle_ref(particle_ref=particle_ref, line=line, p0c=p0c)
        Ekin = this.particle_ref.energy0 - this.particle_ref.mass0
        pdg_id = this.particle_ref.pdg_id
        # TODO: the original Geant4 coupling had -11 for electrons and 11 for positrons. This is wrong. Was this an error in the original code, or is this wrong in BDSIM?
        if abs(pdg_id) == 11:
            pdg_id = -pdg_id

        if seed is None:
            if this._seed is None:
                this._seed = np.int64(np.random.randint(1, int(2^63-1)))
        else:
            this.seed = seed
            # Setting a seed here does not count as manually setting it
            this._seed_set_manually = False
        print(f"Using seed {this.seed}.")
        this.relative_energy_cut = relative_energy_cut

        try:
            import collimasim as cs
        except ImportError as e:
            raise ImportError("Failed to import collimasim. Cannot connect to BDSIM.")
        else:
            this.g4link = cs.XtrackInterface(bdsimConfigFile=bdsim_config_file,
                                             referencePdgId=pdg_id,
                                             referenceEk=Ekin / 1e9, # BDSIM expects GeV
                                             relativeEnergyCut=this.relative_energy_cut,
                                             seed=this.seed, batchMode=batch_mode)

        if line is None:
            if elements is None:
                raise ValueError("Need to provide either `line` or `elements`.")
        elif elements is None:
            elements, _ = line.get_elements_of_type(Geant4Collimator)
        if not hasattr(elements, '__iter__') or isinstance(elements, str):
            elements = [elements]
        elements = [el for el in elements if el.jaw is not None and el.active]
        for el in elements:
            side = 2 if el._side == -1 else el._side
            jaw_L = 0.1 if el.jaw_L is None else el.jaw_L
            jaw_R = -0.1 if el.jaw_R is None else el.jaw_R
            tilt_L = 0.0 if el.tilt_L is None else el.tilt_L
            tilt_R = 0.0 if el.tilt_R is None else el.tilt_R
            this.g4link.addCollimator(el.geant4_id, el.material, el.length,
                                      apertureLeft=jaw_L,
                                      apertureRight=-jaw_R,   # TODO: is this correct?
                                      rotation=np.deg2rad(el.angle),
                                      xOffset=0, yOffset=0, side=side,
                                      jawTiltLeft=tilt_L, jawTiltRight=tilt_R)
        print('Geant4Engine initialised')


    @classmethod
    def stop(cls, clean=False, **kwargs):
        cls(**kwargs)
        this = cls.instance
        this.g4link = None


    @classmethod
    def is_running(cls, **kwargs):
        cls(**kwargs)
        this = cls.instance
        return this.g4link is not None


    @property
    def seed(self):
        return self._seed

    @seed.setter
    def seed(self, val):
        if val is None:
            self._seed_set_manually = False
        else:
            self._seed_set_manually = True
            new_val = np.int64(abs(np.int64(val)))
            if new_val != int(val):
                print(f"Warning: overflow for seed {val}. Using {new_val}.")
        self._seed = val


    @classmethod
    def set_particle_ref(cls, particle_ref=None, line=None, p0c=None, **kwargs):
        cls(**kwargs)
        this = cls.instance
        if this._has_particle_ref():
            print("Reference particle already set!")
            return

        if particle_ref is None:
            if line is None or line.particle_ref is None:
                raise ValueError("Line has no reference particle! "
                               + "Please provide one using `particle_ref`.")
            particle_ref = line.particle_ref
        elif isinstance(particle_ref, xp.Particles):
            if particle_ref._capacity > 1:
                raise ValueError("`particle_ref` has to be a single particle!")
        elif p0c is not None:
                particle_ref = xp.Particles.reference_from_pdg_id(particle_ref, p0c=p0c)
        else:
            raise ValueError("When providing `particle_ref`, it should be an "
                             "xp.Particles object or a PDG ID. In the latter case, "
                             "provide `p0c` as well.")

        if particle_ref.pdg_id == 0:
            particle_ref.pdg_id = xp.pdg.get_pdg_id_from_mass_charge(particle_ref.mass0, particle_ref.q0)
            # TODO: this should be updated in xpart: antiparticle not correctly recognised (missing positron and antimuon etc)
            q0, _, _, _ = xp.pdg.get_properties_from_pdg_id(particle_ref.pdg_id)
            if particle_ref.q0 == -q0:
                pdg_id = -pdg_id

        # TODO: test PDG ID consistent with mass and charge

        this.particle_ref = particle_ref
        if line is not None and line.particle_ref is not None \
        and not xt.line._dicts_equal(line.particle_ref.to_dict(), particle_ref.to_dict()):
            print("Warning: Found different reference particle in line! Overwritten.")
            line.particle_ref = particle_ref


    def _has_particle_ref(self):
        initial = xp.Particles().to_dict()
        current = self.particle_ref.to_dict()
        return not xt.line._dicts_equal(initial, current)

