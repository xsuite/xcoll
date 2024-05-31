# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2023.                 #
# ######################################### #

import numpy as np
import xobjects as xo

_geant4_warning_given = False

class Geant4Engine(xo.HybridClass):

    _xofields = {
        'random_generator_seed':    xo.Int64,
        'reference_pdg_id':         xo.Int64,
        'reference_kinetic_energy': xo.Float64,
        'relative_energy_cut':      xo.Float64,
        'bdsim_config_file':        xo.String
#         'random_freeze_state':    xo.Int64,  # to be implemented; number of randoms already sampled, such that this can be taken up again later
#         'collimators':            Geant4Collimator[:],  # Does not work, need a pointer of sorts
    }

    # The engine is a singleton
    def __new__(cls, *args, **kwargs):
        if not hasattr(cls, 'instance'):
            cls.instance = super().__new__(cls)
            cls.instance._initialised = False
        return cls.instance


    def __init__(self, batch_mode=True, **kwargs):
        if(self._initialised):
            return
        self._initialised = True

        if '_xobject' not in kwargs:
            # Allow seed to be set to None to get default:
            kwargs.setdefault('random_generator_seed', None)
            if kwargs['random_generator_seed'] is None:
                kwargs['random_generator_seed'] = np.random.randint(1, 10000000)
            kwargs.setdefault('reference_pdg_id', -1)
            kwargs.setdefault('reference_kinetic_energy', -1)
            kwargs.setdefault('relative_energy_cut', -1)
            kwargs.setdefault('bdsim_config_file', '')
            self.registered_collimators = {}
            self._geometry_constructed = False
            self._built_collimators = {}
        super().__init__(**kwargs)
        if '_xobject' not in kwargs:
            try:
                import collimasim as cs
            except ImportError:
                global _geant4_warning_given
                if not _geant4_warning_given:
                    print("Warning: Failed to import collimasim. " \
                        + "Geant4Collimators will be installed but are not trackable.")
                    _geant4_warning_given = True
                    self.g4link = None
            else:
                unit_GeV = 1e9 # GeV to eV
                self.g4link = cs.XtrackInterface(bdsimConfigFile=self.bdsim_config_file,
                                                 referencePdgId=self.reference_pdg_id,
                                                 referenceEk=self.reference_kinetic_energy / unit_GeV, # BDSIM expects GeV
                                                 relativeEnergyCut=self.relative_energy_cut,
                                                 seed=self.random_generator_seed, batchMode=batch_mode)
            print('Geant4 engine initialised')

    @property
    def connected(self):
        return self.g4link is not None

    def register_collimator(self, collimator):
        # Register the collimators by reference
        self.registered_collimators[collimator.collimator_id] = collimator

    def deregister_collimator(self, collimator):
        if collimator.collimator_id in self.registered_collimators:
            del self.registered_collimators[collimator.collimator_id]

    def assert_geometry(self):
        if not self._geometry_constructed:
            self._construct_geometry()
        else:
            for collimator_id, collimator in self.registered_collimators.items():
                if not collimator.active:
                    continue
                parameters_for_geometry = self._extract_parameters(collimator)
                if collimator_id not in self._built_collimators \
                or self._built_collimators[collimator_id] != parameters_for_geometry:
                    # TODO: flush and re-start the BDSIM instance in this case (hard, as not re-entry safe)
                    raise Exception("Collimator settings changed after the Geant4 geometry "
                                 + f"has been build for collimator with id {collimator_id}. "
                                 +  "This is not currently supported.")

    def _extract_parameters(self, collimator):
        try:
            iter(collimator.angle)
        except TypeError:
            angle = float(collimator.angle)
        else:
            raise Exception('The Geant4 scattering engine does not '
                          + 'support unequal jaw rotation angles')

        material_rename_map = {
            'c': 'AC150GPH',
            'cu': 'Cu',
            'mogr': 'MG6403Fc',
            'cucd': 'CUDIAM75',
            'iner': 'INERM180',
        }

        material_name = collimator.material
        bdsim_material = material_rename_map.get(material_name.lower(), material_name) 

        parameters_for_geometry = (collimator.collimator_id, 
                                    bdsim_material, 
                                    collimator.length, 
                                    collimator.jaw_L, collimator.jaw_R,
                                    collimator.tilt_L, collimator.tilt_R,
                                    angle, 0, 0, collimator._side, collimator.active)
        return parameters_for_geometry


    def _construct_geometry(self):
        for collimator_id, collimator in self.registered_collimators.items():
            # TODO: impement a more elegant way to enable and disable colliamtors
            if not collimator.active:
                continue

            parameters_for_geometry = self._extract_parameters(collimator)
            self._built_collimators[collimator_id] = parameters_for_geometry
            self.add_collimator(*parameters_for_geometry)

        self._geometry_constructed = True


    def add_collimator(self, element_id, material, length, jaw_L, jaw_R, tilt_L, tilt_R,
                       angle, centre_x, centre_y, side, active):
        if not self.connected:
            raise ValueError("Geant4Engine not linked to BDSIM! Cannot add collimator.")

        self.g4link.addCollimator(element_id, material, length, 
                                  apertureLeft=jaw_L, 
                                  apertureRight=-jaw_R,
                                  rotation=np.deg2rad(angle), 
                                  xOffset=centre_x, yOffset=centre_y,
                                  jawTiltLeft=tilt_L, jawTiltRight=tilt_R, 
                                  side=side)

