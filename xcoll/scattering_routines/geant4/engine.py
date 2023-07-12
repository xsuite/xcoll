import numpy as np
import xobjects as xo

from .interaction import Geant4Interaction

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
        if not kwargs: # If trying to use before intialising
            super().__init__(random_generator_seed = -1,
                             reference_pdg_id = -1,   
                             reference_kinetic_energy = -1,
                             relative_energy_cut = -1,
                             bdsim_config_file = -1)
            return
        
        # Allow seed to be set to None to get default:
        kwargs.setdefault('random_generator_seed', None)
        if kwargs['random_generator_seed'] is None:
            kwargs['random_generator_seed'] = np.random.randint(1, 10000000)

        super().__init__(**kwargs)
        try:
            import collimasim as cs
        except ImportError:
            global _geant4_warning_given
            if not _geant4_warning_given:
                print("Warning: Failed to import collimasim. " \
                      + "Geant4Collimators will be installed but are not trackable.")
                _geant4_warning_given = True
        else:
            unit_GeV = 1e9 # GeV to eV
            self.g4link = cs.XtrackInterface(bdsimConfigFile=self.bdsim_config_file,
                                             referencePdgId=self.reference_pdg_id,
                                             referenceEk=self.reference_kinetic_energy / unit_GeV, # BDSIM expects GeV
                                             relativeEnergyCut=self.relative_energy_cut,
                                             seed=self.random_generator_seed, batchMode=batch_mode)


        self._initialised = True

        self.registered_collimators = {}
        self._geometry_constructed = False
        self._built_collimators = {}


    def register_collimator(self, collimator):
        # Register the collimators by reference
        self.registered_collimators[collimator.collimator_id] = collimator
        interaction_process = Geant4Interaction(element_id=collimator.collimator_id)

        return interaction_process

    def deregister_collimator(self, collimator):
        if collimator.collimator_id in self.registered_collimators:
            del self.registered_collimators[collimator.collimator_id]

    def enable_scattering(self):
        if not self._geometry_constructed:
            self._construct_geometry()
        else:
            for collimator_id, collimator in self.registered_collimators.items():
                if not collimator.active:
                    continue
                parameters_for_geometry = self._extract_parameters(collimator)
                if not self._built_collimators[collimator_id] == parameters_for_geometry:
                    # TODO: flush and re-start the BDSIM instance in this case (hard, as not re-entry safe)
                    raise Exception('Collimator settings changed after the Geant4 geometry'
                                    'has been build. This is not currently supported.')

    def _extract_parameters(self, collimator):
        try:
            iter(collimator.angle)
        except TypeError:
            angle = float(collimator.angle)
        else:
            raise Exception('The Geant4 scattering engine does not'
                            'support unequal jaw rotation angles')

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
                                    collimator.jaw_LU, collimator.jaw_LD, 
                                    collimator.jaw_RU, collimator.jaw_RD,
                                    angle, 
                                    collimator.reference_center[0], 
                                    collimator.reference_center[1],
                                    collimator._side, collimator.active)
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


    def add_collimator(self, element_id, material, length, jaw_LU, jaw_LD, jaw_RU, jaw_RD, 
                       angle, centre_x, centre_y, side, active):
        
        centre_l = jaw_LU + (jaw_LD - jaw_LU) / 2
        centre_r = jaw_RU + (jaw_RD - jaw_RU) / 2

        halfgap = (centre_l - centre_r) / 2
        offs = (centre_l + centre_r) / 2

        tilt_l = np.arctan((jaw_LD - jaw_LU) / length)
        tilt_r = np.arctan((jaw_RD - jaw_RU) / length)

        
        self.g4link.addCollimator(element_id, material, length, 
                                  apertureLeft=halfgap + offs, 
                                  apertureRight=halfgap + offs,
                                  rotation=np.deg2rad(angle), 
                                  xOffset=centre_x, yOffset=centre_y,
                                  jawTiltLeft=tilt_l, jawTiltRight=tilt_r, 
                                  side=side)
        
