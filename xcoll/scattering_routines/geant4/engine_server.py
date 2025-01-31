import subprocess as _subprocess
import rpyc as _rpyc
import collimasim as cs

class BDSIMServer:
    def __init__(self):
        try:
            import collimasim as cs
        except ImportError as e:
            raise ImportError("Failed to import collimasim. Cannot connect to BDSIM.")
        if not hasattr(self, 'g4link'):
            self.g4link = None
        self._initialised = True

    def close(self):
        pass

    def XtrackInterface(self,bdsimConfigFile=None,referencePdgId=None,referenceEk=None,
                        relativeEnergyCut=None,seed=None,batchMode=None):
        self.g4link = cs.XtrackInterface(bdsimConfigFile=bdsimConfigFile,
                                         referencePdgId=referencePdgId,
                                         referenceEk=referenceEk, # BDSIM expects GeV
                                         relativeEnergyCut=relativeEnergyCut,
                                         seed=seed, batchMode=batchMode)

    def addCollimator(self,geant4_id,material,length,apertureLeft=None,apertureRight=None,
                      rotation=None,xOffset=0,yOffset=0,side=None,jawTiltLeft=None,jawTiltRight=None):
        self.g4link.addCollimator(geant4_id, material, length,
                                  apertureLeft=apertureLeft,
                                  apertureRight=apertureRight,   # TODO: is this correct?
                                  rotation=rotation,
                                  xOffset=xOffset, yOffset=yOffset, side=side,
                                  jawTiltLeft=jawTiltLeft, jawTiltRight=jawTiltRight)

    def addParticles(self, particles_x, particles_y, particles_px, particles_py,particles_zeta, delta_temp,
                     particles_chi,particles_charge_ratio, particles_s,particles_pdg_id,particles_particle_id,
                     particles_state,particles_at_element, particles_at_turn):
        coordinates = [particles_x, particles_y, particles_px, particles_py,
                   particles_zeta, delta_temp, particles_chi,
                   particles_charge_ratio, particles_s,
                   particles_pdg_id,particles_particle_id, particles_state,
                   particles_at_element, particles_at_turn]
        self.g4link.addParticles(coordinates)

    def collimateReturn(self, particles_x, particles_y, particles_px, particles_py,particles_zeta, delta_temp,
                        particles_chi,particles_charge_ratio, particles_s,particles_pdg_id,particles_particle_id,
                        particles_state,particles_at_element, particles_at_turn, secondaries_x):
        coordinates = [particles_x, particles_y, particles_px, particles_py,
                       particles_zeta, delta_temp, particles_chi,
                       particles_charge_ratio, particles_s,
                       particles_pdg_id,particles_particle_id, particles_state,
                       particles_at_element, particles_at_turn]
        products = self.g4link.collimateReturn(coordinates)
        for i,x in enumerate(products['x']):
            secondaries_x[i] = x
        return products

    def selectCollimator(self,geant4_id):
        self.g4link.selectCollimator(geant4_id)

    def collimate(self):
        self.g4link.collimate()

    def clearData(self):
        pass