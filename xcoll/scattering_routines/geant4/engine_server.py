import subprocess as _subprocess
import rpyc as _rpyc
import collimasim as cs

class BDSIMServer:
    def __init__(self):
        try:
            import collimasim as cs
        except ImportError as e:
            raise ImportError("Failed to import collimasim. Cannot connect to BDSIM.")
        self.g4link = None

    def close(self):
        pass

    def XtrackInterface(self,bdsimConfigFile=None,referencePdgId=None,referenceEk=None,
                        relativeEnergyCut=None,seed=None,batchMode=None):
        print(f"{bdsimConfigFile=}")
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
    def addParticles2(self, coordinates):
        self.g4link.addParticles(list(coordinates))

    def collimateReturn(self, particles_x, particles_y, particles_px, particles_py,particles_zeta, delta_temp,
                        particles_chi,particles_charge_ratio, particles_s,particles_pdg_id,particles_particle_id,
                        particles_state,particles_at_element, particles_at_turn, secondaries_x,
                        secondaries_y,secondaries_px,secondaries_py,secondaries_zeta,
                        secondaries_delta,secondaries_charge_ratio,secondaries_s,
                        secondaries_pdg_id,secondaries_parent_particle_id,secondaries_at_element,
                        secondaries_at_turn,secondaries_mass_ratio,secondaries_state):
        coordinates = [particles_x, particles_y, particles_px, particles_py,
                       particles_zeta, delta_temp, particles_chi,
                       particles_charge_ratio, particles_s,
                       particles_pdg_id,particles_particle_id, particles_state,
                       particles_at_element, particles_at_turn]
        products = self.g4link.collimateReturn(coordinates)
        for i,x in enumerate(products['x']):
            secondaries_x[i] = x
        for i,x in enumerate(products['y']):
            secondaries_y[i] = x
        for i,x in enumerate(products['px']):
            secondaries_px[i] = x
        for i,x in enumerate(products['py']):
            secondaries_py[i] = x
        for i,x in enumerate(products['zeta']):
            secondaries_zeta[i] = x
        for i,x in enumerate(products['delta']):
            secondaries_delta[i] = x
        for i,x in enumerate(products['charge_ratio']):
            secondaries_charge_ratio[i] = x
        for i,x in enumerate(products['s']):
            secondaries_s[i] = x
        for i,x in enumerate(products['pdg_id']):
            secondaries_pdg_id[i] = x
        for i,x in enumerate(products['parent_particle_id']):
            secondaries_parent_particle_id[i] = x
        for i,x in enumerate(products['at_element']):
            secondaries_at_element[i] = x
        for i,x in enumerate(products['at_turn']):
            secondaries_at_turn[i] = x
        for i,x in enumerate(products['mass_ratio']):
            secondaries_mass_ratio[i] = x
        for i,x in enumerate(products['state']):
            secondaries_state[i] = x
        #return products

    def selectCollimator(self,geant4_id):
        self.g4link.selectCollimator(geant4_id)

    def collimate(self):
        self.g4link.collimate()

    def clearData(self):
        pass