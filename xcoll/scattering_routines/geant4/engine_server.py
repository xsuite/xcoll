# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2025                  #
# ######################################### #

import subprocess as _subprocess
import rpyc as _rpyc
import time
import numpy as np
import io

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
                      rotation=None,xOffset=0,yOffset=0,side=None,jawTiltLeft=None,jawTiltRight=None, isACrystal=False):
        self.g4link.addCollimator(geant4_id, material, length,
                                  apertureLeft=apertureLeft,
                                  apertureRight=apertureRight,   # TODO: is this correct?
                                  rotation=rotation,
                                  xOffset=xOffset, yOffset=yOffset, side=side,
                                  jawTiltLeft=jawTiltLeft, jawTiltRight=jawTiltRight, isACrystal=isACrystal)


    def addParticles(self, particles_x, particles_y, particles_px, particles_py,particles_zeta, delta_temp,
                     particles_chi,particles_charge_ratio, particles_s,particles_pdg_id,particles_particle_id,
                     particles_state,particles_at_element, particles_at_turn):
        coordinates = [particles_x, particles_y, particles_px, particles_py,
                   particles_zeta, delta_temp, particles_chi,
                   particles_charge_ratio, particles_s,
                   particles_pdg_id,particles_particle_id, particles_state,
                   particles_at_element, particles_at_turn]
        self.g4link.addParticles(coordinates)


    def add_particles_and_collimate_return(self, blob, geant4_id):
        buf = io.BytesIO(blob)
        loaded = np.load(buf)
        coords = [np.array(loaded[key]) for key in loaded.files]
        self.g4link.addParticles(coords)
        self.g4link.selectCollimator(geant4_id)
        self.g4link.collimate()
        products = self.g4link.collimateReturn(coords)
        # Pack results as a dict for serialization
        result_dict = {
            'x': products['x'],
            'y': products['y'],
            'px': products['px'],
            'py': products['py'],
            'zeta': products['zeta'],
            'delta': products['delta'],
            'charge_ratio': products['charge_ratio'],
            's': products['s'],
            'pdg_id': products['pdg_id'],
            'parent_particle_id': products['parent_particle_id'],
            'at_element': products['at_element'],
            'at_turn': products['at_turn'],
            'mass_ratio': products['mass_ratio'],
            'state': products['state'],
        }

        # Serialize to binary blob
        out_buf = io.BytesIO()
        np.savez(out_buf, **result_dict)
        return out_buf.getvalue()


    def receive_serialized_particles(self, blob):
        buf = io.BytesIO(blob)
        loaded = np.load(buf)
        coords = [np.array(loaded[key]) for key in loaded.files]
        self.g4link.addParticles(coords)


    def addParticles2(self, coordinates):
        import time
        t0 = time.time()
        t1 = time.time()

        coords = coordinates  # already native
        self.g4link.addParticles([
            coords['x'],
            coords['y'],
            coords['px'],
            coords['py'],
            coords['zeta'],
            coords['delta'],
            coords['chi'],
            coords['charge_ratio'],
            coords['s'],
            coords['pdg_id'],
            coords['particle_id'],
            coords['state'],
            coords['at_element'],
            coords['at_turn'],
        ])
        t2 = time.time()
        print(f"[server] proxy unwrap: {t1 - t0:.6f}s, g4link.addParticles: {t2 - t1:.6f}s")


    def collimateReturn(self, particles_x, particles_y, particles_px, particles_py,particles_zeta, delta_temp,
                        particles_chi,particles_charge_ratio, particles_s,particles_pdg_id,particles_particle_id,
                        particles_state,particles_at_element, particles_at_turn):
        coordinates = [particles_x, particles_y, particles_px, particles_py,
                       particles_zeta, delta_temp, particles_chi,
                       particles_charge_ratio, particles_s,
                       particles_pdg_id,particles_particle_id, particles_state,
                       particles_at_element, particles_at_turn]
        products = self.g4link.collimateReturn(coordinates)

        return products

    def collimateReturnSerialized(self, binary_blob):
        # Deserialize incoming data
        buf = io.BytesIO(binary_blob)
        loaded = np.load(buf)
        coords = [np.array(loaded[key]) for key in loaded.files]
        # Call the Geant4 collimation backend
        products = self.g4link.collimateReturn(coords)
        # Pack results as a dict for serialization
        result_dict = {
            'x': products['x'],
            'y': products['y'],
            'px': products['px'],
            'py': products['py'],
            'zeta': products['zeta'],
            'delta': products['delta'],
            'charge_ratio': products['charge_ratio'],
            's': products['s'],
            'pdg_id': products['pdg_id'],
            'parent_particle_id': products['parent_particle_id'],
            'at_element': products['at_element'],
            'at_turn': products['at_turn'],
            'mass_ratio': products['mass_ratio'],
            'state': products['state'],
        }

        # Serialize to binary blob
        out_buf = io.BytesIO()
        np.savez(out_buf, **result_dict)
        return out_buf.getvalue()


        return products

        return products['x'], products['y'], products['px'], products['py'], \
            products['zeta'], products['delta'], products['charge_ratio'], \
            products['s'], products['pdg_id'], products['parent_particle_id'], \
            products['at_element'], products['at_turn'], products['mass_ratio'], products['state']

        return list(products['x']), list(products['y']), list(products['px']), list(products['py']), \
            list(products['zeta']), list(products['delta']), list(products['charge_ratio']), \
            list(products['s']), np.array(products['pdg_id']), list(products['parent_particle_id']), \
            list(products['at_element']), list(products['at_turn']), list(products['mass_ratio']), list(products['state'])




    def selectCollimator(self,geant4_id):
        self.g4link.selectCollimator(geant4_id)

    def collimate(self):
        self.g4link.collimate()

    def clearData(self):
        self.g4link.clearData()

