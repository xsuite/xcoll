# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2025                  #
# ######################################### #

import io
import numpy as np


class BDSIMServer:

    def __init__(self):
        self.g4link = None


    def XtrackInterface(self, bdsimConfigFile=None, referencePdgId=None, referenceEk=None,
                        relativeEnergyCut=None, seed=None, batchMode=None):
        from g4interface import XtrackInterface
        self.g4link = XtrackInterface(bdsimConfigFile=bdsimConfigFile,
                                      referencePdgId=referencePdgId,
                                      referenceEk=referenceEk,
                                      relativeEnergyCut=relativeEnergyCut,
                                      seed=seed, batchMode=batchMode)


    def addCollimator(self, geant4_id, material, length, apertureLeft=None, apertureRight=None,
                      rotation=None, xOffset=0, yOffset=0, side=None, jawTiltLeft=None,
                      jawTiltRight=None, isACrystal=False):
        self.g4link.addCollimator(geant4_id, material, length,
                                  apertureLeft=apertureLeft,
                                  apertureRight=apertureRight,
                                  rotation=rotation,
                                  xOffset=xOffset, yOffset=yOffset, side=side,
                                  jawTiltLeft=jawTiltLeft, jawTiltRight=jawTiltRight,
                                  isACrystal=isACrystal)


    def add_particles_and_collimate_return(self, blob, geant4_id, num_sent, output_size):
        buf = io.BytesIO(blob)
        loaded = np.load(buf)
        coords = [np.array(loaded[key]) for key in loaded.files]
        self.g4link.addParticles(coords)
        self.g4link.selectCollimator(geant4_id)
        self.g4link.collimate()
        products = self.g4link.collimateReturn(num_sent, output_size)
        # Pack results as a dict for serialization
        result_dict = {
            'x': products['x'],
            'xp': products['xp'],
            'y': products['y'],
            'yp': products['yp'],
            'zeta': products['zeta'],
            'p': products['p'],
            'q': products['q'],
            'm': products['m'],
            'weight': products['weight'],
            'pdg_id': products['pdg_id'],
            'parent_particle_id': products['parent_particle_id'],
            'state': products['state'],
        }

        # Serialize to binary blob
        out_buf = io.BytesIO()
        np.savez(out_buf, **result_dict)
        return out_buf.getvalue()

    def clearData(self):
        self.g4link.clearData()
