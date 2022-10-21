import numpy as np
# from .k2.materials import Material, Crystal
from .collimators import BaseCollimator

from xtrack import BeamInteraction

class Geant4Collimator(BaseCollimator, BeamInteraction):

    iscollective = True

    def __init__(self, **kwargs):
        super().__init__(self, inactive_front=10)
        # We skip the BeamInteraction initialisation (to avoid issues with setting the length)
        # and manually set the necessary flags
        self.name = 'TEST'
        self.interaction_process = None
        self.isthick = True


    def track(self, part):
        BeamInteraction.track(self, part)
    