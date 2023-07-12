import numpy as np
import random
import string

import xobjects as xo
from xtrack import BeamInteraction, Drift

from .base_collimator import BaseCollimator
from ..scattering_routines.geant4 import Geant4Engine


def _new_id64(len=16):
    chars = string.ascii_letters + string.digits + '+/'
    return ''.join(random.choice(chars) for i in range(len))


class Geant4Collimator(BeamInteraction, BaseCollimator ):
    _xofields = BaseCollimator._xofields | {
        'collimator_id': xo.String,
        'material':      xo.String,
        '_tracking':     xo.Int8
    }

    isthick = True
    iscollective = True
    behaves_like_drift = True
    skip_in_loss_location_refinement = True

    _store_in_to_dict      = [ *BaseCollimator._store_in_to_dict, 'material' ]
    _internal_record_class = BaseCollimator._internal_record_class

    _depends_on = [Drift, BeamInteraction, BaseCollimator, Geant4Engine]


    def __init__(self, **kwargs):
        BaseCollimator.__init__(self, **kwargs)
        # We skip the BeamInteraction initialisation (to avoid issues with setting the length)
        # and manually set the necessary flags
        if '_xobject' not in kwargs:
            kwargs.setdefault('collimator_id', _new_id64())
            kwargs.setdefault('_tracking', True)
            if kwargs.get('material') is None:
                raise ValueError("Need to provide a material to the collimator!")

            self.material = kwargs['material']
            self.collimator_id = kwargs['collimator_id']
            self._tracking = kwargs['_tracking']
            
        self.interaction_process = Geant4Engine().register_collimator(self)
        self.equivalent_drift = Drift(length=self.length)
        # TODO: is there a smarter way to hangle the drift tracking when inactive
        # TODO: should the inactive length before/after be handled here?

    # def __del__(self):
    #     # TODO: with a custom destructor, should all base class destructors be called explicitly?
    #     Geant4Engine().deregister_collimator(self)

    def track(self, part):
        if self.active and self._tracking and Geant4Engine().connected:
            Geant4Engine().enable_scattering()
            BeamInteraction.track(self, part)
        else:
            self.equivalent_drift.track(part)
    
    # NEED TO HAVE UPDATE RULES!!! WHEN BeamElement.jaw_F_L etc are updated, this needs to be communicated to Geant4...
    # jaws are currently not set