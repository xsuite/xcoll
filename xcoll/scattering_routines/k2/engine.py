import numpy as np

import xobjects as xo


class K2Engine(xo.HybridClass):

    _xofields = {
        'n_alloc':                xo.Int64,
        'random_generator_seed':  xo.Int64,
#         'random_freeze_state':    xo.Int64,  # to be implemented; number of randoms already sampled, such that this can be taken up again later
#         'collimators':            K2Collimator[:],  # Does not work, need a pointer of sorts
    }

    def __init__(self, **kwargs):
        kwargs.setdefault('_n_alloc', 50000)
        kwargs.setdefault('random_generator_seed', np.random.randint(1, 10000000))
#         kwargs.setdefault('random_freeze_state', -1)
#         kwargs.setdefault('collimators', [])
        super().__init__(**kwargs)
        try:
            from .pyk2f import pyk2_init
        except ImportError:
            print("Warning: Failed importing pyK2 (did you compile?). " \
                  + "K2collimators will be installed but are not trackable.")
        else:
            pyk2_init(n_alloc=self.n_alloc, random_generator_seed=self.random_generator_seed)


