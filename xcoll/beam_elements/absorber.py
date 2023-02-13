from .base_collimator import BaseCollimator
from ..general import _pkg_root

class BlackAbsorber(BaseCollimator):
    _skip_in_to_dict       = BaseCollimator._skip_in_to_dict
    _store_in_to_dict      = BaseCollimator._store_in_to_dict
    _internal_record_class = BaseCollimator._internal_record_class

    skip_in_loss_location_refinement = True

    _extra_c_sources = [
        _pkg_root.joinpath('beam_elements/collimators_src/absorber.h')
    ]

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

