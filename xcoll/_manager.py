# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #


class CollimatorManager:

    def __init__(self, **kwargs):
        raise ValueError("`CollimatorManager` is deprecated. Use `CollimatorDatabase` instead.")

    @classmethod
    def from_yaml(cls, file, **kwargs):
        raise ValueError("`CollimatorManager` is deprecated. Use `CollimatorDatabase` instead.")

    @classmethod
    def from_json(cls, file, **kwargs):
        raise ValueError("`CollimatorManager` is deprecated. Use `CollimatorDatabase` instead.")

    @classmethod
    def from_dict(cls, file, **kwargs):
        raise ValueError("`CollimatorManager` is deprecated. Use `CollimatorDatabase` instead.")
