import json
import numpy as np
from pathlib import Path
import xtrack as xt
import xcoll as xc


def test_create_settings():
    collset = xc.collimator_settings.CollimatorSettings()
    collset.to_dict()
    