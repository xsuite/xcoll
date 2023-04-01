from .beam_elements import K2Collimator, K2Crystal, PyEverestCollimator, PyEverestCrystal
from .scattering_routines.k2 import K2Engine


def install_k2_collimators(self, coll_manager, names=None, *, max_part=50000, seed=None, verbose=False):
    # Check for the existence of a K2Engine; warn if settings are different
    k2engine = K2Engine(_capacity=max_part, random_generator_seed=seed)
    if k2engine._capacity != max_part:
        print(f"Warning: K2 already initiated with a maximum allocation of {k2engine._capacity} particles.\n"
              + f"Ignoring the requested max_part={max_part}.")
    if seed is not None and k2engine.random_generator_seed != seed:
        print(f"Warning: K2 already initiated with seed {k2engine.random_generator_seed}.\n"
              + f"Ignoring the requested seed={seed}.")

    if names is None:
        names = self.collimator_names

    # Do the installation
    def install_func(thiscoll, name):
        return K2Collimator(
                inactive_front=thiscoll['inactive_front'],
                inactive_back=thiscoll['inactive_back'],
                active_length=thiscoll['active_length'],
                angle=thiscoll['angle'],
                material=SixTrack_to_xcoll[thiscoll['material']][0],
                active=False,
                _buffer=self._buffer
               )
    coll_manager._install_collimators(names, install_func=install_func, verbose=verbose, support_legacy_elements=True)
    return k2engine


def install_pyeverest_collimators(self, coll_manager, names=None, *, verbose=False, random_seed=None):
    from .scattering_routines.pyeverest import set_random_seed
    set_random_seed(random_seed)

    if names is None:
        names = self.collimator_names

    # Do the installation
    def install_func(thiscoll, name):
        return PyEverestCollimator(
                inactive_front=thiscoll['inactive_front'],
                inactive_back=thiscoll['inactive_back'],
                active_length=thiscoll['active_length'],
                angle=thiscoll['angle'],
                material=SixTrack_to_xcoll[thiscoll['material']][0],
                active=False,
                _buffer=self._buffer
               )
    coll_manager._install_collimators(names, install_func=install_func, verbose=verbose, support_legacy_elements=True)

