import xobjects as xo
import xpart as xp
from xtrack.base_element import _handle_per_particle_blocks
from ...general import _pkg_root

from functools import partial


# TODO: auto-initialise Rutherford parameters at construction using the C code (make a kernel?)
#       such that it can be set at creation of the collimator, and when a material is updated
class EverestRandom(xo.HybridClass):
    _xofields = {
        'rutherford_iterations':    xo.Int8,
        'rutherford_lower_val':     xo.Float64,
        'rutherford_upper_val':     xo.Float64,
        'rutherford_A':             xo.Float64,
        'rutherford_B':             xo.Float64
    }

    _depends_on = [xp.Particles]

    _extra_c_sources = [
        xp.general._pkg_root.joinpath('random_number_generator/rng_src/base_rng.h'),
        xp.general._pkg_root.joinpath('random_number_generator/rng_src/local_particle_rng.h'),
        _pkg_root.joinpath('scattering_routines','everest','exponential_integral_Ei.h'),
        _pkg_root.joinpath('scattering_routines','everest','random.h')
    ]

    _kernels = {
        'rsample': xo.Kernel(
            c_name='EverestRandomData_get_random',
            args=[
                xo.Arg(xp.Particles, name='particles'),
                xo.Arg(xo.Int32, name='n')
            ],
            ret=xo.Arg(xo.Float64[:]),
            n_threads='n_part')
        }

    def __init__(self, **kwargs):
        kwargs.setdefault('rutherford_iterations', 7)
        kwargs.setdefault('rutherford_lower_val', 0.0009982)
        kwargs.setdefault('rutherford_upper_val', 0.0009982)
        kwargs.setdefault('rutherford_A', 0.)
        kwargs.setdefault('rutherford_B', 0.)
        super().__init__(**kwargs)

    def sample(self, particles, n=1):
        self.compile_kernels(only_if_needed=True, save_source_as="randomtest.c",
                      apply_to_source=[partial(_handle_per_particle_blocks, local_particle_src=xp.gen_local_particle_api())])
        context = self._buffer.context
        return context.kernels.rsample(particles=particles, n=n)
