import xobjects as xo
import xpart as xp
from xtrack.base_element import _handle_per_particle_blocks
from ...general import _pkg_root

import numpy as np
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
        'usample': xo.Kernel(
            c_name='EverestRandomData_get_random',
            args=[
                xo.Arg(xp.Particles, name='particles'),
                xo.Arg(xo.Float64[:], name='samples'),
                xo.Arg(xo.Int32, name='n')
            ]),
        'gsample': xo.Kernel(
            c_name='EverestRandomData_get_random_gauss',
            args=[
                xo.Arg(xp.Particles, name='particles'),
                xo.Arg(xo.Float64[:], name='samples'),
                xo.Arg(xo.Int32, name='n')
            ]),
        'esample': xo.Kernel(
            c_name='EverestRandomData_get_random_exp',
            args=[
                xo.Arg(xp.Particles, name='particles'),
                xo.Arg(xo.Float64[:], name='samples'),
                xo.Arg(xo.Int32, name='n')
            ]),
        'rsample': xo.Kernel(
            c_name='EverestRandomData_get_random_ruth',
            args=[
                xo.Arg(xo.ThisClass, name='ran'),
                xo.Arg(xp.Particles, name='particles'),
                xo.Arg(xo.Float64[:], name='samples'),
                xo.Arg(xo.Int32, name='n')
            ]),
        'set_rutherford': xo.Kernel(
            c_name='EverestRandomData_set_rutherford',
            args=[
                xo.Arg(xo.ThisClass, name='ran'),
                xo.Arg(xo.Float64, name='z'),
                xo.Arg(xo.Float64, name='emr'),
                xo.Arg(xo.Float64, name='upper_val')
            ])
        }


    def __init__(self, **kwargs):
        if '_xobject' not in kwargs:
            kwargs.setdefault('rutherford_iterations', 7)
            kwargs.setdefault('rutherford_lower_val', 0.0009982)
            kwargs.setdefault('rutherford_upper_val', 1.)
            kwargs.setdefault('rutherford_A', 1.)
            kwargs.setdefault('rutherford_B', 1.)
            kwargs.setdefault('_samples', np.array([]))
        super().__init__(**kwargs)


    def set_rutherford_by_material(self, material):
        self.compile_kernels(only_if_needed=True, apply_to_source=[
                            partial(_handle_per_particle_blocks, local_particle_src=xp.gen_local_particle_api())
        ])
        context = self._buffer.context
        context.kernels.set_rutherford(ran=self, z=material.Z, emr=material.nuclear_radius, upper_val=material.hcut)


    def sample(self, n_seeds=1, n_samples=1, dist='uniform'):
        dist_options = ['uniform', 'gaussian', 'normal', 'exponential', 'rutherford']
        if dist not in dist_options:
            raise valueError(f"The variable 'dist' should be one of {dist_options}.")
        self.compile_kernels(only_if_needed=True, apply_to_source=[
                            partial(_handle_per_particle_blocks, local_particle_src=xp.gen_local_particle_api())
        ])
        context = self._buffer.context
        particles = xp.Particles(_capacity=n_seeds)
        particles._init_random_number_generator()
        Array = xo.Float64[:]
        samples = Array(n_seeds*n_samples, _context=context)
        if dist == 'uniform':
            context.kernels.usample(particles=particles, samples=samples, n=n_samples)
        elif dist == 'gaussian' or dist == 'normal':
            context.kernels.gsample(particles=particles, samples=samples, n=n_samples)
        elif dist == 'exponential':
            context.kernels.esample(particles=particles, samples=samples, n=n_samples)
        elif dist == 'rutherford':
            context.kernels.rsample(ran=self, particles=particles, samples=samples, n=n_samples)
        return np.reshape(samples.to_nparray(), (-1, n_samples))

