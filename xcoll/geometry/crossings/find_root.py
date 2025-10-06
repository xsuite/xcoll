import numpy as np

import xobjects as xo

from ..segments import LocalSegment
from ..trajectories import LocalTrajectory
from ..c_init import define_src, PyMethod
from ...general import _pkg_root

# Kernel for root finder
class FindRoot(xo.Struct):
    solution_t    = xo.Float64[:]
    solution_l    = xo.Float64[:]
    converged     = xo.Int8[:]
    num_solutions = xo.Int16
    max_solutions = xo.Int8

    _depends_on = [LocalTrajectory, LocalSegment]
    _kernels = {'newton': xo.Kernel(
                            c_name="FindRoot_newton",
                            args=[
                                xo.Arg(xo.ThisClass, name="finder"),
                                xo.Arg(LocalSegment, name="seg"),
                                xo.Arg(LocalTrajectory, name="traj"),
                                xo.Arg(xo.Float64, pointer=True, name="guess_t"),
                                xo.Arg(xo.Float64, pointer=True, name="guess_l"),
                            ], ret=None),
                'find_crossing': xo.Kernel(
                            c_name="FindRoot_find_crossing",
                            args=[
                                xo.Arg(xo.ThisClass, name="finder"),
                                xo.Arg(LocalSegment, name="seg"),
                                xo.Arg(LocalTrajectory, name="traj"),
                            ], ret=None)}
    _needs_compilation = True
    _extra_c_sources = [
        define_src,
        _pkg_root / 'geometry' / 'crossings' / 'find_root.h']

    def __init__(self, **kwargs):
        kwargs.setdefault('max_solutions', 100)
        kwargs.setdefault('num_solutions', 0)
        kwargs.setdefault('solution_t', np.zeros(kwargs['max_solutions'], dtype=np.float64))
        kwargs.setdefault('solution_l', np.zeros(kwargs['max_solutions'], dtype=np.float64))
        kwargs.setdefault('converged',  np.zeros(kwargs['max_solutions'], dtype=np.int8))
        super().__init__(**kwargs)

    def __getattr__(self, attr):
        kernel_name = attr
        if kernel_name in self._kernels:
            return PyMethod(kernel_name=kernel_name, element=self, element_name='finder')
        raise ValueError(f"Attribute {attr} not found in {self.__class__.__name__}")
