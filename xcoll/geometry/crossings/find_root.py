import xobjects as xo
from ..segments.segment import LocalSegment
from ..trajectories.trajectory import LocalTrajectory
from ...general import _pkg_root
from ..c_init.c_init import define_src, PyMethod

# Kernel for root finder
class find_root(xo.Struct):
    solution_l = xo.Float64
    solution_t = xo.Float64
    _depends_on = [LocalTrajectory, LocalSegment]
    _kernels = {'newton': xo.Kernel(
                            c_name="find_root_newton",
                            args=[
                                xo.Arg(xo.ThisClass, name="finder"),
                                xo.Arg(LocalTrajectory, name="traj"),
                                xo.Arg(LocalSegment, name="seg"),
                                xo.Arg(xo.Float64, pointer=False, name="guess_l"),
                                xo.Arg(xo.Float64, pointer=False, name="guess_t"),
                            ], ret=None)}
    _needs_compilation = True
    _extra_c_sources = [
        define_src,
        _pkg_root / 'geometry' / 'c_init' / 'find_root.h']

    def __getattr__(self, attr):
        kernel_name = attr
        if kernel_name in self._kernels:
            return PyMethod(kernel_name=kernel_name, element=self, element_name='finder')
        raise ValueError(f"Attribute {attr} not found in {self.__class__.__name__}")