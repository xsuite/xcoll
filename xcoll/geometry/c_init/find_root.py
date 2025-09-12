import xobjects as xo
from ..segments.segment import LocalSegment
from ..trajectories.trajectory import LocalTrajectory
from ...general import _pkg_root
from .c_init import define_src, PyMethod

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

    def find_guess_l_and_t(self, traj, seg):
        # we assume we checked beforehand if there is an overlap (wrt analytical solutions)
        og_t1 = seg.t1
        og_t2 = seg.t2
        og_l1 = traj.l1
        og_l2 = traj.l2
        tol = 0.1
        multiple = True
        # atm I will do them one by one. Could be optimized later
        # Find guess l for trajectory 
        #while multiple == True:
        #    multiple = False
        while (seg.t2 - seg.t1) > tol:
            t1_old = seg.t1
            t2_old = seg.t2
            t_marked = (t1_old + t2_old) / 2
            seg.t2 = t_marked
            overlaps_1  = seg.box.overlaps(b2=traj.box)
            seg.t2 = t2_old
            seg.t1 = t_marked
            overlaps_2  = seg.box.overlaps(b2=traj.box)

            # if overlaps_1 and overlaps_2:
                # more_crossings = True
                # t1_2 = t_marked 
                # t2_2 = t2_old
            if overlaps_1 and not overlaps_2:
                seg.t1 = t1_old
                seg.t2 = t_marked
            elif not overlaps_1 and overlaps_2:
                seg.t1 = t_marked
                seg.t2 = t2_old
        while (traj.l2 - traj.l1) > tol:
            l1_old = traj.l1
            l2_old = traj.l2
            t_marked = (l1_old + l2_old) / 2
            traj.l2 = t_marked
            overlaps_1  = traj.box.overlaps(b2=traj.box)
            traj.l2 = l2_old
            traj.l1 = t_marked
            overlaps_2  = traj.box.overlaps(b2=traj.box)

            # if overlaps_1 and overlaps_2:
                # more_crossings = True
                # l1_2 = t_marked 
                # l2_2 = l2_old
            if overlaps_1 and not overlaps_2:
                traj.l1 = l1_old
                traj.l2 = t_marked
            elif not overlaps_1 and overlaps_2:
                traj.l1 = t_marked
                traj.l2 = l2_old
        l1 = traj.l1
        l2 = traj.l2
        t1 = seg.t1
        t2 = seg.t2
        seg.t1 = og_t1
        seg.t2 = og_t2
        traj.l1 = og_l1
        traj.l2 = og_l2
        return (t2-t1)/2, (l2-l1)/2
        
    
