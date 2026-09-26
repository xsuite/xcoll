# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2026.                 #
# ######################################### #

import shutil
import pytest

import xobjects as xo
import xcoll as xc
from xcoll.beam_elements import element_classes


EVEREST_CLASSES = {
    xc.EverestBlock,
    xc.EverestCollimator,
    xc.EverestCrystal,
}

# Present in element_classes, but deliberately not implemented.
NON_INSTANTIABLE_CLASSES = {
    xc.Geant4Crystal,
}


FACTORIES = {

    xc.BlackAbsorber:
        lambda ctx: xc.BlackAbsorber(
            length=1.0,
            jaw=1e-3,
            _context=ctx,
        ),

    xc.BlackCrystal:
        lambda ctx: xc.BlackCrystal(
            length=1.0,
            jaw=1e-3,
            side="+",
            bending_radius=20.0,
            _context=ctx,
        ),

    xc.TransparentCollimator:
        lambda ctx: xc.TransparentCollimator(
            length=1.0,
            jaw=1e-3,
            _context=ctx,
        ),

    xc.TransparentCrystal:
        lambda ctx: xc.TransparentCrystal(
            length=1.0,
            jaw=1e-3,
            side="+",
            bending_radius=20.0,
            _context=ctx,
        ),

    xc.EverestBlock:
        lambda ctx: xc.EverestBlock(
            length=1.0,
            material=xc.materials.CarbonFibreComposite,
            _context=ctx,
        ),

    xc.EverestCollimator:
        lambda ctx: xc.EverestCollimator(
            length=1.0,
            jaw=1e-3,
            material=xc.materials.CarbonFibreComposite,
            _context=ctx,
        ),

    xc.EverestCrystal:
        lambda ctx: xc.EverestCrystal(
            length=1.0,
            jaw=1e-3,
            side="+",
            bending_radius=20.0,
            material=xc.materials.Silicon,
            _context=ctx,
        ),

    xc.FlukaCollimator:
        lambda ctx: xc.FlukaCollimator(
            length=1.0,
            material=xc.materials.CarbonFibreComposite,
            _context=ctx,
        ),

    xc.FlukaCrystal:
        lambda ctx: xc.FlukaCrystal(
            length=1.0,
            jaw=1e-3,
            side="+",
            bending_radius=20.0,
            material=xc.materials.Silicon,
            _context=ctx,
        ),

    xc.Geant4Collimator:
        lambda ctx: xc.Geant4Collimator(
            length=1.0,
            material=xc.materials.CarbonFibreComposite,
            _context=ctx,
        ),

    xc.Geant4CollimatorTip:
        lambda ctx: xc.Geant4CollimatorTip(
            length=1.0,
            material=xc.materials.CarbonFibreComposite,
            tip_thickness=0.05,
            tip_material=xc.materials.Boron,
            _context=ctx,
        ),

    xc.BlowUp:
        lambda ctx: xc.BlowUp(
            plane="H",
            start_at_turn=0,
            stop_at_turn=1,
            _context=ctx,
        ),

    xc.EmittanceMonitor:
        lambda ctx: xc.EmittanceMonitor(
            start_at_turn=0,
            stop_at_turn=1,
            _context=ctx,
        ),
}


def _context_specs():
    specs = ["cpu", "openmp"]
    if xo.ContextCupy in xo.context.available:
        specs.append("cuda-nvrtc")
        clang = xo.settings.cuda_compiler or shutil.which("clang++")
        if clang is not None:
            specs.append("cuda-clang")

    if xo.ContextPyopencl in xo.context.available:
        try:
            specs.extend(
                f"opencl:{device}"
                for device in xo.ContextPyopencl.get_devices()
            )
        except Exception:
            # PyOpenCL importable but no actual OpenCL platform/device.
            pass

    return specs


def _make_context(spec):
    if spec == "cpu":
        return xo.ContextCpu()
    if spec == "openmp":
        return xo.ContextCpu(omp_num_threads='auto')
    if spec == "cuda-nvrtc":
        return xo.ContextCupy(backend="nvrtc")
    if spec == "cuda-clang":
        return xo.ContextCupy(backend="clang")
    if spec.startswith("opencl:"):
        return xo.ContextPyopencl(
            device=spec.split(":", 1)[1]
        )
    raise ValueError(spec)


def test_compile_factory_coverage():
    """Fail when a new element is added without adding it here."""
    concrete = set(element_classes) - NON_INSTANTIABLE_CLASSES
    assert set(FACTORIES) == concrete


@pytest.mark.context_dependent
@pytest.mark.parametrize(
    "prebuilt_kernel",
    [True, False],
    ids=["prebuilt", "compiled"]
)
@pytest.mark.parametrize(
    "context_spec",
    _context_specs(),
)
@pytest.mark.parametrize(
    "element_cls",
    sorted(FACTORIES, key=lambda cls: cls.__name__),
    ids=lambda cls: cls.__name__,
)
def test_compile_element(element_cls, context_spec, prebuilt_kernel):
    # Everest is currently CPU + CUDA only.
    if context_spec.startswith("opencl:") \
    and element_cls in EVEREST_CLASSES:
        pytest.skip("Everest does not currently support OpenCL.")

    xo.settings.force_kernel_compilation = not prebuilt_kernel
    xo.settings.allow_kernel_compilation = not prebuilt_kernel

    context = _make_context(context_spec)
    element = FACTORIES[element_cls](context)

    # Directly compile this element's kernel. In particular, do not rely on
    # Line.build_tracker(), since collective elements can be replaced by
    # drifts in the non-collective tracking kernel.
    element.compile_kernels(
        only_if_needed=False,
    )
