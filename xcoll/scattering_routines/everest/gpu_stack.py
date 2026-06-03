# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

"""Per-thread CUDA device-stack sizing for the Everest scattering kernels.

The Everest crystal/collimator kernels build several sizeable per-track structs
on the (per-thread) CUDA stack -- the geometry descriptor, the collimation
data, and the Everest scattering state -- and recurse a little through the
Amorphous <-> Channelling <-> volume-interaction sub-segments. Together these
exceed the CUDA default per-thread stack frame of 1024 bytes, so a GPU
(``xobjects.ContextCupy``) track of an ``EverestCrystal`` traps with
``cudaErrorIllegalAddress`` unless the device stack limit is raised first.

This module raises ``cudaLimitStackSize`` once per process, only on a Cupy
context, and never below an already-larger limit. It is a *no-op on CPU
contexts*, so CPU numerics are completely untouched.

It is deliberately NOT run as an import side-effect: the limit is only raised
when an Everest element is actually placed on / tracked on a Cupy context.
"""

# cudaLimitStackSize enum value (CUDA runtime). Same on all supported CUDA
# versions; we use the literal so we do not depend on a cupy enum symbol name.
_CUDA_LIMIT_STACK_SIZE = 0x00

# Default per-thread device stack, in bytes. Empirically the standard Everest
# crystal track needs >= 1920 bytes and the worst-case (long crystal + grazing
# pencil beam, deepest Amorphous/Channel/volume recursion) needs >= 2048 bytes;
# both fail at the CUDA default of 1024. 16 KiB gives a comfortable ~8x margin
# over the measured worst case while staying far below the device's stack
# reservation ceiling (a too-large request fails with cudaErrorMemoryAllocation).
CRYSTAL_STACK_LIMIT_BYTES = 16 * 1024

# Track the largest limit we have already requested in this process, so the
# call is idempotent and so we never lower a limit another component raised.
_requested_stack_limit = 0


def _is_cupy_context(context):
    """True iff `context` is an xobjects Cupy (GPU) context. No-op-safe."""
    try:
        import xobjects as xo
    except Exception:  # pragma: no cover - xobjects is always present here
        return False
    cupy_context_cls = getattr(xo, "ContextCupy", None)
    return cupy_context_cls is not None and isinstance(context, cupy_context_cls)


def set_crystal_stack_limit(context, nbytes=CRYSTAL_STACK_LIMIT_BYTES):
    """Ensure the CUDA per-thread stack is large enough for Everest kernels.

    Parameters
    ----------
    context : xobjects context
        The context the Everest element lives on. On a CPU context this is a
        no-op (returns ``None``) so CPU numerics are never affected.
    nbytes : int
        Requested per-thread device stack size in bytes. The effective limit is
        only ever raised, never lowered.

    Returns
    -------
    int or None
        The effective ``cudaLimitStackSize`` after the call (bytes) on a Cupy
        context, or ``None`` on a CPU context / if cupy is unavailable.
    """
    global _requested_stack_limit

    if not _is_cupy_context(context):
        return None

    try:
        import cupy as cp
    except Exception:  # pragma: no cover - cupy is present on a Cupy context
        return None

    current = cp.cuda.runtime.deviceGetLimit(_CUDA_LIMIT_STACK_SIZE)
    target = max(nbytes, current, _requested_stack_limit)

    # Only call the (device-wide) setter if we actually need to raise it; this
    # keeps the guard idempotent and avoids resetting it on every element.
    if target > current:
        cp.cuda.runtime.deviceSetLimit(_CUDA_LIMIT_STACK_SIZE, target)

    _requested_stack_limit = max(_requested_stack_limit, target)
    return cp.cuda.runtime.deviceGetLimit(_CUDA_LIMIT_STACK_SIZE)
