"""Compile selected Xsuite prebuilt CPU kernels, without parallel gcc processes.

The CLI `xsuite-prebuild regenerate` accepts a `--threads` parameter but does
not expose a kernel exclusion flag. The public Python API supports selecting
kernel definitions by name, so this module filters them before compilation.

Usage (repository root)::

    python -m benchmarks.perf.prebuild_cpu --skip-tpsa
"""

from __future__ import annotations

import argparse


def select_kernel_names(definitions, *, skip_tpsa: bool = False) -> list[str]:
    """Preserve Xsuite's order while optionally excluding TPSA configurations."""
    names = []
    for name, metadata in definitions:
        is_tpsa = name.lower().startswith("tpsa_") or bool(
            metadata.get("config", {}).get("XTRACK_TPSA_TRACK", False)
        )
        if not (skip_tpsa and is_tpsa):
            names.append(name)
    return names


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--skip-tpsa", action="store_true",
                        help="Do not build TPSA tracking kernels")
    args = parser.parse_args(argv)

    from xsuite.kernel_definitions import kernel_definitions
    from xsuite.prebuild_kernels import regenerate_kernels

    selected = select_kernel_names(kernel_definitions, skip_tpsa=args.skip_tpsa)
    if not selected:
        raise RuntimeError("No kernel definitions selected")

    print("Rebuilding Xsuite serial prebuilt kernels sequentially:", flush=True)
    for name in selected:
        print(f"  - {name}", flush=True)
    regenerate_kernels(kernels=selected, n_threads=0, context="serial")


if __name__ == "__main__":
    main()
