# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

# # How to use:
# # ===========

# # In master package __init__.py:
# from .aggregate_constants import merge_constants
# merge_constants( # Merge all providers into master namespace
#     entry_point_group="xtrack.constants",
#     target_module=__name__,
#     include_meta=True,      # also attach <group>_<plural>_meta
#     include_groups=True,    # also attach <group>_<plural>_groups (if providers defined any)
#     import_vars=False       # set to True if you also want the constants themselves at root
# )

# # In dependent package pyproject.toml:
# [project.entry-points."xtrack.constants"]
# particle_states = "xcoll.headers.particle_states:XcollParticleState"
# # you can list multiple tables:
# # particle_masses = "xcoll.headers.particle_masses:XcollParticleMass"



from __future__ import annotations

import sys
import numpy as np
from types import MappingProxyType
from typing import Dict, Tuple, Optional, Union

try:
    from importlib.metadata import entry_points  # Py3.10+
except ImportError:  # pragma: no cover
    from importlib_metadata import entry_points  # backport

from .constants import Constants

Number = Union[int, float]

class MergeError(ValueError): ...
class TableKey(tuple):
    # (group, noun); we also keep plural/reverse alongside the data store
    pass

def merge_constants(
    *,
    entry_point_group: str = "constants",
    target_module: Optional[str] = None,
    include_meta: bool = True,
    include_groups: bool = True,
    import_vars: bool = False,          # also re-export the individual constants
    overwrite_vars: bool = False,       # if True, allow overwriting same-valued vars
    order_by: str = "entrypoint",       # "entrypoint" | "alpha"
) -> Dict[Tuple[str, str], Dict[str, object]]:
    """
    Discover provider Constants via entry points and merge tables by (group, noun).

    Returns a dict keyed by (group, noun) with:
      {
        "plural": <plural>,
        "reverse": <reverse>,
        "map": Ordered[name -> number],
        "names": Optional[number -> name],   # only when reverse == "unique"
        "meta": Ordered[name -> {value, c_name, info}]   # when include_meta
        "groups": Ordered[group_name -> frozenset[number]]  # when include_groups
      }
    """
    # discover
    eps = entry_points()
    eps = eps.select(group=entry_point_group) if hasattr(eps, "select") else eps.get(entry_point_group, [])
    eplist = list(eps)
    if order_by == "alpha":
        eplist.sort(key=lambda e: (getattr(e, "name", ""), getattr(getattr(e, "dist", None), "name", "")))

    stores: Dict[Tuple[str, str], Dict[str, object]] = {}

    # track global collisions (in case providers didn’t share the same constants module)
    global_names: Dict[str, str] = {}            # NAME -> "dist/module"
    global_values: Dict[Tuple[str, str], Dict[Number, str]] = {}  # per (group, noun) when unique

    for ep in eplist:
        obj = ep.load()
        if not isinstance(obj, type) or not issubclass(obj, Constants):
            raise MergeError(f"Entry point {ep} does not load a Constants subclass (got {obj!r}).")

        g, n, p, r = obj.__group__, obj.__noun__, obj.__plural__, obj.__reverse__

        key = (g, n)
        store = stores.setdefault(key, {
            "plural": p,
            "reverse": r,
            "map": {},           # ordered
            "meta": {} if include_meta else None,
            "names": {} if r == "unique" else None,
            "groups": {} if include_groups and hasattr(obj, f"{p}_groups") else None,
        })

        # check compatibility
        if store["reverse"] != r:
            raise MergeError(f"Incompatible reverse modes for ({g},{n}): {store['reverse']!r} vs {r!r} from {ep}.")
        if store["plural"] != p:
            # We could normalize to first plural; better to require consistency
            raise MergeError(f"Inconsistent plural for ({g},{n}): {store['plural']!r} vs {p!r} from {ep}.")

        # provider maps
        prov_map: Dict[str, Number] = dict(getattr(obj, p))               # <group>_<plural>
        prov_meta: Optional[Dict[str, dict]] = dict(getattr(obj, f"{p}_meta")) if include_meta else None
        prov_groups: Optional[Dict[str, frozenset]] = None
        if include_groups and hasattr(obj, f"{p}_groups"):
            prov_groups = dict(getattr(obj, f"{p}_groups"))

        # name collisions (global)
        prov_source = f"{getattr(getattr(ep, 'dist', None), 'name', 'unknown')}/{obj.__module__}"
        for name in prov_map:
            if name in global_names and global_names[name] != prov_source:
                raise MergeError(f"Name {name!r} defined by both {global_names[name]} and {prov_source}.")
            global_names[name] = prov_source

        # merge forward map (ordered union)
        merged_map: Dict[str, Number] = store["map"]  # type: ignore
        for k, v in prov_map.items():
            if k in merged_map and merged_map[k] != v:
                raise MergeError(f"Different values for {k!r} in ({g},{n}): {merged_map[k]!r} vs {v!r} (from {prov_source}).")
            merged_map.setdefault(k, v)

        # reverse map (if unique)
        if r == "unique":
            merged_names: Dict[Number, str] = store["names"]  # type: ignore
            gv = global_values.setdefault(key, {})
            for k, v in prov_map.items():
                if v in merged_names and merged_names[v] != k:
                    raise MergeError(f"Value {v!r} used by {merged_names[v]!r} and {k!r} in ({g},{n}).")
                if v in gv and gv[v] != k:
                    raise MergeError(f"Global value collision {v!r} for {gv[v]!r} and {k!r} in ({g},{n}).")
                merged_names.setdefault(v, k)
                gv.setdefault(v, k)

        # meta merge (must not disagree)
        if include_meta and prov_meta is not None:
            merged_meta: Dict[str, dict] = store["meta"]  # type: ignore
            for k, md in prov_meta.items():
                if k in merged_meta:
                    # require same c_name and info if both present
                    prev = merged_meta[k]
                    if any(prev.get(f) != md.get(f) for f in ("value", "c_name", "info")):
                        raise MergeError(f"Metadata conflict for {k!r} in ({g},{n}).")
                merged_meta.setdefault(k, md)

        # groups merge (only if provider classes defined groups)
        if include_groups and hasattr(obj, f"{p}_groups"):
            prov_groups = dict(getattr(obj, f"{p}_groups"))  # name -> tuple
            merged_groups = store["groups"]  # name -> tuple
            for gname, gseq in prov_groups.items():
                if gname in merged_groups and merged_groups[gname] != gseq:
                    raise MergeError(f"Group {gname!r} differs between providers for ({g},{n}).")
                merged_groups.setdefault(gname, gseq)


        # optionally import individual variables into target (numbers/bools only)
        if import_vars and target_module:
            tgt = sys.modules[target_module]
            # NOTE: at master we re-export plain numbers/bools (no per-const docstrings)
            for k, v in prov_map.items():
                if hasattr(tgt, k) and not overwrite_vars:
                    raise MergeError(f"{target_module} already has {k}; set overwrite_vars=True to replace.")
                setattr(tgt, k, v)

    # write merged surfaces into target module
    if target_module:
        tgt = sys.modules[target_module]
        for (g, n), data in stores.items():
            p = data["plural"]
            # forward map
            setattr(tgt, f"{g}_{p}", MappingProxyType(dict(data["map"])))  # type: ignore
            # reverse (if unique)
            if data["reverse"] == "unique":
                setattr(tgt, f"{g}_{n}_names", MappingProxyType(dict(data["names"])))  # type: ignore
            # meta (optional)
            if include_meta and data["meta"] is not None:
                setattr(tgt, f"{g}_{p}_meta", MappingProxyType(dict(data["meta"])))  # type: ignore
            # groups (optional)
            if include_groups and data["groups"] is not None:
                setattr(tgt, f"{g}_{p}_groups", MappingProxyType(dict(data["groups"])))  # type: ignore
            # __all__ housekeeping (minimal)
            all_list = list(getattr(tgt, "__all__", []))
            for nm in (f"{g}_{p}",):
                if nm not in all_list: all_list.append(nm)
            if data["reverse"] == "unique":
                nm = f"{g}_{n}_names"
                if nm not in all_list: all_list.append(nm)
            if include_meta:
                nm = f"{g}_{p}_meta"
                if nm not in all_list: all_list.append(nm)
            if include_groups and data["groups"] is not None:
                nm = f"{g}_{p}_groups"
                if nm not in all_list: all_list.append(nm)
            setattr(tgt, "__all__", all_list)

    return stores
