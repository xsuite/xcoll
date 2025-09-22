# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

from __future__ import annotations

import sys
from enum import Enum
from types import MappingProxyType, ModuleType
from typing import Dict, Mapping, Optional, Union, Tuple, Type, List


ConstantType = Union[int, float, bool]


__all__ = ["constant", "group", "Constants"]


# ---- global registries (cross-module uniqueness) ----------------------------
_XO_CONST_GLOBAL_TYPE_REG: Dict[str, List[Type[Constants]]] = {}    # category -> [defining classes]
_XO_CONST_GLOBAL_NAME_REG: Dict[str, str] = {}                      # NAME -> module
_XO_CONST_GLOBAL_UNIQUE_VALUE_REG: Dict[str, Dict[Union[int, float], Tuple[str, str]]] = {}  # category -> value -> (NAME, module)
_XO_CONST_GLOBAL_COUNT_REG: Dict[str, int] = {}  # For unnamed constants in groups and include guards in C source

# ---- helpers ----------------------------------------------------------------
def _cjoin(prefix: str, name: str, sep: str = "_") -> str:
    if not prefix:
        return name
    return f"{prefix}{'' if prefix.endswith(sep) else sep}{name}"

def _pluralise(category: str, explicit: Optional[str]) -> str:
    if explicit:
        return explicit
    n = category.lower()
    if n.endswith(("s", "x", "z", "ch", "sh")): return category + "es"
    if len(category) >= 2 and n.endswith("y") and category[-2] not in "aeiou": return category[:-1] + "ies"
    return category + "s"

def _to_number(val: ConstantType) -> Union[int, float]:
    if isinstance(val, bool):
        return int(val)           # for C and numeric maps
    if isinstance(val, int):
        return val
    if isinstance(val, float):
        return val
    raise TypeError(f"Unsupported constant type: {type(val).__name__}")

def _to_c_header(meta, include_guard, contents: Optional[List[str]] = None) -> str:
    kk = f"{include_guard}_guard"
    if kk not in _XO_CONST_GLOBAL_COUNT_REG:
        _XO_CONST_GLOBAL_COUNT_REG[kk] = -1
    _XO_CONST_GLOBAL_COUNT_REG[kk] += 1
    include_guard += f"_I{_XO_CONST_GLOBAL_COUNT_REG[kk]}"
    lines = [f"#ifndef {include_guard}",
            f"#define {include_guard}"]
    if contents:
        lines.extend(contents)
    else:
        for vv in meta.values():     # already in definition order
            lines.append(f"  #define {vv['c_name']:35} {vv['value']:25}     // {vv['info']}")
    lines.append(f"#endif /* {include_guard} */")
    lines.append("")
    return "\n".join(lines)

def _merge_map(mod: ModuleType, attr: str, new_items: Mapping, can_extend_vals: bool = False) -> None:
    existing = getattr(mod, attr, None)
    if existing is None:
        merged = dict(new_items)  # preserve insertion
    elif isinstance(existing, ModuleType):
        raise TypeError(f"Cannot merge '{attr}' in {mod.__name__}: existing is a module, not a mapping. "
                         "Do you have a python file named like the attribute?")
    else:
        merged = dict(existing)
        for nk, nv in new_items.items():
            if nk in merged and merged[nk] != nv:
                if can_extend_vals:
                    merged[nk] = tuple(dict.fromkeys(merged[nk] + nv))  # merge groups, preserve order
                else:
                    raise ValueError(f"Duplicate {attr}[{nk!r}] in {mod.__name__}: {merged[nk]} vs {nv}")
            else:
                merged[nk] = nv
    setattr(mod, attr, MappingProxyType(merged))

def _merge_c_header(mod: ModuleType, attr: str, src: str, include_guard: str) -> None:
    if hasattr(mod, attr):
        # Append to existing source
        mod_src = getattr(mod, attr).splitlines()
        new_src = ['  ' + ll for ll in src.splitlines()]
        # We need a wrapper around the individual sources to avoid clashes
        if mod_src[2].startswith(f"  #ifndef"):
            # Wrapper exists already
            new_src = "\n".join([*mod_src[:-1], *new_src, mod_src[-1], ''])
        else:
            mod_src = ['  ' + ll for ll in mod_src]
            new_src = _to_c_header(None, include_guard, contents=[*mod_src, *new_src])
        setattr(mod, attr, new_src)
    else:
        setattr(mod, attr, src)


# ---- spec holders used in class bodies --------------------------------------
class ConstantSpec:
    __slots__ = ("py_value", "info", "c_name", "name")
    def __init__(self, value: ConstantType, info: str = "",
                 c_name: Optional[str] = None, name: Optional[str] = None):
        if not isinstance(value, (int, float, bool)):
            raise TypeError(f"Unsupported constant type: {type(value).__name__}")
        self.py_value = value      # keep original type for Python export
        self.info = info
        self.c_name = c_name
        self.name = name
    def __str__(self):
        return f"ConstantSpec(name={self.name!r}, py_value={self.py_value!r}, info={self.info!r}, c_name={self.c_name!r})"
    def __repr__(self):
        return f"<{self} at {hex(id(self))}>"

def constant(value: ConstantType, info: str = "", *,
             c_name: Optional[str] = None) -> ConstantSpec:
    """Declare a numeric constant with optional metadata inside class bodies."""
    return ConstantSpec(value, info, c_name)

class GroupSpec:
    __slots__ = ("names", "py_values", "info", "name")
    def __init__(self, names: Tuple[Union[ConstantSpec, GroupSpec, Enum, str], ...],
                 info: str = "", py_values: Optional[Tuple[ConstantType, ...]] = None,
                 name: Optional[str] = None):
        if not names:
            raise ValueError("group(...) needs at least one constant name")
        self.names = names
        self.py_values = py_values
        self.info = info
        self.name = name
    def __str__(self):
        return f"GroupSpec(name={self.name!r}, names={self.names!r}, py_values={self.py_values!r}, info={self.info!r})"
    def __repr__(self):
        return f"<{self} at {hex(id(self))}>"

def group(*names: Union[ConstantSpec, GroupSpec, Enum, str], info: str = "") -> GroupSpec:
    """Define a named group of constants by their UPPERCASE names (in definition order)."""
    return GroupSpec(tuple(names), info=info)


# ---- Enum Mixin types, to add metadata to constants -------------------------
class _IntMixin(int):
    def __new__(cls, value, info: str = "", c_name: Optional[str] = None):
        return int.__new__(cls, value)
    def __init__(self, value, info: str = "", c_name: Optional[str] = None):
        self.info = info
        self.c_name = c_name
        try:
            doc = [self.info] if self.info != '' else []
            c_name = [f"Represented in C as {c_name}."] if c_name is not None else []
            self.__doc__ = ' '.join([*doc, *c_name])
        except Exception:
            pass

class _FloatMixin(float):
    def __new__(cls, value, info: str = "", c_name: Optional[str] = None):
        return float.__new__(cls, value)
    def __init__(self, value, info: str = "", c_name: Optional[str] = None):
        self.info = info
        self.c_name = c_name
        try:
            doc = [self.info] if self.info != '' else []
            c_name = [f"Represented in C as {c_name}."] if c_name is not None else []
            self.__doc__ = ' '.join([*doc, *c_name])
        except Exception:
            pass

class _TupleMixin(tuple):
    def __new__(cls, *args):
        if len(args) < 3:
            raise ValueError("_TupleMixin requires *args, info, names.")
        return tuple.__new__(cls, args[:-2])
    def __init__(self, *args):
        self.info = args[-2]
        self.names = args[-1]
        if len(self.names) != len(self):
            raise ValueError("Tuple names length does not match values length.")
        try:
            self.__doc__ = f"{self.info} Represents the following constants: {self.names}."
        except Exception:
            pass


# ---- enum member builder (returns *members*, not necessarily one enum) ------
def _build_members(module_name: str,
                   items: Dict[str, Tuple[ConstantType, str, Optional[str]]],
                   c_prefix: str,
                   allow_mixed: bool,
                   enum_type_name: str) -> Dict[str, Enum]:
    """
    Create enum members for ints and floats separately, keep bools as True/False.
    Returns: mapping name -> member (int-enum, float-enum, or True/False).
    """
    # Partition by Python type (bool stays separate to preserve True/False)
    ints: Dict[str, Tuple[int, str, Optional[str]]] = {}
    floats: Dict[str, Tuple[float, str, Optional[str]]] = {}
    bools: Dict[str, Tuple[bool, str, Optional[str]]] = {}

    for k, spec in items.items():
        if isinstance(spec.py_value, bool):
            bools[k] = (spec.py_value, spec.info, spec.c_name)
        elif isinstance(spec.py_value, int):
            ints[k] = (spec.py_value, spec.info, spec.c_name)
        elif isinstance(spec.py_value, float):
            floats[k] = (spec.py_value, spec.info, spec.c_name)
        else:
            raise TypeError(f"Unsupported type for {k}: {type(spec.py_value).__name__}")

    if ints and floats and not allow_mixed:
        raise TypeError("Mixed int/float values are not allowed when __reverse__ == 'unique'.")

    members: Dict[str, object] = {}

    if ints:
        int_enum = Enum(f"{enum_type_name}Int", {k: (v, i, c or _cjoin(c_prefix, k))
                                              for k, (v, i, c) in ints.items()},
                        module=module_name, type=_IntMixin)
        int_enum.__doc__ = ''
        for k in ints:
            members[k] = getattr(int_enum, k)

    if floats:
        float_enum = Enum(f"{enum_type_name}Float", {k: (v, i, c or _cjoin(c_prefix, k))
                                                  for k, (v, i, c) in floats.items()},
                          module=module_name, type=_FloatMixin)
        float_enum.__doc__ = ''
        for k in floats:
            members[k] = getattr(float_enum, k)

    # Bools exported as plain True/False (metadata lives in *_meta)
    for k, (bv, _, _) in bools.items():
        members[k] = bv

    return members


def _build_group_members(module_name: str,
                   items: Dict[str, Tuple[ConstantType, str, Optional[str]]],
                   enum_type_name: str) -> Dict[str, Enum]:

    members: Dict[str, object] = {}
    tup_enum = Enum(f"{enum_type_name}Tuple", {k: (*spec.py_values, spec.info, spec.names)
                                            for k, spec in items.items()},
                    module=module_name, type=_TupleMixin)
    tup_enum.__doc__ = ''
    for k in items:
        members[k] = getattr(tup_enum, k)

    return members


# ---- metaclass --------------------------------------------------------------
class _ConstantsMeta(type):
    """
    Single-inheritance constants container.
    Subclass Constants and define UPPERCASE = constant(value, info=..., c_name=...) OR plain int/float/bool.

    Enforces:
      - names unique across the process (all modules),
      - values unique across the process if __reverse__ == "unique".

    Module exports (defining module): ordered, as defined
      - <categories>, <categories>_meta, <categories>_src
      - reverse if requested: unique -> <category>_names, multi -> <category>_by_value

    export_to_module (package root): minimal
      - <categories> (+ <category>_names if unique)
    """
    def __new__(mcls, name, bases, ns, **kw):
        cls = super().__new__(mcls, name, bases, dict(ns))
        if name == "Constants":
            return cls  # Nothing to do for the base class

        if not bases or not bases[0] is Constants:
            raise TypeError(f"{name} must inherit from Constants directly; no further subclassing allowed.")

        # Config from this class only (no inheritance)
        category: str      = (ns.get("__category__",   "constant") or "constant").strip()
        plural: str        = _pluralise(category, ns.get("__plural__", None))
        reverse            = ns.get("__reverse__", "unique")  # "unique" | "multi" | None
        c_prefix: str      = (ns.get("__c_prefix__", "") or "").strip()
        include_guard: str = ns.get("__include_guard__", f"{_cjoin(c_prefix, category.upper())}_H")
        if reverse not in (None, "unique", "multi"):
            raise ValueError(f"Invalid __reverse__ {reverse!r}; must be None, 'unique', or 'multi'.")
        cls.__category__ = category
        cls.__plural__ = plural
        cls.__reverse__ = reverse
        cls.__c_prefix__ = c_prefix
        cls.__include_guard__ = include_guard

        # Register type (for cross-module consistency checks)
        if category in _XO_CONST_GLOBAL_TYPE_REG:
            for other_cls in _XO_CONST_GLOBAL_TYPE_REG[category]:
                if other_cls.__plural__ != plural:
                    raise ValueError(f"Category {category!r} already defined with plural "
                                + f"{other_cls.__plural__!r} by {other_cls.__module__}"
                                + f".{other_cls.__name__}. Please make __plural__ consistent.")
                if other_cls.__reverse__ != reverse:
                    raise ValueError(f"Category {category!r} already defined with reverse "
                                + f"{other_cls.__reverse__!r} by {other_cls.__module__}"
                                + f".{other_cls.__name__}. Please make __reverse__ consistent.")
            _XO_CONST_GLOBAL_TYPE_REG[category].append(cls)
        else:
            _XO_CONST_GLOBAL_TYPE_REG[category] = [cls]

        # Collect specifications in definition order
        own_specs: Dict[str, ConstantSpec] = {}
        own_group_specs: Dict[str, GroupSpec] = {}
        for k, v in ns.items():
            if not k.isupper() or k.startswith("_"):
                continue
            if isinstance(v, ConstantSpec):
                v.name = k  # Remember constant name for lookup in groups (potentially in other classes)
                if v.c_name is None:
                    v.c_name = _cjoin(c_prefix, k)
                own_specs[k] = v
            elif isinstance(v, ConstantType):
                own_specs[k] = ConstantSpec(v, info ='', c_name=_cjoin(c_prefix, k), name=k)
            elif isinstance(v, GroupSpec):
                v.name = k  # Remember group name for lookup in other classes
                own_group_specs[k] = v
            else:
                continue

        # Enforce global uniqueness
        mod_name = cls.__module__
        for const_name in {**own_specs, **own_group_specs}:
            if const_name in _XO_CONST_GLOBAL_NAME_REG and _XO_CONST_GLOBAL_NAME_REG[const_name] != mod_name:
                other = _XO_CONST_GLOBAL_NAME_REG[const_name]
                raise ValueError(
                    f"Constant name {const_name!r} already defined in module {other}; global uniqueness required."
                )
        for const_name in {**own_specs, **own_group_specs}:
            _XO_CONST_GLOBAL_NAME_REG[const_name] = mod_name

        if reverse == "unique":
            if category not in _XO_CONST_GLOBAL_UNIQUE_VALUE_REG:
                _XO_CONST_GLOBAL_UNIQUE_VALUE_REG[category] = {}
            for const_name, spec in own_specs.items():
                vnum = _to_number(spec.py_value)
                if vnum in _XO_CONST_GLOBAL_UNIQUE_VALUE_REG[category]:
                    other_name, other_mod = _XO_CONST_GLOBAL_UNIQUE_VALUE_REG[category][vnum]
                    raise ValueError(
                        f"Value {vnum!r} for {const_name!r} in {mod_name} "
                        f"already used by {other_name!r} in {other_mod}; "
                        f"values must be globally unique when __reverse__ == 'unique'."
                    )
            for const_name, spec in own_specs.items():
                _XO_CONST_GLOBAL_UNIQUE_VALUE_REG[category][_to_number(spec.py_value)] = (const_name, mod_name)

        # Store specs in class
        cls._own_specs_ = own_specs
        for kk, vv in own_specs.items():
            setattr(cls, kk, vv)
        name_to_val: Dict[str, Union[int, float]] = {k: _to_number(spec.py_value) for k, spec in own_specs.items()}

        # Resolve groups (ordered, no duplicates, resolved to actual values)
        if own_group_specs:
            for gname, spec in own_group_specs.items():
                seq = []
                names = []
                for nm in spec.names:
                    # References to other groups
                    if isinstance(nm, GroupSpec):       # Reference to group spec
                        if nm.name in own_group_specs:  # Defined in this class, so take resolved values
                            seq.extend(own_group_specs[nm.name].py_values)
                            names.extend(own_group_specs[nm.name].names)
                        else:                           # Defined in other class, take original spec values
                            seq.extend(nm.py_values)
                            names.extend(nm.names)
                        continue
                    elif isinstance(nm, _TupleMixin):   # Reference to group defined in other classes
                        seq.extend(nm.value)
                        names.extend(nm.names)
                        continue
                    elif isinstance(nm, ConstantSpec):  # Reference to constant spec
                        if nm.name in name_to_val:      # Defined in this class, so take resolved value
                            seq.append(name_to_val[nm.name])
                            names.append(nm.name)
                        else:                           # Defined in other class, take original spec value
                            seq.append(nm.py_value)
                            names.append(nm.name)
                        continue
                    elif isinstance(nm, (_IntMixin, _FloatMixin)):  # Reference to constant defined in other classes
                        seq.append(nm.value)
                        names.append(nm.name)
                        continue
                    elif isinstance(nm, str):           # String reference to constant or group defined in this class
                        if nm in own_group_specs:
                            seq.extend(own_group_specs[nm].py_values)
                            names.extend(own_group_specs[nm].names)
                            continue
                        elif nm in name_to_val:
                            seq.append(name_to_val[nm])
                            names.append(nm)
                            continue
                    elif isinstance(nm, ConstantType):  # Reference got evaluated already - named constant lost
                        seq.append(_to_number(nm))
                        kk = f"{category}_unnamed"
                        if kk not in _XO_CONST_GLOBAL_COUNT_REG:
                            _XO_CONST_GLOBAL_COUNT_REG[kk] = -1
                        _XO_CONST_GLOBAL_COUNT_REG[kk] += 1
                        new_name = f"{category.upper()}_UNNAMED_{_XO_CONST_GLOBAL_COUNT_REG[kk]}"
                        names.append(new_name)
                        print(f"Warning: group {gname!r} in {mod_name} contains constant value {nm!r}. "
                            + f"Name of member is lost. Named it {new_name!r}.")
                        continue
                    raise KeyError(f"Group {gname!r} references unknown constant {nm!r}")
                # Ensure no duplicates and sort
                name_to_value = {}
                for vv, nn in zip(seq, names):
                    if nn in name_to_value:
                        if name_to_value[nn] != vv:
                            raise ValueError(f"Constant {nn!r} appears with two different "
                                           + f"values: {name_to_value[nn]} and {vv}")
                    else:
                        name_to_value[nn] = vv
                sorted_items = sorted(name_to_value.items(), key=lambda kv: kv[1])
                names = tuple(ll for ll, _ in sorted_items)
                seq   = tuple(vv for _, vv in sorted_items)
                own_group_specs[gname] = GroupSpec(names=names, py_values=seq, info=spec.info, name=gname)
                setattr(cls, gname, own_group_specs[gname])
            cls._own_group_specs_ = own_group_specs

        # Class-level views (ordered)
        meta_map = {k: {"value": _to_number(spec.py_value), "c_name": spec.c_name, "info": spec.info}
                    for k, spec in own_specs.items()}
        if reverse == "unique":
            names_map = {v: k for k, v in name_to_val.items()}
        elif reverse == "multi":
            multi: Dict[Union[int, float], list[str]] = {}
            for k, v in name_to_val.items():
                multi.setdefault(v, []).append(k)
            names_map = {k: tuple(v) for k, v in multi.items()}

        # Store the maps in the class to allow merging when exporting to module
        cls._map  = name_to_val
        cls._meta = meta_map
        if reverse:
            cls._names = names_map
        if own_group_specs:
            cls._groups = {k: spec.py_values for k, spec in own_group_specs.items()}

        # Store C source in class to allow merging when exporting to module
        cls._src = _to_c_header(cls._meta, cls.__include_guard__)

        # Export to defining module (ordered, no sorting anywhere) - this is the file where the class is defined
        # Multiple classes of the same category can coexist in the same module; the variables will be joined at the module level
        cls.export_to_module(mod_name, include_groups=True, include_src=True, include_meta=True, import_vars=True)

        return cls

    def export_to_module(cls, module_name: str, include_groups: bool = False,
                         include_src: bool = False, include_meta: bool = False,
                         import_vars: bool = False) -> None:
        """Re-export minimal surface to another module (e.g. package __init__)."""
        mod = sys.modules[module_name]
        category = cls.__category__
        plural   = cls.__plural__
        reverse  = cls.__reverse__
        c_prefix = cls.__c_prefix__
        include_guard = cls.__include_guard__

        if import_vars:
            # Build members, preserving Python types
            allow_mixed = (reverse != "unique")
            enum_type_name = ''.join(nn.capitalize() for nn in category.split('_'))
            own_members = _build_members(module_name, cls._own_specs_, c_prefix, allow_mixed, enum_type_name)
            if include_groups and hasattr(cls, "_own_group_specs_"):
                own_group_members = _build_group_members(module_name, cls._own_group_specs_, enum_type_name)
            else:
                own_group_members = {}
            for k, member in {**own_members, **own_group_members}.items():
                if hasattr(mod, k):
                    raise ValueError(f"Module {mod.__name__} already has {k}; cannot redefine.")
                setattr(mod, k, member)

        _merge_map(mod, f"{plural}", cls._map)
        if reverse == "unique":
            _merge_map(mod, f"{category}_names", cls._names)
        elif reverse == "multi":
            _merge_map(mod, f"{category}_names", cls._names, can_extend_vals=True)
        if include_meta:
            _merge_map(mod, f"{plural}_meta", cls._meta)
        if include_src:
            _merge_c_header(mod, f"{plural}_src", cls._src, include_guard)
        if include_groups and hasattr(cls, "_groups"):
            _merge_map(mod, f"{category}_groups", cls._groups)

        # __all__
        all_list = list(getattr(mod, "__all__", []))
        if import_vars:
            for k in {**own_members, **own_group_members}:
                if k not in all_list:
                    all_list.append(k)
        if f"{plural}" not in all_list:
            all_list.append(f"{plural}")
        if reverse and f"{category}_names" not in all_list:
            all_list.append(f"{category}_names")
        if include_meta and f"{plural}_meta" not in all_list:
            all_list.append(f"{plural}_meta")
        if include_src and f"{plural}_src" not in all_list:
            all_list.append(f"{plural}_src")
        if include_groups and hasattr(cls, "_groups") and f"{category}_groups" not in all_list:
            all_list.append(f"{category}_groups")
        setattr(mod, "__all__", all_list)


class Constants(metaclass=_ConstantsMeta):
    """
    Subclass me and define UPPERCASE = constant(value, info=..., c_name=...) OR plain int/float/bool.

    Configure on the subclass:
        __category__       e.g. "particle_state", "interaction_type", ...
        __plural__         e.g. "particle_states" (optional; auto if omitted)
        __reverse__        "unique" | "multi" | None  (default: "unique")
        __c_prefix__       e.g. "XC" (no trailing underscore needed)
        __include_guard__  e.g. "XC_PARTICLE_STATES_H" (default derived)
    """
    pass
