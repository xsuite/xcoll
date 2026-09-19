# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #


class XcollAccessor:
    """A uniform way to access elements from an underlying database.

    The class provides dictionary-like access to elements, and element
    attributes can be retrieved and set directly on the accessor,
    providing a dict(element -> attribute) interface.

    The class expects an attribute '_db' (which can be anything that has
    a 'get' method) that holds the elements. Furthermore, the class needs
    'names', which can be a (potentially dynamic) property or set as an
    attribute in the init, to specify which elements can be addressed by
    the accessor. If not specified, all elements in the db are considered
    (in which case names is equal to self._db.keys()).

    The attributes inside the elements are expected to be accessed as an
    item (if the element is a dict), as an attribute (any other object),
    or via a get method (if the previous two lookups failed). Similarly,
    when setting an attribute, the accessor will first try to set it
    as an item, then as an attribute, and finally via a set method if
    available (otherwise it fails).

    Additionally, each element can have an attribute 'family' that gives
    the name of the family it belongs to, and 'non_family_attributes'
    which lists its attributes that will not be updated when the family
    attribute is set. This will auto-generate a families dict, and
    sub-accessors when accessing by family name.

    Accessing elements:
        Elements of the database can be accessed using the standard
        subscription syntax (e.g. accessor['name'] which returns the
        element itself). If `name` is a family name instead, it will
        return a sub-accessor for that family. The subscription syntax
        cannot be used to generated new elements (i.e.
            accessor['new_element'] = ...
        does not work), but it can be used to set attributes for an
        element with a settings dict:
            acc['el1'] = {'a': -10, 'd': -10, 'f': -999}
        where non-specified attributes will be left unchanged and missing
        ones are generated. When using this approach on a family name,
        all elements in the family will be updated, except for attributes
        that are listed in the non_family_attributes of an element.

    Accessing attributes:
        When accessing an attribute on the accessor, a dict of element
        names and the value of the respective attribute is returned.
        Elements that do not have the requested attribute are omitted. If
        all elements have the same value for the requested attribute,
        this value is returned directly instead of a dict. If no elements
        have the requested attribute, an AttributeError is raised.

    Accessor API:
        - accessor.names: list of element names that are represented by
            the accessor (manually set or dynamic property).
        - accessor.families: dict of family name -> list of element names
            in that family.
        - accessor.family_names: list of family names.
        - accessor.show(): print the content of the accessor.
        - iteration over the accessor iterates over the elements in the
            accessor (unlike a dict which iterates over the keys). The
            elements are returned in the same order as they appear in
            accessor.names.
        - accessor.keys(), accessor.values(), accessor.items() work as
            for a dict.


    Setting an attribute on the accessor will set it for all elements, unless the attribute is listed in an element's non_family_attributes, in which case it will be ignored for that element. Setting an item on the accessor will set attributes for that element, or for all elements in the family if a family name is used as key.
    """

    # These are only used for printing and error messages:
    _typename = 'element'   # What are we accessing? (e.g. elements, collimators, ..)
    _dbtype   = 'line'      # What is the underlying db? (e.g. line, colldb, ..)
    _eltype   = None        # What type of elements are we accessing? (e.g. BeamElement, settings, ..)

    def __init__(self, db, names=None, **kwargs_to_set):
        super().__setattr__('_db', db)
        kwargs_to_set.setdefault('_is_family_sub_accessor', False)
        if names:
            super().__setattr__('names', names)
        elif not 'names' in self.__class__.__dict__:
            super().__setattr__('names', list(db.keys()))
        for key, value in kwargs_to_set.items():
            super().__setattr__(key, value)
        self._check_family_consistency()

    def __repr__(self):
        return f"<{self.__class__.__name__} at {hex(id(self))} (use .show() " \
             + f"to see the content)>"

    def __str__(self):
        if len(self.names) == 0:
            return ''
        res = []
        if len(self.families) > 0:
            res.append('Families:')
            name_len = max(len(name) for name in self.family_names) + 1
            name_len = max(name_len, 10)
            for family, names in self.families.items():
                ff = f'{family}:'
                res.append(f"    {ff:<{name_len}}  {', '.join(names)}")
            fam_names = {vvv for vv in self.families.values() for vvv in vv}
            nofam_names = set(self.names) - fam_names
            if nofam_names:
                head = f"{'no family:':<{name_len}}"
                res.append(f"    {head}  {', '.join(nofam_names)}")
            res.append('')
        res.append(f'{self._typename.capitalize()}s:')
        name_len = max(len(name) for name in self.names)
        for name in self.names:
            if self._eltype:
                cls_name = self._eltype
            else:
                cls_name = self._element_dict[name].__class__.__name__
            res.append(f"    {name:<{name_len}}  ({cls_name})  {self[name]}")
        return "\n".join(res)

    def show(self):
        """Print the content of the accessor."""
        print(self)

    @property
    def families(self):
        families = {}
        try:
            prop_families = self.family
        except AttributeError:
            return families
        else:
            for name in self.names:
                if name in prop_families:
                    if prop_families[name] not in families:
                        families[prop_families[name]] = []
                    families[prop_families[name]].append(name)
            return families

    @property
    def family_names(self):
        return list(self.families.keys())

    @property
    def _element_dict(self):
        # Helper attribute to get correct class instances for keys,
        # values, and items methods
        return {name: self._db.get(name) for name in self.names}

    def keys(self):
        return self._element_dict.keys()

    def values(self):
        return self._element_dict.values()

    def items(self):
        return self._element_dict.items()

    def __iter__(self):
        super().__setattr__('_iter_names', iter(self.names))
        return self

    def __next__(self):
        try:
            name = next(self._iter_names)
        except StopIteration:
            raise StopIteration
        else:
            # Use _get_db_element instead of _db.get to protect against
            # potential changes in between different calls to __next__.
            return self._get_db_element(name)

    def __len__(self):
        return len(self.names)

    def __contains__(self, key):
        return key in self.names

    # @property   # This only works if the underlying db is a dict...
    # def properties(self):
    #     return {attr for ddd in self.values() for attr in ddd.keys()}

    def _get_db_element(self, name):
        if name in self.names:
            return self._db.get(name)
        else:
            raise KeyError(f"{self._typename.capitalize()} `{name}` not "
                           f"found in {self._dbtype}!")

    def _get_element_attr(self, name, attr, default=None,
                          allow_missing=False):
        if default is not None:
            # If a default is provided, we allow missing attributes
            allow_missing = True
        el = self._get_db_element(name)
        if isinstance(el, dict):
            if attr in el:
                return el[attr]
            elif not allow_missing:
                raise AttributeError(f"Attribute `{attr}` not found in "
                                     f"{self._typename} `{name}`!")
            else:
                return default
        elif hasattr(el, attr):
            return getattr(el, attr)
        else:
            try:
                return el.get(attr, default)
            except:
                if not allow_missing:
                    raise AttributeError(f"Attribute `{attr}` not found in "
                                         f"{self._typename} `{name}`!")
                return default

    def _set_element_attr(self, name, attr, value, allow_missing=True,
                          add_new=False):
        el = self._get_db_element(name)
        if not allow_missing:
            if add_new:
                raise ValueError("Cannot use `add_new=True` and "
                                 "`allow_missing=False` at the same time!")
            # Check that the attribute exists, raise error if not
            self._get_element_attr(name, attr, allow_missing=False)
        if isinstance(el, dict):
            if attr in el or add_new:
                el[attr] = value
        elif hasattr(el, attr) or add_new:
            try:
                setattr(el, attr, value)
            except:
                try:
                    el.set(attr, value)
                except:
                    raise AttributeError(f"Cannot set attribute `{attr}` in "
                                         f"{self._typename} `{name}`!")

    def _get_element_non_family_attributes(self, name):
        return self._get_element_attr(name, 'non_family_attributes',
                                      default=[])

    def _element_has_attr(self, name, attr):
        prop = self._get_element_attr(name, attr,
                                      default='__xc_acc_not_found__')
        return prop != '__xc_acc_not_found__'

    def __getattr__(self, attr):
        properties = {}
        for name in self.names:
            if self._element_has_attr(name, attr):
                properties[name] = self._get_element_attr(name, attr)
        if len(properties) == 0:
            raise AttributeError(f"Attribute `{attr}` not found in "
                                 f"{self._dbtype}!")
        # If all values are the same, return a single value
        if len(properties) == len(self.names) \
        and len({tuple(ii) if isinstance(ii, list) else ii
                for ii in properties.values()}) == 1:
            return next(iter(properties.values()))
        # If this is a family accessor, all members should have it
        if self._is_family_sub_accessor:
            missing = set(self.names) - set(properties.keys())
            for name in missing:
                # It is fine if the attribute is declared non-family
                # in this element
                if attr not in self._get_element_non_family_attributes(name):
                    raise AttributeError(f"Attribute `{attr}` is not a family "
                        + f"attribute, as {self._typename} `{name}` does not "
                        + f"have it! If this is intentional, add the attribute"
                        + f" to `{name}`s non_family_attributes list!")
        return properties

    def __setattr__(self, attr, value):
        # First check that the attribute exists 
        getattr(self, attr)
        if isinstance(value, dict):
            if self._is_family_sub_accessor:
                raise ValueError("Cannot set attribute with a dict on a "
                                 "family!")
            for name, val in value.items():
                self._set_element_attr(name, attr, val, allow_missing=False)
        else:
            # If the value is not a dict, it is either a family attribute...
            if self._is_family_sub_accessor:
                for name in self.names:
                    if attr in self._get_element_non_family_attributes(name):
                        # This attribute is non-family, so ignore it when
                        # setting the family
                        continue
                    self._set_element_attr(name, attr, value)
            # ... or a single value to set for all elements
            else:
                for name in self.names:
                    self._set_element_attr(name, attr, value)

    def __getitem__(self, name):
        # We can getitem by name or family, so we overwrite the super method
        if name in self.families:
            return XcollAccessor(db=self._db, names=self.families[name],
                        _typename=self._typename, _dbtype=self._dbtype,
                        _is_family_sub_accessor=True)
        else:
            return self._get_db_element(name)

    def __setitem__(self, name, value):
        # Set attributes of a single element. This is the only way to add new
        # attributes to an element.
        # TODO: allow to define a family this way
        if name in self.families:
            # Set attributes for all family members
            if not isinstance(value, dict):
                raise ValueError(f"Can only set family `{name}` to a "
                               + f"settings dict!")
            for el_name in self.families[name]:
                for attr, val in value.items():
                    if attr in self._get_element_non_family_attributes(el_name):
                        # This attribute is non-family, so ignore it when setting the family
                        continue
                    self._set_element_attr(el_name, attr, val, add_new=True)
        else:
            if not isinstance(value, dict):
                raise ValueError(f"Can only set {self._typename} `{name}` to a "
                               + f"settings dict!")
            # Set attributes for this element
            for attr, val in value.items():
                self._set_element_attr(name, attr, val, add_new=True)

    def _check_family_consistency(self):
        for name in self.names:
            fam = self._get_element_attr(name, 'family', allow_missing=True)
            if fam is not None and not isinstance(fam, str):
                raise ValueError(f"Inconsistent db in XcollAccessor: "
                                 f"{self._typename.capitalize()} `{name}` has "
                                 f"a `family` attribute that is not a string!")
            nf = self._get_element_non_family_attributes(name)
            if nf:
                if fam is None:
                    raise AttributeError(f"Inconsistent db in XcollAccessor: "
                        f"{self._typename.capitalize()} `{name}` has "
                        f"`non_family_attributes` but no `family` attribute!")
                if isinstance(nf, str):
                    self._set_element_attr(name, 'non_family_attributes', [nf])
                elif not hasattr(nf, '__iter__'):
                    raise ValueError(f"Inconsistent db in XcollAccessor: "
                                     f"{self._typename.capitalize()} `{name}` "
                                     f"has `non_family_attributes` but it is "
                                     f"not iterable!")
                for attr in nf:
                    if not self._element_has_attr(name, attr):
                        raise AttributeError(f"Inconsistent db in "
                            f"XcollAccessor: {self._typename.capitalize()} "
                            f"`{name}` has attribute `{attr}` in "
                            f"`non_family_attributes` but it is not present!")
                    if attr in ['family', 'non_family_attributes']:
                        raise ValueError(f"Inconsistent db in XcollAccessor: "
                            f"{self._typename.capitalize()} `{name}` has "
                            f"`{attr}` in `non_family_attributes`!")


class XcollAccessorElementProxy:
    """A view on a single element of an XcollAccessor."""

    def __init__(self, accessor: XcollAccessor, name: str):
        self._accessor = accessor
        self._name = name
        self._element = accessor._get_db_element(name)

    def __str__(self):
        return self._element.__str__()

    def __repr__(self):
        return self._element.__repr__()

    def __getattr__(self, name):
        return getattr(self._target, name)

    def __getattr__(self, attr):
        return self._accessor._get_element_attr(self._name, attr)

    def __setattr__(self, attr, value):
        if attr in ['_accessor', '_name']:
            super().__setattr__(attr, value)
        else:
            self._accessor._set_element_attr(self._name, attr, value)
