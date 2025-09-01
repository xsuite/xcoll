# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #


class XcollAccessor:
    '''This class can be used as a parent to provide a uniform way to access
    elements from an underlying database (like a settings dictionary or an
    xtrack.Line). It provides a dictionary-like access to the elements, and
    furthermore the element attributes can be retrieved and set as attributes
    on the accessor directly (providing a dict(element -> attribute) interface).

    The class expects an attribute '_db' (which can be anything that has a 'get'
    method) that holds the elements. Furthermore, the class needs 'names', which
    can be a (potentially dynamic) property or set as an attribute in the init,
    to specify which elements can be addressed by the accessor. If not specified,
    all elements in the db are considered (names is equal to self._db.keys()).

    Additionally, each element can have an attribute 'family' that gives
    the name of the family it belongs to, and 'non_family_attributes' which
    lists its attributes that will not be updated when the family attribute is
    set. This will auto-generate a families dict, and sub-accessors when
    accessing by family name.
    '''

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
            name_len = max(max(len(name) for name in self.family_names) + 1, 10)
            for family, names in self.families.items():
                ff = f'{family}:'
                res.append(f"    {ff:<{name_len}}  {', '.join(names)}")
            nofam_names = set(self.names) - {vvv for vv in self.families.values() for vvv in vv}
            res.append(f"    {'no family:':<{name_len}}  {', '.join(nofam_names)}")
            res.append('')
        res.append(f'{self._typename.capitalize()}s:')
        name_len = max(len(name) for name in self.names)
        for name in self.names:
            cls_name = self._eltype or self._element_dict[name].__class__.__name__
            res.append(f"    {name:<{name_len}}  ({cls_name})  {self[name]}")
        return "\n".join(res)

    def show(self):
        """Print the content of the accessor."""
        print(self)

    @property
    def _element_dict(self):
        return {name: self._db.get(name) for name in self.names}

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
            return self._element_dict[name]

    def __len__(self):
        return len(self.names)

    def __contains__(self, key):
        return key in self.names

    # @property   # This only works if the underlying db is a dict...
    # def properties(self):
    #     return {attr for ddd in self.values() for attr in ddd.keys()}

    def _get_db_element(self, name):
        if name in self.names:
            return self._element_dict[name]
        else:
            raise KeyError(f"{self._typename.capitalize()} `{name}` not "
                           + f"found in {self._dbtype}!")

    def _get_element_attr(self, name, attr, default=None):
        el = self._get_db_element(name)
        if isinstance(el, dict):
            if attr in el:
                return el[attr]
            else:
                return default
        elif hasattr(el, attr):
            return getattr(el, attr)
        else:
            try:
                return el.get(attr, default)
            except:
                return default

    def _set_element_attr(self, name, attr, value, allow_missing=True):
        el = self._get_db_element(name)
        if isinstance(el, dict):
            if attr in el:
                el[attr] = value
            elif not allow_missing:
                raise AttributeError(f"Attribute `{attr}` not found in "
                                     + f"{self._typename} `{name}`!")
        elif hasattr(el, attr):
            setattr(el, attr, value)
        else:
            try:
                el.set(attr, value)
            except:
                if not allow_missing:
                    raise AttributeError(f"Attribute `{attr}` not found in "
                                         + f"{self._typename} `{name}`!")

    def _get_element_non_family_attributes(self, name):
        return self._get_element_attr(name, 'non_family_attributes', default=[])

    def __getattr__(self, attr):
        properties = {}
        for name, el in self.items():
            prop = self._get_element_attr(name, attr, default='__xc_acc_not_found__')
            if prop != '__xc_acc_not_found__':
                properties[name] = prop
        if len(properties) == 0:
            raise AttributeError(f"Attribute `{attr}` not found in {self._dbtype}!")
        # If all values are the same, return a single value
        if len(properties) == len(self.names) \
        and len({tuple(ii) if isinstance(ii, list) else ii
                for ii in properties.values()}) == 1:
            return next(iter(properties.values()))
        # If this is a family accessor, all members should have the attribute
        if self._is_family_sub_accessor:
            missing = set(self.names) - set(properties.keys())
            for name in missing:
                # It is fine if the attribute is declared non-family
                if attr not in self._get_element_non_family_attributes(name):
                    raise AttributeError(f"Attribute `{attr}` is not a family "
                        + f"attribute, as {self._typename} `{name}` does not "
                        + f"have it! If this is intentional, add the attribute "
                        + f"to `{name}`s non_family_attributes list!")
        return properties

    def __setattr__(self, attr, value):
        this_attr = getattr(self, attr) # Just to check that the attribute exists
        if isinstance(value, dict):
            if self._is_family_sub_accessor:
                raise ValueError("Cannot set attribute with a dict on a family!")
            for name, val in value.items():
                self._set_element_attr(name, attr, val, allow_missing=False)
        else:
            # If the value is not a dict, it is either a family attribute...
            if self._is_family_sub_accessor:
                for name, el in self.items():
                    if attr in self._get_element_non_family_attributes(name):
                        # This attribute is non-family, so ignore it when setting the family
                        continue
                    self._set_element_attr(name, attr, value)
            # ... or a single value to set for all elements
            else:
                for name, el in self.items():
                    self._set_element_attr(name, attr, value)

    def __getitem__(self, name):
        # We can getitem by name or family, so we overwrite the super method
        if name in self.families:
            return XcollAccessor(db=self._db, names=self.families[name], _typename=self._typename,
                                 _dbtype=self._dbtype, _is_family_sub_accessor=True)
        else:
            return self._get_db_element(name)

    # TODO: implement __setitem__ for setting as update (do not delete other attributes)

    def _check_family_consistency(self):
        for name, el in self.items():
            if hasattr(el, 'family') and not isinstance(el.family, str):
                raise TypeError(f"{self._typename.capitalize()} `{name}` has "
                                f"a `family` attribute that is not a string!")
            if hasattr(el, 'non_family_attributes'):
                if not hasattr(el, 'family'):
                    raise AttributeError(f"{self._typename.capitalize()} "
                        + f"`{name}` has `non_family_attributes` but no "
                        + f"`family` attribute!")
                if isinstance(el.non_family_attributes, str):
                    el.non_family_attributes = [el.non_family_attributes]
                elif not hasattr(el.non_family_attributes, '__iter__'):
                    raise TypeError(f"{self._typename.capitalize()} `{name}` "
                                  + f"has `non_family_attributes` but it is "
                                  + f"not iterable!")
                for attr in el.non_family_attributes:
                    if not (isinstance(el, dict) and attr in el) \
                    and not hasattr(el, attr):
                        raise AttributeError(f"{self._typename.capitalize()} "
                          + f"`{name}` has attribute `{attr}` in "
                          + f"`non_family_attributes` but it is not present!")
                    if attr == 'family':
                        raise ValueError(f"{self._typename.capitalize()} "
                          + f"`{name}` has `family` in `non_family_attributes`!")
                    if attr == 'non_family_attributes':
                        raise ValueError(f"{self._typename.capitalize()} "
                          + f"`{name}` has `non_family_attributes` in "
                          + f"`non_family_attributes`!")
