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

    Note that attributes cannot be set in the child class' methods (otherwise
    you'll get a RecursionError). Instead, they should be passed via
    super().__init__(**kwargs_to_set) in the child class __init__ method.
    '''

    # These are only used for printing and error messages:
    _typename = 'element'   # What are we accessing? (e.g. elements, collimators, ..)
    _dbtype = 'line'        # What is the underlying db? (e.g. line, colldb, ..)
    _eltype = None          # What type of elements are we accessing? (e.g. BeamElement, settings, ..)

    def __init__(self, db, names=None, **kwargs_to_set):
        super().__setattr__('_db', db)
        if names:
            super().__setattr__('names', names)
        elif not 'names' in self.__class__.__dict__:
            super().__setattr__('names', list(db.keys()))
        for key, value in kwargs_to_set.items():
            super().__setattr__(key, value)

    @property
    def _element_dict(self):
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
            return self._element_dict[name]

    def __len__(self):
        return len(self.names)

    def __contains__(self, key):
        return key in self.names

    # @property   # This only works if the underlying db is a dict...
    # def properties(self):
    #     return {attr for ddd in self.values() for attr in ddd.keys()}

    def __getattr__(self, attr):
        properties = {}
        for name, el in self.items():
            if hasattr(el, attr):
                properties[name] = getattr(el, attr)
            elif isinstance(el, dict):
                if attr in el:
                    properties[name] = el[attr]
            else:
                try:
                    properties[name] = el.get(attr)
                except:
                    pass
        if len(properties) == 0:
            raise AttributeError(f"Attribute `{attr}` not found in {self._dbtype}!")
        # If all values are the same, return a single value
        if len(properties) == len(self.names) \
        and len({tuple(ii) if isinstance(ii, list) else ii
                for ii in properties.values()}) == 1:
            return next(iter(properties.values()))
        return properties

    def __setattr__(self, attr, value):
        _ = getattr(self, attr) # Just to check that the attribute exists
        if isinstance(value, dict):
            for el, val in value.items():
                if el not in self.names:
                    raise KeyError(f"{self._typename.capitalize()} `{el}` not "
                                   + f"found in {self._dbtype}!")
                elif isinstance(self[el], dict):
                    if attr not in self[el]:
                        raise AttributeError(f"Attribute `{attr}` not found in "
                                            + f"{self._typename} `{el}`!")
                    self[el][attr] = value[el]
                else:
                    if not hasattr(self[el], attr):
                        raise AttributeError(f"Attribute `{attr}` not found in "
                                            + f"{self._typename} `{el}`!")
                    setattr(self[el], attr, value[el])
        else:
            # If value is not a dict, we assume it is a single value to set for all elements
            for name, el in self.items():
                if isinstance(el, dict):
                    if attr in el:
                        el[attr] = value
                elif hasattr(el, attr):
                    setattr(el, attr, value)

    def __getitem__(self, name):
        if name in self.names:
            return self._element_dict[name]
        else:
            raise KeyError(f"{self._typename.capitalize()} `{name}` not "
                           + f"found in {self._dbtype}!")

    def __repr__(self):
        return f"<{self.__class__.__name__} at {hex(id(self))} (use .show() " \
             + f"to see the content)>"

    def __str__(self):
        if len(self.names) == 0:
            return ''
        res = [f'{self._typename.capitalize()}s:']
        name_len = max(len(name) for name in self.names)
        for name in self.names:
            cls_name = self._eltype or self._element_dict[name].__class__.__name__
            res.append(f"    {name:<{name_len}}  ({cls_name})  {self[name]}")
        return "\n".join(res)

    def show(self):
        """Print the content of the accessor."""
        print(self)


class XcollFamilyAccessor(XcollAccessor):
    '''This class can be used as a parent to provide a uniform way to access
    collimators from an underlying database (like a settings dictionary or an
    xtrack.Line). It works similarily to the parent XcollAccessor class, but
    additionally it has an attribute '_family_db' that holds the common
    settings for families of collimators.

    Each collimator can have an attribute 'family' that gives the name of the
    family it belongs to and an attribute 'overwritten_keys' that lists its
    attributes that are not taken from the family (even though the latter
    provides them).
    '''

    # These are only used for printing and error messages:
    _typename = 'collimator'

    def __init__(self, db, names=None, family_db=None, family_names=None, **kwargs_to_set):
        kwargs_to_set['_family_db'] = family_db or {}
        if family_names:
            kwargs_to_set['family_names'] = family_names
        elif not 'family_names' in self.__class__.__dict__:
            kwargs_to_set['family_names'] = list(family_db.keys()) if family_db else []
        super().__init__(db=db, names=names, **kwargs_to_set)

    def __str__(self):
        res = []
        if len(self.families) > 0:
            res.append('Families:')
            name_len = max(len(name) for name in self.family_names)
            for family, names in self.families.items():
                ff = f'{family}:'
                res.append(f"    {ff:<{name_len}}  {', '.join(names)}")
            res.append('')
        res.append(super().__str__())
        return "\n".join(res)

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
                else:
                    if 'no family' not in families:
                        families['no family'] = []
                    families['no family'].append(name)
            return families

    def __getitem__(self, name):
        # We can getitem by name or family, so we overwrite the super method
        if name in self.families:
            return XcollAccessor(db=self._db, names=self.families[name], _typename='collimator',
                                 _dbtype=self._dbtype)
        elif name in self.names:
            return self._db[name]
        else:
            raise KeyError(f"Neither family nor collimator `{name}` found in {self._dbtype}!")
