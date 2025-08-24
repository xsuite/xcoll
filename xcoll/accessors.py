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
    to specify which elements can be addressed by the accessor.

    Note that attributes cannot be set in the child class' methods (otherwise
    you'll get a RecursionError). Instead, they should be passed via
    super().__init__(**kwargs_to_set) in the child class __init__ method.
    '''

    _typename = 'element'
    _dbtype = 'line'
    _eltype = None

    def __init__(self, db, names=None, **kwargs_to_set):
        super().__setattr__('_db', db)
        if names:
            super().__setattr__('names', names)
        elif not 'names' in self.__class__.__dict__:
            super().__setattr__('names', list(db.keys()))
        for key, value in kwargs_to_set.items():
            super().__setattr__(key, value)

    @property
    def names(self):
        # TODO TODO TODO
        # This can be overwritten by the child class for a dynamic definition
        return self.line.get_elements_of_type(element_classes)[1]

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
            return self._db.get(name)

    def __len__(self):
        return len(self.names)

    def __contains__(self, key):
        return key in self.names

    def __getattr__(self, attr):
        properties = {}
        for name, el in self.items():
            if hasattr(el, attr):
                properties[name] = getattr(el, attr)
            elif isinstance(el, dict) and attr in el:
                properties[name] = el.get(attr)
        if len(properties) == 0:
            raise AttributeError(f"Attribute `{attr}` not found.")
        if len({tuple(ii) if isinstance(ii, list) else ii
                for ii in properties.values()}) == 1:
            # If all values are the same, return a single value
            return next(iter(properties.values()))
        return properties

    def __setattr__(self, attr, value):
        if isinstance(value, dict):
            for name, el in self.items():
                if name in value:
                    if not hasattr(el, attr):
                        raise AttributeError(f"Attribute `{attr}` not found in "
                                           + f"{self._typename} `{name}`.")
                    setattr(el, attr, value[name])
        else:
            # If value is not a dict, we assume it is a single value to set for all elements
            for name, el in self.items():
                if hasattr(el, attr):
                    setattr(el, attr, value)

    def __getitem__(self, name):
        if name in self.names:
            return self._db.get(name)
        else:
            raise ValueError(f"{self._typename.capitalize()} `{name}` not "
                           + f"found in {self._dbtype}!")

    def __repr__(self):
        return f"<{self.__class__.__name__} at {hex(id(self))} (use .show() " \
             + f"to see the content)>"

    def __str__(self):
        if len(self.names) == 0:
            return ''
        res = [f'{self._typename.capitalize()}s:']
        for name in self.names:
            cls_name = self._eltype or self._db.get(name).__class__.__name__
            res.append(f"    {name:<16} ({cls_name})")
        return "\n".join(res)

    def show(self):
        """Print the content of the accessor."""
        print(self)


class XcollCollimatorAccessor(XcollAccessor):
    '''This class can be used as a parent to provide a uniform way to access
    collimators from an underlying database (like a settings dictionary or an
    xtrack.Line). It works as the parent XcollAccessor class but additionally it
    has an attribute '_family_db' that holds the settings for a family of
    collimators.

    Each collimator can have an attribute 'family' that gives the name of the
    family it belongs to and an attribute 'overwritten_keys' that lists its
    attributes that are not taken from the family (even though the latter
    provides them).

====================
      It provides a dictionary-like access to the elements, and
    furthermore the element attributes can be retrieved and set as attributes
    on the accessor directly (providing a dict(element -> attribute) interface).

    The class expects an attribute '_db' (which can be anything that has a 'get'
    method) that holds the elements. Furthermore, the class needs 'names', which
    can be a (potentially dynamic) property or set as an attribute in the init,
    to specify which elements can be addressed by the accessor.

    Note that attributes cannot be set in the child class' methods (otherwise
    you'll get a RecursionError). Instead, they should be passed via
    super().__init__(**kwargs_to_set) in the child class __init__ method.
    '''
    _typename = 'collimator'

    def __init__(self, db, names=None, family_db=None, family_names=None, **kwargs_to_set):
        kwargs_to_set['_family_db'] = family_db or {}
        if family_names:
            kwargs_to_set['family_names'] = family_names
        elif not 'family_names' in self.__class__.__dict__:
            super().__setattr__('family_names', list(family_db.keys()))
        super().__init__(_db=db, names=names, **kwargs_to_set)

    def __str__(self):
        res = []
        if len(self.families) > 0:
            res.append('Families:')
            for family, names in self.families.items():
                ff = f'{family}:'
                res.append(f"    {ff:<16} {', '.join(names)}")
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
            # TODO: klopt niet want families kan niet autogegenereerd worden.
            # Family heeft bepaalde settings, en sommige collimators overriden die
            return XcollAccessor(db=self._db, names=self.families[name], _typename='collimator',
                                 _dbtype=self._dbtype)
        elif name in self.names:
            return self._db[name]
        else:
            raise ValueError(f"Neither family nor collimator `{name}` found in {self._dbtype}!")
