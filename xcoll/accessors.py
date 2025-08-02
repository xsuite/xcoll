# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #


class XcollAccessor:
    _typename = 'element'
    _dbtype = 'line'

    def __init__(self, db, names=None, **kwargs_to_set):
        super().__setattr__('_db', db)
        if names:
            super().__setattr__('names', names)
        for key, value in kwargs_to_set.items():
            super().__setattr__(key, value)

    @property
    def _coll_dict(self):
        return {name: self._db.get(name) for name in self.names}

    def keys(self):
        return self._coll_dict.keys()

    def values(self):
        return self._coll_dict.values()

    def items(self):
        return self._coll_dict.items()

    def __iter__(self):
        super().__setattr__('_iter_names', iter(self.names))
        return self

    def __next__(self):
        try:
            name = next(self._iter_names)
        except StopIteration:
            raise StopIteration
        else:
            return self._db[name]

    def __len__(self):
        return len(self.names)

    def __contains__(self, key):
        return key in self.names

    def __getattr__(self, attr):
        properties = {}
        for name, el in self.items():
            if hasattr(el, attr):
                properties[name] = getattr(el, attr)
        if len(properties) == 0:
            raise AttributeError(f"Attribute `{attr}` not found.")
        if len({tuple(ii) if isinstance(ii, list) else ii for ii in properties.values()}) == 1:
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
            # If value is not a dict, we assume it is a single value to set for all collimators
            for name, el in self.items():
                if hasattr(el, attr):
                    setattr(el, attr, value)

    def __getitem__(self, name):
        if name in self.names:
            return self._db[name]
        else:
            raise ValueError(f"{self._typename.capitalize()} `{name}` not found in {self._dbtype}!")

    def __repr__(self):
        return f"<{self.__class__.__name__} at {hex(id(self))}>"

    def __str__(self):
        if len(self.names) == 0:
            return ''
        res = [f'{self._typename.capitalize()}s:']
        for name in self.names:
            res.append(f"    {name:<16} ({self._db[name].__class__.__name__})")
        return "\n".join(res)


class XcollCollimatorAccessor(XcollAccessor):
    _typename = 'collimator'

    def __str__(self):
        res = []
        if len(self.families) > 0:
            res.append('Families:')
            for family, names in self.families.items():
                res.append(f"    {family:8}: {', '.join(names)}")
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
            raise ValueError(f"Neither family nor collimator `{name}` found in {self._dbtype}!")
