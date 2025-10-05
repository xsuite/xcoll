# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

from .material import Material, CompoundMaterial, CrystalMaterial, MixtureMaterial
from ..pretty_print import style, pad_styled


class MaterialsDatabase:
    """Class to handle the materials database for Xcoll."""

    def __init__(self):
        self._materials = {}    # Case-insensitive keys. Each entry is a dict with 'material' and 'name' keys; the latter is meant to keep formatting for easier recognition when printing
        self._aliases = {}      # Case-insensitive keys, aliases to _materials. Each entry is a dict with 'material' and 'name' keys; the latter is meant to keep formatting for easier recognition when printing
        self._fluka_names = {}  # Case-sensitive keys, aliases to _materials
        self._geant4_names = {} # Case-sensitive keys, aliases to _materials

    def __repr__(self):
        return f"<{self.__class__.__name__} at {hex(id(self))} (use .show() " \
             + f"to see the contents)>"

    def __str__(self):
        return self._str(format=False)

    def show(self):
        """Print the content of the accessor."""
        print(self._str(format=True))

    def _str(self, format):
        res = []
        padding = 40
        for cls in [Material, CompoundMaterial, MixtureMaterial, CrystalMaterial]:
            if cls is Material:
                clsname = 'Atomic Elements'
            elif cls is CompoundMaterial:
                clsname = 'Compounds'
            elif cls is MixtureMaterial:
                clsname = 'Mixtures'
            elif cls is CrystalMaterial:
                clsname = 'Crystal Materials'
            header  = style(f'{clsname}', bold=True, colour='navy', enabled=format)
            header += style('  (aliases', italic=True, colour='navy', enabled=format)
            header += style(' <case-insensitive>', dim=True, colour='navy', enabled=format)
            header += style(' | ', italic=True, colour='navy', enabled=format)
            header += style('FLUKA-only aliases', italic=True, colour='forest_green', enabled=format)
            header += style(' <case-sensitive>', dim=True, colour='navy', enabled=format)
            header += style(' | ', italic=True, colour='navy', enabled=format)
            header += style('Geant4-only aliases', italic=True, colour='crimson', enabled=format)
            header += style(' <case-sensitive>', dim=True, colour='navy', enabled=format)
            header += style(')', italic=True, colour='navy', enabled=format)
            res.append(header)
            for name, mat in self._materials.items():
                if cls is Material:
                    if isinstance(mat['material'], CompoundMaterial) or isinstance(mat['material'], CrystalMaterial) \
                    or isinstance(mat['material'], MixtureMaterial):
                        continue
                if cls is CompoundMaterial and not isinstance(mat['material'], CompoundMaterial):
                    continue
                if cls is MixtureMaterial and not isinstance(mat['material'], MixtureMaterial):
                    continue
                if cls is CrystalMaterial and not isinstance(mat['material'], CrystalMaterial):
                    continue
                mess_name = style(mat['name'], bold=True, colour='navy', enabled=format)
                names = []
                fluka_names = []
                geant4_names = []
                for refname in self._aliases.values():
                    if refname['refname'] == name:
                        names.append(style(refname['name'], italic=True, colour='navy', enabled=format))
                for fluka_name, refname in self._fluka_names.items():
                    if refname == name:
                        fluka_names.append(style(fluka_name, italic=True, colour='forest_green', enabled=format))
                for geant4_name, refname in self._geant4_names.items():
                    if refname == name:
                        geant4_names.append(style(geant4_name, italic=True, colour='crimson', enabled=format))
                if names or fluka_names or geant4_names:
                    mess_name = pad_styled(mess_name, 15)
                    mess_name += style(' (', italic=True, colour='navy', enabled=format)
                    mess_name += style(', ', italic=True, colour='navy', enabled=format).join(names)
                    mess_name += style('|', italic=True, colour='navy', enabled=format)
                    mess_name += style(',', italic=True, colour='forest_green', enabled=format).join(fluka_names)
                    mess_name += style('|', italic=True, colour='navy', enabled=format)
                    mess_name += style(', ', italic=True, colour='crimson', enabled=format).join(geant4_names)
                    mess_name += style(')', italic=True, colour='navy', enabled=format)
                mess_name += ':'
                res.append(f"    {pad_styled(mess_name, padding)}  {mat['material']}")
            res.append('')
        return '\n'.join(res[:-1])  # Remove last empty line

    @property
    def fluka(self):
        return MaterialsSubDatabase(self._fluka_names, self, 'FLUKA')

    @property
    def geant4(self):
        return MaterialsSubDatabase(self._geant4_names, self, 'Geant4')

    def __getitem__(self, name):
        if self._strip(name) in self._materials:
            return self._get(name)
        elif self._strip(name) in self._aliases:
            aliased_name = self._get_alias(name)
            if aliased_name not in self._materials:
                raise KeyError(f"Broken database: Material with alias '{name}' "
                               f"not found in the database.")
            return self._get(aliased_name)
        if name in self.fluka and name in self.geant4:
            if self.fluka[name] is self.geant4[name]:
                return self.fluka[name]
            raise KeyError(f"Material with name '{name}' is ambiguous: exists "
                           f"both as FLUKA and Geant4 name but points to different "
                           f"materials. Please use db.fluka['{name}'] or "
                           f"db.geant4['{name}']")
        elif name in self.fluka:
            return self.fluka[name]
        elif name in self.geant4:
            return self.geant4[name]
        raise KeyError(f"Material with name '{name}' not found in the database.")

    def _resolve(self, material):
        if isinstance(material, str):
            material = self[material]
        if not isinstance(material, Material):
            raise ValueError("'material' must be a Material instance or a string")
        return material

    def __setitem__(self, name, material):
        material = self._resolve(material)
        if self._strip(name) in self._materials:
            if self._get(name) is material:
                return  # Already in database
            raise ValueError(f"Material with name '{name}' already exists in the database.")
        elif self._strip(name) in self._aliases:
            if self._get_alias(name) is material:
                return  # Already in database
            raise ValueError(f"Material with name '{name}' already exists in the database.")
        # Check if material already exists (by identity, not by value).
        # If so, just add an alias.
        for kk in self._materials:
            if self._get(kk) is material:
                self._set_alias(name, kk)
                return
        # Otherwise, add new material
        self._set(name, material)

    def __contains__(self, name):
        if self._strip(name) in self._materials:
            return True
        if self._strip(name) in self._aliases:
            return True
        return False

    def items(self):
        for name, mat in self._materials.items():
            yield name, mat['material']

    def keys(self):
        for name in self._materials.keys():
            yield name

    def values(self):
        for mat in self._materials.values():
            yield mat['material']

    def remove_material(self, material):
        material = self._resolve(material)
        to_delete = []
        for name, mat in self.items():
            if mat is material:
                to_delete.append(name)
        if len(to_delete) != 1:
            raise ValueError(f"Broken material database. Material {material.name} "
                             f"not found or found multiple times.")
        for name in to_delete:
            del self._materials[name]
        to_delete = []
        for alias, name in self._aliases.items():
            if name['refname'] not in self._materials:
                to_delete.append(alias)
        for alias in to_delete:
            del self._aliases[alias]
        to_delete = []
        for fluka_name, name in self._fluka_names.items():
            if name not in self._materials:
                to_delete.append(fluka_name)
        for fluka_name in to_delete:
            del self._fluka_names[fluka_name]
        to_delete = []
        for geant4_name, name in self._geant4_names.items():
            if name not in self._materials:
                to_delete.append(geant4_name)
        for geant4_name in to_delete:
            del self._geant4_names[geant4_name]

    def rename_material(self, material, new_name):
        material = self._resolve(material)
        if new_name in self:
            if self[new_name] is material:
                return  # Already has this name or alias
            raise ValueError(f"Material with name '{new_name}' already exists.")
        old_name = self._find_material_by_value(material)
        del self._materials[old_name]
        self._materials[new_name] = material
        to_adapt = []
        for alias, name in self._aliases.items():
            if name['refname'] == old_name:
                to_adapt.append(alias)
        for alias in to_adapt:
            self._aliases[alias]['refname'] = new_name
        to_adapt = []
        to_delete = []
        for fluka_name, name in self._fluka_names.items():
            if fluka_name.lower() in self:
                # This FLUKA name is also a main name or alias, so remove it
                to_delete.append(fluka_name)
            elif name == old_name:
                to_adapt.append(fluka_name)
        for fluka_name in to_adapt:
            self._fluka_names[fluka_name] = new_name
        for fluka_name in to_delete:
            del self._fluka_names[fluka_name]
        to_adapt = []
        to_delete = []
        for geant4_name, name in self._geant4_names.items():
            if geant4_name.lower() in self:
                # This Geant4 name is also a main name or alias, so remove it
                to_delete.append(geant4_name)
            elif name == old_name:
                to_adapt.append(geant4_name)
        for geant4_name in to_adapt:
            self._geant4_names[geant4_name] = new_name
        for geant4_name in to_delete:
            del self._geant4_names[geant4_name]

    # ----- Helpers -----

    def _strip(self, name):
        return name.strip().lower().replace(' ', '').replace('-', '').replace('_', '')

    def _get(self, name):
        return self._materials[self._strip(name)]['material']

    def _set(self, name, material):
        self._materials[self._strip(name)] = {'name': name, 'material': material}

    def _get_alias(self, name):
        return self._aliases[self._strip(name)]['refname']

    def _set_alias(self, name, refname):
        sname = self._strip(name)
        refname = self._strip(refname)
        self._aliases[sname] = {'name': name, 'refname': refname}
        # Remove FLUKA / Geant4 names if they are equal to this alias
        to_delete = []
        for fluka_name, refname in self._fluka_names.items():
            if self._strip(fluka_name) == sname:   # Same alias!
                if refname == refname:              # pointing to the same material
                    to_delete.append(fluka_name)
        for fluka_name in to_delete:
            del self._fluka_names[fluka_name]
        to_delete = []
        for geant4_name, refname in self._geant4_names.items():
            if self._strip(geant4_name) == sname:   # Same alias!
                if refname == refname:           # pointing to the same material
                    to_delete.append(geant4_name)
        for geant4_name in to_delete:
            del self._geant4_names[geant4_name]

    def _find_material_by_value(self, material):
        name = None
        for nn, mat in self.items():
            if mat is material:
                name = nn
        if name is None:
            raise ValueError(f"Broken material database. Material {material.name} not found.")
        return name


db = MaterialsDatabase()


class MaterialsSubDatabase:
    """Class to handle a subset of the materials database for Xcoll."""

    def __init__(self, names, db, typename):
        self.names = names
        self.db = db
        self.typename = typename

    def __getitem__(self, name):
        mat = None
        if name in self.names:
            alias = self.names[name]
            if alias not in self.db:
                raise KeyError(f"Broken database: Material with {self.typename} "
                               f"name '{name}' not found in the database.")
            mat = self.db[alias]
        # If not found in FLUKA/Geant4 names, try main names (case-insensitive),
        # as some FLUKA/Geant4 names are identical to main names (and hence not stored here).
        elif name in self.db:
            mat = self.db[alias]
        if mat is None:
            raise KeyError(f"Material with name '{name}' not found in the database.")
        if self.typename == 'FLUKA' and mat.fluka_name != name:
            raise KeyError(f"Broken database: Material with {self.typename} "
                            f"name '{name}' does not match material's "
                            f"FLUKA name '{mat.fluka_name}'.")
        if self.typename == 'Geant4' and mat.geant4_name != name:
            raise KeyError(f"Broken database: Material with {self.typename} "
                            f"name '{name}' does not match material's "
                            f"Geant4 name '{mat.geant4_name}'.")
        return mat

    def __setitem__(self, name, material):
        material = self.db._resolve(material)
        if name in self.names:
            if self.db[self.names[name]] is material:
                return  # Already in database
            raise ValueError(f"Material with name '{name}' already exists in the {self.typename} database.")
        # Check if material already exists (by identity, not by value).
        for kk, vv in self.db.items():
            if vv is material:
                # If so, just add an alias (keep the name case-sensitive).
                # But only if name is not already an alias
                if name not in self.db:
                    self.names[name] = self.db._strip(kk)
                return  # Return anyway, even if name already existed
        # If not, add new material
        if material.name is None:
            material.name = name
        if self.typename == 'FLUKA' and material.fluka_name is None:
            material.fluka_name = name
        if self.typename == 'Geant4' and material.geant4_name is None:
            material.geant4_name = name
        self.db[material.name] = material # This will add the material to the main database
        self[name] = material                    # This will add the alias

    def __contains__(self, name):
        return name in self.names

    def items(self):
        for name, mat in self.names.items():
            yield name, self.db[mat]

    def keys(self):
        for name in self.names.keys():
            yield name

    def values(self):
        for mat in self.db.values():
            yield self.db[mat]

    def rename_material(self, material, new_name):
        material = self.db._resolve(material)
        if new_name in self:
            if self[new_name] is material:
                return  # Already has this name or alias
            raise ValueError(f"Material with name '{new_name}' already exists.")
        name = self._find_material_by_value(material)
        to_adapt = []
        to_delete = []
        for nn, refname in self._names.items():
            if self.db._strip(nn) in self.db and self.db[self.db._strip(nn)] is material:
                # This name is also a main name or alias, so remove it (no longer needed in the subdatabase)
                to_delete.append(nn)
            elif refname == name:
                to_adapt.append(nn)
        for nn in to_delete:
            del self._names[nn]
        for nn in to_adapt:
            self._names[nn] = new_name
