# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

from .material import Material, CrystalMaterial
from ..pretty_print import style, pad_styled


class MaterialsDatabase:
    """Class to handle the materials database for Xcoll."""

    def __init__(self):
        self._materials = {}    # Case-insensitive keys. Each entry is a dict with 'material' and 'name' keys; the latter is meant to keep formatting for easier recognition when printing
        self._aliases = {}      # Case-insensitive keys, aliases to _materials. Each entry is a dict with 'refname' and 'name' keys; the latter is meant to keep formatting for easier recognition when printing
        self._fluka_names = {}  # Case-sensitive keys, aliases to _materials
        self._geant4_names = {} # Case-sensitive keys, aliases to _materials
        self._iter_index = 0
        self._iter_names = []

    def __repr__(self):
        return f"<{self.__class__.__name__} at {hex(id(self))} (use .show() " \
             + f"to see the contents)>"

    def __str__(self):
        return self._str(format=False)

    def show(self, full=False):
        """Print the content of the accessor."""
        print(self._str(format=True, full=full))

    def _str(self, format, full=True):
        res = []
        padding = 36
        for typ in ["Atomic Elements", "Compounds", "Mixtures", "Crystal Materials",
                    "Old Sixtrack Materials"]:
            header  = style(f'{typ}', bold=True, colour='navy', enabled=format)
            if typ == "Old Sixtrack Materials":
                header += style('   (still default from colldb)', dim=True, colour='navy', enabled=format)
            else:
                header += style('  (aliases', italic=True, colour='navy', enabled=format)
                if full:
                    header += style(' <case-insensitive>', dim=True, colour='navy', enabled=format)
                header += style(' | ', italic=True, colour='navy', enabled=format)
                header += style('FLUKA-only aliases', italic=True, colour='forest_green', enabled=format)
                if full:
                    header += style(' <case-sensitive>', dim=True, colour='navy', enabled=format)
                header += style(' | ', italic=True, colour='navy', enabled=format)
                header += style('Geant4-only aliases', italic=True, colour='crimson', enabled=format)
                if full:
                    header += style(' <case-sensitive>', dim=True, colour='navy', enabled=format)
                header += style(')', italic=True, colour='navy', enabled=format)
            res.append(header)
            for name, mat in self._materials.items():
                this_mat = mat['material']
                if typ == "Atomic Elements" and (not this_mat.is_elemental
                or isinstance(this_mat, CrystalMaterial)
                or this_mat.name.startswith('K2')):
                        continue
                if typ == "Compounds" and (not this_mat.is_compound
                or isinstance(this_mat, CrystalMaterial)
                or this_mat.name.startswith('K2')):
                        continue
                if typ == "Mixtures" and (not this_mat.is_mixture
                or isinstance(this_mat, CrystalMaterial)
                or this_mat.name.startswith('K2')):
                        continue
                if typ == "Crystal Materials" and (this_mat.name.startswith('K2')
                or not isinstance(this_mat, CrystalMaterial)):
                    continue
                if typ == "Old Sixtrack Materials" \
                and not this_mat.name.startswith('K2'):
                    continue
                mess_name = style(mat['name'], bold=True, colour='navy', enabled=format)
                names = []
                names_stripped = [self._strip(mat['name'])]
                fluka_names = []
                geant4_names = []
                for refname in self._aliases.values():
                    if refname['refname'] == name:
                        names.append(style(refname['name'], italic=True, colour='navy', enabled=format))
                        names_stripped.append(self._strip(refname['name']))
                for fluka_name, refname in self._fluka_names.items():
                    if refname == name:
                        if self._strip(fluka_name) not in names_stripped:
                            fluka_names.append(style(fluka_name, italic=True, colour='forest_green', enabled=format))
                for geant4_name, refname in self._geant4_names.items():
                    if refname == name:
                        if self._strip(geant4_name) not in names_stripped:
                            geant4_names.append(style(geant4_name, italic=True, colour='crimson', enabled=format))
                if names or fluka_names or geant4_names:
                    mess_name = pad_styled(mess_name, 16)
                    mess_name += style(' (', italic=True, colour='navy', enabled=format)
                    mess_name += style(', ', italic=True, colour='navy', enabled=format).join(names)
                    mess_name += style('|', italic=True, colour='navy', enabled=format)
                    mess_name += style(',', italic=True, colour='forest_green', enabled=format).join(fluka_names)
                    mess_name += style('|', italic=True, colour='navy', enabled=format)
                    mess_name += style(', ', italic=True, colour='crimson', enabled=format).join(geant4_names)
                    mess_name += style(')', italic=True, colour='navy', enabled=format)
                if full:
                    mess_name += ':'
                    res.append(f"  {pad_styled(mess_name, padding)}  {this_mat}")
                else:
                    res.append(f"  {mess_name}")
            res.append('')
        return '\n'.join(res[:-1])  # Remove last empty line

    @property
    def fluka(self):
        return MaterialsSubDatabase(self._fluka_names, self, 'FLUKA')

    @property
    def geant4(self):
        return MaterialsSubDatabase(self._geant4_names, self, 'Geant4')

    def __getitem__(self, name):
        # Check if in main database
        if self._strip(name) in self._materials:
            return self._get(name)
        # Check if in aliases
        elif self._strip(name) in self._aliases:
            aliased_name = self._get_alias(name)
            if aliased_name not in self._materials:
                raise KeyError(f"Broken database: Material with alias '{name}' "
                               f"not found in the database.")
            return self._get(aliased_name)
        # Check in FLUKA and/or Geant4 names
        elif name in self.fluka and name in self.geant4:
            if self.fluka[name] == self.geant4[name]:
                return self.fluka[name]
            raise KeyError(f"Material with name '{name}' is ambiguous: exists "
                           f"both as FLUKA and Geant4 name but points to different "
                           f"materials. Please use db.fluka['{name}'] or "
                           f"db.geant4['{name}']")
        elif name in self.fluka:
            return self.fluka[name]
        elif name in self.geant4:
            return self.geant4[name]
        # If not found anywhere, raise error
        raise KeyError(f"Material with name '{name}' not found in the database.")

    def __setitem__(self, name, material):
        material = self._resolve(material)
        # Check if name already in main database
        if self._strip(name) in self._materials:
            if self._get(name) == material:
                return  # Already in database
            raise ValueError(f"Material with name '{name}' already exists in the database.")
        # Check if name already in aliases
        elif self._strip(name) in self._aliases:
            if self._get(self._get_alias(name)) == material:
                return  # Already in database
            raise ValueError(f"Material with name '{name}' already exists in the database.")
        # Check if name already in FLUKA names
        elif self._strip(name) in [self._strip(nn) for nn in self.fluka.keys()]:
            for nn in self.fluka.keys():
                if self._strip(nn) == self._strip(name):
                    if self.fluka[nn] == material:
                        return  # Already in database
            raise ValueError(f"Material with name '{name}' already exists in the FLUKA sub-database.")
        # Check if name already in Geant4 names
        elif self._strip(name) in [self._strip(nn) for nn in self.geant4.keys()]:
            for nn in self.geant4.keys():
                if self._strip(nn) == self._strip(name):
                    if self.geant4[nn] == material:
                        return  # Already in database
            raise ValueError(f"Material with name '{name}' already exists in the Geant4 sub-database.")
        # Check if material already exists (by identity, not by value).
        existing_name = self._find_name_by_material(material)
        if existing_name is not None:
            self._set_alias(name, existing_name)
            return
        # Otherwise, add new material
        self._set(name, material)
        if material.name is None:
            material.name = name

    def remove_material(self, material):
        material = self._resolve(material)
        # Delete from the main database
        to_delete = []
        for name, mat in self.items():
            if mat == material:
                to_delete.append(name)
        for name in to_delete:
            del self._materials[name]
        # Delete aliases
        to_delete = []
        for alias, name in self._aliases.items():
            if name['refname'] not in self._materials:
                to_delete.append(alias)
        for alias in to_delete:
            del self._aliases[alias]
        # Delete FLUKA aliases
        to_delete = []
        for fluka_name, refname in self.fluka.names.items():
            if refname not in self._materials:
                to_delete.append(fluka_name)
        for fluka_name in to_delete:
            del self._fluka_names[fluka_name]
        # Delete Geant4 aliases
        to_delete = []
        for geant4_name, refname in self.geant4.names.items():
            if refname not in self._materials:
                to_delete.append(geant4_name)
        for geant4_name in to_delete:
            del self._geant4_names[geant4_name]

    def update_material(self, old_material, new_material):
        old_material = self._resolve(old_material)
        old_name = self._find_name_by_material(old_material)
        if old_name is None:
            raise ValueError("Material not found in database.")
        self._set(old_name, new_material)

    def rename_material(self, material, new_name):
        material = self._resolve(material)
        if new_name in self._materials:
            if self[new_name] == material:
                return  # Already has this name or alias
            raise ValueError(f"Material with name '{new_name}' already exists.")
        old_name = self._find_name_by_material(material)
        del self._materials[old_name]
        self._set(new_name, material)
        # Rename aliases
        to_adapt = []
        for alias, name in self._aliases.items():
            if name['refname'] == old_name:
                to_adapt.append(alias)
        for alias in to_adapt:
            self._aliases[alias]['refname'] = new_name
        # Rename FLUKA aliases
        to_adapt = []
        for fluka_name, refname in self.fluka.items():
            if refname == old_name:
                to_adapt.append(fluka_name)
        for fluka_name in to_adapt:
            self.fluka.names[fluka_name] = new_name
        # Rename Geant4 aliases
        to_adapt = []
        for geant4_name, refname in self.geant4.items():
            if refname == old_name:
                to_adapt.append(geant4_name)
        for geant4_name in to_adapt:
            self.geant4.names[geant4_name] = new_name

    def remove_alias(self, name):
        if self._strip(name) in self._aliases:
            del self._aliases[self._strip(name)]

    def rename_alias(self, name, new_name):
        if self._strip(name) not in self._aliases:
            raise KeyError(f"Alias '{name}' not found.")
        # Check if new_name already exists
        if self._strip(new_name) in self._aliases \
        or self._strip(new_name) in self._materials:
            if self[new_name] == self[name]:
                del self._aliases[self._strip(name)]
                return  # Already has this name
            raise ValueError(f"Material with name '{new_name}' already exists.")
        refname = self._aliases.pop(self._strip(name))['refname']
        self._aliases[self._strip(new_name)] = {'name': new_name, 'refname': refname}

    def __contains__(self, name):
        if name is None:
            return False
        elif isinstance(name, Material):
            return name in self.values()
        elif isinstance(name, str) and self._strip(name) in self._materials:
            return True
        elif isinstance(name, str) and self._strip(name) in self._aliases:
            return True
        elif isinstance(name, str) and self._strip(name) in [self._strip(nn) for nn in self.fluka.keys()]:
            return True
        elif isinstance(name, str) and self._strip(name) in [self._strip(nn) for nn in self.geant4.keys()]:
            return True
        else:
            return False

    def __iter__(self):
        self._iter_index = 0
        self._iter_names = list(self._materials.keys())
        return self

    def __next__(self):
        if self._iter_index >= len(self._iter_names):
            raise StopIteration
        value = self._materials[self._iter_names[self._iter_index]]['material']
        self._iter_index += 1
        return value

    def __len__(self):
        return len(self._materials)

    def items(self):
        for name, mat in self._materials.items():
            yield name, mat['material']

    def keys(self):
        for name in self._materials.keys():
            yield name

    def values(self):
        for mat in self._materials.values():
            yield mat['material']

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

    def _resolve(self, material):
        if isinstance(material, str):
            material = self[material]
        if not isinstance(material, Material):
            raise ValueError("'material' must be a Material instance or a string")
        return material

    def _find_name_by_material(self, material):
        for nn, mat in self.items():
            if mat == material:
                return nn
        # raise ValueError(f"Broken material database. Material {material.name} not found.")


class MaterialsSubDatabase:
    """Class to handle a subset of the materials database for Xcoll."""

    def __init__(self, names, db, typename):
        self.names = names
        self.db = db
        self.typename = typename
        self._iter_index = 0
        self._iter_names = []

    def __getitem__(self, name):
        if name in self.names:
            refname = self.names[name]
            if refname not in self.db:
                raise KeyError(f"Broken database: Material with {self.typename} "
                               f"name '{name}' not found in the database.")
            return self.db[refname]
        raise KeyError(f"Material with name '{name}' not found in the database.")

    def __setitem__(self, name, material):
        if name in self.names:
            return  # Already in database
        material = self.db._resolve(material)
        # Check if material already exists (by identity, not by value).
        existing_name = self.db._find_name_by_material(material)
        if existing_name:
            self.names[name] = self.db._strip(existing_name)
        else:
            # If not, add new material
            if material.name is None:
                material.name = name        # This will add the material to the main database
            if self.typename == 'FLUKA' and material.fluka_name is None:
                material.fluka_name = name  # This will add the alias
            if self.typename == 'Geant4' and material.geant4_name is None:
                material.geant4_name = name # This will add the alias

    def __contains__(self, name):
        return name in self.names

    def __iter__(self):
        self._iter_index = 0
        self._iter_names = list(self.names.keys())
        return self

    def __next__(self):
        if self._iter_index >= len(self._iter_names):
            raise StopIteration
        value = self.names[self._iter_names[self._iter_index]]
        self._iter_index += 1
        return value

    def __len__(self):
        return len(self.names)

    def items(self):
        for name, mat in self.names.items():
            yield name, self.db[mat]

    def keys(self):
        for name in self.names.keys():
            yield name

    def values(self):
        for mat in self.names.values():
            yield self.db[mat]

    def rename_alias(self, name, new_name):
        if name not in self:
            raise KeyError(f"Alias '{name}' not found.")
        if new_name in self:
            if self[new_name] == self[name]:
                del self.names[name]
                return  # Already has this name
            raise ValueError(f"Material with name '{new_name}' already exists.")
        refname = self.names.pop(name)
        self.names[new_name] = refname

    def remove_alias(self, name):
        if name in self.names:
            del self.names[name]


db = MaterialsDatabase()
