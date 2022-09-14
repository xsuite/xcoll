class GeneralMaterial:
    _arguments = ['Z', 'A', 'density', 'excitation_energy', 'nuclear_radius', \
                  'nuclear_elastic_slope', 'cross_section']
    _optional_arguments = ['hcut', 'name']

    def __init_vars__(self, **kwargs):
        self._Z = kwargs.pop('Z')                                   # zatom
        self._A = kwargs.pop('A')                                   # anuc
        self._density = kwargs.pop('density')                       # rho (g cm-3)
        self._excitation_energy = kwargs.pop('excitation_energy')   # exenergy
        self._nuclear_radius = kwargs.pop('nuclear_radius')         # emr
        self._nuclear_elastic_slope = kwargs.pop('nuclear_elastic_slope') # bnref (g cm-2)
            # slope from https://journals.aps.org/prd/pdf/10.1103/PhysRevD.21.3010
        self._cross_section = kwargs.pop('cross_section')           # csref [barn]
        self._hcut = kwargs.pop('hcut', 0.02)
        self._name = kwargs.pop('name','NO NAME')
        return kwargs

    def __check_vars__(self, **kwargs):
        missing = []
        unknown = []
        for arg in self._arguments:
            if arg not in kwargs.keys():
                missing = [*missing, arg]
        if len(missing) > 0:
            raise ValueError(f"Missing required keyword(s) '{', '.join(missing)}' for class {self.__class__.__name__}!")
        for key in kwargs.keys():
            if key not in self._arguments + self._optional_arguments:
                unknown = [*unknown, key]
        if len(unknown) >0:
            raise ValueError(f"Unknown keyword(s) '{', '.join(unknown)}' for class {self.__class__.__name__}!")

    @property
    def Z(self):
        return self._Z

    @property
    def A(self):
        return self._A

    @property
    def density(self):
        return self._density

    @property
    def excitation_energy(self):
        return self._excitation_energy

    @property
    def nuclear_radius(self):
        return self._nuclear_radius

    @property
    def nuclear_elastic_slope(self):
        return self._nuclear_elastic_slope

    @property
    def cross_section(self):
        # Index 0:Total, 1:absorption, 2:nuclear elastic, 3:pp or pn elastic
        #       4:Single Diffractive pp or pn, 5:Coulomb for t above mcs
        return self._cross_section

    @property
    def hcut(self):
        return self._hcut

    @property
    def name(self):
        return self._name


    def to_dict(self):
        thisdict = {'__class__': self.__class__.__name__}
        thisdict |= {key[1:]: val for key, val in self.__dict__.items() }
        return thisdict

    @classmethod
    def from_dict(cls, thisdict):
        thisdict = { key: val for key, val in thisdict.items() if key[0:2] != '__' }
        return cls(**thisdict)






class Material(GeneralMaterial):
    _arguments = GeneralMaterial._arguments + ['radiation_length']
    _optional_arguments = GeneralMaterial._optional_arguments

    def __init__(self, **kwargs):
        self.__check_vars__(**kwargs)
        kwargs = super().__init_vars__(**kwargs)
        kwargs = self.__init_vars__(**kwargs)

    def __init_vars__(self, **kwargs):
        self._radiation_length = kwargs.pop('radiation_length')     # radl
        return kwargs

    @property
    def radiation_length(self):
        return self._radiation_length






class Crystal(GeneralMaterial):
    _arguments = GeneralMaterial._arguments + [
                  'crystal_radiation_length', 'crystal_nuclear_length', 'crystal_plane_distance', \
                  'crystal_potential', 'nuclear_collision_length'
                 ]
    _optional_arguments = GeneralMaterial._optional_arguments

    def __init__(self, **kwargs):
        self.__check_vars__(**kwargs)
        kwargs = super().__init_vars__(**kwargs)
        kwargs = self.__init_vars__(**kwargs)

    def __init_vars__(self, **kwargs):
        self._crystal_radiation_length = kwargs.pop('crystal_radiation_length')   # dlri
        self._crystal_nuclear_length = kwargs.pop('crystal_nuclear_length')       # dlyi
        self._crystal_plane_distance = kwargs.pop('crystal_plane_distance')       # ai
        self._crystal_potential = kwargs.pop('crystal_potential')                 # eUm
        self._nuclear_collision_length = kwargs.pop('nuclear_collision_length')   # collnt [m]
        return kwargs

    @property
    def crystal_radiation_length(self):
        return self._crystal_radiation_length
        
    @property
    def crystal_nuclear_length(self):
        return self._crystal_nuclear_length
        
    @property
    def crystal_plane_distance(self):
        return self._crystal_plane_distance
        
    @property
    def crystal_potential(self):
        return self._crystal_potential
        
    @property
    def nuclear_collision_length(self):
        return self._nuclear_collision_length

    @classmethod
    def from_material(cls, material, **kwargs):
        thisdict = { key: val for key, val in material.to_dict().items() if key[0:2] != '__' }
        thisdict.pop('radiation_length')    # not in crystals
        thisdict.update(kwargs)
        return cls(**thisdict)





# BE
Berylium = Material(
        name = 'Berylium',
        Z = 4.00,
        A = 9.01,
        density = 1.848,
        excitation_energy = 63.7e-9,
        radiation_length = 0.353,
        nuclear_radius = 0.22,
        nuclear_elastic_slope = 74.7,
        cross_section = [0.271, 0.192, 0, 0, 0, 0.0035e-2]
)

# AL
Aluminium = Material(
        name = 'Aluminium',
        Z = 13.00,
        A = 26.98,
        density = 2.70,
        excitation_energy = 166.0e-9,
        radiation_length = 0.089,
        nuclear_radius = 0.302,
        nuclear_elastic_slope = 120.3,
        cross_section = [0.643, 0.418, 0, 0, 0, 0.0340e-2]
)

# CU
Copper = Material(
        name = 'Copper',
        Z = 29.00,
        A = 63.55,
        density = 8.96,
        excitation_energy = 322.0e-9,
        hcut = 0.01,
        nuclear_radius = 0.366,
        radiation_length = 0.0143,
        nuclear_elastic_slope = 217.8,
        cross_section = [1.253, 0.769, 0, 0, 0, 0.1530e-2]
)

# W   (nuclear_radius = 0.520, nuclear_elastic_slope = 440.3)
Tungsten = Material(
        name = 'Tungsten',
        Z = 74.00,
        A = 183.85,
        density = 19.30,
        excitation_energy = 727.0e-9,
        hcut = 0.01,
        nuclear_radius = 0.5208318900309039, 
        radiation_length = 0.0035,
        nuclear_elastic_slope = 420.4304986267596,
        cross_section = [2.765, 1.591, 0, 0, 0, 0.7680e-2]
)
TungstenCrystal = Crystal.from_material( Tungsten,
        crystal_radiation_length = 0.0035,
        crystal_nuclear_length = 0.096,
        crystal_plane_distance = 0.56e-7,
        crystal_potential = 21.0,
        nuclear_collision_length = 0
)

# PB
Lead = Material(
        name = 'Lead',
        Z = 82.00,
        A = 207.19,
        density = 11.35,
        excitation_energy = 823.0e-9,
        hcut = 0.01,
        nuclear_radius = 0.542,
        radiation_length = 0.0056,
        nuclear_elastic_slope = 455.3,
        cross_section = [3.016, 1.724, 0, 0, 0, 0.9070e-2]
)

# C
Carbon = Material(
        name = 'Carbon',
        Z = 6.00,
        A = 12.01,
        density = 1.67,
        excitation_energy = 78.0e-9,
        radiation_length = 0.2557,
        nuclear_radius = 0.25,
        nuclear_elastic_slope = 70.0,
        cross_section = [0.337, 0.232, 0, 0, 0, 0.0076e-2]
)
CarbonCrystal = Crystal.from_material( Carbon,
        crystal_radiation_length = 0.188,
        crystal_nuclear_length = 0.400,
        crystal_plane_distance = 0.63e-7,
        crystal_potential = 21.0,
        nuclear_collision_length = 0
)

# C2
Carbon2 = Material(
        name = 'Carbon2',
        Z = 6.00,
        A = 12.01,
        density = 4.52,
        excitation_energy = 78.0e-9,
        radiation_length = 0.094,
        nuclear_radius = 0.25,
        nuclear_elastic_slope = 70.0,
        cross_section = [0.337, 0.232, 0, 0, 0, 0.0076e-2]
)

# Si
Silicon = Material(
        name = 'Silicon',
        Z = 14.00,
        A = 28.08,
        density = 2.33,
        excitation_energy = 173.0e-9,
        radiation_length = 1,
        nuclear_radius = 0.441,
        nuclear_elastic_slope = 120.14,
        cross_section = [0.664, 0.430, 0, 0, 0, 0.0390e-2]
)
SiliconCrystal = Crystal.from_material( Silicon,
        crystal_radiation_length = 0.0937,
        crystal_nuclear_length = 0.4652,
        crystal_plane_distance = 0.96e-7,
        crystal_potential = 21.34,
        nuclear_collision_length = 0.3016
)

# Ge
Germanium = Material(
        name = 'Germanium',
        Z = 32.00,
        A = 72.63,
        density = 5.323,
        excitation_energy = 350.0e-9,
        radiation_length = 1,
        nuclear_radius = 0.605,
        nuclear_elastic_slope = 226.35,
        cross_section = [1.388, 0.844, 0, 0, 0, 0.1860e-2]
)
GermaniumCrystal = Crystal.from_material( Germanium,
        crystal_radiation_length = 0.02302,
        crystal_nuclear_length = 0.2686,
        crystal_plane_distance = 1.0e-7,
        crystal_potential = 40.0,
        nuclear_collision_length = 0.1632
)

# MoGR
MolybdenumGraphite = Material(
        name = 'MolybdenumGraphite',
        Z = 6.65,
        A = 13.53,
        density = 2.500,
        excitation_energy = 87.1e-9,
        radiation_length = 0.1193,
        nuclear_radius = 0.25,
        nuclear_elastic_slope = 76.7,
        cross_section = [0.362, 0.247, 0, 0, 0, 0.0094e-2]
)

# CuCD
CopperDiamond = Material(
        name = 'CopperDiamond',
        Z = 11.90,
        A = 25.24,
        density = 5.40,
        excitation_energy = 152.9e-9,
        radiation_length = 0.0316,
        nuclear_radius = 0.308,
        nuclear_elastic_slope = 115.0,
        cross_section = [0.572, 0.370, 0, 0, 0, 0.0279e-2]
)

# Mo
Molybdenum = Material(
        name = 'Molybdenum',
        Z = 42.00,
        A = 95.96,
        density = 10.22,
        excitation_energy = 424.0e-9,
        radiation_length = 0.0096,
        nuclear_radius = 0.481,
        nuclear_elastic_slope = 273.9,
        cross_section = [1.713, 1.023, 0, 0, 0, 0.2650e-2]
)

# Glid
Glidcop = Material(
        name = 'Glidcop',
        Z = 28.80,
        A = 63.15,
        density = 8.93,
        excitation_energy = 320.8e-9,
        radiation_length = 0.0144,
        nuclear_radius = 0.418,
        nuclear_elastic_slope = 208.7,
        cross_section = [1.246, 0.765, 0, 0, 0, 0.1390e-2]
)

# Iner
Inermet = Material(
        name = 'Inermet',
        Z = 67.70,
        A = 166.70,
        density = 18.00,
        excitation_energy = 682.2e-9,
        radiation_length = 0.00385,
        nuclear_radius = 0.578,
        nuclear_elastic_slope = 392.1,
        cross_section = [2.548, 1.473, 0, 0, 0, 0.5740e-2]
)


SixTrack_to_xcoll = {
    "BE":   [Berylium],
    "AL":   [Aluminium],
    "CU":   [Copper],
    "W":    [Tungsten, TungstenCrystal],
    "PB":   [Lead],
    "C":    [Carbon, CarbonCrystal],
    "C2":   [Carbon2],
    "Si":   [Silicon, SiliconCrystal],
    "Ge":   [Germanium, GermaniumCrystal],
    "MoGR": [MolybdenumGraphite],
    "CuCD": [CopperDiamond],
    "Mo":   [Molybdenum],
    "Glid": [Glidcop],
    "Iner": [Inermet]
}
