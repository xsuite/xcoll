class Material:
    def __init__(self, **kwargs):
        self._exenergy = kwargs.pop('exenergy')
        self._mass_number = kwargs.pop('mass_number')       # anuc
        self._atomic_number = kwargs.pop('atomic_number')   # zatom
        self._density = kwargs.pop('density')               # rho (g cm-3)
        self._emr = kwargs.pop('emr')
        self._hcut = kwargs.pop('hcut',0.02)
        self._radiation_length = kwargs.pop('radiation_length')  # radl
        self._nuclear_interaction_length = kwargs.pop('nuclear_interaction_length') # bnref (g cm-2)
        self._csref = kwargs.pop('csref')
        self._matid = kwargs.pop('matid')

    @property
    def exenergy(self):
        return self._exenergy

    @property
    def mass_number(self):
        return self._mass_number

    @property
    def atomic_number(self):
        return self._atomic_number

    @property
    def density(self):
        return self._density

    @property
    def emr(self):
        return self._emr

    @property
    def hcut(self):
        return self._hcut

    @property
    def radiation_length(self):
        return self._radiation_length

    @property
    def nuclear_interaction_length(self):
        return self._nuclear_interaction_length

    @property
    def matid(self):
        return self._matid

    @property
    def csref(self):
        return self._csref
        
    @property
    def collnt(self):
        return self._collnt

    def to_dict(self):
        thisdict = {}
        thisdict['__class__'] = 'Material'
        thisdict['exenergy'] = self.exenergy
        thisdict['mass_number'] = self.mass_number
        thisdict['atomic_number'] = self.atomic_number
        thisdict['density'] = self.density
        thisdict['emr'] = self.emr
        thisdict['hcut'] = self.hcut
        thisdict['radiation_length'] = self.radiation_length
        thisdict['nuclear_interaction_length'] = self.nuclear_interaction_length
        thisdict['matid'] = self.matid
        thisdict['csref'] = self.csref
        return thisdict

    @classmethod
    def from_dict(cls, thisdict, *, engine=None):
        return cls(
          exenergy = thisdict['exenergy'],
          mass_number = thisdict['mass_number'],
          atomic_number = thisdict['atomic_number'],
          density = thisdict['density'],
          emr = thisdict['emr'],
          hcut = thisdict['hcut'],
          radiation_length = thisdict['radiation_length'],
          nuclear_interaction_length = thisdict['nuclear_interaction_length'],
          matid = thisdict['matid'],
          csref = thisdict['csref']
        )


class CrystalMaterial(Material):
    def __init__(self, **kwargs): 
                 exenergy, mass_number, atomic_number, density, emr, radiation_length,
                 nuclear_interaction_length, matid , csref, hcut=0.02,
                 crystal_radiation_length=0, crystal_nuclear_length=0,
                 crystal_plane_distance=0, crystal_potential=0, collnt=0):
        super().__init__(id, name)
        self._exenergy = exenergy
        self._mass_number = mass_number       # anuc
        self._atomic_number = atomic_number   # zatom
        self._density = density               # rho (g cm-3)
        self._emr = emr
        self._hcut = hcut
        self._radiation_length = radiation_length  # radl
        self._nuclear_interaction_length = nuclear_interaction_length # bnref (g cm-2)
        self._csref = csref
        self._crystal_radiation_length = crystal_radiation_length   # dlri
        self._crystal_nuclear_length = crystal_nuclear_length       # dlyi
        self._crystal_plane_distance = crystal_plane_distance       # ai
        self._crystal_potential = crystal_potential
        self._collnt = collnt
        self._matid = matid

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
    def collnt(self):
        return self._collnt

    def to_dict(self):
        thisdict = super().to_dict()
        thisdict['__class__'] = 'CrystalMaterial'
        thisdict['crystal_radiation_length'] = self.crystal_radiation_length
        thisdict['crystal_nuclear_length'] = self.crystal_nuclear_length
        thisdict['crystal_plane_distance'] = self.crystal_plane_distance
        thisdict['crystal_potential'] = self.crystal_potential
        thisdict['collnt'] = self.collnt
        return thisdict

    @classmethod
    def from_dict(cls, thisdict, *, engine=None):
        return cls(
          exenergy = thisdict['exenergy'],
          mass_number = thisdict['mass_number'],
          atomic_number = thisdict['atomic_number'],
          density = thisdict['density'],
          emr = thisdict['emr'],
          hcut = thisdict['hcut'],
          radiation_length = thisdict['radiation_length'],
          nuclear_interaction_length = thisdict['nuclear_interaction_length'],
          matid = thisdict['matid'],
          csref = thisdict['csref'],
          crystal_radiation_length = thisdict['crystal_radiation_length'],
          crystal_nuclear_length = thisdict['crystal_nuclear_length'],
          crystal_plane_distance = thisdict['crystal_plane_distance'],
          crystal_potential = thisdict['crystal_potential'],
          collnt = thisdict['collnt']
        )





# BE
Berylium = Material(matid=1,exenergy=63.7e-9,mass_number=9.01,atomic_number=4.00,density=1.848,emr=0.22,radiation_length=0.353,nuclear_interaction_length=74.7,csref=[0.271, 0.192, 0, 0, 0, 0.0035e-2])

# AL
Aluminium = Material(matid=2,exenergy=166.0e-9,mass_number=26.98,atomic_number=13.00,density=2.70,emr=0.302,radiation_length=0.089,nuclear_interaction_length=120.3,csref=[0.643, 0.418, 0, 0, 0, 0.0340e-2])

# CU
Copper = Material(matid=3,exenergy=322.0e-9,mass_number=63.55,atomic_number=29.00,density=8.96,emr=0.366,hcut=0.01,radiation_length=0.0143,nuclear_interaction_length=217.8,csref=[1.253, 0.769, 0, 0, 0, 0.1530e-2])

# W
Tungsten = Material(matid=4,exenergy=727.0e-9,mass_number=183.85,atomic_number=74.00,density=19.30,emr=0.5208318900309039, # 0.520,
            hcut=0.01,radiation_length=0.0035,nuclear_interaction_length=420.4304986267596, # 440.3
        csref=[2.765, 1.591, 0, 0, 0, 0.7680e-2],crystal_radiation_length=0.0035,crystal_nuclear_length=0.096,crystal_plane_distance=0.56e-7,crystal_potential=21.0,collnt=0)

# PB
Lead = Material(matid=5,exenergy=823.0e-9,mass_number=207.19,atomic_number=82.00,density=11.35,emr=0.542,hcut=0.01,radiation_length=0.0056,nuclear_interaction_length=455.3,csref=[3.016, 1.724, 0, 0, 0, 0.9070e-2])

# C
Carbon = Material(matid=6,exenergy=78.0e-9,mass_number=12.01,atomic_number=6.00,density=1.67,emr=0.25,radiation_length=0.2557,nuclear_interaction_length=70.0,csref=[0.337, 0.232, 0, 0, 0, 0.0076e-2],crystal_radiation_length=0.188,crystal_nuclear_length=0.400,crystal_plane_distance=0.63e-7,crystal_potential=21.0,collnt=0)

# C2
Carbon2 = Material(matid=7,exenergy=78.0e-9,mass_number=12.01,atomic_number=6.00,density=4.52,emr=0.25,radiation_length=0.094,nuclear_interaction_length=70.0,csref=[0.337, 0.232, 0, 0, 0, 0.0076e-2])

# Si
Silicon = Material(matid=8,exenergy=173.0e-9,mass_number=28.08,atomic_number=14.00,density=2.33,emr=0.441,radiation_length=1,nuclear_interaction_length=120.14,csref=[0.664, 0.430, 0, 0, 0, 0.0390e-2],crystal_radiation_length=0.0937,crystal_nuclear_length=0.4652,crystal_plane_distance=0.96e-7,crystal_potential=21.34,collnt=0.3016)

# Ge
Germanium = Material(matid=9,exenergy=350.0e-9,mass_number=72.63,atomic_number=32.00,density=5.323,emr=0.605,radiation_length=1,nuclear_interaction_length=226.35,csref=[1.388, 0.844, 0, 0, 0, 0.1860e-2],crystal_radiation_length=0.02302,crystal_nuclear_length=0.2686,crystal_plane_distance=1.0e-7,crystal_potential=40.0,collnt=0.1632)

# MoGR
MolybdenumGraphite = Material(matid=10,exenergy=87.1e-9,mass_number=13.53,atomic_number=6.65,density=2.500,emr=0.25,radiation_length=0.1193,nuclear_interaction_length=76.7,csref=[0.362, 0.247, 0, 0, 0, 0.0094e-2])

# CuCD
CopperDiamond = Material(matid=11,exenergy=152.9e-9,mass_number=25.24,atomic_number=11.90,density=5.40,emr=0.308,radiation_length=0.0316,nuclear_interaction_length=115.0,csref=[0.572, 0.370, 0, 0, 0, 0.0279e-2])

# Mo
Molybdenum = Material(matid=12,exenergy=424.0e-9,mass_number=95.96,atomic_number=42.00,density=10.22,emr=0.481,radiation_length=0.0096,nuclear_interaction_length=273.9,csref=[1.713, 1.023, 0, 0, 0, 0.2650e-2])

# Glid
Glidcop = Material(matid=13,exenergy=320.8e-9,mass_number=63.15,atomic_number=28.80,density=8.93,emr=0.418,radiation_length=0.0144,nuclear_interaction_length=208.7,csref=[1.246, 0.765, 0, 0, 0, 0.1390e-2])

# Iner
Inermet = Material(matid=14,exenergy=682.2e-9,mass_number=166.70,atomic_number=67.70,density=18.00,emr=0.578,radiation_length=0.00385,nuclear_interaction_length=392.1,csref=[2.548, 1.473, 0, 0, 0, 0.5740e-2])