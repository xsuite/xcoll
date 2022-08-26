class Material:
    def __init__(self, *, exenergy, mass_number, atomic_number, density, emr, hcut, radiation_length,
                  nuclear_interaction_length, matid , csref, crystal_radiation_length=0, crystal_nuclear_length=0,
                  crystal_plane_distance=0, crystal_potential=0, collnt=0):
        self._exenergy = exenergy
        self._mass_number = mass_number       # anuc
        self._atomic_number = atomic_number   # zatom
        self._density = density               # rho (g cm-3)
        self._emr = emr
        self._hcut = hcut
        self._radiation_length = radiation_length
        self._nuclear_interaction_length = nuclear_interaction_length # bnref (g cm-2)
        self._csref = csref
        self._crystal_radiation_length = crystal_radiation_length   # dlri
        self._crystal_nuclear_length = crystal_nuclear_length       # dlyi
        self._crystal_plane_distance = crystal_plane_distance       # ai
        self._crystal_potential = crystal_potential
        self._collnt = collnt
        self._matid = matid

    @property
    def can_be_crystal(self):
        return self._crystal_radiation_length != 0 and self._crystal_nuclear_length != 0 and self._crystal_plane_distance != 0 and self._crystal_potential != 0 and self._collnt != 0

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
        # TODO how to save ref to impacts?
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




# Berylium
BE = Material(matid=1,exenergy=63.7e-9,anuc=9.01,zatom=4.00,rho=1.848,emr=0.22,hcut=0.02,radl=0.353,bnref=74.7,csref=[0.271, 0.192, 0, 0, 0, 0.0035e-2],dlri=0,dlyi=0,ai=0,eUm=0,collnt=0)

# Aluminium
AL = Material(matid=2,exenergy=166.0e-9,anuc=26.98,zatom=13.00,rho=2.70,emr=0.302,hcut=0.02,radl=0.089,bnref=120.3,csref=[0.643, 0.418, 0, 0, 0, 0.0340e-2],dlri=0,dlyi=0,ai=0,eUm=0,collnt=0)

# Copper
CU = Material(matid=3,exenergy=322.0e-9,anuc=63.55,zatom=29.00,rho=8.96,emr=0.366,hcut=0.01,radl=0.0143,bnref=217.8,csref=[1.253, 0.769, 0, 0, 0, 0.1530e-2],dlri=0,dlyi=0,ai=0,eUm=0,collnt=0)

# Tungsten
W = Material(matid=4,exenergy=727.0e-9,anuc=183.85,zatom=74.00,rho=19.30,emr=0.5208318900309039, # 0.520,hcut=0.01,radl=0.0035,bnref=420.4304986267596, # 440.3
        csref=[2.765, 1.591, 0, 0, 0, 0.7680e-2],dlri=0.0035,dlyi=0.096,ai=0.56e-7,eUm=21.0,collnt=0)

# Lead
PB = Material(matid=5,exenergy=823.0e-9,anuc=207.19,zatom=82.00,rho=11.35,emr=0.542,hcut=0.01,radl=0.0056,bnref=455.3,csref=[3.016, 1.724, 0, 0, 0, 0.9070e-2],dlri=0,dlyi=0,ai=0,eUm=0,collnt=0)

# Carbon
C = Material(matid=6,exenergy=78.0e-9,anuc=12.01,zatom=6.00,rho=1.67,emr=0.25,hcut=0.02,radl=0.2557,bnref=70.0,csref=[0.337, 0.232, 0, 0, 0, 0.0076e-2],dlri=0.188,dlyi=0.400,ai=0.63e-7,eUm=21.0,collnt=1e-16  #0
)

# Carbon2
C2 = Material(matid=7,exenergy=78.0e-9,anuc=12.01,zatom=6.00,rho=4.52,emr=0.25,hcut=0.02,radl=0.094,bnref=70.0,csref=[0.337, 0.232, 0, 0, 0, 0.0076e-2],dlri=0,dlyi=0,ai=0,eUm=0,collnt=0)

# Silicon
Si = Material(matid=8,exenergy=173.0e-9,anuc=28.08,zatom=14.00,rho=2.33,emr=0.441,hcut=0.02,radl=1,bnref=120.14,csref=[0.664, 0.430, 0, 0, 0, 0.0390e-2],dlri=0.0937,dlyi=0.4652,ai=0.96e-7,eUm=21.34,collnt=0.3016)

# Germanium
Ge = Material(matid=9,exenergy=350.0e-9,anuc=72.63,zatom=32.00,rho=5.323,emr=0.605,hcut=0.02,radl=1,bnref=226.35,csref=[1.388, 0.844, 0, 0, 0, 0.1860e-2],dlri=0.02302,dlyi=0.2686,ai=1.0e-7,eUm=40.0,collnt=0.1632)

# 
MoGR = Material(matid=10,exenergy=87.1e-9,anuc=13.53,zatom=6.65,rho=2.500,emr=0.25,hcut=0.02,radl=0.1193,bnref=76.7,csref=[0.362, 0.247, 0, 0, 0, 0.0094e-2],dlri=0,dlyi=0,ai=0,eUm=0,collnt=0)

# 
CuCD = Material(matid=11,exenergy=152.9e-9,anuc=25.24,zatom=11.90,rho=5.40,emr=0.308,hcut=0.02,radl=0.0316,bnref=115.0,csref=[0.572, 0.370, 0, 0, 0, 0.0279e-2],dlri=0,dlyi=0,ai=0,eUm=0,collnt=0)

# 
Mo = Material(matid=12,exenergy=424.0e-9,anuc=95.96,zatom=42.00,rho=10.22,emr=0.481,hcut=0.02,radl=0.0096,bnref=273.9,csref=[1.713, 1.023, 0, 0, 0, 0.2650e-2],dlri=0,dlyi=0,ai=0,eUm=0,collnt=0)

# 
Glid = Material(matid=13,exenergy=320.8e-9,anuc=63.15,zatom=28.80,rho=8.93,emr=0.418,hcut=0.02,radl=0.0144,bnref=208.7,csref=[1.246, 0.765, 0, 0, 0, 0.1390e-2],dlri=0,dlyi=0,ai=0,eUm=0,collnt=0)

# 
Iner = Material(matid=14,exenergy=682.2e-9,anuc=166.70,zatom=67.70,rho=18.00,emr=0.578,hcut=0.02,radl=0.00385,bnref=392.1,csref=[2.548, 1.473, 0, 0, 0, 0.5740e-2],dlri=0,dlyi=0,ai=0,eUm=0,collnt=0)