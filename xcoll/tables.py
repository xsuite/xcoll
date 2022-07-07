import xobjects as xo

class CollimatorImpactsData(xo.Struct):
    _capacity        = xo.Int64
    _row_id          = xo.Int64
    at_element       = xo.Int64[:]
    s                = xo.Float64[:]
    turn             = xo.Float64[:]
    interaction_type = xo.Int64[:]
    id_parent     = xo.Int64[:]
    x_parent      = xo.Float64[:]
    px_parent     = xo.Float64[:]
    y_parent      = xo.Float64[:]
    py_parent     = xo.Float64[:]
    zeta_parent   = xo.Float64[:]
    delta_parent  = xo.Float64[:]
    energy_parent = xo.Float64[:]
    id_child     = xo.Int64[:]
    x_child      = xo.Float64[:]
    px_child     = xo.Float64[:]
    y_child      = xo.Float64[:]
    py_child     = xo.Float64[:]
    zeta_child   = xo.Float64[:]
    delta_child  = xo.Float64[:]
    energy_child = xo.Float64[:]
    
class CollimatorImpacts(xo.dress(CollimatorImpactsData)):
    def __init__(self, _capacity, **kwargs):
        _capacity = int(_capacity)
        for val in 'id_parent', 'x_parent', 'px_parent', 'y_parent', 'py_parent', 'zeta_parent',\
                'delta_parent', 'energy_parent',\
                'id_child', 'x_child', 'px_child', 'y_child', 'py_child', 'zeta_child',\
                'delta_child', 'energy_child',\
                'at_element', 's', 'turn', 'interaction_type':
            kwargs[val] = _capacity
        self.xoinitialize(**kwargs)
        self._capacity = _capacity
        self._row_id = -1


class ParticleInfoData(xo.Struct):
    _capacity       = xo.Int64
    _row_id         = xo.Int64
    particle_id     = xo.Int64[:]
    particle_m      = xo.Float64[:]
    particle_q      = xo.Int64[:]
    particle_z      = xo.Int64[:]
    particle_a      = xo.Int64[:]
    particle_pdgid  = xo.Int64[:]
    
class ParticleInfo(xo.dress(ParticleInfoData)):
    def __init__(self, _capacity, **kwargs):
        _capacity = int(_capacity)
        for val in  'particle_id', 'particle_m', 'particle_q', 'particle_z', 'particle_a', 'particle_pdgid':
            kwargs[val] = _capacity
        self.xoinitialize(**kwargs)
        self._capacity = _capacity
        self._row_id = -1
        