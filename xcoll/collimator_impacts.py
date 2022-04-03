import xobjects as xo

class CollimatorImpactsData(xo.Struct):
#     collimator_name  = xo.String[:]
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
    m_parent      = xo.Float64[:]
    q_parent      = xo.Int64[:]
    a_parent      = xo.Int64[:]
    pdgid_parent  = xo.Int64[:]
    id_child     = xo.Int64[:]
    x_child      = xo.Float64[:]
    px_child     = xo.Float64[:]
    y_child      = xo.Float64[:]
    py_child     = xo.Float64[:]
    zeta_child   = xo.Float64[:]
    delta_child  = xo.Float64[:]
    energy_child = xo.Float64[:]
    m_child      = xo.Float64[:]
    q_child      = xo.Int64[:]
    a_child      = xo.Int64[:]
    pdgid_child  = xo.Int64[:]
    _capacity = xo.Int64
    _row_id   = xo.Int64
    
class CollimatorImpacts(xo.dress(CollimatorImpactsData)):
    def __init__(self, _capacity, **kwargs):
        _capacity = int(_capacity)
        for val in 'id_parent', 'x_parent', 'px_parent', 'y_parent', 'py_parent', 'zeta_parent',\
                'delta_parent', 'energy_parent', 'm_parent', 'q_parent', 'a_parent', 'pdgid_parent',\
                'id_child', 'x_child', 'px_child', 'y_child', 'py_child', 'zeta_child',\
                'delta_child', 'energy_child', 'm_child', 'q_child', 'a_child', 'pdgid_child',\
                's', 'turn', 'interaction_type':
            kwargs[val] = _capacity
        self.xoinitialize(**kwargs)
        self._capacity = _capacity
        self._row_id = -1
