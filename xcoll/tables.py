import xobjects as xo
import xtrack as xt

class CollimatorImpacts(xo.HybridClass):

    _xofields = {
        '_index':            xt.RecordIndex,
        'at_element':        xo.Int64[:],
        'at_turn':           xo.Int64[:],
        's':                 xo.Float64[:],
        'interaction_type':  xo.Int64[:],
        'parent_id':      xo.Int64[:],
        'parent_x':       xo.Float64[:],
        'parent_px':      xo.Float64[:],
        'parent_y':       xo.Float64[:],
        'parent_py':      xo.Float64[:],
        'parent_zeta':    xo.Float64[:],
        'parent_delta':   xo.Float64[:],
        'parent_energy':  xo.Float64[:],
        'parent_mass':    xo.Float64[:],
        'parent_charge':  xo.Int64[:],
        'parent_z':       xo.Int64[:],
        'parent_a':       xo.Int64[:],
        'parent_pdgid':   xo.Int64[:],
        'child_id':      xo.Int64[:],
        'child_x':       xo.Float64[:],
        'child_px':      xo.Float64[:],
        'child_y':       xo.Float64[:],
        'child_py':      xo.Float64[:],
        'child_zeta':    xo.Float64[:],
        'child_delta':   xo.Float64[:],
        'child_energy':  xo.Float64[:],
        'child_mass':    xo.Float64[:],
        'child_charge':  xo.Int64[:],
        'child_z':       xo.Int64[:],
        'child_a':       xo.Int64[:],
        'child_pdgid':   xo.Int64[:],
    }

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
