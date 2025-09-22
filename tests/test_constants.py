# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

import pytest

from xcoll.xoconstants import (Constants, constant, group, ConstantSpec, GroupSpec, _IntMixin,
                               _FloatMixin, _TupleMixin, _XO_CONST_GLOBAL_NAME_REG,
                               _XO_CONST_GLOBAL_UNIQUE_VALUE_REG, _XO_CONST_GLOBAL_TYPE_REG,
                               _XO_CONST_GLOBAL_COUNT_REG)
import test_constants_pkg as tt
import test_constants_pkg.type1 as tt1
import test_constants_pkg.type4 as tt4
import test_constants_pkg.other_type as tto


def test_registries():
    base_pkg_name = 'test_constants_pkg'
    assert _XO_CONST_GLOBAL_NAME_REG['ENABLED'] == f'{base_pkg_name}.type1'
    assert _XO_CONST_GLOBAL_NAME_REG['BIG'] == f'{base_pkg_name}.type1'
    assert _XO_CONST_GLOBAL_NAME_REG['MASS'] == f'{base_pkg_name}.type1'
    assert _XO_CONST_GLOBAL_NAME_REG['COOL'] == f'{base_pkg_name}.type1'
    assert _XO_CONST_GLOBAL_NAME_REG['VAR1'] == f'{base_pkg_name}.type1'
    assert _XO_CONST_GLOBAL_NAME_REG['VAR2'] == f'{base_pkg_name}.type1'
    assert _XO_CONST_GLOBAL_NAME_REG['VAR3'] == f'{base_pkg_name}.type1'
    assert _XO_CONST_GLOBAL_NAME_REG['GR1'] == f'{base_pkg_name}.type1'
    assert _XO_CONST_GLOBAL_NAME_REG['VAR4'] == f'{base_pkg_name}.type1'
    assert _XO_CONST_GLOBAL_NAME_REG['VAR5'] == f'{base_pkg_name}.type1'
    assert _XO_CONST_GLOBAL_NAME_REG['VAR6'] == f'{base_pkg_name}.type1'
    assert _XO_CONST_GLOBAL_NAME_REG['THING_SOME_VAL'] == f'{base_pkg_name}.type4'
    assert _XO_CONST_GLOBAL_NAME_REG['THING_OTHER'] == f'{base_pkg_name}.type4'
    assert _XO_CONST_GLOBAL_NAME_REG['THING_INT'] == f'{base_pkg_name}.type4'
    assert _XO_CONST_GLOBAL_NAME_REG['YAY_NO_META'] == f'{base_pkg_name}.type4'
    assert _XO_CONST_GLOBAL_NAME_REG['YAY_NO_META_2'] == f'{base_pkg_name}.type4'
    assert _XO_CONST_GLOBAL_NAME_REG['YAY_NO_META_3'] == f'{base_pkg_name}.type4'
    assert _XO_CONST_GLOBAL_NAME_REG['GR2'] == f'{base_pkg_name}.type4'
    assert _XO_CONST_GLOBAL_NAME_REG['VAR7'] == f'{base_pkg_name}.type4'
    assert _XO_CONST_GLOBAL_NAME_REG['VAR8'] == f'{base_pkg_name}.type4'
    assert _XO_CONST_GLOBAL_NAME_REG['VAR9'] == f'{base_pkg_name}.type4'
    assert _XO_CONST_GLOBAL_NAME_REG['VAR10'] == f'{base_pkg_name}.type4'
    assert _XO_CONST_GLOBAL_NAME_REG['VAR11'] == f'{base_pkg_name}.type4'
    assert _XO_CONST_GLOBAL_NAME_REG['VAR12'] == f'{base_pkg_name}.type4'
    assert _XO_CONST_GLOBAL_NAME_REG['GR3'] == f'{base_pkg_name}.type4'
    assert _XO_CONST_GLOBAL_NAME_REG['VAS1'] == f'{base_pkg_name}.other_type'
    assert _XO_CONST_GLOBAL_NAME_REG['VAS2'] == f'{base_pkg_name}.other_type'
    assert _XO_CONST_GLOBAL_NAME_REG['GRO'] == f'{base_pkg_name}.other_type'
    assert _XO_CONST_GLOBAL_NAME_REG['GROO'] == f'{base_pkg_name}.other_type'

    assert _XO_CONST_GLOBAL_UNIQUE_VALUE_REG['unique_type'] == {
        2: ('VAR1', 'test_constants_pkg.type1'),
        -3: ('VAR2', 'test_constants_pkg.type1'),
        9: ('VAR3', 'test_constants_pkg.type1'),
        0.734: ('THING_SOME_VAL', 'test_constants_pkg.type4'),
        1.23: ('THING_OTHER', 'test_constants_pkg.type4'),
        12.0: ('THING_INT', 'test_constants_pkg.type4'),
        9.81: ('YAY_NO_META', 'test_constants_pkg.type4'),
        23.0: ('YAY_NO_META_2', 'test_constants_pkg.type4'),
        1: ('YAY_NO_META_3', 'test_constants_pkg.type4')}
    assert _XO_CONST_GLOBAL_UNIQUE_VALUE_REG['other_type'] == {
        2: ('VAS1', 'test_constants_pkg.other_type'),
        -3: ('VAS2', 'test_constants_pkg.other_type')}

    assert _XO_CONST_GLOBAL_TYPE_REG['type'] == [tt.type1.TestType1]
    assert _XO_CONST_GLOBAL_TYPE_REG['unique_type'] == [tt.type1.TestType2, tt.type4.TestType4]
    assert _XO_CONST_GLOBAL_TYPE_REG['multi_type'] == [tt.type1.TestType3, tt.type4.TestType5, tt.type4.TestType6]
    assert _XO_CONST_GLOBAL_TYPE_REG['other_type'] == [tt.other_type.TestOtherType]

    assert _XO_CONST_GLOBAL_COUNT_REG['XF_TYPE_H_guard'] == 0
    assert _XO_CONST_GLOBAL_COUNT_REG['UU_UNIQUE_TYPE_H_guard'] == 2
    assert _XO_CONST_GLOBAL_COUNT_REG['XF_MULTI_TYPE_H_guard'] == 0
    assert _XO_CONST_GLOBAL_COUNT_REG['unique_type_unnamed'] == 1
    assert _XO_CONST_GLOBAL_COUNT_REG['RR_MULTI_TYPE_H_guard'] == 3
    assert _XO_CONST_GLOBAL_COUNT_REG['other_type_unnamed'] == 0
    assert _XO_CONST_GLOBAL_COUNT_REG['XFREOO_OTHER_TYPE_H_guard'] == 0


def test_class():
    assert isinstance(tt._T1, type) and issubclass(tt._T1, Constants)
    assert isinstance(tt._T2, type) and issubclass(tt._T2, Constants)
    assert isinstance(tt._T3, type) and issubclass(tt._T3, Constants)
    assert isinstance(tt._T4, type) and issubclass(tt._T4, Constants)
    assert isinstance(tt._T5, type) and issubclass(tt._T5, Constants)
    assert isinstance(tt._T6, type) and issubclass(tt._T6, Constants)
    assert isinstance(tt._TO, type) and issubclass(tt._TO, Constants)
    assert tt._T1.__category__ == "type"
    assert tt._T2.__category__ == "unique_type"
    assert tt._T3.__category__ == "multi_type"
    assert tt._T4.__category__ == "unique_type"
    assert tt._T5.__category__ == "multi_type"
    assert tt._T6.__category__ == "multi_type"
    assert tt._TO.__category__ == "other_type"
    assert tt._T1.__reverse__ is None
    assert tt._T2.__reverse__ == "unique"
    assert tt._T3.__reverse__ == "multi"
    assert tt._T4.__reverse__ == "unique"
    assert tt._T5.__reverse__ == "multi"
    assert tt._T6.__reverse__ == "multi"
    assert tt._TO.__reverse__ == "unique"


def test_class_constants():
    assert hasattr(tt._T1, 'ENABLED') and isinstance(tt._T1.ENABLED, ConstantSpec)
    assert hasattr(tt._T1, 'BIG') and isinstance(tt._T1.BIG, ConstantSpec)
    assert hasattr(tt._T1, 'MASS') and isinstance(tt._T1.MASS, ConstantSpec)
    assert hasattr(tt._T2, 'VAR1') and isinstance(tt._T2.VAR1, ConstantSpec)
    assert hasattr(tt._T2, 'VAR2') and isinstance(tt._T2.VAR2, ConstantSpec)
    assert hasattr(tt._T2, 'VAR3') and isinstance(tt._T2.VAR3, ConstantSpec)
    assert hasattr(tt._T3, 'VAR4') and isinstance(tt._T3.VAR4, ConstantSpec)
    assert hasattr(tt._T3, 'VAR5') and isinstance(tt._T3.VAR5, ConstantSpec)
    assert hasattr(tt._T3, 'VAR6') and isinstance(tt._T3.VAR6, ConstantSpec)
    assert hasattr(tt._T4, 'THING_SOME_VAL') and isinstance(tt._T4.THING_SOME_VAL, ConstantSpec)
    assert hasattr(tt._T4, 'THING_OTHER') and isinstance(tt._T4.THING_OTHER, ConstantSpec)
    assert hasattr(tt._T4, 'THING_INT') and isinstance(tt._T4.THING_INT, ConstantSpec)
    assert hasattr(tt._T4, 'YAY_NO_META') and isinstance(tt._T4.YAY_NO_META, ConstantSpec)
    assert hasattr(tt._T4, 'YAY_NO_META_2') and isinstance(tt._T4.YAY_NO_META_2, ConstantSpec)
    assert hasattr(tt._T4, 'YAY_NO_META_3') and isinstance(tt._T4.YAY_NO_META_3, ConstantSpec)
    assert hasattr(tt._T5, 'VAR7') and isinstance(tt._T5.VAR7, ConstantSpec)
    assert hasattr(tt._T5, 'VAR8') and isinstance(tt._T5.VAR8, ConstantSpec)
    assert hasattr(tt._T5, 'VAR9') and isinstance(tt._T5.VAR9, ConstantSpec)
    assert hasattr(tt._T6, 'VAR10') and isinstance(tt._T6.VAR10, ConstantSpec)
    assert hasattr(tt._T6, 'VAR11') and isinstance(tt._T6.VAR11, ConstantSpec)
    assert hasattr(tt._T6, 'VAR12') and isinstance(tt._T6.VAR12, ConstantSpec)
    assert hasattr(tt._TO, 'VAS1') and isinstance(tt._TO.VAS1, ConstantSpec)
    assert hasattr(tt._TO, 'VAS2') and isinstance(tt._TO.VAS2, ConstantSpec)

    assert tt._T1.ENABLED.name == 'ENABLED' and tt._T1.ENABLED.py_value is True
    assert tt._T1.ENABLED.info == 'Important info' and tt._T1.ENABLED.c_name == 'OLA_ENABLED'
    assert tt._T1.BIG.name == 'BIG' and tt._T1.BIG.py_value == 2**64 - 1
    assert tt._T1.BIG.info == 'Big int' and tt._T1.BIG.c_name == 'XF_BIG'
    assert tt._T1.MASS.name == 'MASS' and tt._T1.MASS.py_value == 0.938
    assert tt._T1.MASS.info == 'Mass in GeV/c^2' and tt._T1.MASS.c_name == 'XF_MASS'
    assert tt._T2.VAR1.name == 'VAR1' and tt._T2.VAR1.py_value == 2
    assert tt._T2.VAR1.info == 'Var 1' and tt._T2.VAR1.c_name == 'UU_VAR1'
    assert tt._T2.VAR2.name == 'VAR2' and tt._T2.VAR2.py_value == -3
    assert tt._T2.VAR2.info == 'Var 2' and tt._T2.VAR2.c_name == 'UU_VAR2'
    assert tt._T2.VAR3.name == 'VAR3' and tt._T2.VAR3.py_value == 9
    assert tt._T2.VAR3.info == 'Var 3' and tt._T2.VAR3.c_name == 'UU_VAR3'
    assert tt._T3.VAR4.name == 'VAR4' and tt._T3.VAR4.py_value == 4
    assert tt._T3.VAR4.info == 'Var 4' and tt._T3.VAR4.c_name == 'XF_VAR4'
    assert tt._T3.VAR5.name == 'VAR5' and tt._T3.VAR5.py_value == -5
    assert tt._T3.VAR5.info == 'Var 5' and tt._T3.VAR5.c_name == 'XF_VAR5'
    assert tt._T3.VAR6.name == 'VAR6' and tt._T3.VAR6.py_value == 4
    assert tt._T3.VAR6.info == 'Var 6' and tt._T3.VAR6.c_name == 'XF_VAR6'
    assert tt._T4.THING_SOME_VAL.name == 'THING_SOME_VAL' and tt._T4.THING_SOME_VAL.py_value == 0.734
    assert tt._T4.THING_SOME_VAL.info == 'Type4 thing with some value.' and tt._T4.THING_SOME_VAL.c_name == 'UU_THING_SOME_VAL'
    assert tt._T4.THING_OTHER.name == 'THING_OTHER' and tt._T4.THING_OTHER.py_value == 1.23
    assert tt._T4.THING_OTHER.info == 'Another type4 thing.' and tt._T4.THING_OTHER.c_name == 'UU_THING_OTHER'
    assert tt._T4.THING_INT.name == 'THING_INT' and tt._T4.THING_INT.py_value == 12.0
    assert tt._T4.THING_INT.info == 'Integer type4 thing -> has to be float' and tt._T4.THING_INT.c_name == 'UU_THING_INT'
    assert tt._T4.YAY_NO_META.name == 'YAY_NO_META' and tt._T4.YAY_NO_META.py_value == 9.81
    assert tt._T4.YAY_NO_META.info == '' and tt._T4.YAY_NO_META.c_name == 'UU_YAY_NO_META'
    assert tt._T4.YAY_NO_META_2.name == 'YAY_NO_META_2' and tt._T4.YAY_NO_META_2.py_value == 23.0
    assert tt._T4.YAY_NO_META_2.info == '' and tt._T4.YAY_NO_META_2.c_name == 'UU_YAY_NO_META_2'
    assert tt._T4.YAY_NO_META_3.name == 'YAY_NO_META_3' and tt._T4.YAY_NO_META_3.py_value is True
    assert tt._T4.YAY_NO_META_3.info == '' and tt._T4.YAY_NO_META_3.c_name == 'UU_YAY_NO_META_3'
    assert tt._T5.VAR7.name == 'VAR7' and tt._T5.VAR7.py_value == 4
    assert tt._T5.VAR7.info == 'Var 7' and tt._T5.VAR7.c_name == 'RR_VAR7'
    assert tt._T5.VAR8.name == 'VAR8' and tt._T5.VAR8.py_value == -5
    assert tt._T5.VAR8.info == 'Var 8' and tt._T5.VAR8.c_name == 'RR_VAR8'
    assert tt._T5.VAR9.name == 'VAR9' and tt._T5.VAR9.py_value == 10
    assert tt._T5.VAR9.info == 'Var 9' and tt._T5.VAR9.c_name == 'RR_VAR9'
    assert tt._T6.VAR10.name == 'VAR10' and tt._T6.VAR10.py_value == 4
    assert tt._T6.VAR10.info == 'Var 10' and tt._T6.VAR10.c_name == 'RR_VAR10'
    assert tt._T6.VAR11.name == 'VAR11' and tt._T6.VAR11.py_value == 88
    assert tt._T6.VAR11.info == 'Var 11' and tt._T6.VAR11.c_name == 'RR_VAR11'
    assert tt._T6.VAR12.name == 'VAR12' and tt._T6.VAR12.py_value == 10
    assert tt._T6.VAR12.info == 'Var 12' and tt._T6.VAR12.c_name == 'RR_VAR12'
    assert tt._TO.VAS1.name == 'VAS1' and tt._TO.VAS1.py_value == 2
    assert tt._TO.VAS1.info == 'Vas 1' and tt._TO.VAS1.c_name == 'XFREOO_VAS1'
    assert tt._TO.VAS2.name == 'VAS2' and tt._TO.VAS2.py_value == -3
    assert tt._TO.VAS2.info == 'Vas 2' and tt._TO.VAS2.c_name == 'XFREOO_VAS2'


def test_module_constants():
    assert hasattr(tt1, 'ENABLED') and isinstance(tt1.ENABLED, bool)
    assert hasattr(tt1, 'BIG') and isinstance(tt1.BIG, _IntMixin)
    assert hasattr(tt1, 'MASS') and isinstance(tt1.MASS, _FloatMixin)
    assert hasattr(tt1, 'VAR1') and isinstance(tt1.VAR1, _IntMixin)
    assert hasattr(tt1, 'VAR2') and isinstance(tt1.VAR2, _IntMixin)
    assert hasattr(tt1, 'VAR3') and isinstance(tt1.VAR3, _IntMixin)
    assert hasattr(tt1, 'VAR4') and isinstance(tt1.VAR4, _IntMixin)
    assert hasattr(tt1, 'VAR5') and isinstance(tt1.VAR5, _IntMixin)
    assert hasattr(tt1, 'VAR6') and isinstance(tt1.VAR6, _IntMixin)
    assert hasattr(tt4, 'THING_SOME_VAL') and isinstance(tt4.THING_SOME_VAL, _FloatMixin)
    assert hasattr(tt4, 'THING_OTHER') and isinstance(tt4.THING_OTHER, _FloatMixin)
    assert hasattr(tt4, 'THING_INT') and isinstance(tt4.THING_INT, _FloatMixin)
    assert hasattr(tt4, 'YAY_NO_META') and isinstance(tt4.YAY_NO_META, _FloatMixin)
    assert hasattr(tt4, 'YAY_NO_META_2') and isinstance(tt4.YAY_NO_META_2, _FloatMixin)
    assert hasattr(tt4, 'YAY_NO_META_3') and isinstance(tt4.YAY_NO_META_3, bool)
    assert hasattr(tt4, 'VAR7') and isinstance(tt4.VAR7, _IntMixin)
    assert hasattr(tt4, 'VAR8') and isinstance(tt4.VAR8, _IntMixin)
    assert hasattr(tt4, 'VAR9') and isinstance(tt4.VAR9, _IntMixin)
    assert hasattr(tt4, 'VAR10') and isinstance(tt4.VAR10, _IntMixin)
    assert hasattr(tt4, 'VAR11') and isinstance(tt4.VAR11, _IntMixin)
    assert hasattr(tt4, 'VAR12') and isinstance(tt4.VAR12, _IntMixin)
    assert hasattr(tto, 'VAS1') and isinstance(tto.VAS1, _IntMixin)
    assert hasattr(tto, 'VAS2') and isinstance(tto.VAS2, _IntMixin)

    assert tt1.ENABLED is True
    assert tt1.BIG == 2**64 - 1
    assert tt1.MASS == 0.938
    assert tt1.VAR1 == 2
    assert tt1.VAR2 == -3
    assert tt1.VAR3 == 9
    assert tt1.VAR4 == 4
    assert tt1.VAR5 == -5
    assert tt1.VAR6 == 4
    assert tt4.THING_SOME_VAL == 0.734
    assert tt4.THING_OTHER == 1.23
    assert tt4.THING_INT == 12.0
    assert tt4.YAY_NO_META == 9.81
    assert tt4.YAY_NO_META_2 == 23.0
    assert tt4.YAY_NO_META_3 is True
    assert tt4.VAR7 == 4
    assert tt4.VAR8 == -5
    assert tt4.VAR9 == 10
    assert tt4.VAR10 == 4
    assert tt4.VAR11 == 88
    assert tt4.VAR12 == 10
    assert tto.VAS1 == 2
    assert tto.VAS2 == -3


def test_package_constants():
    # Not imported into package namespace
    assert not hasattr(tt, 'ENABLED')
    assert not hasattr(tt, 'BIG')
    assert not hasattr(tt, 'MASS')
    assert not hasattr(tt, 'VAR1')
    assert not hasattr(tt, 'VAR2')
    assert not hasattr(tt, 'VAR3')
    assert not hasattr(tt, 'VAR4')
    assert not hasattr(tt, 'VAR5')
    assert not hasattr(tt, 'VAR6')
    assert not hasattr(tt, 'THING_SOME_VAL')
    assert not hasattr(tt, 'THING_OTHER')
    assert not hasattr(tt, 'THING_INT')
    assert not hasattr(tt, 'YAY_NO_META')
    assert not hasattr(tt, 'YAY_NO_META_2')
    assert not hasattr(tt, 'YAY_NO_META_3')
    assert not hasattr(tt, 'VAR7')
    assert not hasattr(tt, 'VAR8')
    assert not hasattr(tt, 'VAR9')
    assert not hasattr(tt, 'VAR10')
    assert not hasattr(tt, 'VAR11')
    assert not hasattr(tt, 'VAR12')
    assert not hasattr(tt, 'VAS1')
    assert not hasattr(tt, 'VAS2')


def test_class_constant_lists():
    assert tt._T1._map == {'ENABLED': 1, 'BIG': 18446744073709551615, 'MASS': 0.938}
    assert tt._T2._map == {'VAR1': 2, 'VAR2': -3, 'VAR3': 9}
    assert tt._T3._map == {'VAR4': 4, 'VAR5': -5, 'VAR6': 4}
    assert tt._T4._map == {'THING_SOME_VAL': 0.734, 'THING_OTHER': 1.23, 'THING_INT': 12.0,
                           'YAY_NO_META': 9.81, 'YAY_NO_META_2': 23.0, 'YAY_NO_META_3': 1}
    assert tt._T5._map == {'VAR7': 4, 'VAR8': -5, 'VAR9': 10}
    assert tt._T6._map == {'VAR10': 4, 'VAR11': 88, 'VAR12': 10}
    assert tt._TO._map == {'VAS1': 2, 'VAS2': -3}


def test_module_constant_lists():
    assert not hasattr(tt1, '_map')
    assert not hasattr(tt4, '_map')
    assert not hasattr(tto, '_map')
    assert dict(tt1.types) == {'ENABLED': 1, 'BIG': 18446744073709551615, 'MASS': 0.938}
    assert dict(tt1.unique_types) == {'VAR1': 2, 'VAR2': -3, 'VAR3': 9}
    assert dict(tt1.multi_types) == {'VAR4': 4, 'VAR5': -5, 'VAR6': 4}
    assert dict(tt4.unique_types) == {'THING_SOME_VAL': 0.734, 'THING_OTHER': 1.23, 'THING_INT': 12.0,
                                      'YAY_NO_META': 9.81, 'YAY_NO_META_2': 23.0, 'YAY_NO_META_3': 1}
    assert dict(tt4.multi_types)== {'VAR7': 4, 'VAR8': -5, 'VAR9': 10,
                                    'VAR10': 4, 'VAR11': 88, 'VAR12': 10}
    assert dict(tto.other_types) == {'VAS1': 2, 'VAS2': -3}


def test_package_constant_lists():
    assert not hasattr(tt, '_map')
    assert dict(tt.types) == {'ENABLED': 1, 'BIG': 18446744073709551615, 'MASS': 0.938}
    assert dict(tt.unique_types) == {'VAR1': 2, 'VAR2': -3, 'VAR3': 9, 'THING_SOME_VAL': 0.734,
                                      'THING_OTHER': 1.23, 'THING_INT': 12.0, 'YAY_NO_META': 9.81,
                                      'YAY_NO_META_2': 23.0, 'YAY_NO_META_3': 1}
    assert dict(tt.multi_types) == {'VAR4': 4, 'VAR5': -5, 'VAR6': 4,
                                    'VAR7': 4, 'VAR8': -5, 'VAR9': 10,
                                    'VAR10': 4, 'VAR11': 88, 'VAR12': 10}
    assert dict(tto.other_types) == {'VAS1': 2, 'VAS2': -3}


def test_class_groups():
    assert hasattr(tt._T1, 'COOL') and isinstance(tt._T1.COOL, GroupSpec)
    assert hasattr(tt._T2, 'GR1') and isinstance(tt._T2.GR1, GroupSpec)
    assert hasattr(tt._T4, 'GR2') and isinstance(tt._T4.GR2, GroupSpec)
    assert hasattr(tt._T6, 'GR3') and isinstance(tt._T6.GR3, GroupSpec)
    assert hasattr(tt._TO, 'GRO') and isinstance(tt._TO.GRO, GroupSpec)
    assert hasattr(tt._TO, 'GROO') and isinstance(tt._TO.GROO, GroupSpec)

    assert tt._T1.COOL.name == 'COOL' and tt._T1.COOL.py_values == (0.938, 18446744073709551615)
    assert tt._T1.COOL.info == "Isn't this cool?" and tt._T1.COOL.names == ('MASS', 'BIG')
    assert tt._T2.GR1.name == 'GR1' and tt._T2.GR1.py_values == (-3, 0.938, 2, 18446744073709551615)
    assert tt._T2.GR1.info == '' and tt._T2.GR1.names == ('VAR2', 'MASS', 'VAR1', 'BIG')
    assert tt._T4.GR2.name == 'GR2' and tt._T4.GR2.py_values == (0.734, 9.81, 23.0)
    assert tt._T4.GR2.info == '' and tt._T4.GR2.names == ('THING_SOME_VAL', 'UNIQUE_TYPE_UNNAMED_0', 'UNIQUE_TYPE_UNNAMED_1')
    assert tt._T6.GR3.name == 'GR3' and tt._T6.GR3.py_values == (4, 4, 10)
    assert tt._T6.GR3.info == '' and tt._T6.GR3.names == ('VAR4', 'VAR7', 'VAR12')
    assert tt._TO.GRO.name == 'GRO' and tt._TO.GRO.py_values == (-3, 0.734, 0.938, True, 2, 2, 9.81, 23.0, 18446744073709551615)
    assert tt._TO.GRO.info == 'Fancyyyy' and tt._TO.GRO.names == ('VAR2', 'THING_SOME_VAL', 'MASS', 'YAY_NO_META_3', 'VAR1',
                                                                  'VAS1', 'UNIQUE_TYPE_UNNAMED_0', 'UNIQUE_TYPE_UNNAMED_1',
                                                                  'BIG')
    assert tt._TO.GROO.name == 'GROO' and tt._TO.GROO.py_values == (-3, 0.734, 0.938, 1, 1.23, 2, 2, 9.81, 23.0, 23.0,
                                                                    18446744073709551615)
    assert tt._TO.GROO.info == '' and tt._TO.GROO.names == ('VAR2', 'THING_SOME_VAL', 'MASS', 'OTHER_TYPE_UNNAMED_0',
                                                            'THING_OTHER', 'VAR1', 'VAS1', 'UNIQUE_TYPE_UNNAMED_0',
                                                            'UNIQUE_TYPE_UNNAMED_1', 'YAY_NO_META_2', 'BIG')


def test_module_groups():
    assert hasattr(tt1, 'COOL') and isinstance(tt1.COOL, _TupleMixin)
    assert hasattr(tt1, 'GR1') and isinstance(tt1.GR1, _TupleMixin)
    assert hasattr(tt4, 'GR2') and isinstance(tt4.GR2, _TupleMixin)
    assert hasattr(tt4, 'GR3') and isinstance(tt4.GR3, _TupleMixin)
    assert hasattr(tto, 'GRO') and isinstance(tto.GRO, _TupleMixin)
    assert hasattr(tto, 'GROO') and isinstance(tto.GROO, _TupleMixin)

    assert tt1.COOL == (0.938, 18446744073709551615)
    assert tt1.COOL.info == "Isn't this cool?" and tt1.COOL.names == ('MASS', 'BIG')
    assert tt1.GR1 == (-3, 0.938, 2, 18446744073709551615)
    assert tt1.GR1.info == '' and tt1.GR1.names == ('VAR2', 'MASS', 'VAR1', 'BIG')
    assert tt4.GR2 == (0.734, 9.81, 23.0)
    assert tt4.GR2.info == '' and tt4.GR2.names == ('THING_SOME_VAL', 'UNIQUE_TYPE_UNNAMED_0', 'UNIQUE_TYPE_UNNAMED_1')
    assert tt4.GR3 == (4, 4, 10)
    assert tt4.GR3.info == '' and tt4.GR3.names == ('VAR4', 'VAR7', 'VAR12')
    assert tto.GRO == (-3, 0.734, 0.938, True, 2, 2, 9.81, 23.0, 18446744073709551615)
    assert tto.GRO.info == 'Fancyyyy' and tto.GRO.names == ('VAR2', 'THING_SOME_VAL', 'MASS', 'YAY_NO_META_3', 'VAR1',
                                                             'VAS1', 'UNIQUE_TYPE_UNNAMED_0', 'UNIQUE_TYPE_UNNAMED_1',
                                                             'BIG')
    assert tto.GROO == (-3, 0.734, 0.938, 1, 1.23, 2, 2, 9.81, 23.0, 23.0, 18446744073709551615)
    assert tto.GROO.info == '' and tto.GROO.names == ('VAR2', 'THING_SOME_VAL', 'MASS', 'OTHER_TYPE_UNNAMED_0',
                                                      'THING_OTHER', 'VAR1', 'VAS1', 'UNIQUE_TYPE_UNNAMED_0',
                                                      'UNIQUE_TYPE_UNNAMED_1', 'YAY_NO_META_2', 'BIG')


def test_package_groups():
    # Not imported into package namespace
    assert not hasattr(tt, 'COOL')
    assert not hasattr(tt, 'GR1')
    assert not hasattr(tt, 'GR2')
    assert not hasattr(tt, 'GR3')
    assert not hasattr(tt, 'GRO')
    assert not hasattr(tt, 'GROO')


def test_class_group_lists():
    assert tt._T1._groups == {'COOL': (0.938, 18446744073709551615)}
    assert tt._T2._groups == {'GR1': (-3, 0.938, 2, 18446744073709551615)}
    assert tt._T4._groups == {'GR2': (0.734, 9.81, 23.0)}
    assert tt._T6._groups == {'GR3': (4, 4, 10)}
    assert tt._TO._groups == {'GRO': (-3, 0.734, 0.938, True, 2, 2, 9.81, 23.0, 18446744073709551615),
                              'GROO': (-3, 0.734, 0.938, 1, 1.23, 2, 2, 9.81, 23.0, 23.0, 18446744073709551615)}
    assert not hasattr(tt._T3, '_groups')
    assert not hasattr(tt._T5, '_groups')


def test_module_group_lists():
    assert not hasattr(tt1, '_groups')
    assert not hasattr(tt4, '_groups')
    assert not hasattr(tto, '_groups')
    assert dict(tt1.type_groups) == {'COOL': (0.938, 18446744073709551615)}
    assert dict(tt1.unique_type_groups) == {'GR1': (-3, 0.938, 2, 18446744073709551615)}
    assert dict(tt4.unique_type_groups) == {'GR2': (0.734, 9.81, 23.0)}
    assert dict(tt4.multi_type_groups) == {'GR3': (4, 4, 10)}
    assert dict(tto.other_type_groups) == {'GRO': (-3, 0.734, 0.938, True, 2, 2, 9.81, 23.0, 18446744073709551615),
                                           'GROO': (-3, 0.734, 0.938, 1, 1.23, 2, 2, 9.81, 23.0, 23.0, 18446744073709551615)}


def test_package_group_lists():
    assert not hasattr(tt, '_groups')
    assert dict(tt.type_groups) == {'COOL': (0.938, 18446744073709551615)}
    assert dict(tt.unique_type_groups) == {'GR1': (-3, 0.938, 2, 18446744073709551615),
                                           'GR2': (0.734, 9.81, 23.0)}
    assert dict(tt.multi_type_groups) == {'GR3': (4, 4, 10)}
    assert dict(tt.other_type_groups) == {'GRO': (-3, 0.734, 0.938, True, 2, 2, 9.81, 23.0, 18446744073709551615),
                                          'GROO': (-3, 0.734, 0.938, 1, 1.23, 2, 2, 9.81, 23.0, 23.0, 18446744073709551615)}


def test_class_constant_meta():
    assert tt._T1._meta == {'ENABLED': {'value': 1, 'c_name': 'OLA_ENABLED', 'info': 'Important info'},
                            'BIG': {'value': 18446744073709551615, 'c_name': 'XF_BIG', 'info': 'Big int'},
                            'MASS': {'value': 0.938, 'c_name': 'XF_MASS', 'info': 'Mass in GeV/c^2'}}
    assert tt._T2._meta == {'VAR1': {'value': 2, 'c_name': 'UU_VAR1', 'info': 'Var 1'},
                            'VAR2': {'value': -3, 'c_name': 'UU_VAR2', 'info': 'Var 2'},
                            'VAR3': {'value': 9, 'c_name': 'UU_VAR3', 'info': 'Var 3'}}
    assert tt._T3._meta == {'VAR4': {'value': 4, 'c_name': 'XF_VAR4', 'info': 'Var 4'},
                            'VAR5': {'value': -5, 'c_name': 'XF_VAR5', 'info': 'Var 5'},
                            'VAR6': {'value': 4, 'c_name': 'XF_VAR6', 'info': 'Var 6'}}
    assert tt._T4._meta == {'THING_SOME_VAL': {'value': 0.734, 'c_name': 'UU_THING_SOME_VAL', 'info': 'Type4 thing with some value.'},
                            'THING_OTHER': {'value': 1.23, 'c_name': 'UU_THING_OTHER', 'info': 'Another type4 thing.'},
                            'THING_INT': {'value': 12.0,'c_name': 'UU_THING_INT', 'info': 'Integer type4 thing -> has to be float'},
                            'YAY_NO_META': {'value': 9.81, 'c_name': 'UU_YAY_NO_META', 'info': ''},
                            'YAY_NO_META_2': {'value': 23.0, 'c_name': 'UU_YAY_NO_META_2', 'info': ''},
                            'YAY_NO_META_3': {'value': 1, 'c_name': 'UU_YAY_NO_META_3', 'info': ''}}
    assert tt._T5._meta == {'VAR7': {'value': 4, 'c_name': 'RR_VAR7', 'info': 'Var 7'},
                            'VAR8': {'value': -5, 'c_name': 'RR_VAR8', 'info': 'Var 8'},
                            'VAR9': {'value': 10, 'c_name': 'RR_VAR9', 'info': 'Var 9'}}
    assert tt._T6._meta == {'VAR10': {'value': 4, 'c_name': 'RR_VAR10', 'info': 'Var 10'},
                            'VAR11': {'value': 88, 'c_name': 'RR_VAR11', 'info': 'Var 11'},
                            'VAR12': {'value': 10, 'c_name': 'RR_VAR12', 'info': 'Var 12'}}
    assert tt._TO._meta == {'VAS1': {'value': 2, 'c_name': 'XFREOO_VAS1', 'info': 'Vas 1'},
                            'VAS2': {'value': -3, 'c_name': 'XFREOO_VAS2', 'info': 'Vas 2'}}


def test_module_constant_meta():
    assert not hasattr(tt1, '_meta')
    assert not hasattr(tt4, '_meta')
    assert not hasattr(tto, '_meta')
    assert dict(tt1.types_meta) == {'ENABLED': {'value': 1, 'c_name': 'OLA_ENABLED', 'info': 'Important info'},
                                    'BIG': {'value': 18446744073709551615, 'c_name': 'XF_BIG', 'info': 'Big int'},
                                    'MASS': {'value': 0.938, 'c_name': 'XF_MASS', 'info': 'Mass in GeV/c^2'}}
    assert dict(tt1.unique_types_meta) == {'VAR1': {'value': 2, 'c_name': 'UU_VAR1', 'info': 'Var 1'},
                                           'VAR2': {'value': -3, 'c_name': 'UU_VAR2', 'info': 'Var 2'},
                                           'VAR3': {'value': 9, 'c_name': 'UU_VAR3', 'info': 'Var 3'}}
    assert dict(tt1.multi_types_meta) == {'VAR4': {'value': 4, 'c_name': 'XF_VAR4', 'info': 'Var 4'},
                                          'VAR5': {'value': -5, 'c_name': 'XF_VAR5', 'info': 'Var 5'},
                                          'VAR6': {'value': 4, 'c_name': 'XF_VAR6', 'info': 'Var 6'}}
    assert dict(tt4.unique_types_meta) == {'THING_SOME_VAL': {'value': 0.734, 'c_name': 'UU_THING_SOME_VAL', 'info': 'Type4 thing with some value.'},
                                           'THING_OTHER': {'value': 1.23, 'c_name': 'UU_THING_OTHER', 'info': 'Another type4 thing.'},
                                           'THING_INT': {'value': 12.0,'c_name': 'UU_THING_INT', 'info': 'Integer type4 thing -> has to be float'},
                                           'YAY_NO_META': {'value': 9.81, 'c_name': 'UU_YAY_NO_META', 'info': ''},
                                           'YAY_NO_META_2': {'value': 23.0, 'c_name': 'UU_YAY_NO_META_2', 'info': ''},
                                           'YAY_NO_META_3': {'value': 1, 'c_name': 'UU_YAY_NO_META_3', 'info': ''}}
    assert dict(tt4.multi_types_meta) == {'VAR7': {'value': 4, 'c_name': 'RR_VAR7', 'info': 'Var 7'},
                                          'VAR8': {'value': -5, 'c_name': 'RR_VAR8', 'info': 'Var 8'},
                                          'VAR9': {'value': 10, 'c_name': 'RR_VAR9', 'info': 'Var 9'},
                                          'VAR10': {'value': 4, 'c_name': 'RR_VAR10', 'info': 'Var 10'},
                                          'VAR11': {'value': 88, 'c_name': 'RR_VAR11', 'info': 'Var 11'},
                                          'VAR12': {'value': 10, 'c_name': 'RR_VAR12', 'info': 'Var 12'}}
    assert dict(tto.other_types_meta) == {'VAS1': {'value': 2, 'c_name': 'XFREOO_VAS1', 'info': 'Vas 1'},
                                          'VAS2': {'value': -3, 'c_name': 'XFREOO_VAS2', 'info': 'Vas 2'}}


def test_package_constant_meta():
    assert not hasattr(tt, '_meta')
    assert dict(tt.types_meta) == {'ENABLED': {'value': 1, 'c_name': 'OLA_ENABLED', 'info': 'Important info'},
                                   'BIG': {'value': 18446744073709551615, 'c_name': 'XF_BIG', 'info': 'Big int'},
                                   'MASS': {'value': 0.938, 'c_name': 'XF_MASS', 'info': 'Mass in GeV/c^2'}}
    assert dict(tt.unique_types_meta) == {'VAR1': {'value': 2, 'c_name': 'UU_VAR1', 'info': 'Var 1'},
                                          'VAR2': {'value': -3, 'c_name': 'UU_VAR2', 'info': 'Var 2'},
                                          'VAR3': {'value': 9, 'c_name': 'UU_VAR3', 'info': 'Var 3'},
                                          'THING_SOME_VAL': {'value': 0.734, 'c_name': 'UU_THING_SOME_VAL', 'info': 'Type4 thing with some value.'},
                                          'THING_OTHER': {'value': 1.23, 'c_name': 'UU_THING_OTHER', 'info': 'Another type4 thing.'},
                                          'THING_INT': {'value': 12.0,'c_name': 'UU_THING_INT', 'info': 'Integer type4 thing -> has to be float'},
                                          'YAY_NO_META': {'value': 9.81, 'c_name': 'UU_YAY_NO_META', 'info': ''},
                                          'YAY_NO_META_2': {'value': 23.0, 'c_name': 'UU_YAY_NO_META_2', 'info': ''},
                                          'YAY_NO_META_3': {'value': 1, 'c_name': 'UU_YAY_NO_META_3', 'info': ''}}
    assert dict(tt.multi_types_meta) == {'VAR4': {'value': 4, 'c_name': 'XF_VAR4', 'info': 'Var 4'},
                                         'VAR5': {'value': -5, 'c_name': 'XF_VAR5', 'info': 'Var 5'},
                                         'VAR6': {'value': 4, 'c_name': 'XF_VAR6', 'info': 'Var 6'},
                                         'VAR7': {'value': 4, 'c_name': 'RR_VAR7', 'info': 'Var 7'},
                                         'VAR8': {'value': -5, 'c_name': 'RR_VAR8', 'info': 'Var 8'},
                                         'VAR9': {'value': 10, 'c_name': 'RR_VAR9', 'info': 'Var 9'},
                                         'VAR10': {'value': 4, 'c_name': 'RR_VAR10', 'info': 'Var 10'},
                                         'VAR11': {'value': 88, 'c_name': 'RR_VAR11', 'info': 'Var 11'},
                                         'VAR12': {'value': 10, 'c_name': 'RR_VAR12', 'info': 'Var 12'}}
    assert tto.other_types_meta == {'VAS1': {'value': 2, 'c_name': 'XFREOO_VAS1', 'info': 'Vas 1'},
                                    'VAS2': {'value': -3, 'c_name': 'XFREOO_VAS2', 'info': 'Vas 2'}}


def test_class_constant_names():
    assert tt._T2._names == {2: 'VAR1', -3: 'VAR2', 9: 'VAR3'}
    assert tt._T3._names  == {4: ('VAR4', 'VAR6'), -5: ('VAR5',)}
    assert tt._T4._names == {0.734: 'THING_SOME_VAL', 1.23: 'THING_OTHER', 12.0: 'THING_INT',
                                              9.81: 'YAY_NO_META', 23.0: 'YAY_NO_META_2', 1: 'YAY_NO_META_3'}
    assert tt._T5._names == {4: ('VAR7',), -5: ('VAR8',), 10: ('VAR9',)}
    assert tt._T6._names == {4: ('VAR10',), 88: ('VAR11',), 10: ('VAR12',)}
    assert tt._TO._names  == {2: 'VAS1', -3: 'VAS2'}
    assert not hasattr(tt._T1, '_names')


def test_module_constant_names():
    assert not hasattr(tt1, '_names')
    assert not hasattr(tt4, '_names')
    assert not hasattr(tto, '_names')
    assert dict(tt1.unique_type_names) == {2: 'VAR1', -3: 'VAR2', 9: 'VAR3'}
    assert dict(tt1.multi_type_names)  == {4: ('VAR4', 'VAR6'), -5: ('VAR5',)}
    assert dict(tt4.unique_type_names) == {0.734: 'THING_SOME_VAL', 1.23: 'THING_OTHER', 12.0: 'THING_INT',
                                           9.81: 'YAY_NO_META', 23.0: 'YAY_NO_META_2', 1: 'YAY_NO_META_3'}
    assert dict(tt4.multi_type_names) == {4: ('VAR7', 'VAR10'), -5: ('VAR8',), 10: ('VAR9', 'VAR12'),
                                          88: ('VAR11',)}
    assert dict(tto.other_type_names)  == {2: 'VAS1', -3: 'VAS2'}


def test_package_constant_names():
    assert not hasattr(tt, '_names')
    assert dict(tt.unique_type_names) == {2: 'VAR1', -3: 'VAR2', 9: 'VAR3', 0.734: 'THING_SOME_VAL',
                                          1.23: 'THING_OTHER', 12.0: 'THING_INT', 9.81: 'YAY_NO_META',
                                          23.0: 'YAY_NO_META_2', 1: 'YAY_NO_META_3'}
    assert dict(tt.multi_type_names)  == {4: ('VAR4', 'VAR6', 'VAR7', 'VAR10'), -5: ('VAR5', 'VAR8'),
                                          10: ('VAR9', 'VAR12'), 88: ('VAR11',)}
    assert dict(tt.other_type_names)  == {2: 'VAS1', -3: 'VAS2'}


def test_class_constant_sources():
    assert tt._T1._src == """#ifndef XF_TYPE_H_I0
#define XF_TYPE_H_I0
  #define OLA_ENABLED                                                 1     // Important info
  #define XF_BIG                                   18446744073709551615     // Big int
  #define XF_MASS                                                 0.938     // Mass in GeV/c^2
#endif /* XF_TYPE_H_I0 */
"""
    assert tt._T2._src == """#ifndef UU_UNIQUE_TYPE_H_I0
#define UU_UNIQUE_TYPE_H_I0
  #define UU_VAR1                                                     2     // Var 1
  #define UU_VAR2                                                    -3     // Var 2
  #define UU_VAR3                                                     9     // Var 3
#endif /* UU_UNIQUE_TYPE_H_I0 */
"""
    assert tt._T3._src == """#ifndef XF_MULTI_TYPE_H_I0
#define XF_MULTI_TYPE_H_I0
  #define XF_VAR4                                                     4     // Var 4
  #define XF_VAR5                                                    -5     // Var 5
  #define XF_VAR6                                                     4     // Var 6
#endif /* XF_MULTI_TYPE_H_I0 */
"""
    assert tt._T4._src == """#ifndef UU_UNIQUE_TYPE_H_I1
#define UU_UNIQUE_TYPE_H_I1
  #define UU_THING_SOME_VAL                                       0.734     // Type4 thing with some value.
  #define UU_THING_OTHER                                           1.23     // Another type4 thing.
  #define UU_THING_INT                                             12.0     // Integer type4 thing -> has to be float
  #define UU_YAY_NO_META                                           9.81     // 
  #define UU_YAY_NO_META_2                                         23.0     // 
  #define UU_YAY_NO_META_3                                            1     // 
#endif /* UU_UNIQUE_TYPE_H_I1 */
"""
    assert tt._T5._src == """#ifndef RR_MULTI_TYPE_H_I0
#define RR_MULTI_TYPE_H_I0
  #define RR_VAR7                                                     4     // Var 7
  #define RR_VAR8                                                    -5     // Var 8
  #define RR_VAR9                                                    10     // Var 9
#endif /* RR_MULTI_TYPE_H_I0 */
"""
    assert tt._T6._src == """#ifndef RR_MULTI_TYPE_H_I1
#define RR_MULTI_TYPE_H_I1
  #define RR_VAR10                                                    4     // Var 10
  #define RR_VAR11                                                   88     // Var 11
  #define RR_VAR12                                                   10     // Var 12
#endif /* RR_MULTI_TYPE_H_I1 */
"""
    assert tt._TO._src == """#ifndef XFREOO_OTHER_TYPE_H_I0
#define XFREOO_OTHER_TYPE_H_I0
  #define XFREOO_VAS1                                                 2     // Vas 1
  #define XFREOO_VAS2                                                -3     // Vas 2
#endif /* XFREOO_OTHER_TYPE_H_I0 */
"""


def test_module_constant_sources():
    assert tt1.types_src == """#ifndef XF_TYPE_H_I0
#define XF_TYPE_H_I0
  #define OLA_ENABLED                                                 1     // Important info
  #define XF_BIG                                   18446744073709551615     // Big int
  #define XF_MASS                                                 0.938     // Mass in GeV/c^2
#endif /* XF_TYPE_H_I0 */
"""
    assert tt1.unique_types_src == """#ifndef UU_UNIQUE_TYPE_H_I0
#define UU_UNIQUE_TYPE_H_I0
  #define UU_VAR1                                                     2     // Var 1
  #define UU_VAR2                                                    -3     // Var 2
  #define UU_VAR3                                                     9     // Var 3
#endif /* UU_UNIQUE_TYPE_H_I0 */
"""
    assert tt1.multi_types_src == """#ifndef XF_MULTI_TYPE_H_I0
#define XF_MULTI_TYPE_H_I0
  #define XF_VAR4                                                     4     // Var 4
  #define XF_VAR5                                                    -5     // Var 5
  #define XF_VAR6                                                     4     // Var 6
#endif /* XF_MULTI_TYPE_H_I0 */
"""
    assert tt4.unique_types_src == """#ifndef UU_UNIQUE_TYPE_H_I1
#define UU_UNIQUE_TYPE_H_I1
  #define UU_THING_SOME_VAL                                       0.734     // Type4 thing with some value.
  #define UU_THING_OTHER                                           1.23     // Another type4 thing.
  #define UU_THING_INT                                             12.0     // Integer type4 thing -> has to be float
  #define UU_YAY_NO_META                                           9.81     // 
  #define UU_YAY_NO_META_2                                         23.0     // 
  #define UU_YAY_NO_META_3                                            1     // 
#endif /* UU_UNIQUE_TYPE_H_I1 */
"""
    assert tt4.multi_types_src == """#ifndef RR_MULTI_TYPE_H_I2
#define RR_MULTI_TYPE_H_I2
  #ifndef RR_MULTI_TYPE_H_I0
  #define RR_MULTI_TYPE_H_I0
    #define RR_VAR7                                                     4     // Var 7
    #define RR_VAR8                                                    -5     // Var 8
    #define RR_VAR9                                                    10     // Var 9
  #endif /* RR_MULTI_TYPE_H_I0 */
  #ifndef RR_MULTI_TYPE_H_I1
  #define RR_MULTI_TYPE_H_I1
    #define RR_VAR10                                                    4     // Var 10
    #define RR_VAR11                                                   88     // Var 11
    #define RR_VAR12                                                   10     // Var 12
  #endif /* RR_MULTI_TYPE_H_I1 */
#endif /* RR_MULTI_TYPE_H_I2 */
"""
    assert tto.other_types_src == """#ifndef XFREOO_OTHER_TYPE_H_I0
#define XFREOO_OTHER_TYPE_H_I0
  #define XFREOO_VAS1                                                 2     // Vas 1
  #define XFREOO_VAS2                                                -3     // Vas 2
#endif /* XFREOO_OTHER_TYPE_H_I0 */
"""


def test_package_constant_sources():
    assert tt.types_src == """#ifndef XF_TYPE_H_I0
#define XF_TYPE_H_I0
  #define OLA_ENABLED                                                 1     // Important info
  #define XF_BIG                                   18446744073709551615     // Big int
  #define XF_MASS                                                 0.938     // Mass in GeV/c^2
#endif /* XF_TYPE_H_I0 */
"""
    assert tt.unique_types_src == """#ifndef UU_UNIQUE_TYPE_H_I2
#define UU_UNIQUE_TYPE_H_I2
  #ifndef UU_UNIQUE_TYPE_H_I0
  #define UU_UNIQUE_TYPE_H_I0
    #define UU_VAR1                                                     2     // Var 1
    #define UU_VAR2                                                    -3     // Var 2
    #define UU_VAR3                                                     9     // Var 3
  #endif /* UU_UNIQUE_TYPE_H_I0 */
  #ifndef UU_UNIQUE_TYPE_H_I1
  #define UU_UNIQUE_TYPE_H_I1
    #define UU_THING_SOME_VAL                                       0.734     // Type4 thing with some value.
    #define UU_THING_OTHER                                           1.23     // Another type4 thing.
    #define UU_THING_INT                                             12.0     // Integer type4 thing -> has to be float
    #define UU_YAY_NO_META                                           9.81     // 
    #define UU_YAY_NO_META_2                                         23.0     // 
    #define UU_YAY_NO_META_3                                            1     // 
  #endif /* UU_UNIQUE_TYPE_H_I1 */
#endif /* UU_UNIQUE_TYPE_H_I2 */
"""
    assert tt.multi_types_src == """#ifndef RR_MULTI_TYPE_H_I3
#define RR_MULTI_TYPE_H_I3
  #ifndef XF_MULTI_TYPE_H_I0
  #define XF_MULTI_TYPE_H_I0
    #define XF_VAR4                                                     4     // Var 4
    #define XF_VAR5                                                    -5     // Var 5
    #define XF_VAR6                                                     4     // Var 6
  #endif /* XF_MULTI_TYPE_H_I0 */
  #ifndef RR_MULTI_TYPE_H_I0
  #define RR_MULTI_TYPE_H_I0
    #define RR_VAR7                                                     4     // Var 7
    #define RR_VAR8                                                    -5     // Var 8
    #define RR_VAR9                                                    10     // Var 9
  #endif /* RR_MULTI_TYPE_H_I0 */
  #ifndef RR_MULTI_TYPE_H_I1
  #define RR_MULTI_TYPE_H_I1
    #define RR_VAR10                                                    4     // Var 10
    #define RR_VAR11                                                   88     // Var 11
    #define RR_VAR12                                                   10     // Var 12
  #endif /* RR_MULTI_TYPE_H_I1 */
#endif /* RR_MULTI_TYPE_H_I3 */
"""
    assert tt.other_types_src == """#ifndef XFREOO_OTHER_TYPE_H_I0
#define XFREOO_OTHER_TYPE_H_I0
  #define XFREOO_VAS1                                                 2     // Vas 1
  #define XFREOO_VAS2                                                -3     // Vas 2
#endif /* XFREOO_OTHER_TYPE_H_I0 */
"""


def test_consistency_and_name_clash():
    with pytest.raises(TypeError, match="FailThingType1 must inherit from Constants directly; "
                       "no further subclassing allowed."):
        class FailThingType1(tt._T1):
            __category__ = "type"
            __reverse__  = "unique"

    with pytest.raises(ValueError, match="Category 'type' already defined with reverse None "
                       "by test_constants_pkg.type1.TestType1. Please make __reverse__ consistent."):
        class FailThingType2(Constants):
            __category__ = "type"
            __reverse__  = "unique"

    with pytest.raises(ValueError, match="Category 'type' already defined with plural 'types' "
                       "by test_constants_pkg.type1.TestType1. Please make __plural__ consistent."):
        class FailThingType3(Constants):
            __category__ = "type"
            __reverse__  = None
            __plural__   = "yikes"

    with pytest.raises(ValueError, match="Constant name 'MASS' already defined in module "
                       "test_constants_pkg.type1; global uniqueness required."):
        class FailThingType4(Constants):
            __category__ = "type"
            __reverse__  = None
            __c_prefix__ = "XF"

            VAL1 = constant(5, "Feature enabled")
            VAL2 = constant(2**64 - 1, "Big int")
            MASS = constant(0.938, "GeV/c^2")   # Double definition in other table

    with pytest.raises(ValueError, match="Value 2 for 'VAL17' in test_constants already used by 'VAR1' "
                   "in test_constants_pkg.type1; values must be globally unique when __reverse__ == 'unique'."):
        class FailThingType5(Constants):
            __category__ = "unique_type"
            __reverse__  = "unique"
            __c_prefix__ = "XF"

            VAL17 = constant(2, "Feature enabled") # Val 2 taken by TestType2

    with pytest.raises(TypeError, match="Mixed int/float values are not allowed when __reverse__ == 'unique'."):
        class FailThingType6(Constants):
            __category__ = "unique_type_bis"
            __reverse__  = "unique"
            __c_prefix__ = "BB"

            VAL18 = constant(2, "Int")
            VAL19 = constant(3.14, "Float")

    with pytest.raises(KeyError, match=f"Group 'GR5' references unknown constant 'VAL22'"):
        class FailThingType7(Constants):
            __category__ = "unique_type_tris"
            __reverse__  = "unique"
            __c_prefix__ = "BB"

            VAL20 = constant(2., "Int")
            VAL21 = constant(3.14, "Float")
            GR5   = group("VAL20", "VAL21", "VAL22")
