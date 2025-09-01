# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

import pytest

import xtrack as xt
import xcoll as xc
from xcoll.accessors import XcollAccessor, XcollFamilyAccessor


@pytest.mark.parametrize("outer_is_instance", [False, True
                        ], ids=["outer_dict", "outer_instance"])
@pytest.mark.parametrize("inner_is_instance", [False, True
                        ], ids=["inner_dict", "inner_instance"])
@pytest.mark.parametrize("extra_elements", [False, True
                        ], ids=["no_extra_elements", "with_extra_elements"])
def test_accessor(outer_is_instance, inner_is_instance, extra_elements):
    class InnerClass:
        def __init__(self, **kwargs):
            for k, v in kwargs.items():
                setattr(self, k, v)
        def get(self, kk):
            try:
                return getattr(self, kk)
            except AttributeError:
                raise KeyError
        def __str__(self): # For test debugging
            return f'InnerClass({self.__dict__})'

    class OuterClass:
        def __init__(self, **kwargs):
            for k, v in kwargs.items():
                setattr(self, f'_{k}', v)
        def get(self, kk):
            try:
                return getattr(self, f'_{kk}')
            except AttributeError:
                raise KeyError
        def __str__(self): # For test debugging
            return f'OuterClass({self.__dict__})'

    # Create accessor with different combinations of db types:
    inner = InnerClass if inner_is_instance else dict
    outer = OuterClass if outer_is_instance else dict
    if extra_elements:
        db = outer(
            el1=inner(a=1, b=2, d=10),
            el2=inner(a=3, b=4, d=10),
            el3=inner(a=5, c=6, d=10),
            el4=inner(b=7, e=23)
        )
    else:
        db = outer(
            el1=inner(a=1, b=2, d=10),
            el2=inner(a=3, b=4, d=10),
            el3=inner(a=5, c=6, d=10),
        )
    if outer_is_instance:
        _dbtype = 'test class'
        acc = XcollAccessor(db=db, names=['el1', 'el2', 'el3'], _dbtype=_dbtype, _typename='test element')
    else:
        _dbtype = 'test dict'
        if extra_elements:
            acc = XcollAccessor(db=db, names=['el1', 'el2', 'el3'], _dbtype=_dbtype, _typename='test element')
        else:
            acc = XcollAccessor(db=db, _dbtype=_dbtype, _typename='test element')

    # The accessor should behave like a dict:
    assert len(acc) == 3
    assert set(acc.names) == {'el1', 'el2', 'el3'}
    assert set(acc.keys()) == {'el1', 'el2', 'el3'}
    assert all(isinstance(val, inner) for val in acc.values())
    assert all(isinstance(val, inner) for _, val in acc.items())
    assert {nn for nn, _ in acc.items()} == {'el1', 'el2', 'el3'}
    if not inner_is_instance:
        assert list(acc.values()) == [{'a': 1, 'b': 2, 'd': 10},
                                    {'a': 3, 'b': 4, 'd': 10},
                                    {'a': 5, 'c': 6, 'd': 10}]
        assert list(acc.items()) == [('el1', {'a': 1, 'b': 2, 'd': 10}),
                                    ('el2', {'a': 3, 'b': 4, 'd': 10}),
                                    ('el3', {'a': 5, 'c': 6, 'd': 10})]
    assert acc['el1'] is db.get('el1')
    assert acc['el2'] is db.get('el2')
    assert acc['el3'] is db.get('el3')
    if inner_is_instance:
        assert {nn: acc['el1'].get(nn) for nn in acc['el1'].__dict__.keys()} == {'a': 1, 'b': 2, 'd': 10}
        assert {nn: acc['el2'].get(nn) for nn in acc['el2'].__dict__.keys()} == {'a': 3, 'b': 4, 'd': 10}
        assert {nn: acc['el3'].get(nn) for nn in acc['el3'].__dict__.keys()} == {'a': 5, 'c': 6, 'd': 10}
    else:
        assert acc['el1'] == {'a': 1, 'b': 2, 'd': 10}
        assert acc['el2'] == {'a': 3, 'b': 4, 'd': 10}
        assert acc['el3'] == {'a': 5, 'c': 6, 'd': 10}
    assert 'el1' in acc
    assert 'el2' in acc
    assert 'el3' in acc

    # Accessing a non-existing element or elements that are not represented
    # (i.e. not in names) should raise an error:
    assert 'el4' not in acc
    assert 'non_existing' not in acc
    with pytest.raises(KeyError) as err:
        acc['el4']
    assert str(err.value) == f"'Test element `el4` not found in {_dbtype}!'"
    with pytest.raises(KeyError) as err:
        acc['non_existing']
    assert str(err.value) == f"'Test element `non_existing` not found in {_dbtype}!'"
    with pytest.raises(TypeError) as err:
        acc['el4'] = 8
    assert str(err.value) == "'XcollAccessor' object does not support item assignment"

    if extra_elements:
        assert isinstance(db.get('el4'), inner) # Should not have been modified
    else:
        with pytest.raises(KeyError):           # Should not have been added
            _ = db['el4'] if outer == dict else db.get('el4')
    with pytest.raises(TypeError) as err:
        acc['non_existing'] = 8
    assert str(err.value) == "'XcollAccessor' object does not support item assignment"
    with pytest.raises(KeyError):  # Should not have been added
        _ = db['non_existing'] if outer == dict else db.get('non_existing')

    # Iteration over the accessor should yield the elements (not like a dict, which yields keys):
    for i, nn in enumerate(acc):
        assert nn == db.get(f'el{i+1}')

    # Accessing an attribute should return a dict of element -> attribute value, omitting
    # elements that do not have the attribute:
    assert acc.a == {'el1': 1, 'el2': 3, 'el3': 5}
    assert acc.b == {'el1': 2, 'el2': 4}
    assert acc.c == {'el3': 6}
    with pytest.raises(AttributeError) as err:
        _ = acc.e
    assert str(err.value) == f"Attribute `e` not found in {_dbtype}!"

    # If all elements have the same value for an attribute, return that single value:
    assert acc.d == 10

    # Setting an attribute with a dict should set the attribute for the specified elements:
    acc.a = {'el1': 10, 'el2': 30}
    assert acc.a == {'el1': 10, 'el2': 30, 'el3': 5}
    # And the attribute a in the original db should have been updated:
    assert db.get('el1').get('a') == 10
    assert db.get('el2').get('a') == 30
    assert db.get('el3').get('a') == 5
    # But not the other attributes:
    assert db.get('el1').get('b') == 2
    assert db.get('el2').get('b') == 4
    if extra_elements:
        assert db.get('el4').get('b') == 7
    assert db.get('el3').get('c') == 6
    assert db.get('el1').get('d') == 10
    assert db.get('el2').get('d') == 10
    assert db.get('el3').get('d') == 10
    if extra_elements:
        assert db.get('el4').get('e') == 23

    # Setting an attribute with a single value should set it for all elements that have that attribute:
    acc.b = 5000
    assert acc.b == {'el1': 5000, 'el2': 5000}
    # And the attribute b in the original db should have been updated:
    assert db.get('el1').get('b') == 5000
    assert db.get('el2').get('b') == 5000
    # But not the other attributes nor b in elements that are not represented by the accessor:
    if extra_elements:
        assert db.get('el4').get('b') == 7
    assert db.get('el1').get('a') == 10
    assert db.get('el2').get('a') == 30
    assert db.get('el3').get('a') == 5
    assert db.get('el3').get('c') == 6
    assert db.get('el1').get('d') == 10
    assert db.get('el2').get('d') == 10
    assert db.get('el3').get('d') == 10
    if extra_elements:
        assert db.get('el4').get('e') == 23

    # Testing invalid attribute setting:
    with pytest.raises(AttributeError) as err:
        acc.c = {'el1': -800}
    assert acc.c == {'el3': 6}
    assert str(err.value) == "Attribute `c` not found in test element `el1`!"
    with pytest.raises(KeyError) as err:
        acc.c = {'el7': -800}
    assert acc.c == {'el3': 6}
    assert str(err.value) == f"'Test element `el7` not found in {_dbtype}!'"
    with pytest.raises(AttributeError) as err:
        acc.e = 80
    assert str(err.value) == f"Attribute `e` not found in {_dbtype}!"
    # The db should not have been modified:
    assert db.get('el1').get('a') == 10
    assert db.get('el2').get('a') == 30
    assert db.get('el3').get('a') == 5
    assert db.get('el1').get('b') == 5000
    assert db.get('el2').get('b') == 5000
    if extra_elements:
        assert db.get('el4').get('b') == 7
    assert db.get('el3').get('c') == 6
    assert db.get('el1').get('d') == 10
    assert db.get('el2').get('d') == 10
    assert db.get('el3').get('d') == 10
    if extra_elements:
        assert db.get('el4').get('e') == 23
