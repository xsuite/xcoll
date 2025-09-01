# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

import pytest

import xtrack as xt
import xcoll as xc
from xcoll.accessors import XcollAccessor, XcollFamilyAccessor


def test_dict_of_dicts_accessor():
    db = {
        'el1': {'a': 1, 'b': 2, 'd': 10},
        'el2': {'a': 3, 'b': 4, 'd': 10},
        'el3': {'a': 5, 'c': 6, 'd': 10},
    }
    acc = XcollAccessor(db=db, _dbtype='test dict', _typename='test element')

    # The accessor should behave like a dict:
    assert len(acc) == 3
    assert set(acc.names) == {'el1', 'el2', 'el3'}
    assert set(acc.keys()) == {'el1', 'el2', 'el3'}
    assert list(acc.values()) == [db['el1'], db['el2'], db['el3']]
    assert list(acc.items()) == [('el1', db['el1']), ('el2', db['el2']), ('el3', db['el3'])]
    assert acc['el1'] == db['el1']
    assert acc['el2'] == db['el2']
    assert acc['el3'] == db['el3']
    assert 'el1' in acc
    assert 'el2' in acc
    assert 'el3' in acc
    assert 'non_existing' not in acc
    with pytest.raises(KeyError) as err:
        _ = acc['non_existing']
    assert str(err.value) == "'Test element `non_existing` not found in test dict!'"
    with pytest.raises(TypeError) as err:
        acc['non_existing'] = 8
    assert str(err.value) == "'XcollAccessor' object does not support item assignment"

    # Iteration over the accessor should yield the elements (not like a dict, which yields keys):
    for i, nn in enumerate(acc):
        assert nn == db[f'el{i+1}']

    # The accessor should provide access to the element attributes:
    assert set(acc.properties) == {'a', 'b', 'c', 'd'}

    # Accessing an attribute should return a dict of element -> attribute value, omitting
    # elements that do not have the attribute:
    assert acc.a == {'el1': 1, 'el2': 3, 'el3': 5}
    assert acc.b == {'el1': 2, 'el2': 4}
    assert acc.c == {'el3': 6}
    with pytest.raises(AttributeError) as err:
        _ = acc.e
    assert str(err.value) == "Attribute `e` not found in test dict!"

    # If all elements have the same value for an attribute, return that single value:
    assert acc.d == 10

    # Setting an attribute with a dict should set the attribute for the specified elements:
    acc.a = {'el1': 10, 'el2': 30}
    assert acc.a == {'el1': 10, 'el2': 30, 'el3': 5}

    # Setting an attribute with a single value should set it for all elements that have that attribute:
    acc.b = 5000
    assert acc.b == {'el1': 5000, 'el2': 5000}
    assert db['el1']['b'] is 5000
    assert db['el2']['b'] is 5000

    with pytest.raises(AttributeError) as err:
        acc.c = {'el1': -800}
    assert acc.c == {'el3': 6}
    assert str(err.value) == "Attribute `c` not found in test element `el1`!"

    with pytest.raises(KeyError) as err:
        acc.c = {'el7': -800}
    assert acc.c == {'el3': 6}
    assert str(err.value) == "'Test element `el7` not found in test dict!'"


def test_dict_of_instances_accessor():
    class MyClass:
        def __init__(self, **kwargs):
            for k, v in kwargs.items():
                setattr(self, k, v)
    db = {
        'el1': MyClass(a=1, b=2, d=10),
        'el2': MyClass(a=3, b=4, d=10),
        'el3': MyClass(a=5, c=6, d=10),
    }
    acc = XcollAccessor(db=db, _dbtype='test dict', _typename='test element')

    # The accessor should behave like a dict:
    assert len(acc) == 3
    assert set(acc.names) == {'el1', 'el2', 'el3'}
    assert set(acc.keys()) == {'el1', 'el2', 'el3'}
    assert list(acc.values()) == [db['el1'], db['el2'], db['el3']]
    assert list(acc.items()) == [('el1', db['el1']), ('el2', db['el2']), ('el3', db['el3'])]
    assert acc['el1'] == db['el1']
    assert acc['el2'] == db['el2']
    assert acc['el3'] == db['el3']
    assert 'el1' in acc
    assert 'el2' in acc
    assert 'el3' in acc
    assert 'non_existing' not in acc
    with pytest.raises(KeyError) as err:
        _ = acc['non_existing']
    assert str(err.value) == "'Test element `non_existing` not found in test dict!'"
    with pytest.raises(TypeError) as err:
        acc['non_existing'] = 8
    assert str(err.value) == "'XcollAccessor' object does not support item assignment"

    # Iteration over the accessor should yield the elements (not like a dict, which yields keys):
    for i, nn in enumerate(acc):
        assert nn == db[f'el{i+1}']

    # The accessor should provide access to the element attributes:
    assert set(acc.properties) == {'a', 'b', 'c', 'd'}

    # Accessing an attribute should return a dict of element -> attribute value, omitting
    # elements that do not have the attribute:
    assert acc.a == {'el1': 1, 'el2': 3, 'el3': 5}
    assert acc.b == {'el1': 2, 'el2': 4}
    assert acc.c == {'el3': 6}
    with pytest.raises(AttributeError) as err:
        _ = acc.e
    assert str(err.value) == "Attribute `e` not found in test dict!"

    # If all elements have the same value for an attribute, return that single value:
    assert acc.d == 10

    # Setting an attribute with a dict should set the attribute for the specified elements:
    acc.a = {'el1': 10, 'el2': 30}
    assert acc.a == {'el1': 10, 'el2': 30, 'el3': 5}

    # Setting an attribute with a single value should set it for all elements that have that attribute:
    acc.b = 5000
    assert acc.b == {'el1': 5000, 'el2': 5000}
    assert db['el1']['b'] is 5000
    assert db['el2']['b'] is 5000

    with pytest.raises(AttributeError) as err:
        acc.c = {'el1': -800}
    assert acc.c == {'el3': 6}
    assert str(err.value) == "Attribute `c` not found in test element `el1`!"

    with pytest.raises(KeyError) as err:
        acc.c = {'el7': -800}
    assert acc.c == {'el3': 6}
    assert str(err.value) == "'Test element `el7` not found in test dict!'"


def test_instance_of_dicts_accessor():
    class DBClass:
        _el1 = {'a': 1, 'b': 2, 'd': 10}
        _el2 = {'a': 3, 'b': 4, 'd': 10}
        _el3 = {'a': 5, 'c': 6, 'd': 10}
        def get(self, kk):
            if kk == 'el1':
                return self._el1
            elif kk == 'el2':
                return self._el2
            elif kk == 'el3':
                return self._el3
            else:
                raise KeyError

    acc = XcollAccessor(db=DBClass(), _dbtype='test dict', _typename='test element')

    # The accessor should behave like a dict:
    assert len(acc) == 3
    assert set(acc.names) == {'el1', 'el2', 'el3'}
    assert set(acc.keys()) == {'el1', 'el2', 'el3'}
    assert list(acc.values()) == [db['el1'], db['el2'], db['el3']]
    assert list(acc.items()) == [('el1', db['el1']), ('el2', db['el2']), ('el3', db['el3'])]
    assert acc['el1'] == db['el1']
    assert acc['el2'] == db['el2']
    assert acc['el3'] == db['el3']
    assert 'el1' in acc
    assert 'el2' in acc
    assert 'el3' in acc
    assert 'non_existing' not in acc
    with pytest.raises(KeyError) as err:
        _ = acc['non_existing']
    assert str(err.value) == "'Test element `non_existing` not found in test dict!'"
    with pytest.raises(TypeError) as err:
        acc['non_existing'] = 8
    assert str(err.value) == "'XcollAccessor' object does not support item assignment"

    # Iteration over the accessor should yield the elements (not like a dict, which yields keys):
    for i, nn in enumerate(acc):
        assert nn == db[f'el{i+1}']

    # The accessor should provide access to the element attributes:
    assert set(acc.properties) == {'a', 'b', 'c', 'd'}

    # Accessing an attribute should return a dict of element -> attribute value, omitting
    # elements that do not have the attribute:
    assert acc.a == {'el1': 1, 'el2': 3, 'el3': 5}
    assert acc.b == {'el1': 2, 'el2': 4}
    assert acc.c == {'el3': 6}
    with pytest.raises(AttributeError) as err:
        _ = acc.e
    assert str(err.value) == "Attribute `e` not found in test dict!"

    # If all elements have the same value for an attribute, return that single value:
    assert acc.d == 10

    # Setting an attribute with a dict should set the attribute for the specified elements:
    acc.a = {'el1': 10, 'el2': 30}
    assert acc.a == {'el1': 10, 'el2': 30, 'el3': 5}

    # Setting an attribute with a single value should set it for all elements that have that attribute:
    acc.b = 5000
    assert acc.b == {'el1': 5000, 'el2': 5000}
    assert db['el1']['b'] is 5000
    assert db['el2']['b'] is 5000

    with pytest.raises(AttributeError) as err:
        acc.c = {'el1': -800}
    assert acc.c == {'el3': 6}
    assert str(err.value) == "Attribute `c` not found in test element `el1`!"

    with pytest.raises(KeyError) as err:
        acc.c = {'el7': -800}
    assert acc.c == {'el3': 6}
    assert str(err.value) == "'Test element `el7` not found in test dict!'"


