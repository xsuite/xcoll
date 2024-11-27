#!/usr/bin/env python

import pytest
import contextlib
import os
import sys
from pathlib import Path

from xobjects.context import get_test_contexts

path = Path(__file__).parent

def _collect_tests():
    # Create a list to store the collected test names
    collected_tests = []

    class TestCollector:
        @pytest.hookimpl
        def pytest_collection_modifyitems(self, session, config, items):
            # This hook gets called after test collection
            # Add the test names (node IDs) to the collected_tests list
            for item in items:
                collected_tests.append(item.nodeid)

    # Run pytest collection with our custom plugin
    with contextlib.redirect_stdout(open(os.devnull, 'w')):
        pytest.main(["--collect-only"], plugins=[TestCollector()])

    # do not keep GPU contexts as tests in the listing, as the names can change
    collected_tests = [tst for tst in collected_tests
                       if 'Context' not in tst or 'ContextCpu' in tst]

    return set(collected_tests)

current_tests = _collect_tests()
# available_contexts = [ctx.__class__.__name__ for ctx in get_test_contexts()]
with open(path / 'all_tests.list', 'r') as fid:
    expected_tests = set([line.replace('\n', '') for line in fid.readlines()])
only_current = current_tests.difference(expected_tests)
only_expected = expected_tests.difference(current_tests)
if len(only_current) > 0:
    sys.exit("Please run data/store_all_tests.py as there are some new "
           + "tests that are not logged yet:\n" + '\n'.join(only_current))
if len(only_expected) > 0:
    sys.exit("The following tests were expected but not found:\n" + '\n'.join(only_expected) + "Please run data/store_all_tests.py")
