import pytest
import contextlib
import os


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

    return set(collected_tests)


def test_listing():
    current_tests = _collect_tests()
    with open('all_tests.list', 'r') as fid:
        expected_tests = set([line.replace('\n', '') for line in fid.readlines()])
    only_current = current_tests.difference(expected_tests)
    only_expected = expected_tests.difference(current_tests)
    if len(only_current) > 0:
        pytest.exit(f"Please run store_all_tests.py as there are some new "
                  + f"tests that are not logged yet:{only_current}.")
    if len(only_expected) > 0:
        pytest.exit(f"The following tests were expected but not found: {only_expected}.")
