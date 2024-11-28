#!/usr/bin/env python

import pytest
import contextlib
import os
from pathlib import Path


path = Path(__file__).parent


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
os.chdir(path.parent)
with contextlib.redirect_stdout(open(os.devnull, 'w')):
    pytest.main(["--collect-only"], plugins=[TestCollector()])

# do not keep GPU contexts as tests in the listing, as the names can change
collected_tests = [tst for tst in collected_tests
                    if 'Context' not in tst or 'ContextCpu' in tst]

with open(path / 'all_tests.list', 'w') as fid:
    for line in collected_tests:
        fid.write(f"{line}\n")
