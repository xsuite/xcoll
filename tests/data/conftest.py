# import pytest

# def pytest_collection_modifyitems(items):
#     MODULE_ORDER = ["test__list_tests", "test__regenerate_kernels"]
#     module_mapping = {item: item.module.__name__ for item in items}

#     sorted_items = items[:]

#     for module in reversed(MODULE_ORDER):
#         sorted_items = [it for it in sorted_items if module_mapping[it] == module] + [
#             it for it in sorted_items if module_mapping[it] != module
#         ]
#     items[:] = sorted_items


# def pytest_addoption(parser):
#     parser.addoption("--skip-first-tests", action='store_true', dest='skip_first_tests', default=False)


# @pytest.hookimpl(tryfirst=True)
# def pytest_runtestloop(session):
#     cli_args = session.config.invocation_params.args
#     if '--collect-only' in cli_args:
#         # We are running inside the test__listing test
#         return

#     # Tests to be ran first
#     run_first = ['test__listing', 'test__regenerate_kernels']
#     run_first_tests = [item for item in session.items if item.module.__name__ in run_first]

#     if run_first_tests:
#         for item in run_first_tests:
#             session.items.remove(item)

#         if '--skip-first-tests' not in cli_args:
#             # A separate pytest.main call is needed to run the first tests (to avoid parallelisation etc)
#             exitstatus = pytest.main([str(item.fspath) for item in run_first_tests])

#             # If the first tests fail, stop execution
#             if exitstatus != 0:
#                 return exitstatus

# from xdist.scheduler import LoadScopeScheduling
# # This makes sure only remaining tests run in parallel
# def pytest_xdist_make_scheduler(config, log):
#     return LoadScopeScheduling(config, log)
