def pytest_collection_modifyitems(items):
    MODULE_ORDER = ["test__list_tests", "test__regenerate_kernels"]
    module_mapping = {item: item.module.__name__ for item in items}

    sorted_items = items[:]

    for module in reversed(MODULE_ORDER):
        sorted_items = [it for it in sorted_items if module_mapping[it] == module] + [
            it for it in sorted_items if module_mapping[it] != module
        ]
    items[:] = sorted_items
