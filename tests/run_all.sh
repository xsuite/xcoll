#!/usr/bin/bash

# Verify all tests
echo "Verifying all tests..."
data/assert_listing.py

# If the first step is successful, continue with the remaining tests in parallel
if [ $? -eq 0 ]
then
    echo "done."
    # Prebuild kernels
    echo "Regenerating prebuilt kernels..."
    xsuite-prebuild r
    # "xsuite-prebuild c" to remove kernels
    echo "done."
    pytest -n 10 -vv --html="pytest_results.html" --self-contained-html
else
    echo "Not all tests are verified. Aborting."
    exit 1
fi





