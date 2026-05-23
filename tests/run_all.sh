#!/bin/bash

# Verify all tests
echo "Verifying all tests..."
data/assert_listing.py

# If the first step is successful, continue with the remaining tests in parallel
if [ $? -ne 0 ]
then
    echo "Not all tests are verified. Aborting."
    exit 1
fi

# Prebuild kernels
echo "Regenerating prebuilt kernels..."
xsuite-prebuild r
# "xsuite-prebuild c" to remove kernels
echo "done."

pytest -n 14 -vv
