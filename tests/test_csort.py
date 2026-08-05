# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2026.                 #
# ######################################### #

import pytest
import numpy as np

import xobjects as xo
from xobjects.test_helpers import for_all_test_contexts


class SortTest(xo.HybridClass):
    _xofields = {}
    allow_no_prebuilt_kernel = True
    _extra_c_sources = ["""
#include "xobjects/headers/common.h"
#include "xcoll/headers/sort.h"
GPUKERN void SortTest_sort_rows_asc_double(GPUGLMEM double* values,
                            int64_t row_length, int64_t num_rows) {
    VECTORIZE_OVER(row_id, num_rows);
        GPUGLMEM double* row = values + row_id * row_length;
        XC_SORT_ASC(row, row_length);
    END_VECTORIZE;
}
GPUKERN void SortTest_sort_rows_desc_double(GPUGLMEM double* values,
                                int64_t row_length, int64_t num_rows) {
    VECTORIZE_OVER(row_id, num_rows);
        GPUGLMEM double* row = values + row_id * row_length;
        XC_SORT_DESC(row, row_length);
    END_VECTORIZE;
}
GPUKERN void SortTest_sort_rows_asc_int64_t(GPUGLMEM int64_t* values,
                            int64_t row_length, int64_t num_rows) {
    VECTORIZE_OVER(row_id, num_rows);
        GPUGLMEM int64_t* row = values + row_id * row_length;
        XC_SORT_ASC(row, row_length);
    END_VECTORIZE;
}
GPUKERN void SortTest_sort_rows_desc_int64_t(GPUGLMEM int64_t* values,
                                int64_t row_length, int64_t num_rows) {
    VECTORIZE_OVER(row_id, num_rows);
        GPUGLMEM int64_t* row = values + row_id * row_length;
        XC_SORT_DESC(row, row_length);
    END_VECTORIZE;
}
"""]

    _kernels = {
        "sort_rows_asc_double": xo.Kernel(
            c_name="SortTest_sort_rows_asc_double",
            args=[
                xo.Arg(xo.Float64, pointer=True, name="values"),
                xo.Arg(xo.Int64, name="row_length"),
                xo.Arg(xo.Int64, name="num_rows"),
            ],
            n_threads="num_rows",
        ),
        "sort_rows_desc_double": xo.Kernel(
            c_name="SortTest_sort_rows_desc_double",
            args=[
                xo.Arg(xo.Float64, pointer=True, name="values"),
                xo.Arg(xo.Int64, name="row_length"),
                xo.Arg(xo.Int64, name="num_rows"),
            ],
            n_threads="num_rows",
        ),
        "sort_rows_asc_int64_t": xo.Kernel(
            c_name="SortTest_sort_rows_asc_int64_t",
            args=[
                xo.Arg(xo.Int64, pointer=True, name="values"),
                xo.Arg(xo.Int64, name="row_length"),
                xo.Arg(xo.Int64, name="num_rows"),
            ],
            n_threads="num_rows",
        ),
        "sort_rows_desc_int64_t": xo.Kernel(
            c_name="SortTest_sort_rows_desc_int64_t",
            args=[
                xo.Arg(xo.Int64, pointer=True, name="values"),
                xo.Arg(xo.Int64, name="row_length"),
                xo.Arg(xo.Int64, name="num_rows"),
            ],
            n_threads="num_rows",
        ),
    }


@pytest.mark.xcother
@for_all_test_contexts
@pytest.mark.parametrize("descending", [False, True])
@pytest.mark.parametrize("double", [True, False])
def test_fixed_sorting_networks_exhaustive(test_context, descending, double):
    sorter = SortTest(_context=test_context)
    sorter.compile_kernels()

    for length in range(1, 17):
        print(f"Testing length {length} (descending={descending}, "
              f"double={double})")
        values = _all_binary_rows(length, double)
        result = _run_sort(sorter, values, descending, double)
        expected = np.sort(values, axis=1)
        if descending:
            expected = expected[:, ::-1]
        np.testing.assert_array_equal(result, expected)


@pytest.mark.xcother
@for_all_test_contexts
@pytest.mark.parametrize("descending", [False, True])
@pytest.mark.parametrize("double", [True, False])
@pytest.mark.parametrize("length,num_rows",
    [
        (2, 128),
        (7, 128),
        (16, 128),
        (17, 128),
        (32, 128),
        (33, 137),
        (57, 64),
        (200, 79),
        (2000, 34),
        (37520, 31),
    ],
)
def test_sort_random_lengths(test_context, descending, double, length, num_rows):
    sorter = SortTest(_context=test_context)
    sorter.compile_kernels()
    rng = np.random.default_rng(12345 + length)
    if double:
        values = rng.normal(0.0, 10.0, size=(num_rows, length))
    else:
        values = rng.integers(-100, 100, size=(num_rows, length))

    result = _run_sort(sorter, values, descending, double)
    expected = np.sort(values, axis=1)
    if descending:
        expected = expected[:, ::-1]
    np.testing.assert_array_equal(result, expected)


@pytest.mark.xcother
@for_all_test_contexts
@pytest.mark.parametrize("descending", [False, True])
@pytest.mark.parametrize("double", [True, False])
@pytest.mark.parametrize("length", [16, 17, 32, 33, 200])
def test_sort_duplicates_and_patterns(test_context, descending, double, length):
    sorter = SortTest(_context=test_context)
    sorter.compile_kernels()
    rng = np.random.default_rng(9842 + length)
    if double:
        random_duplicates = rng.integers(-4, 5, size=(32, length)
                                         ).astype(np.float64)
        ascending = np.tile(np.arange(length, dtype=np.float64), (4, 1))
        descending_input = ascending[:, ::-1].copy()
        constant = np.full((4, length), 3.25, dtype=np.float64)
    else:
        random_duplicates = rng.integers(-4, 5, size=(32, length))
        ascending = np.tile(np.arange(length, dtype=np.int64), (4, 1))
        descending_input = ascending[:, ::-1].copy()
        constant = np.full((4, length), 3, dtype=np.int64)

    values = np.concatenate([
            random_duplicates,
            ascending,
            descending_input,
            constant,
        ], axis=0,
    )

    result = _run_sort(sorter, values, descending, double)
    expected = np.sort(values, axis=1)
    if descending:
        expected = expected[:, ::-1]
    np.testing.assert_array_equal(result, expected)


def _all_binary_rows(length, double):
    num_rows = 1 << length
    values = np.arange(num_rows, dtype=np.uint64)[:, None]
    shifts = np.arange(length, dtype=np.uint64)[None, :]
    if double:
        return ((values >> shifts) & 1).astype(np.float64)
    else:
        return ((values >> shifts) & 1).astype(np.int64)

def _run_sort(sorter, values, descending, double):
    test_context = sorter._context
    if double:
        values = np.ascontiguousarray(values, dtype=np.float64)
    else:
        values = np.ascontiguousarray(values, dtype=np.int64)
    if values.ndim != 2:
        raise ValueError("values must be a two-dimensional array")
    num_rows, row_length = values.shape
    flat_host = values.reshape(-1)
    flat_context = test_context.nparray_to_context_array(flat_host)

    if descending:
        if double:
            kernel = test_context.kernels['sort_rows_desc_double']
        else:
            kernel = test_context.kernels['sort_rows_desc_int64_t']
    else:
        if double:
            kernel = test_context.kernels['sort_rows_asc_double']
        else:
            kernel = test_context.kernels['sort_rows_asc_int64_t']
    kernel(
        values=flat_context,
        row_length=row_length,
        num_rows=num_rows,
    )
    result = test_context.nparray_from_context_array(flat_context)
    return result.reshape(num_rows, row_length)
