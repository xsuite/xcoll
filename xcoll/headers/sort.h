// copyright ############################### #
// This file is part of the Xcoll package.   #
// Copyright (c) CERN, 2024.                 #
// ######################################### #

#ifndef XCOLL_GEOM_SORT_H
#define XCOLL_GEOM_SORT_H

#ifdef XO_CONTEXT_CPU
#include <stdint.h>  // for int64_t etc
#endif  // XO_CONTEXT_CPU

#include "xobjects/headers/common.h"
#include "xcoll/headers/helpers.h"

// #pragma GCC diagnostic push
// #pragma GCC diagnostic ignored "-Warray-bounds"


// Fast methods
// ------------

/*
 * Fixed-size sorting networks.
 *
 * OP must have the interface OP(d, i, j), placing d[i] and d[j] in the
 * requested order. XC_ASORT and XC_DSORT therefore use exactly the same
 * networks.
 */
#define XC_SORT_NETWORK_3(OP, d) do {      \
    OP(d, 0, 1); OP(d, 1, 2); OP(d, 0, 1); \
} while (0)

#define XC_SORT_NETWORK_4(OP, d) do {      \
    OP(d, 0, 1); OP(d, 2, 3); OP(d, 0, 2); \
    OP(d, 1, 3); OP(d, 1, 2);              \
} while (0)

#define XC_SORT_NETWORK_5(OP, d) do {      \
    OP(d, 0, 1); OP(d, 3, 4); OP(d, 2, 4); \
    OP(d, 2, 3); OP(d, 1, 4); OP(d, 0, 3); \
    OP(d, 0, 2); OP(d, 1, 3); OP(d, 1, 2); \
} while (0)

#define XC_SORT_NETWORK_6(OP, d) do {      \
    OP(d, 1, 2); OP(d, 4, 5); OP(d, 0, 2); \
    OP(d, 3, 5); OP(d, 0, 1); OP(d, 3, 4); \
    OP(d, 1, 4); OP(d, 0, 3); OP(d, 2, 5); \
    OP(d, 1, 3); OP(d, 2, 4); OP(d, 2, 3); \
} while (0)

#define XC_SORT_NETWORK_7(OP, d) do {      \
    OP(d, 1, 2); OP(d, 3, 4); OP(d, 5, 6); \
    OP(d, 0, 2); OP(d, 3, 5); OP(d, 4, 6); \
    OP(d, 0, 1); OP(d, 4, 5); OP(d, 2, 6); \
    OP(d, 0, 4); OP(d, 1, 5); OP(d, 0, 3); \
    OP(d, 2, 5); OP(d, 1, 3); OP(d, 2, 4); \
    OP(d, 2, 3);                           \
} while (0)

#define XC_SORT_NETWORK_8(OP, d) do {                   \
    OP(d, 0, 1); OP(d, 2, 3); OP(d, 4, 5); OP(d, 6, 7); \
    OP(d, 0, 2); OP(d, 1, 3); OP(d, 4, 6); OP(d, 5, 7); \
    OP(d, 1, 2); OP(d, 5, 6); OP(d, 0, 4); OP(d, 3, 7); \
    OP(d, 1, 5); OP(d, 2, 6); OP(d, 1, 4); OP(d, 3, 6); \
    OP(d, 2, 4); OP(d, 3, 5); OP(d, 3, 4);              \
} while (0)

#define XC_SORT_NETWORK_9(OP, d) do {                                \
    OP(d, 0, 3); OP(d, 1, 7); OP(d, 2, 5); OP(d, 4, 8); OP(d, 0, 7); \
    OP(d, 2, 4); OP(d, 3, 8); OP(d, 5, 6); OP(d, 0, 2); OP(d, 1, 3); \
    OP(d, 4, 5); OP(d, 7, 8); OP(d, 1, 4); OP(d, 3, 6); OP(d, 5, 7); \
    OP(d, 0, 1); OP(d, 2, 4); OP(d, 3, 5); OP(d, 6, 8); OP(d, 2, 3); \
    OP(d, 4, 5); OP(d, 6, 7); OP(d, 1, 2); OP(d, 3, 4); OP(d, 5, 6); \
} while (0)

#define XC_SORT_NETWORK_10(OP, d) do {                               \
    OP(d, 0, 8); OP(d, 1, 9); OP(d, 2, 7); OP(d, 3, 5); OP(d, 4, 6); \
    OP(d, 0, 2); OP(d, 1, 4); OP(d, 5, 8); OP(d, 7, 9); OP(d, 0, 3); \
    OP(d, 2, 4); OP(d, 5, 7); OP(d, 6, 9); OP(d, 0, 1); OP(d, 3, 6); \
    OP(d, 8, 9); OP(d, 1, 5); OP(d, 2, 3); OP(d, 4, 8); OP(d, 6, 7); \
    OP(d, 1, 2); OP(d, 3, 5); OP(d, 4, 6); OP(d, 7, 8); OP(d, 2, 3); \
    OP(d, 4, 5); OP(d, 6, 7); OP(d, 3, 4); OP(d, 5, 6);              \
} while (0)

#define XC_SORT_NETWORK_11(OP, d) do {                                    \
    OP(d, 0, 9);  OP(d, 1, 6);  OP(d, 2, 4);  OP(d, 3, 7);  OP(d, 5, 8);  \
    OP(d, 0, 1);  OP(d, 3, 5);  OP(d, 4, 10); OP(d, 6, 9);  OP(d, 7, 8);  \
    OP(d, 1, 3);  OP(d, 2, 5);  OP(d, 4, 7);  OP(d, 8, 10); OP(d, 0, 4);  \
    OP(d, 1, 2);  OP(d, 3, 7);  OP(d, 5, 9);  OP(d, 6, 8);  OP(d, 0, 1);  \
    OP(d, 2, 6);  OP(d, 4, 5);  OP(d, 7, 8);  OP(d, 9, 10); OP(d, 2, 4);  \
    OP(d, 3, 6);  OP(d, 5, 7);  OP(d, 8, 9);  OP(d, 1, 2);  OP(d, 3, 4);  \
    OP(d, 5, 6);  OP(d, 7, 8);  OP(d, 2, 3);  OP(d, 4, 5);  OP(d, 6, 7);  \
} while (0)

#define XC_SORT_NETWORK_12(OP, d) do {                                    \
    OP(d, 0, 8);  OP(d, 1, 7);  OP(d, 2, 6);  OP(d, 3, 11); OP(d, 4, 10); \
    OP(d, 5, 9);  OP(d, 0, 1);  OP(d, 2, 5);  OP(d, 3, 4);  OP(d, 6, 9);  \
    OP(d, 7, 8);  OP(d, 10, 11);OP(d, 0, 2);  OP(d, 1, 6);  OP(d, 5, 10); \
    OP(d, 9, 11); OP(d, 0, 3);  OP(d, 1, 2);  OP(d, 4, 6);  OP(d, 5, 7);  \
    OP(d, 8, 11); OP(d, 9, 10); OP(d, 1, 4);  OP(d, 3, 5);  OP(d, 6, 8);  \
    OP(d, 7, 10); OP(d, 1, 3);  OP(d, 2, 5);  OP(d, 6, 9);  OP(d, 8, 10); \
    OP(d, 2, 3);  OP(d, 4, 5);  OP(d, 6, 7);  OP(d, 8, 9);  OP(d, 4, 6);  \
    OP(d, 5, 7);  OP(d, 3, 4);  OP(d, 5, 6);  OP(d, 7, 8);                \
} while (0)

#define XC_SORT_NETWORK_13(OP, d) do {                                    \
    OP(d, 0, 12); OP(d, 1, 10); OP(d, 2, 9);  OP(d, 3, 7);  OP(d, 5, 11); \
    OP(d, 6, 8);  OP(d, 1, 6);  OP(d, 2, 3);  OP(d, 4, 11); OP(d, 7, 9);  \
    OP(d, 8, 10); OP(d, 0, 4);  OP(d, 1, 2);  OP(d, 3, 6);  OP(d, 7, 8);  \
    OP(d, 9, 10); OP(d, 11, 12);OP(d, 4, 6);  OP(d, 5, 9);  OP(d, 8, 11); \
    OP(d, 10, 12);OP(d, 0, 5);  OP(d, 3, 8);  OP(d, 4, 7);  OP(d, 6, 11); \
    OP(d, 9, 10); OP(d, 0, 1);  OP(d, 2, 5);  OP(d, 6, 9);  OP(d, 7, 8);  \
    OP(d, 10, 11);OP(d, 1, 3);  OP(d, 2, 4);  OP(d, 5, 6);  OP(d, 9, 10); \
    OP(d, 1, 2);  OP(d, 3, 4);  OP(d, 5, 7);  OP(d, 6, 8);  OP(d, 2, 3);  \
    OP(d, 4, 5);  OP(d, 6, 7);  OP(d, 8, 9);  OP(d, 3, 4);  OP(d, 5, 6);  \
} while (0)

#define XC_SORT_NETWORK_14(OP, d) do {                                     \
    OP(d, 0, 1);  OP(d, 2, 3);  OP(d, 4, 5);  OP(d, 6, 7);  OP(d, 8, 9);   \
    OP(d, 10, 11);OP(d, 12, 13);OP(d, 0, 2);  OP(d, 1, 3);  OP(d, 4, 8);   \
    OP(d, 5, 9);  OP(d, 10, 12);OP(d, 11, 13);OP(d, 0, 4);  OP(d, 1, 2);   \
    OP(d, 3, 7);  OP(d, 5, 8);  OP(d, 6, 10); OP(d, 9, 13); OP(d, 11, 12); \
    OP(d, 0, 6);  OP(d, 1, 5);  OP(d, 3, 9);  OP(d, 4, 10); OP(d, 7, 13);  \
    OP(d, 8, 12); OP(d, 2, 10); OP(d, 3, 11); OP(d, 4, 6);  OP(d, 7, 9);   \
    OP(d, 1, 3);  OP(d, 2, 8);  OP(d, 5, 11); OP(d, 6, 7);  OP(d, 10, 12); \
    OP(d, 1, 4);  OP(d, 2, 6);  OP(d, 3, 5);  OP(d, 7, 11); OP(d, 8, 10);  \
    OP(d, 9, 12); OP(d, 2, 4);  OP(d, 3, 6);  OP(d, 5, 8);  OP(d, 7, 10);  \
    OP(d, 9, 11); OP(d, 3, 4);  OP(d, 5, 6);  OP(d, 7, 8);  OP(d, 9, 10);  \
    OP(d, 6, 7);                                                           \
} while (0)

#define XC_SORT_NETWORK_15(OP, d) do {                                     \
    OP(d, 1, 2);  OP(d, 3, 10); OP(d, 4, 14); OP(d, 5, 8);  OP(d, 6, 13);  \
    OP(d, 7, 12); OP(d, 9, 11); OP(d, 0, 14); OP(d, 1, 5);  OP(d, 2, 8);   \
    OP(d, 3, 7);  OP(d, 6, 9);  OP(d, 10, 12);OP(d, 11, 13);OP(d, 0, 7);   \
    OP(d, 1, 6);  OP(d, 2, 9);  OP(d, 4, 10); OP(d, 5, 11); OP(d, 8, 13);  \
    OP(d, 12, 14);OP(d, 0, 6);  OP(d, 2, 4);  OP(d, 3, 5);  OP(d, 7, 11);  \
    OP(d, 8, 10); OP(d, 9, 12); OP(d, 13, 14);OP(d, 0, 3);  OP(d, 1, 2);   \
    OP(d, 4, 7);  OP(d, 5, 9);  OP(d, 6, 8);  OP(d, 10, 11);OP(d, 12, 13); \
    OP(d, 0, 1);  OP(d, 2, 3);  OP(d, 4, 6);  OP(d, 7, 9);  OP(d, 10, 12); \
    OP(d, 11, 13);OP(d, 1, 2);  OP(d, 3, 5);  OP(d, 8, 10); OP(d, 11, 12); \
    OP(d, 3, 4);  OP(d, 5, 6);  OP(d, 7, 8);  OP(d, 9, 10); OP(d, 2, 3);   \
    OP(d, 4, 5);  OP(d, 6, 7);  OP(d, 8, 9);  OP(d, 10, 11);OP(d, 5, 6);   \
    OP(d, 7, 8);                                                           \
} while (0)

#define XC_SORT_NETWORK_16(OP, d) do {                                     \
    OP(d, 0, 13); OP(d, 1, 12); OP(d, 2, 15); OP(d, 3, 14); OP(d, 4, 8);   \
    OP(d, 5, 6);  OP(d, 7, 11); OP(d, 9, 10); OP(d, 0, 5);  OP(d, 1, 7);   \
    OP(d, 2, 9);  OP(d, 3, 4);  OP(d, 6, 13); OP(d, 8, 14); OP(d, 10, 15); \
    OP(d, 11, 12);OP(d, 0, 1);  OP(d, 2, 3);  OP(d, 4, 5);  OP(d, 6, 8);   \
    OP(d, 7, 9);  OP(d, 10, 11);OP(d, 12, 13);OP(d, 14, 15);OP(d, 0, 2);   \
    OP(d, 1, 3);  OP(d, 4, 10); OP(d, 5, 11); OP(d, 6, 7);  OP(d, 8, 9);   \
    OP(d, 12, 14);OP(d, 13, 15);OP(d, 1, 2);  OP(d, 3, 12); OP(d, 4, 6);   \
    OP(d, 5, 7);  OP(d, 8, 10); OP(d, 9, 11); OP(d, 13, 14);OP(d, 1, 4);   \
    OP(d, 2, 6);  OP(d, 5, 8);  OP(d, 7, 10); OP(d, 9, 13); OP(d, 11, 14); \
    OP(d, 2, 4);  OP(d, 3, 6);  OP(d, 9, 12); OP(d, 11, 13);OP(d, 3, 5);   \
    OP(d, 6, 8);  OP(d, 7, 9);  OP(d, 10, 12);OP(d, 3, 4);  OP(d, 5, 6);   \
    OP(d, 7, 8);  OP(d, 9, 10); OP(d, 11, 12);OP(d, 6, 7);  OP(d, 8, 9);   \
} while (0)

#define XC_SORT_NETWORK_CASE(N, OP, d) \
    case N:                            \
        XC_SORT_NETWORK_##N(OP, d);    \
        return;
#define XC_SORT_NETWORK_CASES(OP, d)   \
    XC_SORT_NETWORK_CASE(3,  OP, d)    \
    XC_SORT_NETWORK_CASE(4,  OP, d)    \
    XC_SORT_NETWORK_CASE(5,  OP, d)    \
    XC_SORT_NETWORK_CASE(6,  OP, d)    \
    XC_SORT_NETWORK_CASE(7,  OP, d)    \
    XC_SORT_NETWORK_CASE(8,  OP, d)    \
    XC_SORT_NETWORK_CASE(9,  OP, d)    \
    XC_SORT_NETWORK_CASE(10, OP, d)    \
    XC_SORT_NETWORK_CASE(11, OP, d)    \
    XC_SORT_NETWORK_CASE(12, OP, d)    \
    XC_SORT_NETWORK_CASE(13, OP, d)    \
    XC_SORT_NETWORK_CASE(14, OP, d)    \
    XC_SORT_NETWORK_CASE(15, OP, d)    \
    XC_SORT_NETWORK_CASE(16, OP, d)


/* ========================================================================= */
/* CUDA & C++                                                                */
/* ========================================================================= */

#if defined(XO_CONTEXT_CUDA) || (defined(XO_CONTEXT_CPU) && defined(__cplusplus))

    template <typename T>
    GPUFUN void xc_insertion_sort_asc_impl(T* arr, int64_t length) {
        for (int64_t i = 1; i < length; ++i) {
            T key = arr[i];
            int64_t j = i;
            while (j > 0 && key < arr[j - 1]) {
                arr[j] = arr[j - 1];
                --j;
            }
            arr[j] = key;
        }
    }

    template <typename T>
    GPUFUN void xc_insertion_sort_desc_impl(T* arr, int64_t length) {
        for (int64_t i = 1; i < length; ++i) {
            T key = arr[i];
            int64_t j = i;
            while (j > 0 && arr[j - 1] < key) {
                arr[j] = arr[j - 1];
                --j;
            }
            arr[j] = key;
        }
    }

    template <typename T>
    GPUFUN void xc_shell_sort_asc_impl(T* arr, int64_t length) {
        const int64_t n_gaps = 7;
        const int64_t gaps[] = {701, 301, 132, 57, 23, 10, 4};
        for (int64_t igap = 0; igap < n_gaps; ++igap) {
            const int64_t gap = gaps[igap];
            if (gap >= length) {
                continue;
            }
            for (int64_t i = gap; i < length; ++i) {
                const T value = arr[i];
                int64_t j = i;
                while (j >= gap && value < arr[j - gap]) {
                    arr[j] = arr[j - gap];
                    j -= gap;
                }
                arr[j] = value;
            }
        }
        xc_insertion_sort_asc_impl(arr, length);
    }

    template <typename T>
    GPUFUN void xc_shell_sort_desc_impl(T* arr, int64_t length){
        const int64_t n_gaps = 7;
        const int64_t gaps[] = {701, 301, 132, 57, 23, 10, 4};
        for (int64_t igap = 0; igap < n_gaps; ++igap) {
            const int64_t gap = gaps[igap];
            if (gap >= length) {
                continue;
            }
            for (int64_t i = gap; i < length; ++i) {
                const T value = arr[i];
                int64_t j = i;
                while (j >= gap && value > arr[j - gap]) {
                    arr[j] = arr[j - gap];
                    j -= gap;
                }
                arr[j] = value;
            }
        }
        xc_insertion_sort_desc_impl(arr, length);
    }

    template <typename T>
    GPUFUN void xc_sort_array_asc_impl(T* arr, int64_t length) {
        switch (length) {
            case 0:
                return;
            case 1:
                return;
            case 2:
                XC_ASORT(arr, 0, 1);
                return;
            XC_SORT_NETWORK_CASES(XC_ASORT, arr)
            default:
                if (length <= 32)
                    xc_insertion_sort_asc_impl(arr, length);
                else
                    xc_shell_sort_asc_impl(arr, length);
                return;
        }
    }

    template <typename T>
    GPUFUN void xc_sort_array_desc_impl(T* arr, int64_t length) {
        switch (length) {
            case 0:
                return;
            case 1:
                return;
            case 2:
                XC_DSORT(arr, 0, 1);
                return;
            XC_SORT_NETWORK_CASES(XC_DSORT, arr)
            default:
                if (length <= 32)
                    xc_insertion_sort_desc_impl(arr, length);
                else
                    xc_shell_sort_desc_impl(arr, length);
                return;
        }
    }

    #define XC_SORT_ASC(arr, length) xc_sort_array_asc_impl((arr), (length))
    #define XC_SORT_DESC(arr, length) xc_sort_array_desc_impl((arr), (length))


/* ========================================================================= */
/* CPU                                                                       */
/* ========================================================================= */

#elif defined(XO_CONTEXT_CPU)

    // Helper macro to define sorting functions for each type: descending and ascending order,
    // using insertion sort for larger arrays and sorting networks for small arrays.
    #define XC_DEFINE_SORT_FUNCTIONS(T, SUFFIX)                               \
        GPUFUN void xc_insertion_sort_asc_##SUFFIX(T* arr, int64_t length) {  \
            for (int64_t i = 1; i < length; ++i) {                            \
                T key = arr[i];                                               \
                int64_t j = i;                                                \
                while (j > 0 && key < arr[j - 1]) {                           \
                    arr[j] = arr[j - 1];                                      \
                    --j;                                                      \
                }                                                             \
                arr[j] = key;                                                 \
            }                                                                 \
        }                                                                     \
        GPUFUN void xc_insertion_sort_desc_##SUFFIX(T* arr, int64_t length) { \
            for (int64_t i = 1; i < length; ++i) {                            \
                T key = arr[i];                                               \
                int64_t j = i;                                                \
                while (j > 0 && arr[j - 1] < key) {                           \
                    arr[j] = arr[j - 1];                                      \
                    --j;                                                      \
                }                                                             \
                arr[j] = key;                                                 \
            }                                                                 \
        }                                                                     \
        GPUFUN void xc_shell_sort_asc_##SUFFIX(T* arr, int64_t length) {      \
            const int64_t n_gaps = 7;                                         \
            const int64_t gaps[] = {701, 301, 132, 57, 23, 10, 4};            \
            for (int64_t igap = 0; igap < n_gaps; ++igap) {                   \
                const int64_t gap = gaps[igap];                               \
                if (gap >= length) {                                          \
                    continue;                                                 \
                }                                                             \
                for (int64_t i = gap; i < length; ++i) {                      \
                    const T value = arr[i];                                   \
                    int64_t j = i;                                            \
                    while (j >= gap && value < arr[j - gap]) {                \
                        arr[j] = arr[j - gap];                                \
                        j -= gap;                                             \
                    }                                                         \
                    arr[j] = value;                                           \
                }                                                             \
            }                                                                 \
            xc_insertion_sort_asc_##SUFFIX(arr, length);                      \
        }                                                                     \
        GPUFUN void xc_shell_sort_desc_##SUFFIX(T* arr, int64_t length){      \
            const int64_t n_gaps = 7;                                         \
            const int64_t gaps[] = {701, 301, 132, 57, 23, 10, 4};            \
            for (int64_t igap = 0; igap < n_gaps; ++igap) {                   \
                const int64_t gap = gaps[igap];                               \
                if (gap >= length) {                                          \
                    continue;                                                 \
                }                                                             \
                for (int64_t i = gap; i < length; ++i) {                      \
                    const T value = arr[i];                                   \
                    int64_t j = i;                                            \
                    while (j >= gap && value > arr[j - gap]) {                \
                        arr[j] = arr[j - gap];                                \
                        j -= gap;                                             \
                    }                                                         \
                    arr[j] = value;                                           \
                }                                                             \
            }                                                                 \
            xc_insertion_sort_desc_##SUFFIX(arr, length);                     \
        }                                                                     \
        GPUFUN void xc_sort_array_asc_##SUFFIX(T* arr, int64_t length) {      \
            switch (length) {                                                 \
                case 0:                                                       \
                    return;                                                   \
                case 1:                                                       \
                    return;                                                   \
                case 2:                                                       \
                    XC_ASORT(arr, 0, 1);                                      \
                    return;                                                   \
                XC_SORT_NETWORK_CASES(XC_ASORT, arr)                          \
                default:                                                      \
                    if (length <= 32)                                         \
                        xc_insertion_sort_asc_##SUFFIX(arr, length);          \
                    else                                                      \
                        xc_shell_sort_asc_##SUFFIX(arr, length);              \
                    return;                                                   \
            }                                                                 \
        }                                                                     \
        GPUFUN void xc_sort_array_desc_##SUFFIX(T* arr, int64_t length) {     \
            switch (length) {                                                 \
                case 0:                                                       \
                    return;                                                   \
                case 1:                                                       \
                    return;                                                   \
                case 2:                                                       \
                    XC_DSORT(arr, 0, 1);                                      \
                    return;                                                   \
                XC_SORT_NETWORK_CASES(XC_DSORT, arr)                          \
                default:                                                      \
                    if (length <= 32)                                         \
                        xc_insertion_sort_desc_##SUFFIX(arr, length);         \
                    else                                                      \
                        xc_shell_sort_desc_##SUFFIX(arr, length);             \
                    return;                                                   \
            }                                                                 \
        }

    // Call sorting function definition macro for each type
    #define XC_CPU_SORT_TYPES(X) \
        X(int8_t,      int8_t)   \
        X(int16_t,     int16_t)  \
        X(int32_t,     int32_t)  \
        X(int64_t,     int64_t)  \
        X(uint8_t,     uint8_t)  \
        X(uint16_t,    uint16_t) \
        X(uint32_t,    uint32_t) \
        X(uint64_t,    uint64_t) \
        X(float,       float)    \
        X(double,      double)   \
        X(long double, ldouble)
    XC_CPU_SORT_TYPES(XC_DEFINE_SORT_FUNCTIONS)

    #define XC_SORT_ARRAY_ASC(arr) _Generic((arr),  \
        int8_t*:      xc_sort_array_asc_int8_t,     \
        int16_t*:     xc_sort_array_asc_int16_t,    \
        int32_t*:     xc_sort_array_asc_int32_t,    \
        int64_t*:     xc_sort_array_asc_int64_t,    \
        uint8_t*:     xc_sort_array_asc_uint8_t,    \
        uint16_t*:    xc_sort_array_asc_uint16_t,   \
        uint32_t*:    xc_sort_array_asc_uint32_t,   \
        uint64_t*:    xc_sort_array_asc_uint64_t,   \
        float*:       xc_sort_array_asc_float,      \
        double*:      xc_sort_array_asc_double,     \
        long double*: xc_sort_array_asc_ldouble     \
    )

    #define XC_SORT_ARRAY_DESC(arr) _Generic((arr), \
        int8_t*:      xc_sort_array_desc_int8_t,    \
        int16_t*:     xc_sort_array_desc_int16_t,   \
        int32_t*:     xc_sort_array_desc_int32_t,   \
        int64_t*:     xc_sort_array_desc_int64_t,   \
        uint8_t*:     xc_sort_array_desc_uint8_t,   \
        uint16_t*:    xc_sort_array_desc_uint16_t,  \
        uint32_t*:    xc_sort_array_desc_uint32_t,  \
        uint64_t*:    xc_sort_array_desc_uint64_t,  \
        float*:       xc_sort_array_desc_float,     \
        double*:      xc_sort_array_desc_double,    \
        long double*: xc_sort_array_desc_ldouble    \
    )

    #undef XC_DEFINE_SORT_FUNCTIONS
    #undef XC_CPU_SORT_TYPES

    #define XC_SORT_ASC(arr, length) XC_SORT_ARRAY_ASC(arr)((arr), (length))
    #define XC_SORT_DESC(arr, length) XC_SORT_ARRAY_DESC(arr)((arr), (length))


/* ========================================================================= */
/* OpenCL                                                                    */
/* ========================================================================= */

#elif defined(XO_CONTEXT_CL)

    // Helper macro to define sorting functions for each type: descending and ascending order,
    // using insertion sort for larger arrays and sorting networks for small arrays.
    #define XC_DEFINE_OCL_SORT_FUNCTIONS(T, ADDRQ)                                           \
        GPUFUN void OCL_OVERLOAD xc_insertion_sort_asc_impl(ADDRQ T* arr, int64_t length) {  \
            for (int64_t i = 1; i < length; ++i) {                                           \
                T key = arr[i];                                                              \
                int64_t j = i;                                                               \
                while (j > 0 && key < arr[j - 1]) {                                          \
                    arr[j] = arr[j - 1];                                                     \
                    --j;                                                                     \
                }                                                                            \
                arr[j] = key;                                                                \
            }                                                                                \
        }                                                                                    \
        GPUFUN void OCL_OVERLOAD xc_insertion_sort_desc_impl(ADDRQ T* arr, int64_t length) { \
            for (int64_t i = 1; i < length; ++i) {                                           \
                T key = arr[i];                                                              \
                int64_t j = i;                                                               \
                while (j > 0 && arr[j - 1] < key) {                                          \
                    arr[j] = arr[j - 1];                                                     \
                    --j;                                                                     \
                }                                                                            \
                arr[j] = key;                                                                \
            }                                                                                \
        }                                                                                    \
        GPUFUN void OCL_OVERLOAD xc_shell_sort_asc_impl(ADDRQ T* arr, int64_t length) {      \
            const int64_t n_gaps = 7;                                                        \
            const int64_t gaps[] = {701, 301, 132, 57, 23, 10, 4};                           \
            for (int64_t igap = 0; igap < n_gaps; ++igap) {                                  \
                const int64_t gap = gaps[igap];                                              \
                if (gap >= length) {                                                         \
                    continue;                                                                \
                }                                                                            \
                for (int64_t i = gap; i < length; ++i) {                                     \
                    const T value = arr[i];                                                  \
                    int64_t j = i;                                                           \
                    while (j >= gap && value < arr[j - gap]) {                               \
                        arr[j] = arr[j - gap];                                               \
                        j -= gap;                                                            \
                    }                                                                        \
                    arr[j] = value;                                                          \
                }                                                                            \
            }                                                                                \
            xc_insertion_sort_asc_impl(arr, length);                                         \
        }                                                                                    \
        GPUFUN void OCL_OVERLOAD xc_shell_sort_desc_impl(ADDRQ T* arr, int64_t length){      \
            const int64_t n_gaps = 7;                                                        \
            const int64_t gaps[] = {701, 301, 132, 57, 23, 10, 4};                           \
            for (int64_t igap = 0; igap < n_gaps; ++igap) {                                  \
                const int64_t gap = gaps[igap];                                              \
                if (gap >= length) {                                                         \
                    continue;                                                                \
                }                                                                            \
                for (int64_t i = gap; i < length; ++i) {                                     \
                    const T value = arr[i];                                                  \
                    int64_t j = i;                                                           \
                    while (j >= gap && value > arr[j - gap]) {                               \
                        arr[j] = arr[j - gap];                                               \
                        j -= gap;                                                            \
                    }                                                                        \
                    arr[j] = value;                                                          \
                }                                                                            \
            }                                                                                \
            xc_insertion_sort_desc_impl(arr, length);                                        \
        }                                                                                    \
        GPUFUN void OCL_OVERLOAD xc_sort_array_asc_impl(ADDRQ T* arr, int64_t length) {      \
            switch (length) {                                                                \
                case 0:                                                                      \
                    return;                                                                  \
                case 1:                                                                      \
                    return;                                                                  \
                case 2:                                                                      \
                    XC_ASORT(arr, 0, 1);                                                     \
                    return;                                                                  \
                XC_SORT_NETWORK_CASES(XC_ASORT, arr)                                         \
                default:                                                                     \
                    if (length <= 32)                                                        \
                        xc_insertion_sort_asc_impl(arr, length);                             \
                    else                                                                     \
                        xc_shell_sort_asc_impl(arr, length);                                 \
                    return;                                                                  \
            }                                                                                \
        }                                                                                    \
        GPUFUN void OCL_OVERLOAD xc_sort_array_desc_impl(ADDRQ T* arr, int64_t length) {     \
            switch (length) {                                                                \
                case 0:                                                                      \
                    return;                                                                  \
                case 1:                                                                      \
                    return;                                                                  \
                case 2:                                                                      \
                    XC_DSORT(arr, 0, 1);                                                     \
                    return;                                                                  \
                XC_SORT_NETWORK_CASES(XC_DSORT, arr)                                         \
                default:                                                                     \
                    if (length <= 32)                                                        \
                        xc_insertion_sort_desc_impl(arr, length);                            \
                    else                                                                     \
                        xc_shell_sort_desc_impl(arr, length);                                \
                    return;                                                                  \
            }                                                                                \
        }

    #define XC_DEFINE_OCL_SORT_TYPE(T)             \
        XC_DEFINE_OCL_SORT_FUNCTIONS(T, __global)  \
        XC_DEFINE_OCL_SORT_FUNCTIONS(T, __local)   \
        XC_DEFINE_OCL_SORT_FUNCTIONS(T, __private)

    XC_DEFINE_OCL_SORT_TYPE(int8_t)
    XC_DEFINE_OCL_SORT_TYPE(int16_t)
    XC_DEFINE_OCL_SORT_TYPE(int32_t)
    XC_DEFINE_OCL_SORT_TYPE(int64_t)
    XC_DEFINE_OCL_SORT_TYPE(uint8_t)
    XC_DEFINE_OCL_SORT_TYPE(uint16_t)
    XC_DEFINE_OCL_SORT_TYPE(uint32_t)
    XC_DEFINE_OCL_SORT_TYPE(uint64_t)
    XC_DEFINE_OCL_SORT_TYPE(float)
    XC_DEFINE_OCL_SORT_TYPE(double)

    #undef XC_DEFINE_OCL_SORT_FUNCTIONS
    #undef XC_DEFINE_OCL_SORT_TYPE

    #define XC_SORT_ASC(arr, length) xc_sort_array_asc_impl((arr), (length))
    #define XC_SORT_DESC(arr, length) xc_sort_array_desc_impl((arr), (length))

#else
#error "Xcoll header: No context defined!"
#endif // context selection


// Clean up macros to avoid polluting the global namespace
#undef XC_SORT_NETWORK_3
#undef XC_SORT_NETWORK_4
#undef XC_SORT_NETWORK_5
#undef XC_SORT_NETWORK_6
#undef XC_SORT_NETWORK_7
#undef XC_SORT_NETWORK_8
#undef XC_SORT_NETWORK_9
#undef XC_SORT_NETWORK_10
#undef XC_SORT_NETWORK_11
#undef XC_SORT_NETWORK_12
#undef XC_SORT_NETWORK_13
#undef XC_SORT_NETWORK_14
#undef XC_SORT_NETWORK_15
#undef XC_SORT_NETWORK_16
#undef XC_SORT_NETWORK_CASE
#undef XC_SORT_NETWORK_CASES

// #pragma GCC diagnostic pop
#endif /* XCOLL_GEOM_SORT_H */
