// copyright ############################### #
// This file is part of the Xcoll package.   #
// Copyright (c) CERN, 2024.                 #
// ######################################### #

#ifndef XCOLL_GEOM_SORT_H
#define XCOLL_GEOM_SORT_H

#ifdef XO_CONTEXT_CPU
#include <stdint.h>  // for int64_t etc
#endif  // XO_CONTEXT_CPU

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Warray-bounds"

// NVRTC rejects GNU statement-expression macros (the `({ ... })` form), so
// MAX/MIN/SWAP are written as plain expression / do-while macros that compile
// on both the CPU and the CUDA device path. Numerics are unchanged: MAX/MIN
// use the same `>`/`<` comparisons as before, and SWAP reads both array slots
// into temporaries (of the array element type) before writing, so it still
// places the smaller value in slot x and the larger in slot y. All call sites
// pass side-effect-free operands (array element reads / locals), so the
// unavoidable double-evaluation in the ternary macros is harmless.
//
// The element-type keyword differs by context: NVRTC compiles as C++ (use
// `decltype`); the CPU path is C (use the GNU `__typeof__`). NVRTC rejects
// `__typeof__`, hence the guard. `decltype((d)[x])` on the subscript lvalue
// yields a reference type, so the GPU branch strips it back to the value type
// with `decltype((d)[x] + 0)` (the `+ 0` makes a prvalue; for the only element
// types used here, double and int64_t, it leaves the type unchanged).
#ifdef XO_CONTEXT_CPU
#define XCOLL_GEOM_TYPEOF(expr) __typeof__(expr)
#else
#define XCOLL_GEOM_TYPEOF(expr) decltype((expr) + 0)
#endif

#ifdef MAX
#undef MAX
#pragma message ("Xcoll geometry: Compiler macro MAX redefined")
#endif
#define MAX(x, y) ((x) > (y) ? (x) : (y))
#ifdef MIN
#undef MIN
#pragma message ("Xcoll geometry: Compiler macro MIN redefined")
#endif
#define MIN(x, y) ((x) < (y) ? (x) : (y))
#ifdef SWAP
#error "Xcoll geometry: Compiler macro SWAP already defined!"
#endif
#define SWAP(d, x, y) do { \
        const XCOLL_GEOM_TYPEOF((d)[x]) _swap_lo = MIN((d)[x], (d)[y]); \
        const XCOLL_GEOM_TYPEOF((d)[x]) _swap_hi = MAX((d)[x], (d)[y]); \
        (d)[x] = _swap_lo; (d)[y] = _swap_hi; \
    } while (0)


// Fast methods
// ------------

static inline void sort_array_of_2_double(double* d){
    SWAP(d, 0, 1);
}

static inline void sort_array_of_3_double(double* d){
    SWAP(d, 0, 1); SWAP(d, 1, 2); SWAP(d, 0, 1);
}

static inline void sort_array_of_4_double(double* d){
    SWAP(d, 0, 1); SWAP(d, 2, 3); SWAP(d, 0, 2); SWAP(d, 1, 3); SWAP(d, 1, 2);
}

static inline void sort_array_of_5_double(double* d){
    SWAP(d, 0, 1); SWAP(d, 3, 4); SWAP(d, 2, 4); SWAP(d, 2, 3); SWAP(d, 1, 4);
    SWAP(d, 0, 3); SWAP(d, 0, 2); SWAP(d, 1, 3); SWAP(d, 1, 2);
}

static inline void sort_array_of_6_double(double* d){
    SWAP(d, 1, 2); SWAP(d, 4, 5); SWAP(d, 0, 2); SWAP(d, 3, 5); SWAP(d, 0, 1);
    SWAP(d, 3, 4); SWAP(d, 1, 4); SWAP(d, 0, 3); SWAP(d, 2, 5); SWAP(d, 1, 3);
    SWAP(d, 2, 4); SWAP(d, 2, 3);
}

static inline void sort_array_of_7_double(double* d){
    SWAP(d, 1, 2); SWAP(d, 3, 4); SWAP(d, 5, 6); SWAP(d, 0, 2); SWAP(d, 3, 5);
    SWAP(d, 4, 6); SWAP(d, 0, 1); SWAP(d, 4, 5); SWAP(d, 2, 6); SWAP(d, 0, 4);
    SWAP(d, 1, 5); SWAP(d, 0, 3); SWAP(d, 2, 5); SWAP(d, 1, 3); SWAP(d, 2, 4);
    SWAP(d, 2, 3);
}

static inline void sort_array_of_8_double(double* d){
    SWAP(d, 0, 1); SWAP(d, 2, 3); SWAP(d, 4, 5); SWAP(d, 6, 7); SWAP(d, 0, 2);
    SWAP(d, 1, 3); SWAP(d, 4, 6); SWAP(d, 5, 7); SWAP(d, 1, 2); SWAP(d, 5, 6);
    SWAP(d, 0, 4); SWAP(d, 3, 7); SWAP(d, 1, 5); SWAP(d, 2, 6); SWAP(d, 1, 4);
    SWAP(d, 3, 6); SWAP(d, 2, 4); SWAP(d, 3, 5); SWAP(d, 3, 4);
}

static inline void sort_array_of_2_int64(int64_t* d){
    SWAP(d, 0, 1);
}

static inline void sort_array_of_3_int64(int64_t* d){
    SWAP(d, 0, 1); SWAP(d, 1, 2); SWAP(d, 0, 1);
}

static inline void sort_array_of_4_int64(int64_t* d){
    SWAP(d, 0, 1); SWAP(d, 2, 3); SWAP(d, 0, 2); SWAP(d, 1, 3); SWAP(d, 1, 2);
}

static inline void sort_array_of_5_int64(int64_t* d){
    SWAP(d, 0, 1); SWAP(d, 3, 4); SWAP(d, 2, 4); SWAP(d, 2, 3); SWAP(d, 1, 4);
    SWAP(d, 0, 3); SWAP(d, 0, 2); SWAP(d, 1, 3); SWAP(d, 1, 2);
}

static inline void sort_array_of_6_int64(int64_t* d){
    SWAP(d, 1, 2); SWAP(d, 4, 5); SWAP(d, 0, 2); SWAP(d, 3, 5); SWAP(d, 0, 1);
    SWAP(d, 3, 4); SWAP(d, 1, 4); SWAP(d, 0, 3); SWAP(d, 2, 5); SWAP(d, 1, 3);
    SWAP(d, 2, 4); SWAP(d, 2, 3);
}

static inline void sort_array_of_7_int64(int64_t* d){
    SWAP(d, 1, 2); SWAP(d, 3, 4); SWAP(d, 5, 6); SWAP(d, 0, 2); SWAP(d, 3, 5);
    SWAP(d, 4, 6); SWAP(d, 0, 1); SWAP(d, 4, 5); SWAP(d, 2, 6); SWAP(d, 0, 4);
    SWAP(d, 1, 5); SWAP(d, 0, 3); SWAP(d, 2, 5); SWAP(d, 1, 3); SWAP(d, 2, 4);
    SWAP(d, 2, 3);
}

static inline void sort_array_of_8_int64(int64_t* d){
    SWAP(d, 0, 1); SWAP(d, 2, 3); SWAP(d, 4, 5); SWAP(d, 6, 7); SWAP(d, 0, 2);
    SWAP(d, 1, 3); SWAP(d, 4, 6); SWAP(d, 5, 7); SWAP(d, 1, 2); SWAP(d, 5, 6);
    SWAP(d, 0, 4); SWAP(d, 3, 7); SWAP(d, 1, 5); SWAP(d, 2, 6); SWAP(d, 1, 4);
    SWAP(d, 3, 6); SWAP(d, 2, 4); SWAP(d, 3, 5); SWAP(d, 3, 4);
}


// Generic methods
// ---------------

// The comparator + qsort default branch below are host-only (qsort and
// function-pointer comparators are unsupported on the CUDA device path).
// The crystal geometry path sorts <=8 crossings, so it always takes a
// hard-coded sort_array_of_N sorting-network branch above and never reaches
// the default branch; the GPU fallback here is a bounded in-place insertion
// sort, which yields an identical ascending order to qsort.
#ifdef XO_CONTEXT_CPU
int cmpfunc_double(const void * a, const void * b) {
   return ( *(double*)a - *(double*)b );
}
#endif  // XO_CONTEXT_CPU

/*gpufun*/
void insertion_sort_double(double* arr, int64_t length){
    for (int64_t i = 1; i < length; i++){
        double key = arr[i];
        int64_t j = i - 1;
        while (j >= 0 && arr[j] > key){
            arr[j + 1] = arr[j];
            j--;
        }
        arr[j + 1] = key;
    }
}

static inline void sort_array_of_double(double* arr, int64_t length){
    switch(length){
        case 2:
            sort_array_of_2_double(arr);
            break;
        case 3:
            sort_array_of_3_double(arr);
            break;
        case 4:
            sort_array_of_4_double(arr);
            break;
        case 5:
            sort_array_of_5_double(arr);
            break;
        case 6:
            sort_array_of_6_double(arr);
            break;
        case 7:
            sort_array_of_7_double(arr);
            break;
        case 8:
            sort_array_of_8_double(arr);
            break;
        default:
#ifdef XO_CONTEXT_CPU
            qsort(arr, length, sizeof(double), cmpfunc_double);
#else
            insertion_sort_double(arr, length);
#endif  // XO_CONTEXT_CPU
    }
}

#ifdef XO_CONTEXT_CPU
int cmpfunc_int64(const void * a, const void * b) {
   return ( *(int64_t*)a - *(int64_t*)b );
}
#endif  // XO_CONTEXT_CPU

/*gpufun*/
void insertion_sort_int64(int64_t* arr, int64_t length){
    for (int64_t i = 1; i < length; i++){
        int64_t key = arr[i];
        int64_t j = i - 1;
        while (j >= 0 && arr[j] > key){
            arr[j + 1] = arr[j];
            j--;
        }
        arr[j + 1] = key;
    }
}

static inline void sort_array_of_int64(int64_t* arr, int64_t length){
    switch(length){
        case 2:
            sort_array_of_2_int64(arr);
            break;
        case 3:
            sort_array_of_3_int64(arr);
            break;
        case 4:
            sort_array_of_4_int64(arr);
            break;
        case 5:
            sort_array_of_5_int64(arr);
            break;
        case 6:
            sort_array_of_6_int64(arr);
            break;
        case 7:
            sort_array_of_7_int64(arr);
            break;
        case 8:
            sort_array_of_8_int64(arr);
            break;
        default:
#ifdef XO_CONTEXT_CPU
            qsort(arr, length, sizeof(int64_t), cmpfunc_int64);
#else
            insertion_sort_int64(arr, length);
#endif  // XO_CONTEXT_CPU
    }
}

#pragma GCC diagnostic pop
#endif /* XCOLL_GEOM_SORT_H */
