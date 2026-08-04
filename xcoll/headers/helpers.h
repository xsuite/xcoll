// copyright ############################### #
// This file is part of the Xcoll package.   #
// Copyright (c) CERN, 2026.                 #
// ######################################### #

#ifndef XCOLL_HELPERS_H
#define XCOLL_HELPERS_H

#ifdef XO_CONTEXT_CPU
#include <stdint.h>  // for int64_t etc
#endif /* XO_CONTEXT_CPU */

#include "xobjects/headers/common.h"
#include "xobjects/headers/atomicadd.h"

/*
 * Type-safe helper macros for computing minima/maxima and sorting pairs.
 *
 * Public interface:
 *   XC_MAX(x, y)      Returns the larger of x and y.
 *   XC_MIN(x, y)      Returns the smaller of x and y.
 *   XC_ASORT(d, i, j) Ensures d[i] <= d[j].
 *   XC_DSORT(d, i, j) Ensures d[i] >= d[j].
 *
 * Design goals:
 *   - Evaluate every argument exactly once.
 *   - Preserve the input type without implicit promotion.
 *   - Behave consistently on CPU, CUDA, and OpenCL.
 *   - Avoid compiler-specific extensions where practical.
 *
 * Assumptions:
 *   - x and y have the same scalar arithmetic type.
 *   - d is an array (or pointer) of a supported scalar type.
 *   - Floating-point behaviour follows the comparison operators
 *         x > y ? x : y
 *     and
 *         x < y ? x : y
 *     (i.e. this intentionally does not implement the NaN semantics
 *     of fmax()/fmin()).
 *   - Mixed-type arguments are not supported; cast explicitly if needed.
 *
 * Implementation notes:
 *   A naive macro such as
 *
 *       #define MAX(x,y) ((x) > (y) ? (x) : (y))
 *
 *   may evaluate its arguments multiple times, leading to incorrect
 *   behaviour for expressions with side effects (e.g. MAX(i++, j)).
 *   These helpers instead dispatch to backend-specific implementations
 *   (C11 _Generic, C++ templates, or OpenCL overloads), ensuring that
 *   each argument is evaluated exactly once while retaining efficient,
 *   fully inlined code.
 */


/* ========================================================================= */
/* CUDA & C++                                                                */
/* ========================================================================= */

#if defined(XO_CONTEXT_CUDA) || (defined(XO_CONTEXT_CPU) && defined(__cplusplus))

    // CUDA source is C++, so use a function template.
    template <typename T>
    GPUFUN T xc_max_impl(T x, T y) {
        return x > y ? x : y;
    }
    #define XC_MAX(x, y) xc_max_impl((x), (y))

    template <typename T>
    GPUFUN T xc_min_impl(T x, T y) {
        return x < y ? x : y;
    }
    #define XC_MIN(x, y) xc_min_impl((x), (y))

    template <typename T>
    GPUFUN void xc_sort_pair_asc(T& a, T& b) {
        if (b < a) {
            T tmp = a;
            a = b;
            b = tmp;
        }
    }
    #define XC_ASORT(d, x, y) xc_sort_pair_asc((d)[(x)], (d)[(y)])

    template <typename T>
    GPUFUN void xc_sort_pair_desc(T& a, T& b) {
        if (a < b) {
            T tmp = a;
            a = b;
            b = tmp;
        }
    }
    #define XC_DSORT(d, x, y) xc_sort_pair_desc((d)[(x)], (d)[(y)])


/* ========================================================================= */
/* CPU                                                                       */
/* ========================================================================= */

#elif defined(XO_CONTEXT_CPU)

    // Typed inline implementations. We prefer these + _Generic dispatch over using
    // a general macro with typeof() because typeof() is a GCC extension and not standard C11.
    GPUFUN int8_t      xc_max_int8_t   (int8_t x, int8_t y)           {return x > y ? x : y;}
    GPUFUN int16_t     xc_max_int16_t  (int16_t x, int16_t y)         {return x > y ? x : y;}
    GPUFUN int32_t     xc_max_int32_t  (int32_t x, int32_t y)         {return x > y ? x : y;}
    GPUFUN int64_t     xc_max_int64_t  (int64_t x, int64_t y)         {return x > y ? x : y;}
    GPUFUN uint8_t     xc_max_uint8_t  (uint8_t x, uint8_t y)         {return x > y ? x : y;}
    GPUFUN uint16_t    xc_max_uint16_t (uint16_t x, uint16_t y)       {return x > y ? x : y;}
    GPUFUN uint32_t    xc_max_uint32_t (uint32_t x, uint32_t y)       {return x > y ? x : y;}
    GPUFUN uint64_t    xc_max_uint64_t (uint64_t x, uint64_t y)       {return x > y ? x : y;}
    GPUFUN float       xc_max_float    (float x, float y)             {return x > y ? x : y;}
    GPUFUN double      xc_max_double   (double x, double y)           {return x > y ? x : y;}

    GPUFUN int8_t      xc_min_int8_t   (int8_t x, int8_t y)           {return x < y ? x : y;}
    GPUFUN int16_t     xc_min_int16_t  (int16_t x, int16_t y)         {return x < y ? x : y;}
    GPUFUN int32_t     xc_min_int32_t  (int32_t x, int32_t y)         {return x < y ? x : y;}
    GPUFUN int64_t     xc_min_int64_t  (int64_t x, int64_t y)         {return x < y ? x : y;}
    GPUFUN uint8_t     xc_min_uint8_t  (uint8_t x, uint8_t y)         {return x < y ? x : y;}
    GPUFUN uint16_t    xc_min_uint16_t (uint16_t x, uint16_t y)       {return x < y ? x : y;}
    GPUFUN uint32_t    xc_min_uint32_t (uint32_t x, uint32_t y)       {return x < y ? x : y;}
    GPUFUN uint64_t    xc_min_uint64_t (uint64_t x, uint64_t y)       {return x < y ? x : y;}
    GPUFUN float       xc_min_float    (float x, float y)             {return x < y ? x : y;}
    GPUFUN double      xc_min_double   (double x, double y)           {return x < y ? x : y;}

    GPUFUN void xc_sort_pair_asc_int8_t   (int8_t* a, int8_t* b)           {if (*b < *a) {int8_t tmp = *a; *a = *b; *b = tmp;}}
    GPUFUN void xc_sort_pair_asc_int16_t  (int16_t* a, int16_t* b)         {if (*b < *a) {int16_t tmp = *a; *a = *b; *b = tmp;}}
    GPUFUN void xc_sort_pair_asc_int32_t  (int32_t* a, int32_t* b)         {if (*b < *a) {int32_t tmp = *a; *a = *b; *b = tmp;}}
    GPUFUN void xc_sort_pair_asc_int64_t  (int64_t* a, int64_t* b)         {if (*b < *a) {int64_t tmp = *a; *a = *b; *b = tmp;}}
    GPUFUN void xc_sort_pair_asc_uint8_t  (uint8_t* a, uint8_t* b)         {if (*b < *a) {uint8_t tmp = *a; *a = *b; *b = tmp;}}
    GPUFUN void xc_sort_pair_asc_uint16_t (uint16_t* a, uint16_t* b)       {if (*b < *a) {uint16_t tmp = *a; *a = *b; *b = tmp;}}
    GPUFUN void xc_sort_pair_asc_uint32_t (uint32_t* a, uint32_t* b)       {if (*b < *a) {uint32_t tmp = *a; *a = *b; *b = tmp;}}
    GPUFUN void xc_sort_pair_asc_uint64_t (uint64_t* a, uint64_t* b)       {if (*b < *a) {uint64_t tmp = *a; *a = *b; *b = tmp;}}
    GPUFUN void xc_sort_pair_asc_float    (float* a, float* b)             {if (*b < *a) {float tmp = *a; *a = *b; *b = tmp;}}
    GPUFUN void xc_sort_pair_asc_double   (double* a, double* b)           {if (*b < *a) {double tmp = *a; *a = *b; *b = tmp;}}

    GPUFUN void xc_sort_pair_desc_int8_t   (int8_t* a, int8_t* b)           {if (*a < *b) {int8_t tmp = *a; *a = *b; *b = tmp;}}
    GPUFUN void xc_sort_pair_desc_int16_t  (int16_t* a, int16_t* b)         {if (*a < *b) {int16_t tmp = *a; *a = *b; *b = tmp;}}
    GPUFUN void xc_sort_pair_desc_int32_t  (int32_t* a, int32_t* b)         {if (*a < *b) {int32_t tmp = *a; *a = *b; *b = tmp;}}
    GPUFUN void xc_sort_pair_desc_int64_t  (int64_t* a, int64_t* b)         {if (*a < *b) {int64_t tmp = *a; *a = *b; *b = tmp;}}
    GPUFUN void xc_sort_pair_desc_uint8_t  (uint8_t* a, uint8_t* b)         {if (*a < *b) {uint8_t tmp = *a; *a = *b; *b = tmp;}}
    GPUFUN void xc_sort_pair_desc_uint16_t (uint16_t* a, uint16_t* b)       {if (*a < *b) {uint16_t tmp = *a; *a = *b; *b = tmp;}}
    GPUFUN void xc_sort_pair_desc_uint32_t (uint32_t* a, uint32_t* b)       {if (*a < *b) {uint32_t tmp = *a; *a = *b; *b = tmp;}}
    GPUFUN void xc_sort_pair_desc_uint64_t (uint64_t* a, uint64_t* b)       {if (*a < *b) {uint64_t tmp = *a; *a = *b; *b = tmp;}}
    GPUFUN void xc_sort_pair_desc_float    (float* a, float* b)             {if (*a < *b) {float tmp = *a; *a = *b; *b = tmp;}}
    GPUFUN void xc_sort_pair_desc_double   (double* a, double* b)           {if (*a < *b) {double tmp = *a; *a = *b; *b = tmp;}}

    // Dispatching on +(x) applies lvalue conversion and removes top-level
    // qualifiers without evaluating x.
    #define XC_MAX_SELECT(x) _Generic(+(x),   \
                int8_t:      xc_max_int8_t,   \
                int16_t:     xc_max_int16_t,  \
                int32_t:     xc_max_int32_t,  \
                int64_t:     xc_max_int64_t,  \
                uint8_t:     xc_max_uint8_t,  \
                uint16_t:    xc_max_uint16_t, \
                uint32_t:    xc_max_uint32_t, \
                uint64_t:    xc_max_uint64_t, \
                float:       xc_max_float,    \
                double:      xc_max_double    \
            )
    #define XC_MAX(x, y) XC_MAX_SELECT(x)((x), (y))

    #define XC_MIN_SELECT(x) _Generic(+(x),   \
                int8_t:      xc_min_int8_t,   \
                int16_t:     xc_min_int16_t,  \
                int32_t:     xc_min_int32_t,  \
                int64_t:     xc_min_int64_t,  \
                uint8_t:     xc_min_uint8_t,  \
                uint16_t:    xc_min_uint16_t, \
                uint32_t:    xc_min_uint32_t, \
                uint64_t:    xc_min_uint64_t, \
                float:       xc_min_float,    \
                double:      xc_min_double    \
            )
    #define XC_MIN(x, y) XC_MIN_SELECT(x)((x), (y))

    #define XC_SORT_PAIR_ASC(ptr) _Generic((ptr),        \
                int8_t*:      xc_sort_pair_asc_int8_t,   \
                int16_t*:     xc_sort_pair_asc_int16_t,  \
                int32_t*:     xc_sort_pair_asc_int32_t,  \
                int64_t*:     xc_sort_pair_asc_int64_t,  \
                uint8_t*:     xc_sort_pair_asc_uint8_t,  \
                uint16_t*:    xc_sort_pair_asc_uint16_t, \
                uint32_t*:    xc_sort_pair_asc_uint32_t, \
                uint64_t*:    xc_sort_pair_asc_uint64_t, \
                float*:       xc_sort_pair_asc_float,    \
                double*:      xc_sort_pair_asc_double    \
            )
    #define XC_ASORT(d, x, y) XC_SORT_PAIR_ASC(&(d)[0])(&(d)[(x)], &(d)[(y)])

    #define XC_SORT_PAIR_DESC(ptr) _Generic((ptr),        \
                int8_t*:      xc_sort_pair_desc_int8_t,   \
                int16_t*:     xc_sort_pair_desc_int16_t,  \
                int32_t*:     xc_sort_pair_desc_int32_t,  \
                int64_t*:     xc_sort_pair_desc_int64_t,  \
                uint8_t*:     xc_sort_pair_desc_uint8_t,  \
                uint16_t*:    xc_sort_pair_desc_uint16_t, \
                uint32_t*:    xc_sort_pair_desc_uint32_t, \
                uint64_t*:    xc_sort_pair_desc_uint64_t, \
                float*:       xc_sort_pair_desc_float,    \
                double*:      xc_sort_pair_desc_double    \
            )
    #define XC_DSORT(d, x, y) XC_SORT_PAIR_DESC(&(d)[0])(&(d)[(x)], &(d)[(y)])


/* ========================================================================= */
/* OpenCL                                                                    */
/* ========================================================================= */

#elif defined(XO_CONTEXT_CL)

    /*
    * In in xobjects/headers/atomicadd.h:
    *   - Overloading by __attribute__((overloadable)) (checked and assigned to OCL_OVERLOAD)
    *   - Double precision (cl_khr_fp64) checked and activated
    */
    GPUFUN int8_t   OCL_OVERLOAD xc_max_impl (int8_t x, int8_t y)     {return x > y ? x : y;}
    GPUFUN int16_t  OCL_OVERLOAD xc_max_impl (int16_t x, int16_t y)   {return x > y ? x : y;}
    GPUFUN int32_t  OCL_OVERLOAD xc_max_impl (int32_t x, int32_t y)   {return x > y ? x : y;}
    GPUFUN int64_t  OCL_OVERLOAD xc_max_impl (int64_t x, int64_t y)   {return x > y ? x : y;}
    GPUFUN uint8_t  OCL_OVERLOAD xc_max_impl (uint8_t x, uint8_t y)   {return x > y ? x : y;}
    GPUFUN uint16_t OCL_OVERLOAD xc_max_impl (uint16_t x, uint16_t y) {return x > y ? x : y;}
    GPUFUN uint32_t OCL_OVERLOAD xc_max_impl (uint32_t x, uint32_t y) {return x > y ? x : y;}
    GPUFUN uint64_t OCL_OVERLOAD xc_max_impl (uint64_t x, uint64_t y) {return x > y ? x : y;}
    GPUFUN float    OCL_OVERLOAD xc_max_impl (float x, float y)       {return x > y ? x : y;}
    GPUFUN double   OCL_OVERLOAD xc_max_impl (double x, double y)     {return x > y ? x : y;}
    #define XC_MAX(x, y) xc_max_impl((x), (y))

    GPUFUN int8_t   OCL_OVERLOAD xc_min_impl (int8_t x, int8_t y)     {return x < y ? x : y;}
    GPUFUN int16_t  OCL_OVERLOAD xc_min_impl (int16_t x, int16_t y)   {return x < y ? x : y;}
    GPUFUN int32_t  OCL_OVERLOAD xc_min_impl (int32_t x, int32_t y)   {return x < y ? x : y;}
    GPUFUN int64_t  OCL_OVERLOAD xc_min_impl (int64_t x, int64_t y)   {return x < y ? x : y;}
    GPUFUN uint8_t  OCL_OVERLOAD xc_min_impl (uint8_t x, uint8_t y)   {return x < y ? x : y;}
    GPUFUN uint16_t OCL_OVERLOAD xc_min_impl (uint16_t x, uint16_t y) {return x < y ? x : y;}
    GPUFUN uint32_t OCL_OVERLOAD xc_min_impl (uint32_t x, uint32_t y) {return x < y ? x : y;}
    GPUFUN uint64_t OCL_OVERLOAD xc_min_impl (uint64_t x, uint64_t y) {return x < y ? x : y;}
    GPUFUN float    OCL_OVERLOAD xc_min_impl (float x, float y)       {return x < y ? x : y;}
    GPUFUN double   OCL_OVERLOAD xc_min_impl (double x, double y)     {return x < y ? x : y;}
    #define XC_MIN(x, y) xc_min_impl((x), (y))

    // For sort macros, need to overload for each memory space (global, local, private)
    // because OpenCL does not allow overloading on pointer types with different address spaces.
    #define XC_DEFINE_OCL_ASC_SORT(T, suffix) do {                                  \
        GPUFUN void OCL_OVERLOAD xc_sort_pair_asc(__global T* a, __global T* b) {   \
            if (*b < *a) {T tmp = *a; *a = *b; *b = tmp;}}                          \
        GPUFUN void OCL_OVERLOAD xc_sort_pair_asc(__local T* a, __local T* b) {     \
            if (*b < *a) {T tmp = *a; *a = *b; *b = tmp;}}                          \
        GPUFUN void OCL_OVERLOAD xc_sort_pair_asc(__private T* a, __private T* b) { \
            if (*b < *a) {T tmp = *a; *a = *b; *b = tmp;}}                          \
    } while(0)  // Do-while for scoping; allows trailing semicolon after macro usage
    XC_DEFINE_OCL_ASC_SORT(int8_t, int8_t)
    XC_DEFINE_OCL_ASC_SORT(int16_t, int16_t)
    XC_DEFINE_OCL_ASC_SORT(int32_t, int32_t)
    XC_DEFINE_OCL_ASC_SORT(int64_t, int64_t)
    XC_DEFINE_OCL_ASC_SORT(uint8_t, uint8_t)
    XC_DEFINE_OCL_ASC_SORT(uint16_t, uint16_t)
    XC_DEFINE_OCL_ASC_SORT(uint32_t, uint32_t)
    XC_DEFINE_OCL_ASC_SORT(uint64_t, uint64_t)
    XC_DEFINE_OCL_ASC_SORT(float, float)
    XC_DEFINE_OCL_ASC_SORT(double, double)
    #define XC_ASORT(d, x, y) xc_sort_pair_asc(&(d)[(x)], &(d)[(y)])

    #define XC_DEFINE_OCL_DESC_SORT(T, suffix) do {                                  \
        GPUFUN void OCL_OVERLOAD xc_sort_pair_desc(__global T* a, __global T* b) {   \
            if (*a < *b) {T tmp = *a; *a = *b; *b = tmp;}}                           \
        GPUFUN void OCL_OVERLOAD xc_sort_pair_desc(__local T* a, __local T* b) {     \
            if (*a < *b) {T tmp = *a; *a = *b; *b = tmp;}}                           \
        GPUFUN void OCL_OVERLOAD xc_sort_pair_desc(__private T* a, __private T* b) { \
            if (*a < *b) {T tmp = *a; *a = *b; *b = tmp;}}                           \
    } while(0)  // Do-while for scoping; allows trailing semicolon after macro usage
    XC_DEFINE_OCL_DESC_SORT(int8_t, int8_t)
    XC_DEFINE_OCL_DESC_SORT(int16_t, int16_t)
    XC_DEFINE_OCL_DESC_SORT(int32_t, int32_t)
    XC_DEFINE_OCL_DESC_SORT(int64_t, int64_t)
    XC_DEFINE_OCL_DESC_SORT(uint8_t, uint8_t)
    XC_DEFINE_OCL_DESC_SORT(uint16_t, uint16_t)
    XC_DEFINE_OCL_DESC_SORT(uint32_t, uint32_t)
    XC_DEFINE_OCL_DESC_SORT(uint64_t, uint64_t)
    XC_DEFINE_OCL_DESC_SORT(float, float)
    XC_DEFINE_OCL_DESC_SORT(double, double)
    #define XC_DSORT(d, x, y) xc_sort_pair_desc(&(d)[(x)], &(d)[(y)])

    #undef XC_DEFINE_OCL_ASC_SORT
    #undef XC_DEFINE_OCL_DESC_SORT

#else
#error "Xcoll header: No context defined!"
#endif /* context selection */


#endif /* XCOLL_HELPERS_H */
