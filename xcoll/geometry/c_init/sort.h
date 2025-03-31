// copyright ############################### #
// This file is part of the Xcoll package.   #
// Copyright (c) CERN, 2025.                 #
// ######################################### #

#ifndef XCOLL_GEOM_SORT_H
#define XCOLL_GEOM_SORT_H

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Warray-bounds"

#ifdef MAX
#undef MAX
#pragma message ("Xcoll geometry: Compiler macro MAX redefined")
#endif
#define MAX(x, y) ({const __typeof__ (x) _x = (x); \
                    const __typeof__ (y) _y = (y); \
                    _x > _y ? _x : _y; })
#ifdef MIN
#undef MIN
#pragma message ("Xcoll geometry: Compiler macro MIN redefined")
#endif
#define MIN(x, y) ({const __typeof__ (x) _x = (x); \
                    const __typeof__ (y) _y = (y); \
                    _x < _y ? _x : _y; })
#ifdef SWAP
#error "Xcoll geometry: Compiler macro SWAP already defined!"
#endif
#define SWAP(d,x,y) ({const __typeof__(*d) _x = MIN(d[x], d[y]); \
                      const __typeof__(*d) _y = MAX(d[x], d[y]); \
                      d[x] = _x; d[y] = _y; })
#ifdef INTERVALS_OVERLAP
#undef INTERVALS_OVERLAP
#pragma message ("Xcoll geometry: Compiler macro INTERVALS_OVERLAP redefined")
#endif
#define INTERVALS_OVERLAP(minA, maxA, minB, maxB) (((maxA) >= (minB)) && ((maxB) >= (minA)))


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

int cmpfunc_double(const void * a, const void * b) {
   return ( *(double*)a - *(double*)b );
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
            qsort(arr, length, sizeof(double), cmpfunc_double);
    }
}

int cmpfunc_int64(const void * a, const void * b) {
   return ( *(int64_t*)a - *(int64_t*)b );
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
            qsort(arr, length, sizeof(int64_t), cmpfunc_int64);
    }
}

#pragma GCC diagnostic pop
#endif /* XCOLL_GEOM_SORT_H */