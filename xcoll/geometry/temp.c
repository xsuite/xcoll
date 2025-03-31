#include <stdint.h>
#ifndef XOBJ_TYPEDEF_BoundingBox
#define XOBJ_TYPEDEF_BoundingBox
typedef   struct BoundingBox_s * BoundingBox;
 static inline BoundingBox BoundingBox_getp(BoundingBox restrict  obj){
  int64_t offset=0;
  return (BoundingBox)(( char*) obj+offset);
}
 static inline double BoundingBox_get_rC(const BoundingBox restrict  obj){
  int64_t offset=0;
  return *( double*)(( char*) obj+offset);
}
 static inline void BoundingBox_set_rC(BoundingBox restrict  obj, double value){
  int64_t offset=0;
  *( double*)(( char*) obj+offset)=value;
}
 static inline  double* BoundingBox_getp_rC(BoundingBox restrict  obj){
  int64_t offset=0;
  return ( double*)(( char*) obj+offset);
}
 static inline double BoundingBox_get_sin_tC(const BoundingBox restrict  obj){
  int64_t offset=0;
  offset+=8;
  return *( double*)(( char*) obj+offset);
}
 static inline void BoundingBox_set_sin_tC(BoundingBox restrict  obj, double value){
  int64_t offset=0;
  offset+=8;
  *( double*)(( char*) obj+offset)=value;
}
 static inline  double* BoundingBox_getp_sin_tC(BoundingBox restrict  obj){
  int64_t offset=0;
  offset+=8;
  return ( double*)(( char*) obj+offset);
}
 static inline double BoundingBox_get_cos_tC(const BoundingBox restrict  obj){
  int64_t offset=0;
  offset+=16;
  return *( double*)(( char*) obj+offset);
}
 static inline void BoundingBox_set_cos_tC(BoundingBox restrict  obj, double value){
  int64_t offset=0;
  offset+=16;
  *( double*)(( char*) obj+offset)=value;
}
 static inline  double* BoundingBox_getp_cos_tC(BoundingBox restrict  obj){
  int64_t offset=0;
  offset+=16;
  return ( double*)(( char*) obj+offset);
}
 static inline double BoundingBox_get_proj_l(const BoundingBox restrict  obj){
  int64_t offset=0;
  offset+=24;
  return *( double*)(( char*) obj+offset);
}
 static inline void BoundingBox_set_proj_l(BoundingBox restrict  obj, double value){
  int64_t offset=0;
  offset+=24;
  *( double*)(( char*) obj+offset)=value;
}
 static inline  double* BoundingBox_getp_proj_l(BoundingBox restrict  obj){
  int64_t offset=0;
  offset+=24;
  return ( double*)(( char*) obj+offset);
}
 static inline double BoundingBox_get_proj_w(const BoundingBox restrict  obj){
  int64_t offset=0;
  offset+=32;
  return *( double*)(( char*) obj+offset);
}
 static inline void BoundingBox_set_proj_w(BoundingBox restrict  obj, double value){
  int64_t offset=0;
  offset+=32;
  *( double*)(( char*) obj+offset)=value;
}
 static inline  double* BoundingBox_getp_proj_w(BoundingBox restrict  obj){
  int64_t offset=0;
  offset+=32;
  return ( double*)(( char*) obj+offset);
}
 static inline double BoundingBox_get_l(const BoundingBox restrict  obj){
  int64_t offset=0;
  offset+=40;
  return *( double*)(( char*) obj+offset);
}
 static inline void BoundingBox_set_l(BoundingBox restrict  obj, double value){
  int64_t offset=0;
  offset+=40;
  *( double*)(( char*) obj+offset)=value;
}
 static inline  double* BoundingBox_getp_l(BoundingBox restrict  obj){
  int64_t offset=0;
  offset+=40;
  return ( double*)(( char*) obj+offset);
}
 static inline double BoundingBox_get_w(const BoundingBox restrict  obj){
  int64_t offset=0;
  offset+=48;
  return *( double*)(( char*) obj+offset);
}
 static inline void BoundingBox_set_w(BoundingBox restrict  obj, double value){
  int64_t offset=0;
  offset+=48;
  *( double*)(( char*) obj+offset)=value;
}
 static inline  double* BoundingBox_getp_w(BoundingBox restrict  obj){
  int64_t offset=0;
  offset+=48;
  return ( double*)(( char*) obj+offset);
}
 static inline double BoundingBox_get_sin_tb(const BoundingBox restrict  obj){
  int64_t offset=0;
  offset+=56;
  return *( double*)(( char*) obj+offset);
}
 static inline void BoundingBox_set_sin_tb(BoundingBox restrict  obj, double value){
  int64_t offset=0;
  offset+=56;
  *( double*)(( char*) obj+offset)=value;
}
 static inline  double* BoundingBox_getp_sin_tb(BoundingBox restrict  obj){
  int64_t offset=0;
  offset+=56;
  return ( double*)(( char*) obj+offset);
}
 static inline double BoundingBox_get_cos_tb(const BoundingBox restrict  obj){
  int64_t offset=0;
  offset+=64;
  return *( double*)(( char*) obj+offset);
}
 static inline void BoundingBox_set_cos_tb(BoundingBox restrict  obj, double value){
  int64_t offset=0;
  offset+=64;
  *( double*)(( char*) obj+offset)=value;
}
 static inline  double* BoundingBox_getp_cos_tb(BoundingBox restrict  obj){
  int64_t offset=0;
  offset+=64;
  return ( double*)(( char*) obj+offset);
}
#endif

#ifndef XCOLL_GEOM_DEFINES_H
#define XCOLL_GEOM_DEFINES_H
#include <math.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>

#ifndef XC_GEOM_EPSILON
#define XC_GEOM_EPSILON 1e-15
#endif

#ifndef XC_GEOM_S_MAX
#define XC_GEOM_S_MAX 1e+21
#endif

#ifndef XC_GEOM_ROOT_NEWTON_EPSILON
#define XC_GEOM_ROOT_NEWTON_EPSILON 1e-10
#endif

#ifndef XC_GEOM_ROOT_NEWTON_MAX_ITER
#define XC_GEOM_ROOT_NEWTON_MAX_ITER 100
#endif

#ifndef XC_GEOM_ROOT_NEWTON_DERIVATIVE_TOL
#define XC_GEOM_ROOT_NEWTON_DERIVATIVE_TOL 1e-10
#endif

#ifndef XC_GEOM_ROOT_GRID_MAX_INTER
#define XC_GEOM_ROOT_GRID_MAX_INTER 10
#endif

#ifndef XC_GEOM_ROOT_GRID_POINTS
#define XC_GEOM_ROOT_GRID_POINTS 1000
#endif


#endif /* XCOLL_GEOM_DEFINES_H */

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
// copyright ############################### #
// This file is part of the Xcoll package.   #
// Copyright (c) CERN, 2025.                 #
// ######################################### #

#ifndef XCOLL_GEOM_METHODS_H
#define XCOLL_GEOM_METHODS_H


// This is a quick macro to use inside a function body on a parameter that is not
// used inside the function (this avoids throwing warnings at compilation time).
#ifndef UNUSED
#define UNUSED(expr) (void)(expr)
#endif


// This function calculates the overlap between an array and a given interval.
// The array comes in pairs of points, e.g. in-out-in-out... or out-in-out-in...
// IMPORTANT:
// The array and interval are assumed to be sorted!
// Furthermore, the array should have one extra slot allocated at the end, in case it needs to be expanded..
// This is always true for the arrays created by get_s, as we create them with 2*n_segments slots.
 static inline
void calculate_overlap_array_interval(double* arr, int8_t* length, double* interval){
    if (arr[0] > interval[1]){
        // No overlap
        *length = 0;
    }
    if ((*length)%2 == 1){
        // Special case: last interval of array is open until infinity,
        // so we add an extra point at the end to represent infinity.
        arr[*length] = XC_GEOM_S_MAX*0.1;
        (*length)++;
    } else if (arr[*length-1] < interval[0]){
        // No overlap
        *length = 0;
    }
    int8_t i_start = 0;
    // Find the start of overlap
    for (int8_t i=0; i<*length; i++){
        if (arr[i] >= interval[0]){
            if (i%2 == 0){
                // This is the first point of overlap
                i_start = i;
            } else {
                // The vertical restriction is the first point of overlap
                i_start = i-1;
                arr[i_start] = interval[0];
            }
            break;
        }
    }
    // Find the end of overlap
    int8_t i_stop = *length - 1;
    for (int8_t i=0; i<*length; i++){
        if (arr[i] >= interval[1]){
            if (i%2 == 0){
                // The previous point is the last point of overlap
                i_stop = i-1;
            } else {
                // The vertical restriction is the first point of overlap
                i_stop = i;
                arr[i_stop] = interval[1];
            }
            break;
        }
    }
    *length = i_stop - i_start + 1;
    if (i_start > 0){
        for (int8_t i=0; i<*length; i++){
            arr[i] = arr[i+i_start];
        }
    }
}


#endif /* XCOLL_GEOM_METHODS_H */
// copyright ############################### #
// This file is part of the Xcoll package.   #
// Copyright (c) CERN, 2025.                 #
// ######################################### #

#ifndef XCOLL_GEOM_FIND_ROOT_H
#define XCOLL_GEOM_FIND_ROOT_H

#include <math.h>
#include <stdio.h>


 static inline
int8_t BoundingBox_overlaps(BoundingBox b1, BoundingBox b2){
    // v1-v4 are the four vertices of the first box in counterclockwise order
    // w1-w4 are the four vertices of the second box in counterclockwise order
    // e1-e2 are the two axes of the first box
    // f1-f2 are the two axes of the second box
    double cos_tb_b1 = BoundingBox_get_cos_tb(b1);
    double sin_tb_b1 = BoundingBox_get_sin_tb(b1);
    double cos_tC_b1 = BoundingBox_get_cos_tC(b1);
    double sin_tC_b1 = BoundingBox_get_sin_tC(b1);
    double l_b1      = BoundingBox_get_l(b1);
    double w_b1      = BoundingBox_get_w(b1);
    double rC_b1     = BoundingBox_get_rC(b1);
    double proj_l_b1 = BoundingBox_get_proj_l(b1);
    double proj_w_b1 = BoundingBox_get_proj_w(b1);

    double cos_tb_b2 = BoundingBox_get_cos_tb(b2);
    double sin_tb_b2 = BoundingBox_get_sin_tb(b2);
    double cos_tC_b2 = BoundingBox_get_cos_tC(b2);
    double sin_tC_b2 = BoundingBox_get_sin_tC(b2);
    double l_b2      = BoundingBox_get_l(b2);
    double w_b2      = BoundingBox_get_w(b2);
    double rC_b2     = BoundingBox_get_rC(b2);
    double proj_l_b2 = BoundingBox_get_proj_l(b2);
    double proj_w_b2 = BoundingBox_get_proj_w(b2);

    double projs[4];
    double cos_tb1_tb2 = cos_tb_b1 * cos_tb_b2 + sin_tb_b1 * sin_tb_b2;  // cos (tb1 - tb2)
    double sin_tb1_tb2 = sin_tb_b1 * cos_tb_b2 - cos_tb_b1 * sin_tb_b2;  // sin (tb1 - tb2)
    // The length of the projection of vector W on E is given by: |W| cos(theta) where theta is the angle between W and E
    // projection of vertices of box 2 on the length axis of box 1 (e1)
    projs[0] = rC_b2 * (cos_tb_b1 * cos_tC_b2 + sin_tb_b1 * sin_tC_b2);  // first vertex w1:  |w1| cos (tb1 - tc2)
    projs[1] = projs[0] + l_b2 * cos_tb1_tb2; // second vertex w2 = w1 + f1
    projs[2] = projs[1] + w_b2 * sin_tb1_tb2; // third vertex w3 = w1 + f1 + f2
    projs[3] = projs[0] + w_b2 * sin_tb1_tb2; // fourth vertex w4 = w1 + f2
    sort_array_of_4_double(projs);
    if (!INTERVALS_OVERLAP(proj_l_b1, proj_l_b1 + l_b1, projs[0], projs[3])){
        return 0;
    }

    // length of projection of vertices of box 2 on the width axis of box 1 (e2)
    projs[0] = rC_b2 * (cos_tb_b1 * sin_tC_b2 - sin_tb_b1 * cos_tC_b2);  // first vertex w1:  |w1| cos(tb1 + pi/2 - tc2)
    projs[1] = projs[0] + l_b2 * sin_tb1_tb2; // second vertex w2 = w1 + f1
    projs[2] = projs[1] + w_b2 * cos_tb1_tb2; // third vertex w3 = w1 + f1 + f2
    projs[3] = projs[0] + w_b2 * cos_tb1_tb2; // fourth vertex w4 = w1 + f2
    sort_array_of_4_double(projs);
    if (!INTERVALS_OVERLAP(proj_w_b1, proj_w_b1 + w_b1, projs[0], projs[3])){ 
        return 0;
    }

    // length of projection of vertices of box 1 on the length axis of box 2 (f1)
    sin_tb1_tb2 = - sin_tb1_tb2; // due to sin being asymmetric
    projs[0] = rC_b1 * (cos_tb_b2 * cos_tC_b1 + sin_tb_b2 * sin_tC_b1);  // first vertex v1:  |v1| cos (tb1 - tc2)
    projs[1] = projs[0] + l_b1 * cos_tb1_tb2; // second vertex v2 = v1 + e1
    projs[2] = projs[1] + w_b1 * sin_tb1_tb2; // third vertex v3 = v1 + e1 + e2
    projs[3] = projs[0] + w_b1 * sin_tb1_tb2; // fourth vertex v4 = v1 + e2
    sort_array_of_4_double(projs);
    if (!INTERVALS_OVERLAP(proj_l_b2, proj_l_b2 + l_b2, projs[0], projs[3])){ 
        return 0;
    }

    // length of projection of vertices of box 1 on the width axis of box 2 (f2)
    projs[0] = rC_b1 * (cos_tb_b2 * sin_tC_b1 - sin_tb_b2 * cos_tC_b1);  // first vertex v1:  |v1| cos(tb1 + pi/2 - tc2)
    projs[1] = projs[0] + l_b1 * sin_tb1_tb2; // second vertex v2 = v1 + e1
    projs[2] = projs[1] + w_b1 * cos_tb1_tb2; // third vertex v3 = v1 + e1 + e2
    projs[3] = projs[0] + w_b1 * cos_tb1_tb2; // fourth vertex v4 = v1 + e2 
    sort_array_of_4_double(projs);
    if (!INTERVALS_OVERLAP(proj_w_b2, proj_w_b2 + w_b2, projs[0], projs[3])){ 
        return 0;
    }
    return 1;
}

//  static inline
// void grid_search_and_newton(LocalTrajectory traj, LocalSegment seg,
//                             double s_min, double s_max, double* roots, double max_crossings,
//                             int8_t* number_of_roots) {
//     /// Find the intervals where the function changes sign within the range [s_min, s_max]
//     //  in which later Newton's method can be applied to find the root(s) for each interval
//     double grid_step = (s_max - s_min) / XC_GRID_POINTS;
//     int interval_count = 0;

//     double prev_s   = s_min;
//     double prev_val = LocalTrajectory_func(traj, prev_s) - LocalSegment_func(seg, prev_s);

//     for (int i = 1; i <= XC_GRID_POINTS - 1; i++) {
//         if (interval_count >= max_crossings) break; // you cannot have more intervals than roots
//         double curr_s = s_min + i * grid_step;
//         double curr_val = LocalTrajectory_func(traj, curr_s) - LocalSegment_func(seg, curr_s);
//         if (prev_val * curr_val < 0) {
//             double initial_guess = 0.5 * (prev_s + curr_s);      // initial guess is midpoint
//             roots[interval_count] = newton(traj, seg, initial_guess);
//             interval_count++;
//         }
//         prev_s = curr_s;
//         prev_val = curr_val;
//     }
//     *number_of_roots = interval_count;
// }

//  static inline
// double newton(LocalTrajectory traj, LocalSegment seg, double guess_s) {
//     for (int i = 0; i < XC_NEWTON_MAX_ITER; i++) {
//         double f = LocalTrajectory_func(traj, guess_s) - LocalSegment_func(seg, guess_s);
//         double f_prime = LocalTrajectory_deriv(traj, guess_s) - LocalSegment_deriv(seg, guess_s);
//         if (fabs(f_prime) < XC_NEWTON_DERIVATIVE_TOL) return guess_s;

//         double guess_new = guess_s - f / f_prime;
//         if (fabs(guess_new - guess_s) < XC_NEWTON_EPSILON) return guess_new;
//         guess_s = guess_new;
//     }
//     return guess_s;
// }

// // --------------------------------------------------------------------------------------------

// void LocalCrossing_func(LocalTrajectory traj, LocalSegment seg, double TS[2], double l, double t){
//     // Here we get the expr from segment and traj, and we connect them to create TSs and TSx
//     // TSs = Ts(l) - Ss(t)
//     // TSx = Tx(l) - Sx(t)
//     TS[0] = LocalTrajectory_func_s(traj, l) - LocalSegment_func_s(seg, t);
//     TS[1] = LocalTrajectory_func_x(traj, l) - LocalSegment_func_x(seg, t);
// }

// void LocalCrossing_inv_J(LocalTrajectory traj, LocalSegment seg, double J_inv[2][2], int8_t* no_crossing, double l, double t){
//     double J[2][2];
//     // get derivatives
//     J[0][0] = LocalTrajectory_deriv_s(traj,l); // get deriv. dsT/dl LocalTrajectory_deriv_s
//     J[0][1] = -LocalSegment_deriv_s(seg, t);     // get deriv. dsS/dt LocalSegment_deriv_s
//     J[1][0] = LocalTrajectory_deriv_x(traj,l); // get deriv. dxT/dl LocalTrajectory_deriv_x
//     J[1][1] = -LocalSegment_deriv_x(seg, t);     // get deriv. dxS/dt LocalSegment_deriv_x

//     double det = J[0][0] * J[1][1] - J[1][0]*J[0][1];
//     if (fabs(det) < XC_GEOM_ROOT_NEWTON_DERIVATIVE_TOL){
//         printf("There is no crossing. \n");
//         *no_crossing = 1;
//         return;
//     }
//     J_inv[0][0] =  J[1][1] / det;
//     J_inv[0][1] = -J[0][1] / det;
//     J_inv[1][0] = -J[1][0] / det;
//     J_inv[1][1] =  J[0][0] / det;
// }

// void newton(LocalTrajectory traj, LocalSegment seg, double* guess_l, double* guess_t, int8_t* no_crossing){
//     double J_inv[2][2];
//     double TS[2];

//     for (int i = 0; i < XC_GEOM_ROOT_NEWTON_MAX_ITER; i++){
//         F_G(traj, seg, TS, *guess_l, *guess_t);
//         get_inv_J(traj, seg, J_inv, no_crossing, *guess_l, *guess_t);
//         if (*no_crossing){
//             return;
//         }
//         double new_t = *guess_t - (J_inv[0][0]*TS[0] + J_inv[0][1]*TS[1]);
//         double new_l = *guess_l - (J_inv[1][0]*TS[0] + J_inv[1][1]*TS[1]);

//         // Check for convergence
//         if ((fabs(new_t -  *guess_t) < XC_GEOM_ROOT_NEWTON_EPSILON) && (fabs(new_l - *guess_l) < XC_GEOM_ROOT_NEWTON_EPSILON)){
//             return;
//         }
//         // Update the guesses for the next iteration
//         *guess_t = new_t;  // Keep *t updated
//         *guess_l = new_l;  // Keep *l updated
//     }
//     *no_crossing = 1;
// }

// int8_t LocalCrossing_box_has_root(double TS_UL[2], double TS_UR[2], double TS_DL[2], double TS_DR[2]){

//     // Evaluate F at corners
//     evaluate_traj_seg(traj, seg, TS_DL, l1, t1);
//     evaluate_traj_seg(traj, seg, TS_DR, l2, t1);
//     evaluate_traj_seg(traj, seg, TS_UL, l1, t2);
//     evaluate_traj_seg(traj, seg, TS_UR, l2, t2);

//     // Check for sign changes for each component across the corners.
//     // (For simplicity, here we check if the min and max differ in sign.)
//     int8_t f1_has_sign_change = (min(TS_DL[0], TS_DR[0], TS_UL[0], TS_UR[0]) < 0 &&
//                                 max(TS_DL[0], TS_DR[0], TS_UL[0], TS_UR[0]) > 0);
//     int8_t f2_has_sign_change = (min(TS_DL[1], TS_DR[1], TS_UL[1], TS_UR[1]) < 0 &&
//                                 max(TS_DL[1], TS_DR[1], TS_UL[1], TS_UR[1]) > 0);
//     return (l1 < 0 && l2 > 0) || (l1 > 0 && l2 < 0) || (t1 < 0 && t2 > 0) || (t1 > 0 && t2 < 0);
// }


// void grid_search_and_newton(LocalTrajectory traj, LocalSegment seg, double* l, double* t, double l_min,
//                             double l_max, double t_min, double t_max, double* roots_l, double* roots_t,
//                             int max_crossings, int* number_of_roots){
//     int N_l = 100;
//     int N_t = 100;
//     double grid_step_l = (l_max - l_min) / N_l; // this is just for now. We need to get interval for t and l.
//     double grid_step_t = (t_max - t_min) / N_t;
//     int n_roots = 0;
//     double TS_prev, TS_curr;
//     double prev_t = t_min;
//     double prev_l = l_min;
//     double curr_t;
//     double curr_l;
//     int8_t no_crossing = 0;

//     LocalCrossing_func(traj, seg, TS_prev, prev_t, prev_l);

//     // Stage 1: Coarse grid search over (l, t)
//     for (int i=0; i < N_l-1; i++){
//         double l1 = l_min + i * grid_step_l;
//         double l2 = l_min + (i+1) * grid_step_l;
//         for (int j=0; j < N_t-1; j++){
//             double t1 = t_min + j * grid_step_t;
//             double t2 = t_min + (j+1) * grid_step_t;
//             F_G(traj, seg, FG_curr, curr_t, curr_l);





//             if (n_roots >= max_crossings) {
//                 return;  // all possible roots have been found
//             }
//             curr_t = t_min + t_step;
//             curr_l = l_min + l_step;
//             F_G(traj, seg, FG_curr, curr_t, curr_l);

//             if ((FG_prev[0] * FG_curr[0] < 0) || (FG_prev[1] * FG_curr[1] < 0)) {
//                 double initial_guess_t = 0.5*(curr_t - prev_t);
//                 double initial_guess_l = 0.5*(curr_l - prev_l);
//                 newton(traj, seg, &initial_guess_l, &initial_guess_t, &no_crossing);
//                 if (no_crossing){
//                     printf("No crossing found with Newton's method. \n");
//                 } else{
//                     roots_t[n_roots] = initial_guess_t;
//                     roots_l[n_roots] = initial_guess_l;
//                     n_roots++;
//                 }
//             }
//         }
//     }
// }

//  double prev_s   = s_min;
//     double prev_val = LocalTrajectory_func(traj, prev_s) - LocalSegment_func(seg, prev_s);

//     for (int i = 1; i <= XC_GRID_POINTS - 1; i++) {
//         if (interval_count >= max_crossings) break; // you cannot have more intervals than roots
//         double curr_s = s_min + i * grid_step;
//         double curr_val = LocalTrajectory_func(traj, curr_s) - LocalSegment_func(seg, curr_s);
//         if (prev_val * curr_val < 0) {
//             double initial_guess = 0.5 * (prev_s + curr_s);      // initial guess is midpoint
//             roots[interval_count] = newton(traj, seg, initial_guess);
//             interval_count++;
//         }
//         prev_s = curr_s;
//         prev_val = curr_val;


// int main() {
//     double s0 = 1.5;
//     double x0 = 0.0;
//     double theta = 20.0*M_PI/180.0;
//     double s1 = 1.0;
//     double x1 =-1.0;
//     double s2 =2.0;
//     double x2 =2.0;
//     double l = -1.0;
//     double t = -1.0;
//     Params p = {s0, x0, theta, s1, x1, s2, x2, l, t};
//     double guess_t = 0.5;
//     double guess_l = 0.5;
//     newton(&guess_t, &guess_l, &p);
//     return 0;
// }

#endif /* XCOLL_GEOM_FIND_ROOT_H */