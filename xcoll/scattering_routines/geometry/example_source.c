#include <stdint.h>
#ifndef XOBJ_TYPEDEF_Arr2Float64
#define XOBJ_TYPEDEF_Arr2Float64
typedef   struct Arr2Float64_s * Arr2Float64;
 static inline Arr2Float64 Arr2Float64_getp(Arr2Float64 restrict  obj){
  int64_t offset=0;
  return (Arr2Float64)(( char*) obj+offset);
}
 static inline int64_t Arr2Float64_len(Arr2Float64 restrict  obj){
  return 2;
}
 static inline double Arr2Float64_get(const Arr2Float64 restrict  obj, int64_t i0){
  int64_t offset=0;
  offset+=i0*8;
  return *( double*)(( char*) obj+offset);
}
 static inline void Arr2Float64_set(Arr2Float64 restrict  obj, int64_t i0, double value){
  int64_t offset=0;
  offset+=i0*8;
  *( double*)(( char*) obj+offset)=value;
}
 static inline  double* Arr2Float64_getp1(Arr2Float64 restrict  obj, int64_t i0){
  int64_t offset=0;
  offset+=i0*8;
  return ( double*)(( char*) obj+offset);
}
#endif
#ifndef XOBJ_TYPEDEF_GeomCInit
#define XOBJ_TYPEDEF_GeomCInit
typedef   struct GeomCInit_s * GeomCInit;
 static inline GeomCInit GeomCInit_getp(GeomCInit restrict  obj){
  int64_t offset=0;
  return (GeomCInit)(( char*) obj+offset);
}
#endif

#ifndef XCOLL_GEOM_DEFINES_H
#define XCOLL_GEOM_DEFINES_H
#include <math.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>

#ifndef XC_EPSILON
#define XC_EPSILON 1e-15
#endif

#ifndef XC_S_MAX
#define XC_S_MAX 1e+21
#endif

#endif /* XCOLL_GEOM_DEFINES_H */

// copyright ############################### #
// This file is part of the Xcoll package.   #
// Copyright (c) CERN, 2024.                 #
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
// Copyright (c) CERN, 2024.                 #
// ######################################### #

#ifndef XCOLL_GEOM_METHODS_H
#define XCOLL_GEOM_METHODS_H

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
        arr[*length] = XC_S_MAX*0.1;
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
// Copyright (c) CERN, 2024.                 #
// ######################################### #

#ifndef XCOLL_GEOM_FIND_ROOT_H
#define XCOLL_GEOM_FIND_ROOT_H

#define XC_NEWTON_EPSILON 1.e-10
#define XC_NEWTON_MAX_ITER 100  // Maximum number of iterations
#define XC_NEWTON_DERIVATIVE_TOL 1e-10  // Threshold for small derivative


// Function to compute the value of the cubic equation at t
 static inline
double _bezier_cubic(double t, double a, double b, double c, double d) {
    return a*t*t*t + b*t*t + c*t + d;
}

// Derivative of the cubic function
 static inline
double _bezier_cubic_prime(double t, double a, double b, double c) {
    return 3*a*t*t + 2*b*t + c;
}

// Newton's method to solve the cubic equation with safety checks
// TODO: This finds one root, but we want all. Need to find a smart combination of newton and bisection
// TODO: Need to find a way to catch oscillatory behaviour and start from another guess..
//       Example: f(x) = x3 − 2x + 2 will oscillate between 0 and 1
// TODO: Need a catch for undefined values (e.g. when derivative has a 1/x).
//       Similarily, the iteration can send itself beyond the valid domain
//       Example: f(x) = ln x  => x_(n+1) = xn(1- ln xn). If one starts at x >= e, it will get NaN
 static inline
double newton_method(double a, double b, double c, double d, double initial_guess, int *status) {
    double t = initial_guess;
    double t_new;
    int iter = 0;

    while (fabs(_bezier_cubic(t, a, b, c, d)) > XC_NEWTON_EPSILON) {
        double derivative = _bezier_cubic_prime(t, a, b, c);

        // Check if the derivative is too small
        if (fabs(derivative) < XC_NEWTON_DERIVATIVE_TOL) {
            *status = -1;  // Indicate failure due to small derivative
            return t;       // Return the best guess we have
        }

        t_new = t - _bezier_cubic(t, a, b, c, d) / derivative;

        // Check if the change in t is very small
        if (fabs(t_new - t) < XC_NEWTON_EPSILON) break;

        t = t_new;
        iter++;

        // Check if we exceed maximum iterations
        if (iter >= XC_NEWTON_MAX_ITER) {
            *status = -2;  // Indicate failure due to max iterations
            return t;
        }
    }

    *status = 0;  // Indicate success
    return t;
}

// Bisection method to find root in an interval [t1, t2]
 static inline
double bisection_method(double a, double b, double c, double d, double t1, double t2) {
    double mid;

    while ((t2 - t1) >= XC_NEWTON_EPSILON) {
        mid = (t1 + t2) / 2;

        // Check if the middle point is a root
        if (fabs(_bezier_cubic(mid, a, b, c, d)) < XC_NEWTON_EPSILON) return mid;

        // Decide the side to repeat the steps
        if (_bezier_cubic(mid, a, b, c, d) * _bezier_cubic(t1, a, b, c, d) < 0)
            t2 = mid;
        else
            t1 = mid;
    }

    return mid;
}

#endif /* XCOLL_GEOM_FIND_ROOT_H */
#ifndef XOBJ_TYPEDEF_DriftTrajectory
#define XOBJ_TYPEDEF_DriftTrajectory
typedef   struct DriftTrajectory_s * DriftTrajectory;
 static inline DriftTrajectory DriftTrajectory_getp(DriftTrajectory restrict  obj){
  int64_t offset=0;
  return (DriftTrajectory)(( char*) obj+offset);
}
#endif
// copyright ############################### #
// This file is part of the Xcoll package.   #
// Copyright (c) CERN, 2024.                 #
// ######################################### #

#ifndef XCOLL_COLL_GEOM_DRIFT_H
#define XCOLL_COLL_GEOM_DRIFT_H


 static inline
int8_t DriftTrajectory_vlimit(double* restrict_s, double s0, double y0, double ym, double ymin, double ymax){
    if (fabs(ym) < XC_EPSILON){
        // Trajectory parallel to s axis
        if (y0 < ymin || y0 > ymax){
            return 0;  // Completely outside - no crossing possible
        } else {
            return -1; // Completely inside - no vertical check needed
        }
    } else {
        restrict_s[0] = (ymin - y0)/ym + s0;
        restrict_s[1] = (ymax - y0)/ym + s0;
        SWAP(restrict_s, 0, 1);   // To make sure these are sorted
        return 1;  // Default behavior: check overlap with horizontal crossings
    }
}


// Curve length is int_s1^s2 sqrt(1 + (dx/ds)^2 + (dy/ds)^2) ds

 static inline
double DriftTrajectory_length(double s0, double x0, double xm, double y0, double ym, double s1, double s2){
    (void) s0;  // Avoid unused parameter warning
    (void) x0;  // Avoid unused parameter warning
    (void) y0;  // Avoid unused parameter warning
    return (s2-s1)*sqrt(1 + xm*xm + ym*ym);
}


#endif /* XCOLL_COLL_GEOM_DRIFT_H */
#ifndef XOBJ_TYPEDEF_LineSegment
#define XOBJ_TYPEDEF_LineSegment
typedef   struct LineSegment_s * LineSegment;
 static inline LineSegment LineSegment_getp(LineSegment restrict  obj){
  int64_t offset=0;
  return (LineSegment)(( char*) obj+offset);
}
 static inline double LineSegment_get_s1(const LineSegment restrict  obj){
  int64_t offset=0;
  return *( double*)(( char*) obj+offset);
}
 static inline void LineSegment_set_s1(LineSegment restrict  obj, double value){
  int64_t offset=0;
  *( double*)(( char*) obj+offset)=value;
}
 static inline  double* LineSegment_getp_s1(LineSegment restrict  obj){
  int64_t offset=0;
  return ( double*)(( char*) obj+offset);
}
 static inline double LineSegment_get_x1(const LineSegment restrict  obj){
  int64_t offset=0;
  offset+=8;
  return *( double*)(( char*) obj+offset);
}
 static inline void LineSegment_set_x1(LineSegment restrict  obj, double value){
  int64_t offset=0;
  offset+=8;
  *( double*)(( char*) obj+offset)=value;
}
 static inline  double* LineSegment_getp_x1(LineSegment restrict  obj){
  int64_t offset=0;
  offset+=8;
  return ( double*)(( char*) obj+offset);
}
 static inline double LineSegment_get_s2(const LineSegment restrict  obj){
  int64_t offset=0;
  offset+=16;
  return *( double*)(( char*) obj+offset);
}
 static inline void LineSegment_set_s2(LineSegment restrict  obj, double value){
  int64_t offset=0;
  offset+=16;
  *( double*)(( char*) obj+offset)=value;
}
 static inline  double* LineSegment_getp_s2(LineSegment restrict  obj){
  int64_t offset=0;
  offset+=16;
  return ( double*)(( char*) obj+offset);
}
 static inline double LineSegment_get_x2(const LineSegment restrict  obj){
  int64_t offset=0;
  offset+=24;
  return *( double*)(( char*) obj+offset);
}
 static inline void LineSegment_set_x2(LineSegment restrict  obj, double value){
  int64_t offset=0;
  offset+=24;
  *( double*)(( char*) obj+offset)=value;
}
 static inline  double* LineSegment_getp_x2(LineSegment restrict  obj){
  int64_t offset=0;
  offset+=24;
  return ( double*)(( char*) obj+offset);
}
#endif
// copyright ############################### #
// This file is part of the Xcoll package.   #
// Copyright (c) CERN, 2024.                 #
// ######################################### #

#ifndef XCOLL_COLL_GEOM_LINESEG_H
#define XCOLL_COLL_GEOM_LINESEG_H

 static inline
void LineSegment_crossing_drift(LineSegment seg, int8_t* n_hit, double* s, double s0, double x0, double xm){
    // Get segment data
    double s1 = LineSegment_get_s1(seg);
    double x1 = LineSegment_get_x1(seg);
    double s2 = LineSegment_get_s2(seg);
    double x2 = LineSegment_get_x2(seg);
    double denom = (x2 - x1) - (s2 - s1)*xm;
    if (fabs(denom) < XC_EPSILON){
        // Trajectory is parallel to the segment
        if (fabs((x0 - x1)/(s0 - s1) - xm) < XC_EPSILON){
            // Trajectory overlaps with the segment
            // TODO: This is situational; we should return s1 if get_s_first and current_s if after current_s
            //       For now we hit twice (because we go nor IN nor OUT)
            s[*n_hit] = s1;
            (*n_hit)++;
            s[*n_hit] = s2;
            (*n_hit)++;
        } else {
            // No crossing
            return;
        }
    } else {
        double t = (x0 - x1 - (s0 - s1)*xm) / denom;
        if (t >= 0 && t <= 1){
            s[*n_hit] = s1*(1-t) + s2*t;
            (*n_hit)++;
        }
    }
}

#endif /* XCOLL_COLL_GEOM_LINESEG_H */
#ifndef XOBJ_TYPEDEF_HalfOpenLineSegment
#define XOBJ_TYPEDEF_HalfOpenLineSegment
typedef   struct HalfOpenLineSegment_s * HalfOpenLineSegment;
 static inline HalfOpenLineSegment HalfOpenLineSegment_getp(HalfOpenLineSegment restrict  obj){
  int64_t offset=0;
  return (HalfOpenLineSegment)(( char*) obj+offset);
}
 static inline double HalfOpenLineSegment_get_s(const HalfOpenLineSegment restrict  obj){
  int64_t offset=0;
  return *( double*)(( char*) obj+offset);
}
 static inline void HalfOpenLineSegment_set_s(HalfOpenLineSegment restrict  obj, double value){
  int64_t offset=0;
  *( double*)(( char*) obj+offset)=value;
}
 static inline  double* HalfOpenLineSegment_getp_s(HalfOpenLineSegment restrict  obj){
  int64_t offset=0;
  return ( double*)(( char*) obj+offset);
}
 static inline double HalfOpenLineSegment_get_x(const HalfOpenLineSegment restrict  obj){
  int64_t offset=0;
  offset+=8;
  return *( double*)(( char*) obj+offset);
}
 static inline void HalfOpenLineSegment_set_x(HalfOpenLineSegment restrict  obj, double value){
  int64_t offset=0;
  offset+=8;
  *( double*)(( char*) obj+offset)=value;
}
 static inline  double* HalfOpenLineSegment_getp_x(HalfOpenLineSegment restrict  obj){
  int64_t offset=0;
  offset+=8;
  return ( double*)(( char*) obj+offset);
}
 static inline double HalfOpenLineSegment_get_t(const HalfOpenLineSegment restrict  obj){
  int64_t offset=0;
  offset+=16;
  return *( double*)(( char*) obj+offset);
}
 static inline void HalfOpenLineSegment_set_t(HalfOpenLineSegment restrict  obj, double value){
  int64_t offset=0;
  offset+=16;
  *( double*)(( char*) obj+offset)=value;
}
 static inline  double* HalfOpenLineSegment_getp_t(HalfOpenLineSegment restrict  obj){
  int64_t offset=0;
  offset+=16;
  return ( double*)(( char*) obj+offset);
}
#endif
// copyright ############################### #
// This file is part of the Xcoll package.   #
// Copyright (c) CERN, 2024.                 #
// ######################################### #

#ifndef XCOLL_COLL_GEOM_HALFOPENLINESEG_H
#define XCOLL_COLL_GEOM_HALFOPENLINESEG_H

 static inline
void HalfOpenLineSegment_crossing_drift(HalfOpenLineSegment seg, int8_t* n_hit, double* s, double s0, double x0, double xm){
    // Get segment data
    double s1 = HalfOpenLineSegment_get_s(seg);
    double x1 = HalfOpenLineSegment_get_x(seg);
    double s2 = s1 + cos(HalfOpenLineSegment_get_t(seg));
    double x2 = x1 + sin(HalfOpenLineSegment_get_t(seg));
    double denom = (x2 - x1) - (s2 - s1)*xm;
    if (fabs(denom) < XC_EPSILON){
        // Trajectory is parallel to the segment
        if (fabs((x0 - x1)/(s0 - s1) - xm) < XC_EPSILON){
            // Trajectory overlaps with the segment
            // TODO: This is situational; we should return s1 if get_s_first and current_s if after current_s
            //       For now we hit twice (because we go nor IN nor OUT)
            s[*n_hit] = s1;
            (*n_hit)++;
            s[*n_hit] = s1;
            (*n_hit)++;
        } else {
            // No hit
            return;
        }
    } else {
        double t = (x0 - x1 - (s0 - s1)*xm) / denom;
        if (t >= 0){  // We do not check for t<=1 as it is a half-open segment
            s[*n_hit] = s1*(1-t) + s2*t;
            (*n_hit)++;
        }
    }
}

#endif /* XCOLL_COLL_GEOM_HALFOPENLINESEG_H */
#ifndef XOBJ_TYPEDEF_CircularSegment
#define XOBJ_TYPEDEF_CircularSegment
typedef   struct CircularSegment_s * CircularSegment;
 static inline CircularSegment CircularSegment_getp(CircularSegment restrict  obj){
  int64_t offset=0;
  return (CircularSegment)(( char*) obj+offset);
}
 static inline double CircularSegment_get_R(const CircularSegment restrict  obj){
  int64_t offset=0;
  return *( double*)(( char*) obj+offset);
}
 static inline void CircularSegment_set_R(CircularSegment restrict  obj, double value){
  int64_t offset=0;
  *( double*)(( char*) obj+offset)=value;
}
 static inline  double* CircularSegment_getp_R(CircularSegment restrict  obj){
  int64_t offset=0;
  return ( double*)(( char*) obj+offset);
}
 static inline double CircularSegment_get_s(const CircularSegment restrict  obj){
  int64_t offset=0;
  offset+=8;
  return *( double*)(( char*) obj+offset);
}
 static inline void CircularSegment_set_s(CircularSegment restrict  obj, double value){
  int64_t offset=0;
  offset+=8;
  *( double*)(( char*) obj+offset)=value;
}
 static inline  double* CircularSegment_getp_s(CircularSegment restrict  obj){
  int64_t offset=0;
  offset+=8;
  return ( double*)(( char*) obj+offset);
}
 static inline double CircularSegment_get_x(const CircularSegment restrict  obj){
  int64_t offset=0;
  offset+=16;
  return *( double*)(( char*) obj+offset);
}
 static inline void CircularSegment_set_x(CircularSegment restrict  obj, double value){
  int64_t offset=0;
  offset+=16;
  *( double*)(( char*) obj+offset)=value;
}
 static inline  double* CircularSegment_getp_x(CircularSegment restrict  obj){
  int64_t offset=0;
  offset+=16;
  return ( double*)(( char*) obj+offset);
}
 static inline double CircularSegment_get_t1(const CircularSegment restrict  obj){
  int64_t offset=0;
  offset+=24;
  return *( double*)(( char*) obj+offset);
}
 static inline void CircularSegment_set_t1(CircularSegment restrict  obj, double value){
  int64_t offset=0;
  offset+=24;
  *( double*)(( char*) obj+offset)=value;
}
 static inline  double* CircularSegment_getp_t1(CircularSegment restrict  obj){
  int64_t offset=0;
  offset+=24;
  return ( double*)(( char*) obj+offset);
}
 static inline double CircularSegment_get_t2(const CircularSegment restrict  obj){
  int64_t offset=0;
  offset+=32;
  return *( double*)(( char*) obj+offset);
}
 static inline void CircularSegment_set_t2(CircularSegment restrict  obj, double value){
  int64_t offset=0;
  offset+=32;
  *( double*)(( char*) obj+offset)=value;
}
 static inline  double* CircularSegment_getp_t2(CircularSegment restrict  obj){
  int64_t offset=0;
  offset+=32;
  return ( double*)(( char*) obj+offset);
}
#endif
// copyright ############################### #
// This file is part of the Xcoll package.   #
// Copyright (c) CERN, 2024.                 #
// ######################################### #

#ifndef XCOLL_COLL_GEOM_CIRCULARSEG_H
#define XCOLL_COLL_GEOM_CIRCULARSEG_H

 static inline
void CircularSegment_crossing_drift(CircularSegment seg, int8_t* n_hit, double* s, double s0, double x0, double xm){
    // Get segment data
    double R  = CircularSegment_get_R(seg);
    double sC = CircularSegment_get_s(seg);
    double xC = CircularSegment_get_x(seg);
    double t1 = CircularSegment_get_t1(seg);
    double t2 = CircularSegment_get_t2(seg);
    // Move the angles to [-pi, pi]
    int8_t reversed = 0, full_circle = 0;
    if (fabs(fabs(t2 - t1) - 2*M_PI) < XC_EPSILON){
        full_circle = 1;
    }
    while (t1 < -M_PI){
        t1 += 2*M_PI;
    }
    while (t1 > M_PI){
        t1 -= 2*M_PI;
    }
    while (t2 < -M_PI){
        t2 += 2*M_PI;
    }
    while (t2 > M_PI){
        t2 -= 2*M_PI;
    }
    if (t2 < t1){
        reversed = 1;
    }
    // Calculate crossings
    double a = 1 + xm*xm;
    double bb = sC - xm*(x0 - xC - xm*s0); // This is -b/2 with b from the quadratic formula
    double c = sC*sC + (x0 - xC - xm*s0)*(x0 - xC - xm*s0) - R*R;
    double disc = bb*bb - a*c; // This is  2*discriminant**2
    if (disc < 0){
        // No crossing
        return;
    }
    for (int8_t i = 0; i < 2; i++) {
        double sgnD = i*2-1; // negative and positive solutions; if multiplicity 2, we add the same solution twice
        double new_s = (bb + sgnD*sqrt(fabs(disc)))/a;
        double new_x = x0 + (new_s - s0)*xm;
        double t = atan2(new_x - xC, new_s - sC);
        if (full_circle){
            // Full circle, so always hit
            s[*n_hit] = new_s;
            (*n_hit)++;
        } else if (reversed){
            // t2 < t1, so we are looking at the inverted region of angles
            if (t1 <= t || t <= t2){
                s[*n_hit] = new_s;
                (*n_hit)++;
            }
        } else {
            if (t1 <= t && t <= t2){
                s[*n_hit] = new_s;
                (*n_hit)++;
            }
        }
    }
}

#endif /* XCOLL_COLL_GEOM_CIRCULARSEG_H */
#ifndef XOBJ_TYPEDEF_BezierSegment
#define XOBJ_TYPEDEF_BezierSegment
typedef   struct BezierSegment_s * BezierSegment;
 static inline BezierSegment BezierSegment_getp(BezierSegment restrict  obj){
  int64_t offset=0;
  return (BezierSegment)(( char*) obj+offset);
}
 static inline double BezierSegment_get_s1(const BezierSegment restrict  obj){
  int64_t offset=0;
  return *( double*)(( char*) obj+offset);
}
 static inline void BezierSegment_set_s1(BezierSegment restrict  obj, double value){
  int64_t offset=0;
  *( double*)(( char*) obj+offset)=value;
}
 static inline  double* BezierSegment_getp_s1(BezierSegment restrict  obj){
  int64_t offset=0;
  return ( double*)(( char*) obj+offset);
}
 static inline double BezierSegment_get_x1(const BezierSegment restrict  obj){
  int64_t offset=0;
  offset+=8;
  return *( double*)(( char*) obj+offset);
}
 static inline void BezierSegment_set_x1(BezierSegment restrict  obj, double value){
  int64_t offset=0;
  offset+=8;
  *( double*)(( char*) obj+offset)=value;
}
 static inline  double* BezierSegment_getp_x1(BezierSegment restrict  obj){
  int64_t offset=0;
  offset+=8;
  return ( double*)(( char*) obj+offset);
}
 static inline double BezierSegment_get_s2(const BezierSegment restrict  obj){
  int64_t offset=0;
  offset+=16;
  return *( double*)(( char*) obj+offset);
}
 static inline void BezierSegment_set_s2(BezierSegment restrict  obj, double value){
  int64_t offset=0;
  offset+=16;
  *( double*)(( char*) obj+offset)=value;
}
 static inline  double* BezierSegment_getp_s2(BezierSegment restrict  obj){
  int64_t offset=0;
  offset+=16;
  return ( double*)(( char*) obj+offset);
}
 static inline double BezierSegment_get_x2(const BezierSegment restrict  obj){
  int64_t offset=0;
  offset+=24;
  return *( double*)(( char*) obj+offset);
}
 static inline void BezierSegment_set_x2(BezierSegment restrict  obj, double value){
  int64_t offset=0;
  offset+=24;
  *( double*)(( char*) obj+offset)=value;
}
 static inline  double* BezierSegment_getp_x2(BezierSegment restrict  obj){
  int64_t offset=0;
  offset+=24;
  return ( double*)(( char*) obj+offset);
}
 static inline double BezierSegment_get_cs1(const BezierSegment restrict  obj){
  int64_t offset=0;
  offset+=32;
  return *( double*)(( char*) obj+offset);
}
 static inline void BezierSegment_set_cs1(BezierSegment restrict  obj, double value){
  int64_t offset=0;
  offset+=32;
  *( double*)(( char*) obj+offset)=value;
}
 static inline  double* BezierSegment_getp_cs1(BezierSegment restrict  obj){
  int64_t offset=0;
  offset+=32;
  return ( double*)(( char*) obj+offset);
}
 static inline double BezierSegment_get_cx1(const BezierSegment restrict  obj){
  int64_t offset=0;
  offset+=40;
  return *( double*)(( char*) obj+offset);
}
 static inline void BezierSegment_set_cx1(BezierSegment restrict  obj, double value){
  int64_t offset=0;
  offset+=40;
  *( double*)(( char*) obj+offset)=value;
}
 static inline  double* BezierSegment_getp_cx1(BezierSegment restrict  obj){
  int64_t offset=0;
  offset+=40;
  return ( double*)(( char*) obj+offset);
}
 static inline double BezierSegment_get_cs2(const BezierSegment restrict  obj){
  int64_t offset=0;
  offset+=48;
  return *( double*)(( char*) obj+offset);
}
 static inline void BezierSegment_set_cs2(BezierSegment restrict  obj, double value){
  int64_t offset=0;
  offset+=48;
  *( double*)(( char*) obj+offset)=value;
}
 static inline  double* BezierSegment_getp_cs2(BezierSegment restrict  obj){
  int64_t offset=0;
  offset+=48;
  return ( double*)(( char*) obj+offset);
}
 static inline double BezierSegment_get_cx2(const BezierSegment restrict  obj){
  int64_t offset=0;
  offset+=56;
  return *( double*)(( char*) obj+offset);
}
 static inline void BezierSegment_set_cx2(BezierSegment restrict  obj, double value){
  int64_t offset=0;
  offset+=56;
  *( double*)(( char*) obj+offset)=value;
}
 static inline  double* BezierSegment_getp_cx2(BezierSegment restrict  obj){
  int64_t offset=0;
  offset+=56;
  return ( double*)(( char*) obj+offset);
}
#endif
// copyright ############################### #
// This file is part of the Xcoll package.   #
// Copyright (c) CERN, 2024.                 #
// ######################################### #

#ifndef XCOLL_COLL_GEOM_BEZIERSEG_H
#define XCOLL_COLL_GEOM_BEZIERSEG_H

 static inline
void _hit_s_bezier(BezierSegment seg, double t, double multiplicity, int8_t* n_hit, double* s){
    double s1  = BezierSegment_get_s1(seg);
    double s2  = BezierSegment_get_s2(seg);
    double cs1 = BezierSegment_get_cs1(seg);
    double cs2 = BezierSegment_get_cs2(seg);
    double new_s = (1-t)*(1-t)*(1-t)*s1 + 3*(1-t)*(1-t)*t*cs1 + 3*(1-t)*t*t*cs2 + t*t*t*s2;
    for (int8_t i = 0; i < multiplicity; i++) {
        s[*n_hit] = new_s;
        (*n_hit)++;
    }
}

 static inline
void BezierSegment_crossing_drift(BezierSegment seg, int8_t* n_hit, double* s, double s0, double x0, double xm){
    // Get segment data
    double s1  = BezierSegment_get_s1(seg);
    double x1  = BezierSegment_get_x1(seg);
    double s2  = BezierSegment_get_s2(seg);
    double x2  = BezierSegment_get_x2(seg);
    double cs1 = BezierSegment_get_cs1(seg);
    double cx1 = BezierSegment_get_cx1(seg);
    double cs2 = BezierSegment_get_cs2(seg);
    double cx2 = BezierSegment_get_cx2(seg);
    // The Bézier curve is defined by the parametric equations (with t in [0, 1]):
    // s(t) = (1-t)^3*s1 + 3(1-t)^2*t*cs1 + 3(1-t)*t^2*cs2 + t^3*s2
    // x(t) = (1-t)^3*x1 + 3(1-t)^2*t*cx1 + 3(1-t)*t^2*cx2 + t^3*x2
    // Plug the parametric eqs into the drift trajectory x(t) = m*(s(t) - s0) + x0 and solve for t
    // The solutions for t (which we get by Cardano's method) are valid if in [0, 1]
    double a = (xm*s1 - x1) - (xm*s2 - x2) - 3*(xm*cs1 - cx1) + 3*(xm*cs2 - cx2);
    double b = 6*(xm*cs1 - cx1) - 3*(xm*cs2 - cx2) - 3*(xm*s1 - x1);
    double c = 3*(xm*s1 - x1) - 3*(xm*cs1 - cx1);
    double d = (xm*s0 - x0) - (xm*s1 - x1);
    double t;
    // Edge cases
    if (fabs(a) < XC_EPSILON){
        if (fabs(b) < XC_EPSILON){
            if (fabs(c) < XC_EPSILON){
                if (fabs(d) < XC_EPSILON){
                    // The trajectory is on the Bézier curve
                    // TODO: This cannot happen because we don't have a cubic trajectory.
                    //       Adapt if these ever would be implemented.
                    return;
                } else {
                    // No solutions
                    return;
                }
            } else {
                // This is a linear equation
                t = -d/c;
                if (0 <= t && t <= 1){
                    _hit_s_bezier(seg, t, 1, n_hit, s);
                }
            }
        } else {
            // This is a quadratic equation
            double disc = c*c - 4*b*d;
            if (disc < 0){
                // No solutions
                return;
            }
            for (int8_t i = 0; i < 2; i++) {
                double sgnD = i*2-1; // negative and positive solutions; if multiplicity 2, we add the same solution twice
                t = (-c + sgnD*sqrt(fabs(disc)))/(2*b);
                if (0 <= t && t <= 1){
                    _hit_s_bezier(seg, t, 1, n_hit, s);
                }
            }
        }
    } else {
        // Full cubic equation. Coefficients for the depressed cubic t^3 + p*t + q = 0:
        double p = (3*a*c - b*b)/(3*a*a);
        double q = (2*b*b*b - 9*a*b*c + 27*a*a*d)/(27*a*a*a);
        double disc = -p*p*p/27 - q*q/4;  // This is the discriminant of the depressed cubic but divided by (4*27)
        if (fabs(disc) < XC_EPSILON){
            if (fabs(p) < XC_EPSILON){
                // One real root with multiplicity 3
                t = -b/(3*a);
                if (0 <= t && t <= 1){
                    _hit_s_bezier(seg, t, 3, n_hit, s);
                }
            } else {
                // Two real roots (one simple and one with multiplicity 2)
                t = 3*q/p - b/(3*a);
                if (0 <= t && t <= 1){
                    _hit_s_bezier(seg, t, 1, n_hit, s);
                }
                t = -3*q/(2*p) - b/(3*a);
                if (0 <= t && t <= 1){
                    _hit_s_bezier(seg, t, 2, n_hit, s);
                }
            }
        } else if (disc < 0){
            // One real root
            t = cbrt(-q/2 + sqrt(fabs(disc))) + cbrt(-q/2 - sqrt(fabs(disc))) - b/(3*a);
            if (0 <= t && t <= 1){
                _hit_s_bezier(seg, t, 1, n_hit, s);
            }
        } else {
            // Three real roots
            double phi = acos(3*q/(2*p)*sqrt(fabs(3/p)));
            t = 2*sqrt(fabs(p/3))*cos(phi/3) - b/(3*a);
            if (0 <= t && t <= 1){
                _hit_s_bezier(seg, t, 1, n_hit, s);
            }
            t = 2*sqrt(fabs(p/3))*cos((phi + 2*M_PI)/3) - b/(3*a);
            if (0 <= t && t <= 1){
                _hit_s_bezier(seg, t, 1, n_hit, s);
            }
            t = 2*sqrt(fabs(p/3))*cos((phi + 4*M_PI)/3) - b/(3*a);
            if (0 <= t && t <= 1){
                _hit_s_bezier(seg, t, 1, n_hit, s);
            }
        }
    }
}

#endif /* XCOLL_COLL_GEOM_BEZIERSEG_H */
#ifndef XOBJ_TYPEDEF_LocalSegment
#define XOBJ_TYPEDEF_LocalSegment
typedef   struct LocalSegment_s * LocalSegment;
enum LocalSegment_e{LocalSegment_LineSegment_t,LocalSegment_HalfOpenLineSegment_t,LocalSegment_CircularSegment_t,LocalSegment_BezierSegment_t};
 static inline LocalSegment LocalSegment_getp(LocalSegment restrict  obj){
  int64_t offset=0;
  return (LocalSegment)(( char*) obj+offset);
}
 static inline int64_t LocalSegment_typeid(const LocalSegment restrict  obj){
  int64_t offset=0;
  offset+=8;
  return *( int64_t*)(( char*) obj+offset);
}
 static inline  void* LocalSegment_member(const LocalSegment restrict  obj){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset);
 return ( void*)(( char*) obj+offset);
}
 static inline void LocalSegment_crossing_drift(const LocalSegment restrict  obj,  int8_t* restrict  n_hit,  double* restrict  s, double s0, double x0, double xm){
   void* member = LocalSegment_member(obj);
  switch (LocalSegment_typeid(obj)){
        #ifndef LOCALSEGMENT_SKIP_LINESEGMENT
        case LocalSegment_LineSegment_t:
            return LineSegment_crossing_drift((LineSegment) member,n_hit,s,s0,x0,xm);
            break;
        #endif
        #ifndef LOCALSEGMENT_SKIP_HALFOPENLINESEGMENT
        case LocalSegment_HalfOpenLineSegment_t:
            return HalfOpenLineSegment_crossing_drift((HalfOpenLineSegment) member,n_hit,s,s0,x0,xm);
            break;
        #endif
        #ifndef LOCALSEGMENT_SKIP_CIRCULARSEGMENT
        case LocalSegment_CircularSegment_t:
            return CircularSegment_crossing_drift((CircularSegment) member,n_hit,s,s0,x0,xm);
            break;
        #endif
        #ifndef LOCALSEGMENT_SKIP_BEZIERSEGMENT
        case LocalSegment_BezierSegment_t:
            return BezierSegment_crossing_drift((BezierSegment) member,n_hit,s,s0,x0,xm);
            break;
        #endif
  }
  return;
}
#endif
#ifndef XOBJ_TYPEDEF_ArrNLocalSegment
#define XOBJ_TYPEDEF_ArrNLocalSegment
typedef   struct ArrNLocalSegment_s * ArrNLocalSegment;
 static inline ArrNLocalSegment ArrNLocalSegment_getp(ArrNLocalSegment restrict  obj){
  int64_t offset=0;
  return (ArrNLocalSegment)(( char*) obj+offset);
}
 static inline int64_t ArrNLocalSegment_len(ArrNLocalSegment restrict  obj){
  int64_t offset=0;
   int64_t* arr = ( int64_t*)(( char*) obj+offset);
  return arr[1];
}
 static inline LocalSegment ArrNLocalSegment_getp1(ArrNLocalSegment restrict  obj, int64_t i0){
  int64_t offset=0;
  offset+=16+i0*16;
  return (LocalSegment)(( char*) obj+offset);
}
 static inline int64_t ArrNLocalSegment_typeid(const ArrNLocalSegment restrict  obj, int64_t i0){
  int64_t offset=0;
  offset+=16+i0*16;
  offset+=8;
  return *( int64_t*)(( char*) obj+offset);
}
 static inline  void* ArrNLocalSegment_member(const ArrNLocalSegment restrict  obj, int64_t i0){
  int64_t offset=0;
  offset+=16+i0*16;
  offset+=*( int64_t*)(( char*) obj+offset);
 return ( void*)(( char*) obj+offset);
}
#endif
#ifndef XOBJ_TYPEDEF_Shape2D
#define XOBJ_TYPEDEF_Shape2D
typedef   struct Shape2D_s * Shape2D;
 static inline Shape2D Shape2D_getp(Shape2D restrict  obj){
  int64_t offset=0;
  return (Shape2D)(( char*) obj+offset);
}
 static inline ArrNLocalSegment Shape2D_getp_segments(Shape2D restrict  obj){
  int64_t offset=0;
  offset+=16;
  return (ArrNLocalSegment)(( char*) obj+offset);
}
 static inline int64_t Shape2D_len_segments(Shape2D restrict  obj){
  int64_t offset=0;
  offset+=16;
   int64_t* arr = ( int64_t*)(( char*) obj+offset);
  return arr[1];
}
 static inline LocalSegment Shape2D_getp1_segments(Shape2D restrict  obj, int64_t i0){
  int64_t offset=0;
  offset+=16;
  offset+=16+i0*16;
  return (LocalSegment)(( char*) obj+offset);
}
 static inline int64_t Shape2D_typeid_segments(const Shape2D restrict  obj, int64_t i0){
  int64_t offset=0;
  offset+=16;
  offset+=16+i0*16;
  offset+=8;
  return *( int64_t*)(( char*) obj+offset);
}
 static inline  void* Shape2D_member_segments(const Shape2D restrict  obj, int64_t i0){
  int64_t offset=0;
  offset+=16;
  offset+=16+i0*16;
  offset+=*( int64_t*)(( char*) obj+offset);
 return ( void*)(( char*) obj+offset);
}
 static inline int64_t Shape2D_get__seg_id(const Shape2D restrict  obj){
  int64_t offset=0;
  offset+=8;
  return *( int64_t*)(( char*) obj+offset);
}
 static inline void Shape2D_set__seg_id(Shape2D restrict  obj, int64_t value){
  int64_t offset=0;
  offset+=8;
  *( int64_t*)(( char*) obj+offset)=value;
}
 static inline  int64_t* Shape2D_getp__seg_id(Shape2D restrict  obj){
  int64_t offset=0;
  offset+=8;
  return ( int64_t*)(( char*) obj+offset);
}
#endif

 static inline
void Shape2D_crossing_drift(Shape2D shape, int8_t* n_hit, double* s, double s0, double x0, double xm){
    int64_t n_segments = Shape2D_len_segments(shape);
    for (int8_t i=0; i<n_segments;i++) {
        LocalSegment seg = Shape2D_getp1_segments(shape, i);
        LocalSegment_crossing_drift(seg, n_hit, s, s0, x0, xm);
    }
    sort_array_of_double(s, (int64_t) *n_hit);
}

 static inline
double Shape2D_crossing_drift_first(Shape2D shape, double s0, double x0, double xm){
    int8_t n_hit = 0;
    int64_t seg_id = Shape2D_get__seg_id(shape);
    switch (seg_id){
        /*START_SEG_ID_CASES_Shape2D_DriftTrajectory_first*/
        case 0: {  // SIZE: 8
            double s[8];
            Shape2D_crossing_drift(shape, &n_hit, s, s0, x0, xm);
            if (n_hit>0){
                return s[0];
            }
            return XC_S_MAX;
        }
        case 1: {  // SIZE: 15
            double s[15];
            Shape2D_crossing_drift(shape, &n_hit, s, s0, x0, xm);
            if (n_hit>0){
                return s[0];
            }
            return XC_S_MAX;
        }
        case 2: {  // SIZE: 2
            double s[2];
            Shape2D_crossing_drift(shape, &n_hit, s, s0, x0, xm);
            if (n_hit>0){
                return s[0];
            }
            return XC_S_MAX;
        }
        /*END_SEG_ID_CASES_Shape2D_DriftTrajectory_first*/
        default:
            printf("Unknown seg_id for Shape2D: %ld\nPlease recompile.", seg_id); fflush(stdout);  //only_for_context cpu_serial
            (void) n_hit;  // Avoid warning
            return XC_S_MAX;
    }
}

 static inline
double Shape2D_crossing_drift_before_s(Shape2D shape, double s0, double x0, double xm, double before_s){
    int8_t n_hit = 0;
    int64_t seg_id = Shape2D_get__seg_id(shape);
    switch (seg_id){
        /*START_SEG_ID_CASES_Shape2D_DriftTrajectory_before_s*/
        case 0: {  // SIZE: 8
            double s[8];
            Shape2D_crossing_drift(shape, &n_hit, s, s0, x0, xm);
            for (int8_t i=n_hit-1; i>=0; i--){
                if (s[i] <= before_s){
                    return s[i];
                }
            }
            return XC_S_MAX;
        }
        case 1: {  // SIZE: 15
            double s[15];
            Shape2D_crossing_drift(shape, &n_hit, s, s0, x0, xm);
            for (int8_t i=n_hit-1; i>=0; i--){
                if (s[i] <= before_s){
                    return s[i];
                }
            }
            return XC_S_MAX;
        }
        case 2: {  // SIZE: 2
            double s[2];
            Shape2D_crossing_drift(shape, &n_hit, s, s0, x0, xm);
            for (int8_t i=n_hit-1; i>=0; i--){
                if (s[i] <= before_s){
                    return s[i];
                }
            }
            return XC_S_MAX;
        }
        /*END_SEG_ID_CASES_Shape2D_DriftTrajectory_before_s*/
        default:
            printf("Unknown seg_id for Shape2D: %ld\nPlease recompile.", seg_id); fflush(stdout);  //only_for_context cpu_serial
            (void) n_hit;  // Avoid warning
            return XC_S_MAX;
    }
}

 static inline
double Shape2D_crossing_drift_after_s(Shape2D shape, double s0, double x0, double xm, double after_s){
    int8_t n_hit = 0;
    int64_t seg_id = Shape2D_get__seg_id(shape);
    switch (seg_id){
        /*START_SEG_ID_CASES_Shape2D_DriftTrajectory_after_s*/
        case 0: {  // SIZE: 8
            double s[8];
            Shape2D_crossing_drift(shape, &n_hit, s, s0, x0, xm);
            for (int8_t i=0; i<n_hit; i++){
                if (s[i] >= after_s){
                    return s[i];
                }
            }
            return XC_S_MAX;
        }
        case 1: {  // SIZE: 15
            double s[15];
            Shape2D_crossing_drift(shape, &n_hit, s, s0, x0, xm);
            for (int8_t i=0; i<n_hit; i++){
                if (s[i] >= after_s){
                    return s[i];
                }
            }
            return XC_S_MAX;
        }
        case 2: {  // SIZE: 2
            double s[2];
            Shape2D_crossing_drift(shape, &n_hit, s, s0, x0, xm);
            for (int8_t i=0; i<n_hit; i++){
                if (s[i] >= after_s){
                    return s[i];
                }
            }
            return XC_S_MAX;
        }
        /*END_SEG_ID_CASES_Shape2D_DriftTrajectory_after_s*/
        default:
            printf("Unknown seg_id for Shape2D: %ld\nPlease recompile.", seg_id); fflush(stdout);  //only_for_context cpu_serial
            (void) n_hit;  // Avoid warning
            return XC_S_MAX;
    }
}

 static inline
double Shape2D_crossing_drift_last(Shape2D shape, double s0, double x0, double xm){
    int8_t n_hit = 0;
    int64_t seg_id = Shape2D_get__seg_id(shape);
    switch (seg_id){
        /*START_SEG_ID_CASES_Shape2D_DriftTrajectory_last*/
        case 0: {  // SIZE: 8
            double s[8];
            Shape2D_crossing_drift(shape, &n_hit, s, s0, x0, xm);
            if (n_hit>0){
                return s[n_hit-1];
            }
            return XC_S_MAX;
        }
        case 1: {  // SIZE: 15
            double s[15];
            Shape2D_crossing_drift(shape, &n_hit, s, s0, x0, xm);
            if (n_hit>0){
                return s[n_hit-1];
            }
            return XC_S_MAX;
        }
        case 2: {  // SIZE: 2
            double s[2];
            Shape2D_crossing_drift(shape, &n_hit, s, s0, x0, xm);
            if (n_hit>0){
                return s[n_hit-1];
            }
            return XC_S_MAX;
        }
        /*END_SEG_ID_CASES_Shape2D_DriftTrajectory_last*/
        default:
            printf("Unknown seg_id for Shape2D: %ld\nPlease recompile.", seg_id); fflush(stdout);  //only_for_context cpu_serial
            (void) n_hit;  // Avoid warning
            return XC_S_MAX;
    }
}
#ifndef XOBJ_TYPEDEF_Shape2DV
#define XOBJ_TYPEDEF_Shape2DV
typedef   struct Shape2DV_s * Shape2DV;
 static inline Shape2DV Shape2DV_getp(Shape2DV restrict  obj){
  int64_t offset=0;
  return (Shape2DV)(( char*) obj+offset);
}
 static inline ArrNLocalSegment Shape2DV_getp_segments(Shape2DV restrict  obj){
  int64_t offset=0;
  offset+=32;
  return (ArrNLocalSegment)(( char*) obj+offset);
}
 static inline int64_t Shape2DV_len_segments(Shape2DV restrict  obj){
  int64_t offset=0;
  offset+=32;
   int64_t* arr = ( int64_t*)(( char*) obj+offset);
  return arr[1];
}
 static inline LocalSegment Shape2DV_getp1_segments(Shape2DV restrict  obj, int64_t i0){
  int64_t offset=0;
  offset+=32;
  offset+=16+i0*16;
  return (LocalSegment)(( char*) obj+offset);
}
 static inline int64_t Shape2DV_typeid_segments(const Shape2DV restrict  obj, int64_t i0){
  int64_t offset=0;
  offset+=32;
  offset+=16+i0*16;
  offset+=8;
  return *( int64_t*)(( char*) obj+offset);
}
 static inline  void* Shape2DV_member_segments(const Shape2DV restrict  obj, int64_t i0){
  int64_t offset=0;
  offset+=32;
  offset+=16+i0*16;
  offset+=*( int64_t*)(( char*) obj+offset);
 return ( void*)(( char*) obj+offset);
}
 static inline Arr2Float64 Shape2DV_getp_vlimit(Shape2DV restrict  obj){
  int64_t offset=0;
  offset+=8;
  return (Arr2Float64)(( char*) obj+offset);
}
 static inline int64_t Shape2DV_len_vlimit(Shape2DV restrict  obj){
  return 2;
}
 static inline double Shape2DV_get_vlimit(const Shape2DV restrict  obj, int64_t i0){
  int64_t offset=0;
  offset+=8;
  offset+=i0*8;
  return *( double*)(( char*) obj+offset);
}
 static inline void Shape2DV_set_vlimit(Shape2DV restrict  obj, int64_t i0, double value){
  int64_t offset=0;
  offset+=8;
  offset+=i0*8;
  *( double*)(( char*) obj+offset)=value;
}
 static inline  double* Shape2DV_getp1_vlimit(Shape2DV restrict  obj, int64_t i0){
  int64_t offset=0;
  offset+=8;
  offset+=i0*8;
  return ( double*)(( char*) obj+offset);
}
 static inline int64_t Shape2DV_get__seg_id(const Shape2DV restrict  obj){
  int64_t offset=0;
  offset+=24;
  return *( int64_t*)(( char*) obj+offset);
}
 static inline void Shape2DV_set__seg_id(Shape2DV restrict  obj, int64_t value){
  int64_t offset=0;
  offset+=24;
  *( int64_t*)(( char*) obj+offset)=value;
}
 static inline  int64_t* Shape2DV_getp__seg_id(Shape2DV restrict  obj){
  int64_t offset=0;
  offset+=24;
  return ( int64_t*)(( char*) obj+offset);
}
#endif

 static inline
void Shape2DV_crossing_drift(Shape2DV shape, int8_t* n_hit, double* s, double s0, double x0, double xm, double y0, double ym){
    double restrict_s[2];
    double ymin = Shape2DV_get_vlimit(shape, 0);
    double ymax = Shape2DV_get_vlimit(shape, 1);
    int8_t v_result = DriftTrajectory_vlimit(restrict_s, s0, y0, ym, ymin, ymax);
    if (v_result == 0){
        return;  // Completely outside - no crossing possible
    } else {
        int64_t n_segments = Shape2DV_len_segments(shape);
        for (int8_t i=0; i<n_segments;i++) {
            LocalSegment seg = Shape2DV_getp1_segments(shape, i);
            LocalSegment_crossing_drift(seg, n_hit, s, s0, x0, xm);
        }
        sort_array_of_double(s, (int64_t) *n_hit);
        if (v_result == 1){
            calculate_overlap_array_interval(s, n_hit, restrict_s);
        }
    }
}

 static inline
double Shape2DV_crossing_drift_first(Shape2DV shape, double s0, double x0, double xm, double y0, double ym){
    int8_t n_hit = 0;
    int64_t seg_id = Shape2DV_get__seg_id(shape);
    switch (seg_id){
        /*START_SEG_ID_CASES_Shape2DV_DriftTrajectory_first*/
        case 0: {  // SIZE: 15
            double s[15];
            Shape2DV_crossing_drift(shape, &n_hit, s, s0, x0, xm, y0, ym);
            if (n_hit>0){
                return s[0];
            }
            return XC_S_MAX;
        }
        /*END_SEG_ID_CASES_Shape2DV_DriftTrajectory_first*/
        default:
            printf("Unknown seg_id for Shape2DV: %ld\nPlease recompile.", seg_id); fflush(stdout);  //only_for_context cpu_serial
            (void) n_hit;  // Avoid warning
            return XC_S_MAX;
    }
}

 static inline
double Shape2DV_crossing_drift_before_s(Shape2DV shape, double s0, double x0, double xm, double y0, double ym, double before_s){
    int8_t n_hit = 0;
    int64_t seg_id = Shape2DV_get__seg_id(shape);
    switch (seg_id){
        /*START_SEG_ID_CASES_Shape2DV_DriftTrajectory_before_s*/
        case 0: {  // SIZE: 15
            double s[15];
            Shape2DV_crossing_drift(shape, &n_hit, s, s0, x0, xm, y0, ym);
            for (int8_t i=n_hit-1; i>=0; i--){
                if (s[i] <= before_s){
                    return s[i];
                }
            }
            return XC_S_MAX;
        }
        /*END_SEG_ID_CASES_Shape2DV_DriftTrajectory_before_s*/
        default:
            printf("Unknown seg_id for Shape2DV: %ld\nPlease recompile.", seg_id); fflush(stdout);  //only_for_context cpu_serial
            (void) n_hit;  // Avoid warning
            return XC_S_MAX;
    }
}

 static inline
double Shape2DV_crossing_drift_after_s(Shape2DV shape, double s0, double x0, double xm, double y0, double ym, double after_s){
    int8_t n_hit = 0;
    int64_t seg_id = Shape2DV_get__seg_id(shape);
    switch (seg_id){
        /*START_SEG_ID_CASES_Shape2DV_DriftTrajectory_after_s*/
        case 0: {  // SIZE: 15
            double s[15];
            Shape2DV_crossing_drift(shape, &n_hit, s, s0, x0, xm, y0, ym);
            for (int8_t i=0; i<n_hit; i++){
                if (s[i] >= after_s){
                    return s[i];
                }
            }
            return XC_S_MAX;
        }
        /*END_SEG_ID_CASES_Shape2DV_DriftTrajectory_after_s*/
        default:
            printf("Unknown seg_id for Shape2DV: %ld\nPlease recompile.", seg_id); fflush(stdout);  //only_for_context cpu_serial
            (void) n_hit;  // Avoid warning
            return XC_S_MAX;
    }
}

 static inline
double Shape2DV_crossing_drift_last(Shape2DV shape, double s0, double x0, double xm, double y0, double ym){
    int8_t n_hit = 0;
    int64_t seg_id = Shape2DV_get__seg_id(shape);
    switch (seg_id){
        /*START_SEG_ID_CASES_Shape2DV_DriftTrajectory_last*/
        case 0: {  // SIZE: 15
            double s[15];
            Shape2DV_crossing_drift(shape, &n_hit, s, s0, x0, xm, y0, ym);
            if (n_hit>0){
                return s[n_hit-1];
            }
            return XC_S_MAX;
        }
        /*END_SEG_ID_CASES_Shape2DV_DriftTrajectory_last*/
        default:
            printf("Unknown seg_id for Shape2DV: %ld\nPlease recompile.", seg_id); fflush(stdout);  //only_for_context cpu_serial
            (void) n_hit;  // Avoid warning
            return XC_S_MAX;
    }
}