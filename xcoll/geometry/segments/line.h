// copyright ############################### #
// This file is part of the Xcoll package.   #
// Copyright (c) CERN, 2024.                 #
// ######################################### #

#ifndef XCOLL_COLL_GEOM_LINESEG_H
#define XCOLL_COLL_GEOM_LINESEG_H
#define XC_LINE_CROSSINGS 2
// #include "../trajectories/mcs.h"
// #include "../c_init/simpson.h" // this is just for me to avoid squiggles


/*gpufun*/
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

/*gpufun*/
double MultipleCoulomb_Line(double s, McsLineParams params){
    // MCS trajectory form PDG rewritted in terms of A, B and s/Xo. 
    // Params_Line* p = (Params_Line*)params;
    const double Ax = McsLineParams_get_Ax(params);
    const double Xo = McsLineParams_get_Xo(params);
    double x2 = McsLineParams_get_x2(params);
    double s2 = McsLineParams_get_s2(params);
    double x1 = McsLineParams_get_x1(params);
    double s1 = McsLineParams_get_s1(params);

    double mcs  = Ax * pow(sqrt(s/Xo),3.0) * (1.0/0.038 + log(s/Xo));
    double line = x1 + (x2 - x1) / (s2 - s1) * (s - s1);
    return mcs - line;
}

/*gpufun*/
double MultipleCoulombDeriv_Line(double s, McsLineParams params){
    // MCS trajectory derivative wrt s
    const double Ax = McsLineParams_get_Ax(params);
    const double Xo = McsLineParams_get_Xo(params);
    double x2 = McsLineParams_get_x2(params);
    double s2 = McsLineParams_get_s2(params);
    double x1 = McsLineParams_get_x1(params);
    double s1 = McsLineParams_get_s1(params);

    double mcs_deriv  = Ax/Xo * (sqrt(s/Xo)*3.0/2.0*log(s/Xo)+1.0/0.038 + sqrt(s/Xo));
    double line_deriv = (x2 - x1) / (s2 - s1);
    return mcs_deriv - line_deriv;
}
// TODO: FIX
/*gpufun*/ 
void LineSegment_crossing_mcs(LineSegment seg, int8_t* n_hit, double* s, const double* Ax, const double Xo, void* params){
//     // Get segment data
//     double s1 = LineSegment_get_s1(seg);
//     double x1 = LineSegment_get_x1(seg);
//     double s2 = LineSegment_get_s2(seg);
//     double x2 = LineSegment_get_x2(seg);
    
//     // Define roots array and parameters
//     Params_Line params = {Xo, *Ax, x2, s2, x1, s1};
//     double roots[XC_LINE_CROSSINGS];
//     int8_t number_of_roots = 0;

//     grid_search_and_newton(MultipleCoulomb_Line, MultipleCoulombDeriv_Line, s1, s2, roots, XC_LINE_CROSSINGS, &params, &number_of_roots);
//     // now we have the roots, but that doesnt mean anything yet

//     // should this function be about finding the crossings, and THEN checking in jaw if nucl or not, and then updating nhit?
//     // Process the roots
//     for (int i = 0; i < XC_LINE_CROSSINGS; ++i) {
//         if (roots[i] >= s1 && roots[i] <= s2 && i < number_of_roots){ 
//             s[*n_hit] = roots[i];
//             (*n_hit)++;
//         }
//     }
}
#endif /* XCOLL_COLL_GEOM_LINESEG_H */