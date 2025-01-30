// copyright ############################### #
// This file is part of the Xcoll package.   #
// Copyright (c) CERN, 2024.                 #
// ######################################### #

#ifndef XCOLL_COLL_GEOM_HALFOPENLINESEG_H
#define XCOLL_COLL_GEOM_HALFOPENLINESEG_H
#define XC_HALFOPENLINE_CROSSINGS 2


/*gpufun*/
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

/*gpufun*/
double MultipleCoulomb_HalfOpenLine(double s, McsHalfOpenLineParams params){
    // MCS trajectory form PDG rewritted in terms of A, B and s/Xo. 
    const double Ax = McsHalfOpenLineParams_get_Ax(params);
    const double Xo = McsHalfOpenLineParams_get_Xo(params);
    double s2       = McsHalfOpenLineParams_get_s2(params);
    double x2       = McsHalfOpenLineParams_get_x2(params);
    double s1       = McsHalfOpenLineParams_get_s1(params);
    double x1       = McsHalfOpenLineParams_get_x1(params);

    double mcs  = Ax * pow(sqrt(s/Xo),3.0) * (1.0/0.038 + log(s/Xo));
    double half_open_line = x2 + (x2 - x1) / (s2 - s1) * (s - s1);
    return mcs - half_open_line;
}

/*gpufun*/
double MultipleCoulombDeriv_HalfOpenLine(double s,  McsHalfOpenLineParams params){
    // MCS trajectory derivative wrt s
    const double Ax = McsHalfOpenLineParams_get_Ax(params);
    const double Xo = McsHalfOpenLineParams_get_Xo(params);
    double s2       = McsHalfOpenLineParams_get_s2(params);
    double x2       = McsHalfOpenLineParams_get_x2(params);
    double s1       = McsHalfOpenLineParams_get_s1(params);
    double x1       = McsHalfOpenLineParams_get_x1(params);

    double mcs_deriv  = Ax/Xo * (sqrt(s/Xo)*3.0/2.0*log(s/Xo)+1.0/0.038 + sqrt(s/Xo));
    double half_open_line_deriv = (x2 - x1) / (s2 - s1);
    return mcs_deriv - half_open_line_deriv;
}

/*gpufun*/
void HalfOpenLineSegment_crossing_mcs(HalfOpenLineSegment seg, int8_t* n_hit, double* s, const double* Ax, const double Xo){
    // // Get segment data
    // double s1 = HalfOpenLineSegment_get_s(seg);
    // double x1 = HalfOpenLineSegment_get_x(seg);
    // double s2 = s1 + cos(HalfOpenLineSegment_get_t(seg));
    // double x2 = x1 + sin(HalfOpenLineSegment_get_t(seg));

    // // Define roots array and parameters
    // double roots[XC_HALFOPENLINE_CROSSINGS];
    // int number_of_roots = 0;
    // Params_HalfOpenLine params = {Xo, *Ax, x2, s2, x1, s1};

    // grid_search_and_newton(MultipleCoulomb_HalfOpenLine, MultipleCoulombDeriv_HalfOpenLine, s1, s2, roots, XC_HALFOPENLINE_CROSSINGS, &params, &number_of_roots);
    // // now we have the roots, but that doesnt mean anything yet
    
    // // should this function be about finding the crossings, and THEN checking in jaw if nucl or not, and then updating nhit?
    // // Process the roots
    // for (int i = 0; i < XC_HALFOPENLINE_CROSSINGS; ++i) {
    //     if (roots[i] >= s1 && roots[i] <= s2 && i < number_of_roots){
    //         s[*n_hit] = roots[i];
    //         (*n_hit)++;
    //     }
    // }
}  
#endif /* XCOLL_COLL_GEOM_HALFOPENLINESEG_H */