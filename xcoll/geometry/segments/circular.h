// copyright ############################### #
// This file is part of the Xcoll package.   #
// Copyright (c) CERN, 2024.                 #
// ######################################### #

#ifndef XCOLL_COLL_GEOM_CIRCULARSEG_H
#define XCOLL_COLL_GEOM_CIRCULARSEG_H
#define XC_CIRCULARSEG_CROSSINGS 2

// im i allowed to do this 
typedef struct {
    double Xo;
    double Ax;
    double R;
    double sC;
    double xC;
    double x;
} Params_Circular;

/*gpufun*/
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


/*gpufun*/
double MultipleCoulomb_Circular(double s, McsCircularParams params){//CircularSegment sehh, double x, const double Xo, const double Ax){
    // MCS trajectory form PDG rewritted in terms of A, B and s/Xo and with circular equation. Note: s is particle s.
    double R        = McsCircularParams_get_R(params);
    double sC       = McsCircularParams_get_sC(params);
    double xC       = McsCircularParams_get_xC(params);
    const double Ax = McsCircularParams_get_Ax(params);
    const double Xo = McsCircularParams_get_Xo(params);
    double x        = McsCircularParams_get_x(params);

    double mcs  = Ax * pow(sqrt(s/Xo),3.0) * (1.0/0.038 + log(s/Xo));
    double circle = pow(x - xC, 2.0) + pow(s - sC, 2.0) - pow(R, 2.0);
    return mcs - circle;
}

/*gpufun*/
double MultipleCoulombDeriv_Circular(double s, McsCircularParams params){
    // MCS trajectory derivative wrt s. Note: s is particle s.
    double R        = McsCircularParams_get_R(params);
    double sC       = McsCircularParams_get_sC(params);
    const double Ax = McsCircularParams_get_Ax(params);
    const double Xo = McsCircularParams_get_Xo(params);

    double mcs_deriv  = Ax/Xo * (sqrt(s/Xo)*3.0/2.0*log(s/Xo)+1.0/0.038 + sqrt(s/Xo));
    double circle_deriv = 2*(s - sC);
    return mcs_deriv - circle_deriv;
}

/*gpufun*/
void CircularSegment_crossing_mcs(CircularSegment seg, int8_t* n_hit, double* s, double x, const double* Ax, const double Xo){
    // Get segment data
    // double R  = CircularSegment_get_R(seg);
    // double sC = CircularSegment_get_s(seg);
    // double xC = CircularSegment_get_x(seg);
    // double t1 = CircularSegment_get_t1(seg);
    // double t2 = CircularSegment_get_t2(seg);
    // double s_min = sC - R;
    // double s_max = sC + R;

    // // define parameters 
    // Params_Circular params = {Xo, *Ax, R, sC, xC, x};
    // int number_of_roots = 0;
    // double roots[XC_CIRCULARSEG_CROSSINGS];

    // grid_search_and_newton(MultipleCoulomb_Circular, MultipleCoulombDeriv_Circular, s_min, s_max, roots, XC_CIRCULARSEG_CROSSINGS, &params, &number_of_roots);
    // // now we have the roots, but that doesnt mean anything yet

    // // should this function be about finding the crossings, and THEN checking in jaw if nucl or not, and then updating nhit?
    // // Process the roots
    // for (int i = 0; i < XC_CIRCULARSEG_CROSSINGS; ++i) {
    //     if (roots[i] >= sC - R && roots[i] <= sC + R && i < number_of_roots) {
    //         s[*n_hit] = roots[i];
    //         (*n_hit)++;                                  // is this correct..? 
    //     }
    // }
}
#endif /* XCOLL_COLL_GEOM_CIRCULARSEG_H */