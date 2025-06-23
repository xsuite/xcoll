// copyright ############################### #
// This file is part of the Xcoll package.   #
// Copyright (c) CERN, 2023.                 #
// ######################################### #

#ifndef XCOLL_EVEREST_CRYSTAL_PARAMETERS_H
#define XCOLL_EVEREST_CRYSTAL_PARAMETERS_H
#include <math.h>
#include <stdio.h>
#include <stdlib.h>


/*gpufun*/
void calculate_initial_angle(EverestData restrict everest, LocalParticle* part, CrystalGeometry restrict cg){
    double R = cg->bending_radius;
    double s = LocalParticle_get_s(part);
    double x = LocalParticle_get_x(part);
    double s_P = cg->s_P;
    double x_P = cg->x_P;
    double r   = sqrt((s-s_P)*(s-s_P) + (x-x_P)*(x-x_P));
    everest->r = r;
    everest->t_I = R/fabs(R)*asin( (s-s_P)/r); // Tangent angle of the channelling planes (not necessarily the same as xp)
}


/*gpufun*/
void calculate_opening_angle(EverestData restrict everest, LocalParticle* part, CrystalGeometry restrict cg){
    double t    = cg->bending_angle;
    double xd   = cg->width;
    double R    = cg->bending_radius;
    double sinp = sin(cg->miscut_angle);
    double cosp = cos(cg->miscut_angle);

    // Radius from starting point to miscut centre
    double r = everest->r;
    // Intersection with exit face
    double bb = R/2. * sin(2.*t) * (1 - cosp - tan(t)*sinp);
    double s_F = bb + sqrt(bb*bb - tan(t)*sin(2.*t)*(R*R*(1-cosp) - r*r/2.));
    double x_F = R - s_F/tan(t);

    // We could intersect with the upper (positive miscut) or lower bend (negative miscut) before the exit face
    // Distance between bending centre and miscut centre:
    double d = R*sqrt(2*(1-cosp));
    if (cg->miscut_angle > 0){
        // Check if intersection with upper bend R-xd is possible
        double Rb = R-xd;
        if (d > fabs(r - Rb)){
            double st_UB = (r*r - Rb*Rb)/(2*d);
            double xt_UB = -sqrt(r*r/2. + Rb*Rb/2. - st_UB*st_UB - d*d/4.);
            double s_UB  = (st_UB/d - 0.5)*R*sinp - xt_UB/d*R*(1-cosp);
            if (s_F > s_UB){
                // Upper bend encountered before exit face
                s_F = s_UB;
                x_F = (st_UB/d + 0.5)*R*(1-cosp) + xt_UB/d*R*sinp + R*cosp;
            }
        }
    } else {
        // Check if intersection with lower bend R is possible
        if (d > fabs(R - r)){
            double st_LB = (R*R - r*r)/(2*d);
            double xt_LB = -sqrt(r*r/2. + R*R/2. - st_LB*st_LB - d*d/4.);
            double s_LB  = -(st_LB/d + 0.5)*R*sinp + xt_LB/d*R*(1-cosp);
            if (s_F > s_LB){
                // Lower bend encountered before exit face
                s_F = s_LB;
                x_F = -(st_LB/d - 0.5)*R*(1-cosp) - xt_LB/d*R*sinp + R*cosp;
            }
        }
    }

    // Opening angle of channelling trajectory
    double s = LocalParticle_get_s(part);
    double x = LocalParticle_get_x(part);
    everest->t_P = acos(1 - ( (s_F-s)*(s_F-s) - (x_F-x)*(x_F-x) ) / (2*r*r) );
}


/*gpufun*/
double _critical_angle0(EverestCollData restrict coll, double pc){
    // Define typical angles/probabilities for orientation 110
    double eum = coll->eum;
    double eta = coll->eta;
    return sqrt(2.e-9*eta*eum/pc); // Critical angle (rad) for straight crystals    // pc is actually beta pc
}

/*gpufun*/
double _critical_radius(EverestCollData restrict coll, double pc){
    // Define typical angles/probabilities for orientation 110
    double eum = coll->eum;
    double ai  = coll->ai;
    double eta = coll->eta;
    return pc/(2.e-6*sqrt(eta)*eum)*ai;  // Critical curvature radius [m]   // pc is actually beta pc
}

/*gpufun*/
double _critical_angle(EverestCollData restrict coll, double t_c0, double Rc_over_R){
    double t_c = 0;
    if (Rc_over_R <= 1.) {
        // Otherwise no channelling possible
        t_c = t_c0*(1 - Rc_over_R); // Critical angle for curved crystal
        if (coll->orient == 2) {
            t_c *= 0.98;
        }
    }
    return t_c;
}

/*gpufun*/
void calculate_critical_angle(EverestData restrict everest, LocalParticle* part, CrystalGeometry restrict cg, double pc){
    // Define typical angles/probabilities for orientation 110
    everest->t_c0  = _critical_angle0(everest->coll, pc);
    double Rcrit = _critical_radius(everest->coll, pc);
    everest->Rc_over_R = Rcrit / fabs(cg->bending_radius);
    everest->t_c = _critical_angle(everest->coll, everest->t_c0, everest->Rc_over_R);
}


/*gpufun*/
void calculate_VI_parameters(EverestData restrict everest, LocalParticle* part, double pc){

    double ratio = everest->Rc_over_R;
    double t_c0  = everest->t_c0;

    // Correction by sasha drozdin/armen    // TODO: pc = energy
    // K=0.0007 is taken based on simulations using CATCH.f (V.Biryukov)
    // TODO: From the paper (arXiv:0808.1486, eq (18), (34)) it seems the power of E should be 0.25
    // Typo in Daniele's thesis in 3.25 (E^0.2 instead of 1/E^0.2)
    everest->Vcapt = 7.e-4*(1./ratio - 0.7)/pow(pc, 0.2);

    double Ang_rms, Ang_avr;
    double c1 = -3./2.;   // Fitting coefficient
    double c2 = 5./3.;  // Fitting coefficient
    double c3 =  1.7;   // Fitting coefficient
    if (ratio > 1.) {
        // no channelling possibile
        Ang_avr = c1*t_c0*5.e-2/ratio;         // Average angle reflection
        Ang_rms = c3*0.42*t_c0*sin(1.4/ratio); // RMS scattering
        everest->Vcapt = 0.;                   // Probability of VC is zero
    } else if (ratio > 1./3.) {
        // Strongly bent crystal
        Ang_avr = c1*t_c0*(0.1972/ratio - 0.1472);        // Average angle reflection
        Ang_rms = c3*0.42*t_c0*sin(0.4713/ratio + 0.85);  // RMS scattering
    } else {
        // Rcrit << R
        Ang_avr = c1*t_c0*(1. - c2*ratio);     // Average angle for VR
        Ang_rms = c3*t_c0*ratio;               // RMS scattering
    }
    if (everest->coll->orient == 2) {
        Ang_avr = Ang_avr*0.93;
        Ang_rms = Ang_rms*1.05;
    }

    everest->Ang_rms = Ang_rms;
    everest->Ang_avr = Ang_avr;
}

#endif /* XCOLL_EVEREST_CRYSTAL_PARAMETERS_H */