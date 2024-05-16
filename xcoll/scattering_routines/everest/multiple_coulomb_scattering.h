// copyright ############################### #
// This file is part of the Xcoll Package.   #
// Copyright (c) CERN, 2024.                 #
// ######################################### #

#ifndef XCOLL_EVEREST_MCS_H
#define XCOLL_EVEREST_MCS_H
#include <math.h>
#include <stdio.h>


/*gpufun*/
double iterat(double a, double b, double dh, double s) {
    double ds = s;
    while (1) {
        ds = ds*0.5;
        if (pow(s,3) < pow((a+b*s),2)) {
            s = s+ds;
        } else {
            s = s-ds;
        }
        if (ds < dh) {
            break;
        } else {
            continue;
        }
    }
    return s;
}


/*gpufun*/
double soln3(double a, double b, double dh, double smax) {
    double s;
    if (b == 0) {
        s = pow(a,0.6666666666666667);
        if (s > smax) {
            s = smax;
        }
        return s;
    }
    if (a == 0) {    
        if (b > 0) {
            s = pow(b,2);
        } else {
            s = 0;
        }
        if (s > smax) {
            s = smax;
        }
        return s;
    }
    if (b > 0) {
        if (pow(smax,3) <= pow((a + b*smax),2)) {
            s = smax;
            return s;
        } else {
            s = smax*0.5;
            s = iterat(a,b,dh,s);
        }
    } else {
        double c = (-1*a)/b;
        if (smax < c) {
            if ((pow(smax,3)) <= pow((a + b*smax),2)) {
                s = smax;
                return s;
            } else {
                s = smax*0.5;
                s = iterat(a,b,dh,s);
            }
        } else {
            s = c*0.5;
            s = iterat(a,b,dh,s);
        }
    }
    return s;
}


/*gpufun*/
double* scamcs(LocalParticle* part, double x0, double xp0, double s) {
    double* result = (double*)malloc(2 * sizeof(double));

    // Generate two Gaussian random numbers z1 and z2
    double r2 = 0;
    double v1 = 0;
    double v2 = 0;
    while (1) {
        v1 = 2*RandomUniform_generate(part) - 1;
        v2 = 2*RandomUniform_generate(part) - 1;
        r2 = pow(v1,2) + pow(v2,2);
        if(r2 < 1) {
            break;
        }
    }
    double a   = sqrt((-2*log(r2))/r2);
    double z1  = v1*a;
    double z2  = v2*a;

    // MCS scaling by length in units of radiation length
    double ss  = sqrt(s) * (1 + 0.038*log(s));

    result[0] = x0  + s*(xp0 + 0.5*ss*(z2 + z1*0.577350269));
    result[1] = xp0 + ss*z2;
    return result;
}


/*gpufun*/
void mcs(EverestData restrict everest, LocalParticle* part, double length, double p, int edge_check){
    InteractionRecordData record = everest->coll->record;
    RecordIndex record_index     = everest->coll->record_index;
    int8_t sc = everest->coll->record_scatterings;

    // First log particle at start of multiple coulomb scattering
    int64_t i_slot;
    if (sc) i_slot = InteractionRecordData_log(record, record_index, part, XC_MULTIPLE_COULOMB_SCATTERING);

    double const radl = everest->coll->radl;
    double s;
    double theta = 13.6e-3/p;
    double h   = 0.001;
    double dh  = 0.0001;
    double bn0 = 0.4330127019;
    double rlen0 = length/radl;
    double rlen  = rlen0;

    double x  = LocalParticle_get_x(part);
    double z  = LocalParticle_get_y(part);
#ifdef XCOLL_USE_EXACT
    double xp = LocalParticle_get_exact_xp(part);  // This is tangent
    double zp = LocalParticle_get_exact_yp(part);  // This is tangent
#else
    double xp = LocalParticle_get_xp(part);
    double zp = LocalParticle_get_yp(part);
#endif

    x  = (x/theta)/radl;
    xp = xp/theta;
    z  = (z/theta)/radl;
    zp = zp/theta;

    if (edge_check){
        while (1) {
            double ae = bn0 * x;
            double be = bn0 * xp;
            // #######################################
            // ae = np.array(ae, dtype=np.double64)
            // be = np.array(be, dtype=np.double64)
            // dh = np.array(dh, dtype=np.double64)
            // rlen = np.array(rlen, dtype=np.double64)
            // s = np.array(s, dtype=np.double64)
            // #######################################
            s = soln3(ae,be,dh,rlen);
            if (s < h) {
                s = h;
            }
            // TODO: should cap s whenever we are out (the two if cases below), because now the scamcs is applied over an s that might be (slightly) too long
            double* res = scamcs(part, x, xp, s);
            x  = res[0];
            xp = res[1];
            free(res);
            if (x <= 0) {
                s = rlen0 - rlen + s;
                break; // go to 20
            }
            if (s + dh >= rlen) {
                s = rlen0;
                break; // go to 20
            }
            // go to 10
            rlen = rlen - s;
        }

    } else {
        double* res = scamcs(part, x, xp, rlen0);
        x  = res[0];
        xp = res[1];
        free(res);
        s = rlen0;
    }

    double* res = scamcs(part, z, zp, s);
    z  = res[0];
    zp = res[1];
    free(res);

    LocalParticle_set_x(part, x*theta*radl);
    LocalParticle_set_y(part, z*theta*radl);
#ifdef XCOLL_USE_EXACT
    LocalParticle_set_exact_xp_yp(part, xp*theta, zp*theta);
#else
    LocalParticle_set_xp_yp(part, xp*theta, zp*theta);
#endif
    LocalParticle_add_to_s(part, s*radl);

    // Finally log particle at end of multiple coulomb scattering
    if (sc) InteractionRecordData_log_child(record, i_slot, part, length);
}

#endif /* XCOLL_EVEREST_MCS_H */