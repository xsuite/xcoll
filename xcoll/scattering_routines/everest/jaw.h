// copyright ############################### #
// This file is part of the Xcoll Package.   #
// Copyright (c) CERN, 2023.                 #
// ######################################### #

#ifndef XCOLL_EVEREST_JAW_H
#define XCOLL_EVEREST_JAW_H
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
    // double x0  = *xx;
    // double xp0 = *xxp;
    double r2 = 0;
    double v1 = 0;
    double v2 = 0;
    double* result = (double*)malloc(2 * sizeof(double));

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
    double ss  = sqrt(s);
    double sss = 1 + 0.038*log(s);

    result[0] = x0  + s*(xp0 + ((0.5*ss)*sss)*(z2 + z1*0.577350269));
    result[1] = xp0 + (ss*z2)*sss;
    return result;
}


/*gpufun*/
double* mcs(LocalParticle* part, MaterialData material, double zlm1, double p, double x, double xp, double z, double zp) {

    double const radl = MaterialData_get_radiation_length(material);
    double s;
    double theta = 13.6e-3/p;
    double h   = 0.001;
    double dh  = 0.0001;
    double bn0 = 0.4330127019;
    double rlen0 = zlm1/radl;
    double rlen  = rlen0;
    double* result = (double*)malloc(5 * sizeof(double));

    x  = (x/theta)/radl;
    xp = xp/theta;
    z  = (z/theta)/radl;
    zp = zp/theta;

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
            s = (rlen0-rlen)+ s;
            break; // go to 20
        }
        if ((s + dh) >= rlen) {
            s = rlen0;
            break; // go to 20
        }
        // go to 10
        rlen = rlen - s;
    }

    double* res = scamcs(part, z,zp,s);
    z  = res[0];
    zp = res[1];
    free(res);

    result[0] = s*radl;
    result[1] = (x*theta)*radl;
    result[2] = xp*theta;
    result[3] = (z*theta)*radl;
    result[4] = zp*theta;
    return result;
}


/*gpufun*/
double* tetat(LocalParticle* part, double t, double p) {
    double teta = sqrt(t)/p;
    double va = 0;
    double vb  = 0;
    double va2 = 0;
    double vb2 = 0;
    double r2  = 0;
    double* result = (double*)malloc(2 * sizeof(double));
    
    while (1) {
        va  = 2*RandomUniform_generate(part) - 1;
        vb  = RandomUniform_generate(part);
        va2 = pow(va,2);
        vb2 = pow(vb,2);
        r2  = va2 + vb2;
        if(r2 < 1) {
            break;
        }
    }
        
    result[0] = (teta*((2*va)*vb))/r2;
    result[1]  = (teta*(va2 - vb2))/r2;
    return result;
}


/*gpufun*/
double* gettran(RandomRutherfordData rng, LocalParticle* part, double inter, double p, struct ScatteringParameters scat) {

    double* res = (double*)malloc(2 * sizeof(double));

    // Neither if-statements below have an else, so defaulting function return to zero.
    double result = 0;

    if (inter==2) { // Nuclear Elastic
        result = RandomExponential_generate(part)/scat.bn;
    } else if (inter==3) { // pp Elastic
        result = RandomExponential_generate(part)/scat.bpp;
    } else if (inter==4) { // Single Diffractive
        double xm2 = exp(RandomUniform_generate(part) * scat.xln15s);
        double bsd = 0;
        p = p * (1 - xm2/scat.ecmsq);
    
        if (xm2 < 2) {
            bsd = 2 * scat.bpp;
        } else if ((xm2 >= 2) & (xm2 <= 5)) {
            bsd = ((106.0 - 17.0*xm2)*scat.bpp)/36.0;
        } else {
            bsd = (7*scat.bpp)/12.0;
        }
        result = RandomExponential_generate(part)/bsd;
    } else if (inter==5) { // Coulomb
        result = RandomRutherford_generate(rng, part);
    }

    res[0] = result;
    res[1] = p;
    return res;
}


/*gpufun*/
int ichoix(LocalParticle* part, struct ScatteringParameters scat) {

    double aran = RandomUniform_generate(part);
    int i;
    for (i = 0; i < 5; ++i) {
        if (aran <= scat.cprob[i]) {
            break;
        }
    }
    return i;
}


/*gpufun*/
double* jaw(LocalParticle* part, MaterialData material, RandomRutherfordData rng, struct ScatteringParameters scat,
            double p, double zlm, CollimatorImpactsData record, RecordIndex record_index) {

    double* result = (double*)malloc(3 * sizeof(double));
    double s;
    double nabs = 0;
    double rlen = zlm;
    double m_dpodx = 0.;
    double t;
    double tx; 
    double tz;

    double rpp_in  = LocalParticle_get_rpp(part);
    double x  = LocalParticle_get_x(part);
    double xp = LocalParticle_get_px(part)*rpp_in;
    double z  = LocalParticle_get_y(part);
    double zp = LocalParticle_get_py(part)*rpp_in;

    // Do a step for a point-like interaction.
    // Get monte-carlo interaction length.
    while (1) {

        struct IonisationProperties properties = calculate_ionisation_properties(p, (GeneralMaterialData) material);
        double zlm1 = scat.xintl*RandomExponential_generate(part);
                        
        // If the monte-carlo interaction length is longer than the
        // remaining collimator length, then put it to the remaining
        // length, do multiple coulomb scattering and return.
        // LAST STEP IN ITERATION LOOP
        if (zlm1 > rlen) { 
            zlm1 = rlen;
            double* res = mcs(part, material, zlm1, p, x, xp, z, zp);
            s = res[0];
            x = res[1];
            xp = res[2];
            z = res[3];
            zp = res[4];
            free(res);

            s = zlm - rlen + s;
            m_dpodx = calcionloss(part, rlen, properties);  // DM routine to include tail // TODO: should not be rlen but s after updating
            p = p-m_dpodx*s; // This is correct: ionisation loss is only calculated and applied at end of while (break)
            break;
        }
        // Otherwise do multi-coulomb scattering.
        // REGULAR STEP IN ITERATION LOOP
        double* res1 = mcs(part, material, zlm1, p, x, xp, z, zp);
        s = res1[0];
        x = res1[1];
        xp = res1[2];
        z = res1[3];
        zp = res1[4];
        free(res1);

        // Check if particle is outside of collimator (X.LT.0) after
        // MCS. If yes, calculate output longitudinal position (s),
        // reduce momentum (output as dpop) and return.
        // PARTICLE LEFT COLLIMATOR BEFORE ITS END.

        if(x <= 0) {
            s = zlm - rlen + s;
            m_dpodx = calcionloss(part, rlen, properties);  // TODO: should not be rlen but s after updating
            p = p-m_dpodx*s;  // correct
            break;
        }

        // Check whether particle is absorbed. If yes, calculate output
        // longitudinal position (s), reduce momentum (output as dpop)
        // and return.
        // PARTICLE WAS ABSORBED INSIDE COLLIMATOR DURING MCS.

        int inter = ichoix(part, scat);
        nabs = inter;
        if (inter == 1) {
            s = (zlm-rlen)+zlm1;
            m_dpodx = calcionloss(part, rlen, properties);
            p = p-m_dpodx*s;
            break;
        }


        // Now treat the other types of interaction, as determined by ICHOIX:

        // Nuclear-Elastic:          inter = 2
        // pp Elastic:               inter = 3
        // Single-Diffractive:       inter = 4    (changes momentum p)
        // Coulomb:                  inter = 5

        // Gettran returns some monte carlo number, that, as I believe, gives the rms transverse momentum transfer.

        double* res2 = gettran(rng, part, inter, p, scat);
        t = res2[0];
        p = res2[1];
        free(res2);

        // Tetat calculates from the rms transverse momentum transfer in
        // monte-carlo fashion the angle changes for x and z planes. The
        // angle change is proportional to SQRT(t) and 1/p, as expected.

        double* res3 = tetat(part,t,p);
        tx = res3[0]; 
        tz = res3[1];
        free(res3);

        // Apply angle changes
        xp = xp + tx;
        zp = zp + tz;

        // Treat single-diffractive scattering.
        if(inter == 4) {
            // added update for s
            s    = (zlm-rlen)+zlm1;
        }

        // Calculate the remaining interaction length and close the iteration loop.
        rlen = rlen-zlm1;
    }

    LocalParticle_set_x(part, x);
    LocalParticle_set_px(part, xp/rpp_in);
    LocalParticle_set_y(part, z);
    LocalParticle_set_py(part, zp/rpp_in);

    result[0] = p;
    result[1] = nabs;
    result[2] = s;
    return result;
}  
  

#endif /* XCOLL_EVEREST_JAW_H */
