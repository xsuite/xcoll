// copyright ############################### #
// This file is part of the Xcoll Package.   #
// Copyright (c) CERN, 2023.                 #
// ######################################### #

#ifndef XCOLL_EVEREST_JAW_H
#define XCOLL_EVEREST_JAW_H
#include <math.h>
#include <stdio.h>



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
double* gettran(EverestData restrict everest, LocalParticle* part, double inter, double p) {

    double* res = (double*)malloc(2 * sizeof(double));

    // Neither if-statements below have an else, so defaulting function return to zero.
    double result = 0;

    if (inter==2) { // Nuclear Elastic
        result = RandomExponential_generate(part)/everest->bn;

    } else if (inter==3) { // pp Elastic
        result = RandomExponential_generate(part)/everest->bpp;

    } else if (inter==4) { // Single Diffractive
        double xm2 = exp(RandomUniform_generate(part) * everest->xln15s);
        double bsd = 0;
        p = p * (1 - xm2/everest->ecmsq);
    
        if (xm2 < 2) {
            bsd = 2 * everest->bpp;
        } else if ((xm2 >= 2) & (xm2 <= 5)) {
            bsd = ((106.0 - 17.0*xm2)*everest->bpp)/36.0;
        } else {
            bsd = (7*everest->bpp)/12.0;
        }
        result = RandomExponential_generate(part)/bsd;

    } else if (inter==5) { // Coulomb
        result = RandomRutherford_generate(everest->coll->rng, part);
    }

    res[0] = result;
    res[1] = p;
    return res;
}


/*gpufun*/
int ichoix(EverestData restrict everest, LocalParticle* part) {

    double aran = RandomUniform_generate(part);
    int i;
    for (i = 0; i < 5; ++i) {
        if (aran < everest->cprob[i]) {
            break;
        }
    }
    return i;
}


/*gpufun*/
double* jaw(EverestData restrict everest, LocalParticle* part, double p, double zlm, int only_mcs, int edge_check) {

    double* result = (double*)malloc(3 * sizeof(double));

    double s;
    double nabs = 0;
    double rlen = zlm;
    double m_dpodx = 0.;
    double t;
    double tx; 
    double tz;

    double rpp_in = LocalParticle_get_rpp(part);
    double x  = LocalParticle_get_x(part);
    double xp = LocalParticle_get_px(part)*rpp_in;
    double z  = LocalParticle_get_y(part);
    double zp = LocalParticle_get_py(part)*rpp_in;

    if (only_mcs) {
        double* res = mcs(everest, part, zlm, p, x, xp, z, zp, edge_check);
        s = res[0];
        x = res[1];
        xp = res[2];
        z = res[3];
        zp = res[4];
        free(res);

    } else {
        // Do a step for a point-like interaction.
        // Get monte-carlo interaction length.
        while (1) {

            calculate_ionisation_properties(everest, p);
            double zlm1 = everest->xintl*RandomExponential_generate(part);

            // If the monte-carlo interaction length is longer than the
            // remaining collimator length, then put it to the remaining
            // length, do multiple coulomb scattering and return.
            // LAST STEP IN ITERATION LOOP
            if (zlm1 > rlen) {
                zlm1 = rlen;
                double* res = mcs(everest, part, zlm1, p, x, xp, z, zp, edge_check);
                s = res[0];
                x = res[1];
                xp = res[2];
                z = res[3];
                zp = res[4];
                free(res);

                s = zlm - rlen + s;
                m_dpodx = calcionloss(everest, part, rlen);  // DM routine to include tail // TODO: should not be rlen but s after updating
                p = p-m_dpodx*s; // This is correct: ionisation loss is only calculated and applied at end of while (break)
                break;
            }
            // Otherwise do multi-coulomb scattering.
            // REGULAR STEP IN ITERATION LOOP
            double* res1 = mcs(everest, part, zlm1, p, x, xp, z, zp, edge_check);
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
                m_dpodx = calcionloss(everest, part, rlen);  // TODO: should not be rlen but s after updating
                p = p-m_dpodx*s;  // correct
                break;
            }

            // Check whether particle is absorbed. If yes, calculate output
            // longitudinal position (s), reduce momentum (output as dpop)
            // and return.
            // PARTICLE WAS ABSORBED INSIDE COLLIMATOR DURING MCS.

            int inter = ichoix(everest, part);
            nabs = inter;
            if (inter == 1) {
                s = (zlm-rlen)+zlm1;
                m_dpodx = calcionloss(everest, part, rlen);
                p = p-m_dpodx*s;
                break;
            }


            // Now treat the other types of interaction, as determined by ICHOIX:

            // Nuclear-Elastic:          inter = 2
            // pp Elastic:               inter = 3
            // Single-Diffractive:       inter = 4    (changes momentum p)
            // Coulomb:                  inter = 5

            // Gettran returns some monte carlo number, that, as I believe, gives the rms transverse momentum transfer.

            double* res2 = gettran(everest, part, inter, p);
            t = res2[0];
            p = res2[1];
            free(res2);

            // Tetat calculates from the rms transverse momentum transfer in
            // monte-carlo fashion the angle changes for x and z planes. The
            // angle change is proportional to SQRT(t) and 1/p, as expected.

            double* res3 = tetat(part, t, p);
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
