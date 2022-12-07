#include <math.h>

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
  


double soln3(double a, double b, double dh, double smax, double s) {


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



double scamcs(double xx, double xxp, double s) {

    double x0  = xx;
    double xp0 = xxp;
    double r2 = 0;
    double v1 = 0;
    double v2 = 0;

    while (1) {
        double v1 = 2*get_random() - 1;
        double v2 = 2*get_random() - 1;
        double r2 = pow(v1,2) + pow(v2,2);

        if(r2 < 1) {
            break;
        }
    }

    double a   = sqrt((-2*log(r2))/r2);
    double z1  = v1*a;
    double z2  = v2*a;
    double ss  = sqrt(s);
    double sss = 1 + 0.038*log(s);
    xx  = x0  + s*(xp0 + ((0.5*ss)*sss)*(z2 + z1*0.577350269));
    xxp = xp0 + (ss*z2)*sss;

    return xx, xxp;

}



double mcs(double s, double mc_radl, double mc_zlm1, double mc_p0, double mc_x, double mc_xp, double mc_z, double mc_zp, double mc_dpop) {

    double theta    = 13.6e-3/(mc_p0*(1+mc_dpop)); // dpop   = (p - p0)/p0;
    double h   = 0.001;
    double dh  = 0.0001;
    double bn0 = 0.4330127019;
    double rlen0 = mc_zlm1/mc_radl;
    double rlen  = rlen0;

    mc_x     = (mc_x/theta)/mc_radl;
    mc_xp    = mc_xp/theta;
    mc_z     = (mc_z/theta)/mc_radl;
    mc_zp    = mc_zp/theta;


    while (1) {
        
        double ae = bn0*mc_x;
        double be = bn0*mc_xp;

        // #######################################
        // ae = np.array(ae, dtype=np.double64)
        // be = np.array(be, dtype=np.double64)
        // dh = np.array(dh, dtype=np.double64)
        // rlen = np.array(rlen, dtype=np.double64)
        // s = np.array(s, dtype=np.double64)
        // #######################################
        s = soln3(ae,be,dh,rlen,s);

        if (s < h) {
            s = h;
        }

        mc_x, mc_xp = scamcs(mc_x,mc_xp,s);

        if (mc_x <= 0) {
            s = (rlen0-rlen)+s;
            break; // go to 20
        }

        if ((s+dh) >= rlen) {
            s = rlen0;
            break; // go to 20
        }
        // go to 10
        rlen = rlen-s;
    }

    mc_z, mc_zp = scamcs(mc_z,mc_zp,s);

    s  = s*mc_radl;
    mc_x  = (mc_x*theta)*mc_radl;
    mc_xp = mc_xp*theta;
    mc_z  = (mc_z*theta)*mc_radl;
    mc_zp = mc_zp*theta;

    return s, mc_x, mc_xp, mc_z, mc_zp, mc_dpop ;


}



double tetat(double t, double p, double tx, double tz) {


    double teta = sqrt(t)/p;
    double va = 0;
    double vb  = 0;
    double va2 = 0;
    double vb2 = 0;
    double r2  = 0;
    
    while (1) {
        va  = 2*get_random() - 1;
        vb  = get_random();
        va2 = pow(va,2);
        vb2 = pow(vb,2);
        r2  = va2 + vb2;

        if(r2 < 1) {
            break;
        }
    }
        
    tx  = (teta*((2*va)*vb))/r2;
    tz  = (teta*(va2 - vb2))/r2;

    return tx, tz;
                
}
                
                
                
double gettran(double inter, double p, double tt_bn, double tt_ecmsq, double tt_xln15s, double tt_bpp, double cgen) {

 
    // Neither if-statements below have an else, so defaulting function return to zero.
    double result = 0;

    if (inter==2) { // Nuclear Elastic
        result = (-1*log(get_random()))/tt_bn;
    }
    
    else if (inter==3) { // pp Elastic
        result = (-1*log(get_random()))/tt_bpp;
    }

    else if (inter==4) { // Single Diffractive
        double xm2 = exp(get_random() * tt_xln15s);
        double bsd = 0;
        p   = p * (1 - xm2/tt_ecmsq);
    
        if (xm2 < 2) {
            bsd = 2 * tt_bpp;
        }

        else if ((xm2 >= 2) & (xm2 <= 5)) {
            bsd = ((106.0 - 17.0*xm2)*tt_bpp)/36.0;
        }

        else {
            bsd = (7*tt_bpp)/12.0;
        }
   
        result = (-1*log(get_random()))/bsd;
    }

    else if (inter==5) { // Coulomb
        result = get_random_ruth();
    }

    return result, p;
                
}
                

double calcionloss(double p, double rlen, double il_exenergy, double il_anuc, double il_zatom, double il_rho, double enlo) {

    double mom    = p*1.0e3; //[GeV/c] -> [MeV/c]
    double enr    = pow((mom*mom + 938.271998*938.271998),0.5); //[MeV]
    double gammar = enr/938.271998;
    double betar  = mom/enr;
    double bgr    = betar*gammar;
    double kine   = ((2*0.510998902)*bgr)*bgr;
    double k = 0.307075;

    // Mean excitation energy
    double exEn = il_exenergy*1.0e3; // [MeV]

    // tmax is max energy loss from kinematics
    double tmax = kine/(1 + (2*gammar)*(0.510998902/938.271998) + pow((0.510998902/938.271998),2)); // [MeV]

    // Plasma energy - see PDG 2010 table 27.1
    double plen = pow(((il_rho*il_zatom)/il_anuc),0.5)*28.816e-6; // [MeV]

    // Calculate threshold energy
    // Above this threshold, the cross section for high energy loss is calculated and then
    // a random number is generated to determine if tail energy loss should be applied, or only mean from Bethe-Bloch
    // below threshold, only the standard Bethe-Bloch is used (all particles get average energy loss)

    // thl is 2*width of Landau distribution (as in fig 27.7 in PDG 2010). See Alfredo's presentation for derivation
    double thl = ((((4*(k*il_zatom))*rlen)*1.0e2)*il_rho)/(il_anuc*pow(betar,2)); // [MeV]

    // Bethe-Bloch mean energy loss
    enlo = ((k*il_zatom)/(il_anuc*pow(betar,2))) * (0.5*log((kine*tmax)/(exEn*exEn)) - pow(betar,2) - log(plen/exEn) - log(bgr) + 0.5);
    enlo = ((enlo*il_rho)*1.0e-1)*rlen; // [GeV]

    // Threshold Tt is Bethe-Bloch + 2*width of Landau distribution
    double Tt = enlo*1.0e3 + thl; // [MeV]

    // Cross section - see Alfredo's presentation for derivation
    double cs_tail = ((k*il_zatom)/(il_anuc*pow(betar,2))) * (0.5*((1/Tt)-(1/tmax)) - (log(tmax/Tt)*pow(betar,2))/(2*tmax) + (tmax-Tt)/((4*pow(gammar,2))*pow(938.271998,2)));

    // Probability of being in tail: cross section * density * path length
    double prob_tail = ((cs_tail*il_rho)*rlen)*1.0e2;

    // Determine based on random number if tail energy loss occurs.
    if (get_random() < prob_tail) {
        enlo = ((k*il_zatom)/(il_anuc*pow(betar,2))) * (0.5*log((kine*tmax)/(exEn*exEn)) - pow(betar,2) - log(plen/exEn) - log(bgr) + 0.5 + pow(tmax,2)/((8*pow(gammar,2))*pow(938.271998,2)));
        enlo = (enlo*il_rho)*1.0e-1; // [GeV/m]
    }
    else {
        // If tail energy loss does not occur, just use the standard Bethe-Bloch
        enlo = enlo/rlen;  // [GeV/m]
    }
        
    return il_exenergy, enlo;
}



int ichoix(int ich_cprob[]) {

    double aran = get_random();
    int i;
    
    for (i = 1; i < 5; ++i) {
        i += 1;
        if (aran <= ich_cprob[i]) {
            break;
        }
    }
    return i;
}

// float calculateSum(float num[]) {
//   float sum = 0.0;

//   for (int i = 0; i < 6; ++i) {
//     sum += num[i];
//   }

//   return sum;
// }

                
double jaw(double run_exenergy, double run_anuc, double run_zatom, double run_rho, double run_radl, int run_cprob[], double run_xintl, double run_bn, double run_ecmsq, double run_xln15s, double run_bpp, double p0, double nabs, double s, double zlm, double x, double xp, double z, double zp, double dpop, double cgen) {
    

    // Note that the input parameter is dpop. Here the momentum p is constructed out of this input.
    double p = p0*(1+dpop);
    nabs = 0;
      
    // Initialize the interaction length to input interaction length
    double rlen = zlm;
    double m_dpodx = 0.;
    double tx = 0.;
    double tz = 0.;

    // Do a step for a point-like interaction.
    // Get monte-carlo interaction length.
    while (1) {

        double run_zlm1 = (-1*run_xintl)*log(get_random());
                        
        // If the monte-carlo interaction length is longer than the
        // remaining collimator length, then put it to the remaining
        // length, do multiple coulomb scattering and return.
        // LAST STEP IN ITERATION LOOP
        if (run_zlm1 > rlen) {
            
            run_zlm1 = rlen;
            s, x, xp, z, zp, dpop = mcs(s,run_radl,run_zlm1,p0,x,xp,z,zp,dpop);

            s = (zlm-rlen)+s;
            run_exenergy, m_dpodx = calcionloss(p,rlen,run_exenergy,run_anuc,run_zatom,run_rho,m_dpodx);  // DM routine to include tail
            p = p-m_dpodx*s;
                    
            dpop = (p-p0)/p0;
            break;
        }
        // Otherwise do multi-coulomb scattering.
        // REGULAR STEP IN ITERATION LOOP
        s, x, xp, z, zp, dpop = mcs(s,run_radl,run_zlm1,p0,x,xp,z,zp,dpop);

        // Check if particle is outside of collimator (X.LT.0) after
        // MCS. If yes, calculate output longitudinal position (s),
        // reduce momentum (output as dpop) and return.
        // PARTICLE LEFT COLLIMATOR BEFORE ITS END.

        if(x <= 0) {

            s = (zlm-rlen)+s;
            run_exenergy, m_dpodx = calcionloss(p,rlen,run_exenergy,run_anuc,run_zatom,run_rho,m_dpodx);

            p = p-m_dpodx*s;
            dpop = (p-p0)/p0;
            break;
        }

        // Check whether particle is absorbed. If yes, calculate output
        // longitudinal position (s), reduce momentum (output as dpop)
        // and return.
        // PARTICLE WAS ABSORBED INSIDE COLLIMATOR DURING MCS.

        int inter = ichoix(run_cprob);
        nabs = inter;

        if (inter == 1) {

            s = (zlm-rlen)+run_zlm1;
            run_exenergy, m_dpodx = calcionloss(p,rlen,run_exenergy,run_anuc,run_zatom,run_rho,m_dpodx);

            p = p-m_dpodx*s;
            dpop = (p-p0)/p0;

            break;
        }


        // Now treat the other types of interaction, as determined by ICHOIX:

        // Nuclear-Elastic:          inter = 2
        // pp Elastic:               inter = 3
        // Single-Diffractive:       inter = 4    (changes momentum p)
        // Coulomb:                  inter = 5

        // Gettran returns some monte carlo number, that, as I believe, gives the rms transverse momentum transfer.

        double t, p = gettran(inter,p,run_bn,run_ecmsq,run_xln15s,run_bpp,cgen);

        // Tetat calculates from the rms transverse momentum transfer in
        // monte-carlo fashion the angle changes for x and z planes. The
        // angle change is proportional to SQRT(t) and 1/p, as expected.

        tx, tz = tetat(t,p,tx,tz);

        // Apply angle changes
        xp = xp+tx;
        zp = zp+tz;

        // Treat single-diffractive scattering.
        if(inter == 4) {
            // added update for s
            s    = (zlm-rlen)+run_zlm1;

            // Add this code to get the momentum transfer also in the calling routine
            dpop = (p-p0)/p0;
        }

        // Calculate the remaining interaction length and close the iteration loop.
        rlen = rlen-run_zlm1;
    }
                  
    return run_exenergy, run_bn, p0, nabs, s, zlm, x, xp, z, zp, dpop;

}  
  