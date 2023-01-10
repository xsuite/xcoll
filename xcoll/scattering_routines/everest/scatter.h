#include <math.h>
#include <stdlib.h>

double rutherford(double t, double zatom, double emr) {
    double cnorm  = 2.607e-5;
    double cnform = 0.8561e3;
    double result = (cnorm*exp(((-1*t)*cnform)*pow(emr,2))) * pow((zatom/t),2);
    return result;
}

// isimp and linside are booleans

double* scatter(double x_in, double  xp_in, double  y_in, double  yp_in, double  s_in, double  p_in, double  val_part_hit, double 
                val_part_abs, double  val_part_impact, double  val_part_indiv, double  val_part_linteract, double  val_nabs_type, double  val_linside, double  run_exenergy, double  
                run_anuc, double  run_zatom, double  run_emr, double  run_rho, double   run_hcut, double  run_bnref, double  run_csref0, double  run_csref1, double  run_csref5, double run_radl, double  
                run_dlri, double  run_dlyi, double  run_eum, double  run_ai, double  run_collnt, double  cprob0, double cprob1, double cprob2, double cprob3, double cprob4, double cprob5, double  run_xintl, double  run_bn, double  
                run_ecmsq, double  run_xln15s, double  run_bpp, double  is_crystal, double  length, double  c_rotation, double  c_aperture, double  c_offset, double  c_tilt0, double c_tilt1, double  
                onesided, double  cry_tilt, double  cry_rcurv, double  cry_bend, double  cry_alayer, double  cry_xmax, double  cry_ymax, double  cry_orient, double  cry_miscut, double  cry_cBend, double  
                cry_sBend, double  cry_cpTilt, double  cry_spTilt, double  cry_cnTilt, double  cry_snTilt, double  p0, double  x0, double  xp0, double  nhit, double  nabs, double  fracab, double  nnuc0, double  ien0, double  nnuc1, double  
                ien1, double  iProc, double  n_chan, double  n_VR, double  n_amorphous, double  s_imp) {

    static double result[28];
    
    val_part_impact = -1.;
    val_part_linteract = -1.;
    val_part_indiv = -1.;

    double x = x_in;
    double xp = xp_in;
    double z = y_in;
    double zp = yp_in;
    double p = p_in;
    double sp = 0;
    double dpop = (p - p0)/p0;
    double tiltangle = 0.;

    double mirror = 1.;

    // Compute rotation factors for collimator rotation
    double cRot   = cos(c_rotation);
    double sRot   = sin(c_rotation);
    double cRRot  = cos(-c_rotation);
    double sRRot  = sin(-c_rotation);

    // Transform particle coordinates to get into collimator coordinate  system
    // First do rotation into collimator frame
    x  =  x_in*cRot + sRot*y_in;
    z  =  y_in*cRot - sRot*x_in;
    xp = xp_in*cRot + sRot*yp_in;
    zp = yp_in*cRot - sRot*xp_in;


    // For one-sided collimators consider only positive X. For negative X jump to the next particle
    if (onesided && (x < 0.)) {
        result[0] = x_in;
        result[1] = xp_in;
        result[2] = y_in;
        result[3] = yp_in;
        result[4] = s_in;
        result[5] = p_in;
        result[6] = val_part_hit;
        result[7] = val_part_abs;
        result[8] = val_part_impact;
        result[9] = val_part_indiv;
        result[10] = val_part_linteract;
        result[11] = val_nabs_type;
        result[12] = val_linside;
        result[13] = p0;
        result[14] = x0;
        result[15] = xp0;
        result[16] = nhit;
        result[17] = nabs;
        result[18] = fracab;
        result[19] = nnuc0;
        result[20] = ien0;
        result[21] = nnuc1;
        result[22] = ien1;
        result[23] = iProc;
        result[24] = n_chan;
        result[25] = n_VR;
        result[26] = n_amorphous;
        result[27] = s_imp;

        return result;
    }
    // Log input energy + nucleons as per the FLUKA coupling
    nnuc0 = nnuc0 + 1.;
    ien0 = ien0 + p_in * 1.0e3;

    // Now mirror at the horizontal axis for negative X offset
    if (x < 0) {
        mirror    = -1;
        tiltangle = -1*c_tilt1;
    }
    else {
        mirror    = 1;
        tiltangle = c_tilt0;
    }
    x  = mirror*x;
    xp = mirror*xp;

    // Shift with opening and offset
    x = (x - c_aperture/2.) - mirror*c_offset;

    // Include collimator tilt
    if (tiltangle > 0.) {
        xp = xp - tiltangle;
    }
    if (tiltangle < 0.) {
        x  = x + sin(tiltangle) * length;
        xp = xp - tiltangle;
    }

    // particle passing above the jaw are discarded => take new event
    // entering by the face, shorten the length (zlm) and keep track of
    // entrance longitudinal coordinate (keeps) for histograms

    // The definition is that the collimator jaw is at x>=0.

    // 1) Check whether particle hits the collimator
    int isimp = 0;
    double s     = 0.;
    double zlm = -1*length;


    if (is_crystal) {

        double* crystal_result = crystal(x,
                                xp,
                                z,
                                zp,
                                s,
                                p,
                                x0,
                                xp0,
                                zlm,
                                s_imp,
                                isimp,
                                val_part_hit, 
                                val_part_abs, 
                                val_part_impact, 
                                val_part_indiv, 
                                length, 
                                run_exenergy, 
                                run_rho, 
                                run_anuc, 
                                run_zatom, 
                                run_emr, 
                                run_dlri, 
                                run_dlyi, 
                                run_ai, 
                                run_eum, 
                                run_collnt,                                                                                                                             
                                run_hcut, 
                                run_bnref, 
                                run_csref0, 
                                run_csref1, 
                                run_csref5,                                                                                                                             
                                nhit, 
                                nabs,
                                cry_tilt,
                                cry_rcurv,
                                cry_bend,
                                cry_alayer,
                                cry_xmax,
                                cry_ymax,
                                cry_orient,
                                cry_miscut,
                                cry_cBend,
                                cry_sBend,
                                cry_cpTilt,
                                cry_spTilt,
                                cry_cnTilt,
                                cry_snTilt,
                                iProc,
                                n_chan,
                                n_VR,
                                n_amorphous
                                );

        val_part_hit = crystal_result[0];
        val_part_abs = crystal_result[1];
        val_part_impact = crystal_result[2];
        val_part_indiv = crystal_result[3];
        nhit = crystal_result[4];
        nabs = crystal_result[5];
        s_imp = crystal_result[6];
        isimp = crystal_result[7];
        s = crystal_result[8];
        zlm = crystal_result[9];
        x = crystal_result[10];
        xp = crystal_result[11];
        x0 = crystal_result[12];
        xp0 = crystal_result[13];
        z = crystal_result[14];
        zp = crystal_result[15];
        p = crystal_result[16];
        iProc = crystal_result[17];
        n_chan = crystal_result[18];
        n_VR = crystal_result[19];
        n_amorphous = crystal_result[20];

        if (nabs != 0.) {
            val_part_abs = 1.;
            val_part_linteract = zlm;
        }
        s_imp  = (s - length) + s_imp;
    }

    else {

        if (x >= 0.) {
        // Particle hits collimator and we assume interaction length ZLM equal
        // to collimator length (what if it would leave collimator after
        // small length due to angle???)
            zlm  = length;
            val_part_impact = x;
            val_part_indiv  = xp;
        }
        else if (xp <= 0.) {
        // Particle does not hit collimator. Interaction length ZLM is zero.
            zlm = 0.;
        }
        else {
        // Calculate s-coordinate of interaction point
            s = -x/xp;
        
            if (s < length) {
                zlm = length - s;
                val_part_impact = 0.;
                val_part_indiv  = xp;
            }
            else {
                zlm = 0.;
            }
        
        }
        // First do the drift part
        // DRIFT PART
        double drift_length = length - zlm;
        if (drift_length > 0) {
            x  = x  + xp * drift_length;
            z  = z  + zp * drift_length;
            sp = sp + drift_length;
        }
        // Now do the scattering part
        if (zlm > 0.) {
            if (!val_linside) {
            // first time particle hits collimator: entering jaw
                val_linside = 1;
            }
            nhit = nhit + 1;

            
            double* jaw_result = jaw(run_exenergy,
                            run_anuc,
                            run_zatom,
                            run_rho,
                            run_radl,
                            cprob0,
                            cprob1,
                            cprob2,
                            cprob3,
                            cprob4,
                            cprob5,
                            run_xintl,
                            run_bn,
                            run_ecmsq,
                            run_xln15s,
                            run_bpp,
                            p0,
                            nabs,
                            s,
                            zlm,
                            x,
                            xp,
                            z,
                            zp,
                            dpop);

            run_exenergy = jaw_result[0];
            run_bn = jaw_result[1];
            p0 = jaw_result[2];
            nabs = jaw_result[3];
            s = jaw_result[4];
            zlm = jaw_result[5];
            x = jaw_result[6];
            xp = jaw_result[7];
            z = jaw_result[8];
            zp = jaw_result[9];
            dpop = jaw_result[10];

            val_nabs_type = nabs;
            val_part_hit  = 1;

            isimp = 1;
            // simp  = s_impact

            // Writeout should be done for both inelastic and single diffractive. doing all transformations
            // in x_flk and making the set to 99.99 mm conditional for nabs=1
            if (nabs == 1 || nabs == 4) {
            // Transform back to lab system for writeout.
            // keep x,y,xp,yp unchanged for continued tracking, store lab system variables in x_flk etc

            // Finally, the actual coordinate change to 99 mm
                if (nabs == 1) {
                    fracab = fracab + 1;
                    x = 99.99e-3;
                    z = 99.99e-3;
                    val_part_linteract = zlm;
                    val_part_abs = 1;
                // Collimator jaw interaction
                }
            }

            if (nabs != 1. && zlm > 0.) {
            // Do the rest drift, if particle left collimator early
                drift_length = (length-(s+sp));

                if (drift_length > 1.0e-15) {
                    val_linside = 0;
                    x  = x  + xp * drift_length;
                    z  = z  + zp * drift_length;
                    sp = sp + drift_length;
                }
                val_part_linteract = zlm - drift_length;
            }

        }
    }

    // Transform back to particle coordinates with opening and offset
    if (x < 99.0e-3) {
        // Include collimator tilt
        if (tiltangle > 0) {
            x  = x  + tiltangle*length;
            xp = xp + tiltangle;
        }
        else if (tiltangle < 0) {
            x  = x  + tiltangle*length;
            xp = xp + tiltangle;
            x  = x  - sin(tiltangle) * length;
        }

        // Transform back to particle coordinates with opening and offset
        x   = (x + c_aperture/2) + mirror*c_offset;

        // Now mirror at the horizontal axis for negative X offset
        x  = mirror * x;
        xp = mirror * xp;

        // Last do rotation into collimator frame
        x_in  =  x*cRRot +  z*sRRot;
        y_in  =  z*cRRot -  x*sRRot;
        xp_in = xp*cRRot + zp*sRRot;
        yp_in = zp*cRRot - xp*sRRot;

        // Log output energy + nucleons as per the FLUKA coupling
        // Do not log dead particles
        nnuc1       = nnuc1 + 1;                           // outcoming nucleons
        ien1        = ien1  + p_in * 1e3;                 // outcoming energy

        if (is_crystal) {
            p_in = p;
            s_in = s_in + s;
        }
        else {
            p_in = (1 + dpop) * p0;
            s_in = sp;
        }
    }

    else {
        x_in = x;
        y_in = z;
    }

    result[0] = x_in;
    result[1] = xp_in;
    result[2] = y_in;
    result[3] = yp_in;
    result[4] = s_in;
    result[5] = p_in;
    result[6] = val_part_hit;
    result[7] = val_part_abs;
    result[8] = val_part_impact;
    result[9] = val_part_indiv;
    result[10] = val_part_linteract;
    result[11] = val_nabs_type;
    result[12] = val_linside;
    result[13] = p0;
    result[14] = x0;
    result[15] = xp0;
    result[16] = nhit;
    result[17] = nabs;
    result[18] = fracab;
    result[19] = nnuc0;
    result[20] = ien0;
    result[21] = nnuc1;
    result[22] = ien1;
    result[23] = iProc;
    result[24] = n_chan;
    result[25] = n_VR;
    result[26] = n_amorphous;
    result[27] = s_imp;

    return result;

}