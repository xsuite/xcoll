// copyright ############################### #
// This file is part of the Xcoll Package.   #
// Copyright (c) CERN, 2023.                 #
// ######################################### #

#ifndef XCOLL_EVEREST_CRYSTAL_INTERACT_H
#define XCOLL_EVEREST_CRYSTAL_INTERACT_H
#include <math.h>
#include <stdio.h>
#include <stdlib.h>


// Processes
int const proc_out         =  -1;     // Crystal not hit
int const proc_AM          =   1;     // Amorphous
int const proc_VR          =   2;     // Volume reflection
int const proc_CH          =   3;     // Channeling
int const proc_VC          =   4;     // Volume capture
int const proc_absorbed    =   5;     // Absorption
int const proc_DC          =   6;     // Dechanneling
int const proc_pne         =   7;     // Proton-neutron elastic interaction
int const proc_ppe         =   8;     // Proton-proton elastic interaction
int const proc_diff        =   9;     // Single diffractive
int const proc_ruth        =  10;     // Rutherford scattering
int const proc_ch_absorbed =  15;     // Channeling followed by absorption
int const proc_ch_pne      =  17;     // Channeling followed by proton-neutron elastic interaction
int const proc_ch_ppe      =  18;     // Channeling followed by proton-proton elastic interaction
int const proc_ch_diff     =  19;     // Channeling followed by single diffractive
int const proc_ch_ruth     =  20;     // Channeling followed by Rutherford scattering
int const proc_TRVR        = 100;     // Volume reflection in VR-AM transition region
int const proc_TRAM        = 101;     // Amorphous in VR-AM transition region


/*gpufun*/
double* movech(RandomRutherfordData rng, LocalParticle* part, double dz, double pc, double r, double rc, 
               CrystalMaterialData material, double iProc) {

    double* result = (double*)malloc(2 * sizeof(double));

    // Material properties
    double const anuc   = CrystalMaterialData_get_A(material);
    double const rho    = CrystalMaterialData_get_density(material);
    double const eum    = CrystalMaterialData_get_crystal_potential(material);
    double       collnt = CrystalMaterialData_get_nuclear_collision_length(material);

    double pc_in = pc;

    //Rescale the total and inelastic cross-section accordigly to the average density seen
    double rpp = LocalParticle_get_rpp(part);
    double x_i = LocalParticle_get_x(part);
    double xp  = LocalParticle_get_px(part)*rpp;
    int np  = x_i/XC_PLANE_DISTANCE;    //Calculate in which crystalline plane the particle enters
    x_i = x_i - (np + 0.5)*XC_PLANE_DISTANCE;   //Rescale the incoming x to the middle of the crystalline plane

    double pv   = pow(pc, 2.)/sqrt(pow(pc, 2.) + pow(XC_PROTON_MASS*1.0e-3, 2.))*1.0e9;          //Calculate pv=P/E   TODO: this is beta?
    double Ueff = eum*4.*pow(x_i/XC_PLANE_DISTANCE, 2.) + pv*x_i/r; //Calculate effective potential
    double Et   = pv*pow(xp, 2.)/2. + Ueff;          //Calculate transverse energy
    double Ec   = eum*pow(1. - rc/r, 2.);              //Calculate critical energy in bent crystals

    //To avoid negative Et
    double xminU = -pow(XC_PLANE_DISTANCE, 2.)*pc*1.0e9/(8.*eum*r);
    double Umin  = fabs(eum*4.*pow(xminU/XC_PLANE_DISTANCE, 2.) + pv*xminU/r);
    Et    = Et + Umin;
    Ec    = Ec + Umin;

    //Calculate min e max of the trajectory between crystalline planes
    double x_min = (-(XC_PLANE_DISTANCE/2.)*rc)/r - (XC_PLANE_DISTANCE/2.)*sqrt(Et/Ec);
    double x_max = (-(XC_PLANE_DISTANCE/2.)*rc)/r + (XC_PLANE_DISTANCE/2.)*sqrt(Et/Ec);

    //Change ref. frame and go back with 0 on the crystalline plane on the left
    x_min = x_min - XC_PLANE_DISTANCE/2.;
    x_max = x_max - XC_PLANE_DISTANCE/2.;

    //Calculate the "normal density" in m^-3
    double N_am  = rho*6.022e23*1.0e6/anuc;

    //Calculate atomic density at min and max of the trajectory oscillation
    // erf returns the error function of complex argument
    double rho_max = N_am*XC_PLANE_DISTANCE/2.*(erf(x_max/sqrt(2*pow(XC_THERMAL_VIBRATIONS, 2.))) - erf((XC_PLANE_DISTANCE-x_max)/sqrt(2*pow(XC_THERMAL_VIBRATIONS, 2.))));
    double rho_min = N_am*XC_PLANE_DISTANCE/2.*(erf(x_min/sqrt(2*pow(XC_THERMAL_VIBRATIONS, 2.))) - erf((XC_PLANE_DISTANCE-x_min)/sqrt(2*pow(XC_THERMAL_VIBRATIONS, 2.))));

    //"zero-approximation" of average nuclear density seen along the trajectory
    double avrrho  = (rho_max - rho_min)/(x_max - x_min);
    avrrho  = 2.*avrrho/N_am;

    struct ScatteringParameters scat = calculate_scattering(pc, (GeneralMaterialData) material, avrrho);

    //Multiple Coulomb Scattering
    //Can nuclear interaction happen?
    //Rescaled nuclear collision length
    if (avrrho == 0) {
        collnt = 1.0e6;
    } else {
        collnt = collnt/avrrho;
    }

    double zlm = collnt*RandomExponential_generate(part);

    if (zlm < dz) {

        //Choose nuclear interaction
        double aran = RandomUniform_generate(part);
        int i = 1;

        while (aran > scat.cprob[i]) {
            i=i+1;
        }
        
        int ichoix = i;

        //Do the interaction
        double t = 0 ; //default value to cover ichoix=1
        
        if (ichoix==1) {
            iProc = proc_ch_absorbed; //deep inelastic, impinging p disappeared
        } else if (ichoix==2) { //p-n elastic
            iProc = proc_ch_pne;
            t     = RandomExponential_generate(part)/scat.bn;
        } else if (ichoix==3) { //p-p elastic
            iProc = proc_ch_ppe;
            t     = RandomExponential_generate(part)/scat.bpp;
        } else if (ichoix==4) { //Single diffractive
            iProc = proc_ch_diff;
            double xm2 = exp(RandomUniform_generate(part)*scat.xln15s);
            pc = pc*(1 - xm2/scat.ecmsq);

            double bsd = 0.0;
            if (xm2 < 2.) {
                bsd = 2*scat.bpp;
            } else if (xm2 >= 2. && xm2 <= 5.) {
                bsd = ((106.0 - 17.0*xm2)*scat.bpp)/36.0;
            } else if (xm2 > 5.) {
                bsd = (7*scat.bpp)/12.0;
            }
            //end if
            t = RandomExponential_generate(part)/bsd;
        } else { //(ichoix==5)
            iProc = proc_ch_ruth;
            t = RandomRutherford_generate(rng, part);
        }

        double teta;
        //Calculate the related kick -----------
        if (ichoix == 4) {
            teta = sqrt(t)/pc_in; //DIFF has changed PC!!!
        } else {
            teta = sqrt(t)/pc;
        }

        double tx = teta*RandomNormal_generate(part);
        double tz = teta*RandomNormal_generate(part);

        //Change p angle
        LocalParticle_add_to_px(part, tx/rpp);
        LocalParticle_add_to_py(part, tz/rpp);
    }

    result[0] = pc;
    result[1] = iProc;
    return result;
}


/*gpufun*/
double* moveam(RandomRutherfordData rng, LocalParticle* part, double dz, double dei, double pc,
              CrystalMaterialData material, double iProc) {

//     double iProc = proc_out;
    double* result = (double*)malloc(2 * sizeof(double));

    struct ScatteringParameters scat = calculate_scattering(pc, (GeneralMaterialData) material, 1.);

    // Material properties
    double const dlr      = CrystalMaterialData_get_crystal_radiation_length(material);
    double const collnt   = CrystalMaterialData_get_nuclear_collision_length(material);

    double pc_in = pc;

    // Multiple Coulomb Scattering
    pc  = pc - dei*dz; // Energy lost because of ionization process[GeV]

    double dya   = (13.6/pc)*sqrt(dz/dlr); // RMS of coloumb scattering MCS (mrad)
    double kxmcs = dya*RandomNormal_generate(part)*1.0e-3;
    double kymcs = dya*RandomNormal_generate(part)*1.0e-3;

    double rpp = LocalParticle_get_rpp(part);
    LocalParticle_add_to_px(part, kxmcs/rpp);
    LocalParticle_add_to_py(part, kymcs/rpp);

    // Can nuclear interaction happen?
    double zlm = collnt*RandomExponential_generate(part);

    if (zlm < dz) {
        // Choose nuclear interaction
        double aran = RandomUniform_generate(part);
        int i=1;

        while (aran > scat.cprob[i]) {
            i = i+1;
            //goto 10
        }

        int ichoix = i;

        // Do the interaction
        double t = 0 ;// default value to cover ichoix=1
        if (ichoix==1) {
        //case(1) // Deep inelastic, impinging p disappeared
            iProc = proc_absorbed;
        } else if (ichoix==2) { // p-n elastic
            iProc = proc_pne;
            t     = RandomExponential_generate(part)/scat.bn;
        } else if (ichoix==3) { // p-p elastic
            iProc = proc_ppe;
            t     = RandomExponential_generate(part)/scat.bpp;
        } else if (ichoix==4) { // Single diffractive
            iProc = proc_diff;
            double xm2 = exp(RandomUniform_generate(part)*scat.xln15s);
            double bsd = 0.0;
            pc    = pc*(1 - xm2/scat.ecmsq);

            if(xm2 < 2) {
                bsd = 2*scat.bpp;
            } else if(xm2 >= 2 && xm2 <= 5) {
                bsd = ((106.0 - 17.0*xm2)*scat.bpp)/36.0;
            } else if(xm2 > 5) {
                bsd = 7.0*scat.bpp/12.0;
            }
        
            t = RandomExponential_generate(part)/bsd;
        } else { //(ichoix==5)
            iProc = proc_ruth;
            t = RandomRutherford_generate(rng, part);
        }

        // end select

        double teta;

        // Calculate the related kick
        if(ichoix == 4) {
            teta = sqrt(t)/pc_in ;// DIFF has changed PC
        } else {
            teta = sqrt(t)/pc;
        }

        double tx = teta*RandomNormal_generate(part);
        double tz = teta*RandomNormal_generate(part);

        // Change p angle
        LocalParticle_add_to_px(part, tx/rpp);
        LocalParticle_add_to_py(part, tz/rpp);
        
    }

    result[0] = pc;
    result[1] = iProc;
    return result; // Turn on/off nuclear interactions
}


/*gpufun*/
double calcionloss_cry(LocalParticle* part, double length, struct CrystalProperties prop) {

    double prob_tail = prop.prob_tail_c1 + prop.prob_tail_c2 * length
                     + prop.prob_tail_c3 * length * log(length) + prop.prob_tail_c4 * length * length;

    if (RandomUniform_generate(part) < prob_tail) {
        return prop.energy_loss_tail;
    } else {
        return prop.energy_loss;
    }
}

/*gpufun*/
double Channel_single_particle_4d(LocalParticle* part, double arc_length, double bend_ang, double miscut, double sigma_ran) {
    // Channeling: happens over an arc length L_chan (potentially less if dechanneling)
    //             This equates to an opening angle tP wrt. to the point P (center of miscut)
    //             The angle xp at the start of channeling (I) is tP/2 + miscut
    //             The angle xp at the end of channeling (F) is tP + miscut
    //             In practice: we drift from start to end, but also update the angle afterwards

    // TODO: random angle distribution: only at exit?

    // The distance from I to F is the chord length of the angle tP: d = 2 r sin(tP/2)
    // Hence the longitudinal distance (the length to be drifted) is the projection of this using the
    // xp at the start of channeling: s = 2 r sin(tP/2)cos(tP/2 + miscut)
    //                                  = 2 r sin(tP/2)cos(tP/2)cos(miscut) - 2 r sin(tP/2)sin(tP/2)sin(miscut)
    //                                  = r sin(tP)cos(miscut) - 2 r sin(tP/2)^2 sin(miscut)
    //                                  = L_chan/tP ( sin(tP)cos(miscut) - 2 sin(tP/2)^2 sin(miscut) )

    double drift_length = arc_length/bend_ang * (
                            sin(bend_ang)*cos(miscut)
                            - 2.*sin(bend_ang/2.)*sin(bend_ang/2.)*sin(miscut)
                          );

    double rpp = LocalParticle_get_rpp(part);
    LocalParticle_set_px(part, (bend_ang/2. + miscut)/rpp);          // Angle at start of channeling
    Drift_single_particle_4d(part, drift_length);
    double ran_angle = 0.5*RandomNormal_generate(part)*sigma_ran;    // Extra random angle spread at exit
    LocalParticle_set_px(part, (bend_ang + miscut + ran_angle)/rpp); // Angle at end of channeling
    return drift_length;
}


/*gpufun*/
double* interact(RandomRutherfordData rng, LocalParticle* part, double pc,
                 double length, double s_P, double x_P, CrystalMaterialData material, double cry_tilt, double bend_r,
                 double cry_alayer, double xdim, double ydim, double cry_orient, double cry_miscut, double bend_ang,
                 double cry_cBend, double cry_sBend, double cry_cpTilt, double cry_spTilt, double cry_cnTilt,
                 double cry_snTilt, double iProc, CollimatorImpactsData record, RecordIndex record_index) {

    double* result = (double*)malloc(6 * sizeof(double));

    // Material properties
    double const exenergy = CrystalMaterialData_get_excitation_energy(material);
    double const rho      = CrystalMaterialData_get_density(material);
    double const anuc     = CrystalMaterialData_get_A(material);
    double const zatom    = CrystalMaterialData_get_Z(material);
    double const emr      = CrystalMaterialData_get_nuclear_radius(material);
    double const dlri     = CrystalMaterialData_get_crystal_radiation_length(material);
//     double const dlyi     = CrystalMaterialData_get_crystal_nuclear_length(material);
    double const ai       = CrystalMaterialData_get_crystal_plane_distance(material);
    double const eum      = CrystalMaterialData_get_crystal_potential(material);
    double const collnt   = CrystalMaterialData_get_nuclear_collision_length(material);
    double const hcut     = CrystalMaterialData_get_hcut(material);
    double const bnref    = CrystalMaterialData_get_nuclear_elastic_slope(material);
    double const csref0   = CrystalMaterialData_get_cross_section(material, 0);
    double const csref1   = CrystalMaterialData_get_cross_section(material, 1);
    double const csref5   = CrystalMaterialData_get_cross_section(material, 5);

    double const rpp  = LocalParticle_get_rpp(part);
    double const x  = LocalParticle_get_x(part);
    double const xp = LocalParticle_get_px(part)*rpp;
    double const y  = LocalParticle_get_y(part);
    double const yp = LocalParticle_get_py(part)*rpp;

    double energy_loss = 0.;

    double c_v1 =  1.7;   // Fitting coefficient
    double c_v2 = -1.5;   // Fitting coefficient

    int zn  = 1;

    struct CrystalProperties properties = calculate_crystal_properties(pc, material);

    double const_dech = properties.const_dech;

    // MISCUT second step: fundamental coordinates (crystal edges and plane curvature radius)
    double s_K = length;
    double x_K = bend_r * (1.-cos(bend_ang));
    double s_M = (bend_r-xdim) * sin(bend_ang);
    double x_M = xdim + (bend_r-xdim)*(1.-cos(bend_ang));
    double r   = sqrt(pow(s_P,2.) + pow((x-x_P),2.));

    // MISCUT third step: F coordinates (exit point) on crystal exit face
    double A_F = pow((tan(bend_ang)),2.) + 1.;
    double B_F = ((-2)*pow((tan(bend_ang)),2.))*bend_r + (2.*tan(bend_ang))*s_P - 2.*x_P;
    double C_F = pow((tan(bend_ang)),2.)*(pow(bend_r,2.)) - ((2.*tan(bend_ang))*s_P)*bend_r + pow(s_P,2.) + pow(x_P,2.) - pow(r,2.);

    // Coordinates of exit point F
    double x_F = (-B_F-sqrt(pow(B_F,2.)-4.*(A_F*C_F)))/(2.*A_F);
    double s_F = (-tan(bend_ang))*(x_F-bend_r);

    if (x_F < x_K || x_F > x_M || s_F < s_M || s_F > s_K) {

        double alpha_F;
        double beta_F;
        
        if (cry_miscut == 0 && fabs(x_F-x_K) <= 1.0e-13 && fabs(s_F-s_K) <= 1.0e3) {
        // no miscut, entrance from below: exit point is K (lower edge)
            x_F = x_K;
            s_F = s_K;
        } else if (cry_miscut == 0 && fabs(x_F-x_M) <= 1.0e3 && fabs(s_F-s_M) <= 1.0e3) {
        // no miscut, entrance from above: exit point is M (upper edge)
            x_F = x_M;
            s_F = s_M;
        } else {
        // MISCUT Third step (bis): F coordinates (exit point)  on bent side
            if (cry_miscut < 0) {
            // Intersect with bottom side
                alpha_F = (bend_r-x_P)/x_P;
                beta_F = -(pow(r,2.)-pow(s_P,2.)-pow(x_P,2.))/(2*s_P);
                A_F = pow(alpha_F,2.) + 1.;
                B_F = 2.*(alpha_F*beta_F) - 2.*bend_r;
                C_F = pow(beta_F,2.);
            } else {
            // Intersect with top side
                alpha_F = (bend_r-x_P)/s_P;
                beta_F = -(pow(r,2.)+xdim*(xdim-(2.*bend_r))-pow(s_P,2.)-pow(x_P,2.))/(2.*s_P);
                A_F = pow(alpha_F,2.) + 1.;
                B_F = 2.*(alpha_F*beta_F) - 2.*bend_r;
                C_F = pow(beta_F,2.) - xdim*(xdim-2.*bend_r);
            }
            
            x_F = (-B_F-sqrt(pow(B_F,2.)-4.*(A_F*C_F)))/(2.*A_F);
            s_F = alpha_F*x_F + beta_F;
        }
    }

    // MISCUT 4th step: deflection and length calculation
    double a = sqrt(pow(s_F,2.) + pow((x-x_F),2.));  // Line from entrance point I to exit point F
    double tP = acos((2.*pow(r,2.) - pow(a,2.))/(2.*pow(r,2.)));  // Angle in P that spans I to F
    double tdefl = asin((s_F-s_P)/r); // final total deflection in crystal reference frame (exit angle)
    double L_chan = r*tP;   // channeling length (if no miscut, this would be r*bend_ang): along the trajectory

    // Channeling: happens over a length L_chan (potentially less if dechanneling)
    //             This equates to an opening angle to the point P (center of miscut) tP
    //             The angle xp at the start of channeling is tP/2 + miscut
    //             The angle xp at the end of channeling is tP + miscut

    double xp_rel = xp - cry_miscut;
    double ymin = -ydim/2.;
    double ymax =  ydim/2.;

    // FIRST CASE: p don't interact with crystal
    if (y < ymin || y > ymax || x > xdim) {
        Drift_single_particle_4d(part, length);
        iProc = proc_out;
        result[0] = pc;
        result[1] = iProc;
        return result;

    } else if (x < cry_alayer || y-ymin < cry_alayer || ymax-y < cry_alayer) {
    // SECOND CASE: p hits the amorphous layer
    // TODO: NOT SUPPORTED
        LocalParticle_kill_particle(part, XC_ERR_NOT_IMPLEMENTED);
        return result;
//         double x0    = x;
//         double y0    = y;
//         double a_eq  = 1. + pow(xp,2.);
//         double b_eq  = (2.*x)*xp - (2.*xp)*bend_r;
//         double c_eq  = pow(x,2.) - (2.*x)*bend_r;
//         double delta = pow(b_eq,2.) - (4.*a_eq)*c_eq;
//         s = (-b_eq+sqrt(delta))/(2.*a_eq);
//         if (s >= length) {
//             s = length;
//         }
//         x   =  xp*s + x0;
//         double len_xs = sqrt(pow((x-x0),2.) + pow(s,2.));
//         double len_ys;
//         if (yp >= 0 && y + yp*s <= ymax) {
//             len_ys = yp*len_xs;
//         } else if (yp < 0 && y + yp*s >= ymin) {
//             len_ys = yp*len_xs;
//         } else {
//             s      = (ymax-y)/yp;
//             len_ys = sqrt(pow((ymax-y),2.) + pow(s,2.));
//             x   = x0 + xp*s;
//             len_xs = sqrt(pow((x-x0),2.) + pow(s,2.));
//         }
        
//         double am_len = sqrt(pow(len_xs,2.) + pow(len_ys,2.));
//         s     = s/2;
//         x  = x0 + xp*s;
//         y     = y0 + yp*s;
//         iProc = proc_AM;

//         energy_loss = calcionloss_cry(part, length, properties);

//         double* result_am = moveam(rng, part, am_len, energy_loss, pc, material);
//         pc = result_am[0];
//         iProc = result_am[1];
//         free(result_am);

//         x = x + xp*(length-s);
//         y = y + yp*(length-s);

//         result[0] = x;
//         result[1] = xp;
//         result[2] = y;
//         result[3] = yp;
//         result[4] = pc;
//         result[5] = iProc;
//         return result;

    } else if (x > xdim-cry_alayer && x < xdim) {
    // TODO: NOT SUPPORTED
        LocalParticle_kill_particle(part, XC_ERR_NOT_IMPLEMENTED);
        return result;
//         iProc = proc_AM;
        
//         energy_loss = calcionloss_cry(part, length, properties);

//         double* result_am = moveam(rng, part, length, energy_loss, pc, material);
//         pc = result_am[0];
//         iProc = result_am[1];
//         free(result_am);

//         result[0] = x;
//         result[1] = xp;
//         result[2] = y;
//         result[3] = yp;
//         result[4] = pc;
//         result[5] = iProc;
//         return result;
    }

    //THIRD CASE: the p interacts with the crystal.
    //Define typical angles/probabilities for orientation 110
    double xpcrit0 = sqrt((2.0e-9*eum)/pc);   //Critical angle (rad) for straight crystals
    double Rcrit   = (pc/(2.0e-6*eum))*ai; //Critical curvature radius [m]

    //If R>Rcritical=>no channeling is possible (ratio<1)
    double ratio  = bend_r/Rcrit;
    double xpcrit = (xpcrit0*(bend_r-Rcrit))/bend_r; //Critical angle for curved crystal

    double Ang_rms;
    double Ang_avr;
    double Vcapt;
    if (ratio <= 1.) { //no possibile channeling
        Ang_rms = ((c_v1*0.42)*xpcrit0)*sin(1.4*ratio); //RMS scattering
        Ang_avr = ((c_v2*xpcrit0)*5.0e-2)*ratio;                         //Average angle reflection
        Vcapt   = 0.;                                                //Probability of VC
    } else if (ratio <= 3.) { //Strongly bent crystal
        Ang_rms = ((c_v1*0.42)*xpcrit0)*sin(0.4713*ratio + 0.85); //RMS scattering
        Ang_avr = (c_v2*xpcrit0)*(0.1972*ratio - 0.1472);                  //Average angle reflection
        Vcapt   = 7.0e-4*(ratio - 0.7)/pow(pc,2.0e-1);                           //Correction by sasha drozdin/armen
        //K=0.0007 is taken based on simulations using CATCH.f (V.Biryukov)
    } else { //Rcry >> Rcrit
        Ang_rms = (c_v1*xpcrit0)*(1./ratio);                //RMS scattering
        Ang_avr = (c_v2*xpcrit0)*(1. - 1.6667/ratio); //Average angle for VR
        Vcapt   = 7.0e-4*(ratio - 0.7)/pow(pc,2.0e-1); //Probability for VC correction by sasha drozdin/armen
        //K=0.0007 is taken based on simulations using CATCH.f (V.Biryukov)
    }
    if (cry_orient == 2) {
        Ang_avr = Ang_avr*0.93;
        Ang_rms = Ang_rms*1.05;
        xpcrit  = xpcrit*0.98;
    }

    if (fabs(xp_rel) < xpcrit) {
    // Channeling is possible but need to check

        double alpha  = xp_rel/xpcrit;
        double Chann  = sqrt(0.9*(1 - pow(alpha,2.)))*sqrt(1.-(1./ratio)); //Saturation at 95%
        double N_atom = 1.0e-1;

        //if they can channel: 2 options
        if (RandomUniform_generate(part) <= Chann) {
        //option 1:channeling

            double TLdech1 = (const_dech*pc)*pow((1.-1./ratio),2.); //Updated calculate typical dech. length(m)

            if(RandomUniform_generate(part) <= N_atom) {
                TLdech1 = ((const_dech/2.0e2)*pc)*pow((1.-1./ratio),2.);  //Updated dechanneling length (m)      
            }
            double Dechan = RandomExponential_generate(part); //Probability of dechanneling
            double Ldech  = TLdech1*Dechan;   //Actual dechan. length

            //careful: the dechanneling lentgh is along the trajectory
            //of the particle -not along the longitudinal coordinate...
            if (Ldech < L_chan) {
            // We are dechanneling: channeling until Ldech, then amorphous (CH -> MCS -> ..)
                iProc = proc_DC;
                double channeled_length = Channel_single_particle_4d(part, Ldech, Ldech/r, cry_miscut, xpcrit);
                energy_loss = calcionloss_cry(part, channeled_length, properties);
                pc = pc - 0.5*energy_loss*channeled_length; //Energy loss to ionization while in CH [GeV]  // TODO: why 0.5 ?

                // Remaining part is amorphous
                Drift_single_particle_4d(part, 0.5*(length-channeled_length));
                energy_loss = calcionloss_cry(part, length-channeled_length, properties);
                double* result_am = moveam(rng, part, length-channeled_length, energy_loss, pc, material, iProc);

                pc = result_am[0];
                iProc = result_am[1];
                free(result_am);

                Drift_single_particle_4d(part, 0.5*(length-channeled_length));

                CollimatorImpactsData_set_interaction(record, record_index, part, 0, channeled_length, XC_CRYSTAL_CHANNELING);
                CollimatorImpactsData_set_interaction(record, record_index, part, channeled_length, channeled_length,
                                                      XC_CRYSTAL_DECHANNELING);
                CollimatorImpactsData_set_interaction(record, record_index, part, channeled_length, length,
                                                      XC_CRYSTAL_AMORPHOUS);

            } else {
            // We are channeling (never dechannel spontaneously)   (CH  or  CH  -> PP/..)
                iProc = proc_CH;
                double xpin = xp;
                double ypin = yp;

                //check if a nuclear interaction happen while in CH
                double* result_ch = movech(rng, part, L_chan, pc, bend_r, Rcrit, material, iProc);
                pc = result_ch[0];
                iProc = result_ch[1];
                free(result_ch);

                if (iProc != proc_CH) {
                    //if an nuclear interaction happened, move until the middle with initial xp,yp:
                    //propagate until the "crystal exit" with the new xp,yp accordingly with the rest
                    //of the code in "thin lens approx"
                    LocalParticle_add_to_x(part, (0.5*L_chan)*xpin);
                    LocalParticle_add_to_y(part, (0.5*L_chan)*ypin);
                    Drift_single_particle_4d(part, 0.5*L_chan);

                    energy_loss = calcionloss_cry(part, length, properties);
                    pc = pc - energy_loss*length; //energy loss to ionization [GeV]

                    CollimatorImpactsData_set_interaction(record, record_index, part, 0, 0.5*L_chan, XC_CRYSTAL_CHANNELING);
                    CollimatorImpactsData_set_interaction(record, record_index, part, 0.5*L_chan, L_chan,
                                                          XC_CRYSTAL_AMORPHOUS);

                } else {
                    // Channel
                    double channeled_length = Channel_single_particle_4d(part, L_chan, tP, cry_miscut, xpcrit);
                    energy_loss = calcionloss_cry(part, channeled_length, properties);
                    pc = pc - 0.5*energy_loss*channeled_length; //energy loss to ionization [GeV]  // TODO: why 0.5 ?
                    // TODO: the px at exit of channeling should be calculated from xp with updated energy

                    // Drift remaining length
                    Drift_single_particle_4d(part, length-channeled_length);

                    CollimatorImpactsData_set_interaction(record, record_index, part, 0, channeled_length,
                                                          XC_CRYSTAL_CHANNELING);
                } 
            }

        } else { //Option 2: VR
            //good for channeling but don't channel (1-2)
            iProc = proc_VR;
            LocalParticle_add_to_px(part, 0.45*(xp_rel/xpcrit + 1)*Ang_avr/rpp);
            Drift_single_particle_4d(part, 0.5*length);

            energy_loss = calcionloss_cry(part, length, properties);

            double* result_am = moveam(rng, part, length, energy_loss, pc, material, iProc);

            pc = result_am[0];
            iProc = result_am[1];
            free(result_am);

            Drift_single_particle_4d(part, 0.5*length);

            CollimatorImpactsData_set_interaction(record, record_index, part, 0, 0.5*L_chan, XC_CRYSTAL_DRIFT);
            CollimatorImpactsData_set_interaction(record, record_index, part, 0.5*L_chan, 0.5*L_chan,
                                                  XC_CRYSTAL_VOLUME_REFLECTION);
            CollimatorImpactsData_set_interaction(record, record_index, part, 0.5*L_chan, L_chan,
                                                  XC_CRYSTAL_AMORPHOUS);
        }

    } else { //case 3-2: no good for channeling. check if the can VR
        double Lrefl = xp_rel*r; //Distance of refl. point [m]
        double Srefl = sin(xp_rel/2. + cry_miscut)*Lrefl;

        if (Lrefl > 0 && Lrefl < L_chan) { //VR point inside

        //2 options: volume capture and volume reflection

            if (RandomUniform_generate(part) > Vcapt || zn == 0) { //Option 1: VR
                iProc = proc_VR;

                Drift_single_particle_4d(part, Srefl);

                double Dxp = Ang_avr;
                LocalParticle_add_to_px(part, Dxp/rpp + Ang_rms*RandomNormal_generate(part)/rpp);
                Drift_single_particle_4d(part, 0.5*(length - Srefl));

                energy_loss = calcionloss_cry(part, length-Srefl, properties);

                double* result_am = moveam(rng, part, length-Srefl, energy_loss, pc, material, iProc);

                pc = result_am[0];
                iProc = result_am[1];
                free(result_am);

                Drift_single_particle_4d(part, 0.5*(length - Srefl));

            } else { //Option 2: VC
                Drift_single_particle_4d(part, Srefl);

                double TLdech2 = (const_dech/1.0e1)*pc*pow((1-1/ratio),2.) ;         //Updated typical dechanneling length(m)
                double Ldech   = TLdech2 * pow((sqrt(1.0e-2 + RandomExponential_generate(part)) - 1.0e-1),2.); //Updated DC length
                double tdech   = Ldech/bend_r;
                double Sdech   = Ldech*cos(xp + 0.5*tdech);

                if (Ldech < length-Lrefl) {
                    iProc = proc_DC;
                    double Dxp = Ldech/bend_r + (0.5*RandomNormal_generate(part))*xpcrit;
                    LocalParticle_add_to_x(part, Ldech*(sin(0.5*Dxp+xp))); //Trajectory at channeling exit
                    LocalParticle_add_to_y(part, Sdech*yp);
                    LocalParticle_set_px(part, Dxp/rpp);
                    double Red_S = (length - Srefl) - Sdech;
                    Drift_single_particle_4d(part, 0.5*Red_S);

                    energy_loss = calcionloss_cry(part, Srefl, properties);

                    pc = pc - energy_loss*Srefl; //"added" energy loss before capture

                    energy_loss = calcionloss_cry(part, Sdech, properties);
                    pc = pc - (0.5*energy_loss)*Sdech; //"added" energy loss while captured

                    energy_loss = calcionloss_cry(part, Red_S, properties);

                    double* result_am = moveam(rng, part, Red_S, energy_loss, pc, material, iProc);

                    pc = result_am[0];
                    iProc = result_am[1];
                    free(result_am);

                    Drift_single_particle_4d(part, 0.5*Red_S);

                } else {
                    iProc   = proc_VC;
                    double Rlength = length - Lrefl;
                    double tchan   = Rlength/bend_r;
                    double Red_S   = Rlength*cos(xp + 0.5*tchan);

                    energy_loss = calcionloss_cry(part, Lrefl, properties);
                    pc   = pc - energy_loss*Lrefl; //"added" energy loss before capture
                    double xpin = xp;
                    double ypin = yp;

                    //Check if a nuclear interaction happen while in ch
                    double* result_ch = movech(rng, part, Rlength, pc, bend_r, Rcrit, material, iProc);
                    pc = result_ch[0];
                    iProc = result_ch[1];
                    free(result_ch);
                                    
                    if (iProc != proc_VC) {
                        //if an nuclear interaction happened, move until the middle with initial xp,yp: propagate until
                        //the "crystal exit" with the new xp,yp aciordingly with the rest of the code in "thin lens approx"
                        LocalParticle_add_to_x(part, (0.5*Rlength)*xpin);
                        LocalParticle_add_to_y(part, (0.5*Rlength)*ypin);
                        Drift_single_particle_4d(part, 0.5*Rlength);

                        energy_loss = calcionloss_cry(part, Rlength, properties);
                        pc = pc - energy_loss*Rlength;

                    } else {
                        double Dxp = (length-Lrefl)/bend_r;
                        LocalParticle_add_to_x(part, sin(0.5*Dxp+xp)*Rlength); //Trajectory at channeling exit
                        LocalParticle_add_to_y(part, Red_S*yp);
                        LocalParticle_set_px(part, tdefl/rpp + (0.5*RandomNormal_generate(part))*xpcrit/rpp); //[mrad]

                        energy_loss = calcionloss_cry(part, Rlength, properties);
                        pc = pc - (0.5*energy_loss)*Rlength;  //"added" energy loss once captured

                    }
                }
            }

        } else {
            //Case 3-3: move in amorphous substance (big input angles)
            //Modified for transition vram daniele
            if (xp_rel > tdefl-cry_miscut + 2*xpcrit || xp_rel < -xpcrit) {
                iProc = proc_AM;
                Drift_single_particle_4d(part, 0.5*length);
                if(zn > 0) {
                    energy_loss = calcionloss_cry(part, length, properties);

                    double* result_am = moveam(rng, part, length, energy_loss, pc, material, iProc);
                    pc = result_am[0];
                    iProc = result_am[1];
                    free(result_am);
                }

                Drift_single_particle_4d(part, 0.5*length);

            } else {
                double Pvr = (xp_rel-(tdefl-cry_miscut))/(2.*xpcrit);
                if(RandomUniform_generate(part) > Pvr) {
                    iProc = proc_TRVR;
                    Drift_single_particle_4d(part, Srefl);

                    double Dxp = (((-3.*Ang_rms)*xp_rel)/(2.*xpcrit) + Ang_avr) + ((3.*Ang_rms)*(tdefl-cry_miscut))/(2.*xpcrit);
                    LocalParticle_add_to_px(part, Dxp/rpp);
                    Drift_single_particle_4d(part, 0.5*(length-Srefl));

                    energy_loss = calcionloss_cry(part, length-Srefl, properties);

                    double* result_am = moveam(rng, part, length-Srefl, energy_loss, pc, material, iProc);
                    pc = result_am[0];
                    iProc = result_am[1];
                    free(result_am);

                    Drift_single_particle_4d(part, 0.5*(length - Srefl));

                } else {
                    iProc = proc_TRAM;
                    Drift_single_particle_4d(part, Srefl);
                    double Dxp = ((((-1.*(13.6/pc))*sqrt(length/dlri))*1.0e-3)*xp_rel)/(2.*xpcrit) + (((13.6/pc)*sqrt(length/dlri))*1.0e-3)*(1.+(tdefl-cry_miscut)/(2.*xpcrit));
                    LocalParticle_add_to_px(part, Dxp/rpp);
                    Drift_single_particle_4d(part, 0.5*(length-Srefl));

                    energy_loss = calcionloss_cry(part, length-Srefl, properties);
                    double* result_am = moveam(rng, part, length-Srefl, energy_loss, pc, material, iProc);
                    pc = result_am[0];
                    iProc = result_am[1];
                    free(result_am);

                    Drift_single_particle_4d(part, 0.5*(length - Srefl));
                }
            }
        }            
    }

    result[0] = pc;
    result[1] = iProc;
    return result;
}


/*gpufun*/
double* crystal(RandomRutherfordData rng, LocalParticle* part, double p,
                double length, CrystalMaterialData material, double cry_tilt, double bend_r, double bend_ang,
                double cry_alayer, double xdim, double ydim, double cry_orient, double cry_miscut,
                CollimatorImpactsData record, RecordIndex record_index) {

    double* crystal_result = (double*)malloc(3 * sizeof(double));
    double iProc = proc_out;

    double const cry_cBend  = cos(bend_ang);
    double const cry_sBend  = sin(bend_ang);
    double const cry_cpTilt = cos(cry_tilt);
    double const cry_spTilt = sin(cry_tilt);

    // Move origin of x to inner front corner (transformation 4 in Figure 3.3 of thesis Valentina Previtali)
    double shift = 0;
    if (cry_tilt < 0) {
        shift = bend_r*(1 - cry_cpTilt);
        if (cry_tilt < -bend_ang) {
            shift = bend_r*(cry_cpTilt - cos(bend_ang - cry_tilt));
        }
        LocalParticle_add_to_x(part, -shift);
    } 

    // Rotate tilt (transformation 5 in Figure 3.3 of thesis Valentina Previtali)
    double s = YRotation_single_particle_rotate_only(part, 0., cry_tilt);

    // 3rd transformation: drift to the new coordinate s=0
    Drift_single_particle_4d(part, -s);

    // Check that particle hit the crystal
    double x = LocalParticle_get_x(part);
    double rpp_in  = LocalParticle_get_rpp(part);
    double xp = LocalParticle_get_px(part)*rpp_in;
    if (x >= 0. && x < xdim) {
        // MISCUT first step: P coordinates (center of curvature of crystalline planes)
        double s_P = (bend_r-xdim)*sin(-cry_miscut);
        double x_P = xdim + (bend_r-xdim)*cos(-cry_miscut);
        double* result = interact(rng, part, p, length, s_P, x_P, material, cry_tilt,
                                  bend_r, cry_alayer, xdim, ydim, cry_orient, cry_miscut,
                                  bend_ang, cry_cBend, cry_sBend, cry_cpTilt, cry_spTilt, cry_cpTilt,
                                  -cry_spTilt, iProc, record, record_index);

        p = result[0];
        iProc = result[1];
        free(result);

        s = bend_r*cry_sBend;

    } else {
        double xp_tangent=0;
        if (x < 0) { // Crystal can be hit from below
            xp_tangent = sqrt((-(2.*x)*bend_r + pow(x,2.))/(pow(bend_r,2.)));
        } else {             // Crystal can be hit from above
            xp_tangent = asin((bend_r*(1. - cry_cBend) - x)/sqrt(((2.*bend_r)*(bend_r - x))*(1 - cry_cBend) + pow(x,2.)));
        }
        // If the hit is below, the angle must be greater or equal than the tangent,
        // or if the hit is above, the angle must be smaller or equal than the tangent
        if ((x < 0. && xp >= xp_tangent) || (x >= 0. && xp <= xp_tangent)) {

            // If it hits the crystal, calculate in which point and apply the transformation and drift to that point
            double a_eq  = 1 + pow(xp,2.);
            double b_eq  = (2.*xp)*(x - bend_r);
            double c_eq  = -(2.*x)*bend_r + pow(x,2.);
            double delta = pow(b_eq,2.) - 4*(a_eq*c_eq);
            double s_int = (-b_eq - sqrt(delta))/(2*a_eq);

            // MISCUT first step: P coordinates (center of curvature of crystalline planes)
            double s_P_tmp = (bend_r-xdim)*sin(-cry_miscut);
            double x_P_tmp = xdim + (bend_r-xdim)*cos(-cry_miscut);

            if (s_int < bend_r*cry_sBend) {
                // Transform to a new reference system: shift and rotate
                double tilt_int = s_int/bend_r;
                double x_int  = xp*s_int + x;
                LocalParticle_add_to_y(part, LocalParticle_get_py(part)*rpp_in*s_int);
                LocalParticle_set_x(part, 0.);
                LocalParticle_add_to_px(part, -tilt_int/rpp_in);

                // MISCUT first step (bis): transform P in new reference system
                // Translation
                s_P_tmp = s_P_tmp - s_int;
                x_P_tmp = x_P_tmp - x_int;
                // Rotation
                double s_P = s_P_tmp*cos(tilt_int) + x_P_tmp*sin(tilt_int);
                double x_P = -s_P_tmp*sin(tilt_int) + x_P_tmp*cos(tilt_int);

                double* result = interact(rng, part, p, length-(tilt_int*bend_r), s_P, x_P,
                                          material, cry_tilt, bend_r, cry_alayer, xdim, ydim, cry_orient, 
                                          cry_miscut, bend_ang, cry_cBend, cry_sBend, cry_cpTilt, cry_spTilt, cry_cpTilt,
                                          -cry_spTilt, iProc, record, record_index);

                p = result[0];
                iProc = result[1];
                free(result);

                s = bend_r*sin(bend_ang - tilt_int);

                // un-rotate
                s = YRotation_single_particle_rotate_only(part, s, -tilt_int);

                // 2nd: shift back the 2 axis
                LocalParticle_add_to_x(part, x_int);
                s = s + s_int;
            } else {
                // Drift
                s = bend_r*sin(length/bend_r);
                Drift_single_particle_4d(part, s);
            }
        } else {
            // Drift
            s = bend_r*sin(length/bend_r);
            Drift_single_particle_4d(part, s);
        }
    }

    // transform back from the crystal to the collimator reference system
    // 1st: un-rotate the coordinates
    s = YRotation_single_particle_rotate_only(part, s, -cry_tilt);

    // 2nd: shift back the reference frame
    if (cry_tilt < 0) {
        LocalParticle_add_to_px(part, shift);
    }

    // 3rd: shift to new S=Length position
    Drift_single_particle_4d(part, length-s);

    int is_hit = 0;
    int is_abs = 0;
    if (iProc != proc_out) {
        is_hit = 1;
    }
    if (iProc == proc_absorbed || iProc == proc_ch_absorbed) {
        is_abs = 1;
    }

    crystal_result[0] = is_hit;
    crystal_result[1] = is_abs;
    crystal_result[2] = p;
    return crystal_result;
}

#endif /* XCOLL_EVEREST_CRYSTAL_INTERACT_H */
