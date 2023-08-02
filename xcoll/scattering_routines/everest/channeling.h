// copyright ############################### #
// This file is part of the Xcoll Package.   #
// Copyright (c) CERN, 2023.                 #
// ######################################### #

#ifndef XCOLL_EVEREST_CHANNELING_H
#define XCOLL_EVEREST_CHANNELING_H
#include <math.h>
#include <stdio.h>
#include <stdlib.h>


// Convention:
//       L: length of curved trajectory
//       s: longitudinal coordinate (in collimator frame)  -  projection of curved trajectory
//       t: opening angle of curve ('t' for theta)
//       r: radius from position to bending centre (including t_in):  L = t*r


/*gpufun*/
double channeling_average_density(EverestData restrict coll, LocalParticle* part, double pc, double r, double ratio) {

    // Material properties
    double const anuc = coll->anuc;
    double const rho  = coll->rho;
    double const eum  = coll->eum;

    //Rescale the total and inelastic cross-section accordigly to the average density seen
    double rpp = LocalParticle_get_rpp(part);
    double x_i = LocalParticle_get_x(part);
    double xp  = LocalParticle_get_px(part)*rpp;
    int np     = x_i/XC_PLANE_DISTANCE;          //Calculate in which crystalline plane the particle enters
    x_i       -= (np + 0.5)*XC_PLANE_DISTANCE;   //Rescale the incoming x to the middle of the crystalline plane

    double pv   = pow(pc, 2.)/sqrt(pow(pc, 2.) + pow(XC_PROTON_MASS*1.0e-3, 2.))*1.0e9; //Calculate pv=P/E   TODO: this is beta?
    double Ueff = eum*4.*pow(x_i/XC_PLANE_DISTANCE, 2.) + pv*x_i/r; //Calculate effective potential
    double Et   = pv*pow(xp, 2.)/2. + Ueff;       //Calculate transverse energy
    double Ec   = eum*pow(1. - ratio, 2.);         //Calculate critical energy in bent crystals

    //To avoid negative Et
    double xminU = -pow(XC_PLANE_DISTANCE, 2.)*pc*1.0e9/(8.*eum*r);
    double Umin  = fabs(eum*4.*pow(xminU/XC_PLANE_DISTANCE, 2.) + pv*xminU/r);
    Et    = Et + Umin;
    Ec    = Ec + Umin;

    //Calculate min e max of the trajectory between crystalline planes
    double x_min = -0.5*XC_PLANE_DISTANCE*ratio - 0.5*XC_PLANE_DISTANCE*sqrt(Et/Ec);
    double x_max = -0.5*XC_PLANE_DISTANCE*ratio + 0.5*XC_PLANE_DISTANCE*sqrt(Et/Ec);

    //Change ref. frame and go back with 0 on the crystalline plane on the left
    x_min = x_min - XC_PLANE_DISTANCE/2.;
    x_max = x_max - XC_PLANE_DISTANCE/2.;

    //Calculate the "normal density" in m^-3
    double N_am  = rho*XC_AVOGADRO*1.0e6/anuc;

    //Calculate atomic density at min and max of the trajectory oscillation
    // erf returns the error function of complex argument
    double rho_max = erf(x_max/sqrt(2*pow(XC_THERMAL_VIBRATIONS, 2.)));
    rho_max       -= erf((XC_PLANE_DISTANCE-x_max)/sqrt(2*pow(XC_THERMAL_VIBRATIONS, 2.)));
    rho_max       *= N_am*XC_PLANE_DISTANCE/2.;
    double rho_min = erf(x_min/sqrt(2*pow(XC_THERMAL_VIBRATIONS, 2.)));
    rho_min       -= erf((XC_PLANE_DISTANCE-x_min)/sqrt(2*pow(XC_THERMAL_VIBRATIONS, 2.)));
    rho_min       *= N_am*XC_PLANE_DISTANCE/2.;

    //"zero-approximation" of average nuclear density seen along the trajectory
    double avrrho  = (rho_max - rho_min)/(x_max - x_min);
    avrrho  = 2.*avrrho/N_am;

    return avrrho;
}


/*gpufun*/
double* channel_transport(EverestData restrict coll, LocalParticle* part, double pc, double L_chan, double t_P, double t_in) {
    // Channeling: happens over an arc length L_chan (potentially less if dechanneling)
    //             This equates to an opening angle t_P wrt. to the point P (center of miscut if at start of crystal)
    //             The angle xp at the start of channeling (I) is t_P/2 + t_in
    //             The angle xp at the end of channeling (F) is t_P + t_in
    //             In practice: we drift from start to end, but overwrite the angle afterwards
    // TODO: why does channeling only have 50% energy loss?

    double* result = (double*)malloc(2 * sizeof(double));

    CollimatorImpactsData record = coll->record;
    RecordIndex record_index     = coll->record_index;

    // First log particle at start of channeling
    int64_t i_slot = CollimatorImpactsData_log(record, record_index, part, XC_CRYSTAL_CHANNELING);

    // Do channeling.
    // The distance from I to F is the chord length of the angle t_P: d = 2 r sin(t_P/2)
    // Hence the longitudinal distance (the length to be drifted) is the projection of this using the
    // xp at the start of channeling: s = 2 r sin(t_P/2)cos(t_P/2 + t_in)
    //                                  = 2 r sin(t_P/2)cos(t_P/2)cos(t_in) - 2 r sin(t_P/2)sin(t_P/2)sin(t_in)
    //                                  = r sin(t_P)cos(t_in) - 2 r sin(t_P/2)^2 sin(t_in)
    //                                  = L_chan/t_P ( sin(t_P)cos(t_in) - 2 sin(t_P/2)^2 sin(t_in) )
    double drift_length = L_chan/t_P * (
                            sin(t_P)*cos(t_in)
                            - 2.*sin(t_P/2.)*sin(t_P/2.)*sin(t_in)
                          );
    double rpp = LocalParticle_get_rpp(part);
    LocalParticle_set_px(part, (t_P/2. + t_in)/rpp);          // Angle at start of channeling
    Drift_single_particle_4d(part, drift_length);
    double sigma_ran = 0.5*coll->xpcrit;                      // Extra random angle spread at exit
    double ran_angle = RandomNormal_generate(part)*sigma_ran;
    LocalParticle_set_px(part, (t_P + t_in + ran_angle)/rpp); // Angle at end of channeling

    // Apply energy loss along trajectory
    double energy_loss = 0.5*calcionloss(coll, part, L_chan);
    // TODO: LocalParticle_add_to_energy(part, - energy_loss*L_chan*1.e9, change_angle);
    // if change_angle = 0  => LocalParticle_scale_px(part, old_rpp / new_rpp) such that xp remains the same
    // this is done in K2, though, this is no longer correct with exact drifts...?
    pc = pc - energy_loss*L_chan; //energy loss to ionization [GeV]

    // Finally log particle at end of channeling
    CollimatorImpactsData_log_child(record, i_slot, part, drift_length);

    result[0] = drift_length;
    result[1] = pc;
    return result;
}


/*gpufun*/
double* Channel(EverestData restrict coll, LocalParticle* part, double pc, double L_chan, double t_chan, double t_in,
                double remaining_length) {

    double* result = (double*)malloc(3 * sizeof(double));
    double iProc = proc_CH;
    double crit_r = coll->Rcrit;
    double bend_r = coll->bend_r;
    double ratio = crit_r/bend_r;

    CollimatorImpactsData record = coll->record;
    RecordIndex record_index     = coll->record_index;

    // TODO: compiler flag to update at every energy change if high precision
    calculate_ionisation_properties(coll, pc);
    double const_dech = calculate_dechanneling_length(coll, pc);

    // Calculate curved position L_dechan of dechanneling
    double TLdech1 = const_dech*pc*pow((1. - ratio), 2.); //Updated calculate typical dech. length(m)
    double N_atom = 1.0e-1;
    if(RandomUniform_generate(part) <= N_atom) {
        TLdech1 /= 200.;   // Updated dechanneling length (m)      
    }
    double L_dechan = TLdech1*RandomExponential_generate(part);   // Actual dechan. length
    
    // Calculate curved position L_nucl of nuclear interaction
    double avrrho = channeling_average_density(coll, part, pc, bend_r, ratio);
    calculate_scattering(coll, pc, avrrho);
    double collnt = coll->collnt;
    if (avrrho == 0) {
        collnt = 1.e10;  // very large because essentially 1/0
    } else {
        collnt = collnt/avrrho;
    }
    double L_nucl = collnt*RandomExponential_generate(part);

    // Do the channeling
    // We compare 3 lengths: L_chan vs L_dechan vs L_nucl
    if (L_chan <= fmin(L_dechan, L_nucl)){
        // Channel full length
        double* result_chan = channel_transport(coll, part, pc, L_chan, t_chan, t_in);
        remaining_length -= result_chan[0];
        pc                = result_chan[1];
        free(result_chan);

    } else if (L_dechan < L_nucl) {
        // Channel up to L_dechan, then amorphous
        double* result_chan = channel_transport(coll, part, pc, L_dechan, t_chan*L_dechan/L_chan, t_in);
        remaining_length -= result_chan[0];
        pc                = result_chan[1];
        free(result_chan);
        CollimatorImpactsData_log(record, record_index, part, XC_CRYSTAL_DECHANNELING);

        double* result_am = Amorphous(coll, part, pc, remaining_length, 1);
        remaining_length  = result_am[0];
        pc                = result_am[1];

    } else {
        // Channel up to L_nucl, then scatter, then amorphous
        double* result_chan = channel_transport(coll, part, pc, L_nucl, t_chan*L_nucl/L_chan, t_in);
        remaining_length -= result_chan[0];
        pc                = result_chan[1];
        free(result_chan);

        double* result_ni = nuclear_interaction(coll, part, pc, iProc);
        pc    = result_ni[0];
        iProc = result_ni[1];
        free(result_ni);

        double* result_am = Amorphous(coll, part, pc, remaining_length, 0);
        remaining_length  = result_am[0];
        pc                = result_am[1];
    }

    result[0] = remaining_length;
    result[1] = pc;
    result[2] = iProc;
    return result;
}

#endif /* XCOLL_EVEREST_CHANNELING_H */
