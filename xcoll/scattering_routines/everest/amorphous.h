// copyright ############################### #
// This file is part of the Xcoll Package.   #
// Copyright (c) CERN, 2023.                 #
// ######################################### #

#ifndef XCOLL_EVEREST_AMORPHOUS_H
#define XCOLL_EVEREST_AMORPHOUS_H
#include <math.h>
#include <stdio.h>
#include <stdlib.h>


/*gpufun*/
double nuclear_interaction(EverestData restrict everest, LocalParticle* part, double pc) {

    CollimatorImpactsData record = everest->coll->record;
    RecordIndex record_index     = everest->coll->record_index;

#ifdef XCOLL_REFINE_ENERGY
    calculate_scattering(everest, pc);
#endif

    //Choose nuclear interaction
    double aran = RandomUniform_generate(part);
    int ichoix = 1;

    while (aran > everest->cprob[ichoix]) {
        ichoix += 1;
    }

    //Do the interaction
    double teta = 0 ; //default value to cover ichoix=1
    int64_t i_slot;
    if (ichoix==1) {
        i_slot = CollimatorImpactsData_log(record, record_index, part, XC_ABSORBED);
        LocalParticle_set_state(part, XC_LOST_ON_EVEREST_CRYSTAL);

    } else {
        if (ichoix==2) { //p-n elastic
            i_slot = CollimatorImpactsData_log(record, record_index, part, XC_PN_ELASTIC);
            teta  = sqrt(RandomExponential_generate(part)/everest->bn)/pc;

        } else if (ichoix==3) { //p-p elastic
            i_slot = CollimatorImpactsData_log(record, record_index, part, XC_PP_ELASTIC);
            teta  = sqrt(RandomExponential_generate(part)/everest->bpp)/pc;

        } else if (ichoix==4) { //Single diffractive
            i_slot = CollimatorImpactsData_log(record, record_index, part, XC_SINGLE_DIFFRACTIVE);
            double xm2 = exp(RandomUniform_generate(part)*everest->xln15s);
            double bsd = 0.0;
            if (xm2 < 2.) {
                bsd = 2*everest->bpp;
            } else if (xm2 >= 2. && xm2 <= 5.) {
                bsd = ((106.0 - 17.0*xm2)*everest->bpp)/36.0;
            } else if (xm2 > 5.) {
                bsd = (7*everest->bpp)/12.0;
            }
            teta = sqrt(RandomExponential_generate(part)/bsd)/pc;
            pc = pc*(1 - xm2/everest->ecmsq);

        } else { //(ichoix==5)
            i_slot = CollimatorImpactsData_log(record, record_index, part, XC_RUTHERFORD);
            teta  = sqrt(RandomRutherford_generate(everest->coll->rng, part))/pc;
        }

        double tx = teta*RandomNormal_generate(part);
        double tz = teta*RandomNormal_generate(part);

        //Change the angles
        LocalParticle_add_to_xp_yp(part, tx, tz);

        CollimatorImpactsData_log_child(record, i_slot, part, 0);
    }

    return pc;
}


/*gpufun*/
void volume_reflection(EverestData restrict everest, LocalParticle* part, int8_t transition) {

    CollimatorImpactsData record = everest->coll->record;
    RecordIndex record_index     = everest->coll->record_index;
    int64_t i_slot;

    double Ang_avr = everest->Ang_avr;
    double Ang_rms = everest->Ang_rms;
    if (transition == XC_VOLUME_REFLECTION_TRANS_CH){
        // We are in transition from CH to VR
        double xp = LocalParticle_get_xp(part);
        // TODO: we believe the original 0.45 comes from the 0.9 saturation factor, so we changed it to 0.5
        Ang_avr *= 0.5*((xp - everest->t_I)/everest->t_c + 1);
        Ang_rms  = 0;   // TODO: why does transition to CH not use any spread?
        i_slot = CollimatorImpactsData_log(record, record_index, part, XC_VOLUME_REFLECTION_TRANS_CH);

    } else if (transition == XC_VOLUME_REFLECTION_TRANS_MCS){
        // We are in transition from VR to MCS
        double t_c = everest->t_c;
        double t_P = everest->t_P;
        double xp_rel = LocalParticle_get_xp(part) - everest->t_I;
        // TODO: where does 3 come from
        Ang_rms *= -3.*(xp_rel-t_P)/(2.*t_c); // TODO: why no random number?
        i_slot = CollimatorImpactsData_log(record, record_index, part, XC_VOLUME_REFLECTION_TRANS_MCS);

    } else {
        Ang_rms *= RandomNormal_generate(part);
        i_slot = CollimatorImpactsData_log(record, record_index, part, XC_VOLUME_REFLECTION);

    }
    double t_VR = Ang_avr + Ang_rms;
    LocalParticle_add_to_xp(part, t_VR);
    CollimatorImpactsData_log_child(record, i_slot, part, 0);
}


// Amorphous transport is just Multiple Coulomb scattering
/*gpufun*/
double amorphous_transport(EverestData restrict everest, LocalParticle* part, double pc, double length) {

    CollimatorImpactsData record = everest->coll->record;
    RecordIndex record_index     = everest->coll->record_index;
    int64_t i_slot;

    // Accumulated effect of mcs on the angles (with initial energy)
    // TODO: pc is energy
    // TODO: missing factor (1 + 0.038*log( L / dlri) )  ( *Z if not protons)
    double dya = (13.6/pc)*sqrt(length/everest->coll->dlri)*1.0e-3; // RMS of coloumb scattering MCS (rad)
    double kxmcs, kymcs;
    
    // Are we in a transition region?
    calculate_opening_angle(everest, part);
#ifdef XCOLL_REFINE_ENERGY
    calculate_critical_angle(everest, part, pc);
#endif
    double t_c = everest->t_c;
    double t_P = everest->t_P;
    double xp_rel = LocalParticle_get_xp(part) - everest->t_I;
    if (xp_rel > 0 && xp_rel <= t_P + 2*t_c){
       // We are in a transition region
        double prob_MCS = (xp_rel - t_P) / (2*t_c);
        if (RandomUniform_generate(part) > prob_MCS){
            // We are on the VR side
            volume_reflection(everest, part, XC_VOLUME_REFLECTION_TRANS_MCS);
            // Now normal MCS
            i_slot = CollimatorImpactsData_log(record, record_index, part, XC_MULTIPLE_COULOMB_SCATTERING);
            kxmcs = dya*RandomNormal_generate(part);
            kymcs = dya*RandomNormal_generate(part);

        } else {
            // We are on the MCS side: modified angles
            i_slot = CollimatorImpactsData_log(record, record_index, part, XC_MULTIPLE_COULOMB_TRANS_VR);
            dya *= 1 - (xp_rel-t_P)/(2.*t_c);
            kxmcs = dya*RandomNormal_generate(part);
            kymcs = dya*RandomNormal_generate(part);
        }

    } else {
        // normal MCS
        i_slot = CollimatorImpactsData_log(record, record_index, part, XC_MULTIPLE_COULOMB_SCATTERING);
        kxmcs = dya*RandomNormal_generate(part);
        kymcs = dya*RandomNormal_generate(part);
    }

    // Transport of Multiple Coulomb Scattering
    Drift_single_particle_4d(part, length);

    // Energy lost because of ionisation process[GeV]
    double energy_loss = calcionloss(everest, part, length);
    pc  = pc - energy_loss*length;

    // Store new angles
    LocalParticle_add_to_xp_yp(part, kxmcs, kymcs);

    // Finally log particle at end of mcs
    CollimatorImpactsData_log_child(record, i_slot, part, length);

    return pc;
}

double Channel(EverestData restrict everest, LocalParticle* part, double pc, double length);

/*gpufun*/
double Amorphous(EverestData restrict everest, LocalParticle* part, double pc, double length) {

    if (LocalParticle_get_state(part) < 1){
        // Do nothing if already absorbed
        return pc;
    }

    CollimatorImpactsData record = everest->coll->record;
    RecordIndex record_index     = everest->coll->record_index;

    calculate_initial_angle(everest, part);

    // -----------------------------------------------
    // Calculate longitudinal position where we go out
    // -----------------------------------------------
    // Assumption: crystal bends away from beam (this does not work when crystal bends towards beam)
    // All s are absolute in the collimator frame
    double R = everest->coll->bend_r;
    double t = everest->coll->bend_ang;
    double d = everest->coll->xdim;
    double s  = LocalParticle_get_s(part);
    double x  = LocalParticle_get_x(part);
    double xp = LocalParticle_get_xp(part);
    // We either exit the crystal at the exit face
    double s1 = (R - x + s*xp) / (xp + cos(t)/sin(t) );
    if (s1 < s) {s1 = 1e10;}
    // or the larger bend
    double dd = (R - x + s*xp) / (1 + xp*xp);
    double s2 = dd*xp + sqrt(R*R / (1 + xp*xp) * dd*dd);
    if (s2 < s) {s2 = 1e10;}
    // or the smaller bend (which potentially has two solutions, in case the particle enters again before the exit face)
    double s3 = dd*xp - sqrt( (R-d)*(R-d) / (1 + xp*xp) * dd*dd);
    if (s3 < s) {s3 = 1e10;}
    // whichever comes first after the current position.
    // This gives us an "exit length" \with respect to the current position
    double exit_point  = fmin(fmin(s1, s2), s3);
    double length_exit = fmin(exit_point - s, length);

    // ----------------------------------------------------
    // Calculate longitudinal length of nuclear interaction
    // ----------------------------------------------------
    double length_nucl = everest->coll->collnt*RandomExponential_generate(part);

    // --------------------------------------------------------------
    // Calculate longitudinal length of volume interaction (VR or VC)
    // --------------------------------------------------------------
    // This happens when the particle is tangential to the crystal planes.
    // The (straight) trajectory until the point of reflection is r sin(xp - t_I),
    // hence the longitudinal length is r sin(xp - miscut) cos xp
    double length_VI = everest->r * sin(xp - everest->t_I) * cos(xp);
    // If length_VI is negative, VI is not possible, so set to a big value
    if (length_VI < 0) {length_VI = 1e10;}

// printf("Amorphous:  L: %f,  L_E: %f,  L_n: %f,  L_VI: %f\n", length, length_exit, length_nucl, length_VI);
    // ------------------------------------------------------------------------
    // Compare the 3 lengths: the first one encountered is what will be applied
    // ------------------------------------------------------------------------
    if (length_exit <= fmin(length_nucl, length_VI)){
        // MCS to exit point
        pc = amorphous_transport(everest, part, pc, length_exit);
        // Now drift the remaining
        // However, if we have exited at s3, and we encounter s4 before s2, we reenter:
        double s4 = dd*xp + sqrt( (R-d)*(R-d) / (1 + xp*xp) * dd*dd);  // second solution for smaller bend
        if (s3 < fmin(s1, s2) && s4 < s2){
            // We drift until re-entry
            Drift_single_particle_4d(part, s4 - exit_point);
            // We call the main Amorphous function for the leftover
            pc = Amorphous(everest, part, pc, length - length_exit - s4 + exit_point);
        } else {
            // Drift leftover out of the crystal
            CollimatorImpactsData_log(record, record_index, part, XC_EXIT_JAW);
            Drift_single_particle_4d(part, length - length_exit);
        }

    } else if (length_nucl < length_VI) {
        // MCS to nuclear interaction
        pc = amorphous_transport(everest, part, pc, length_nucl);
        // interact
        pc = nuclear_interaction(everest, part, pc);
        // We call the main Amorphous function for the leftover
        pc = Amorphous(everest, part, pc, length - length_nucl);

    } else {
        // TODO: we believe we should MCS to the VI point. However, this changes the angles, such that
        // by the time we arrive, the angle is certainly not correct anymore for VR or VC.
        // This is why in the original code the proton is just drifted to the VI point (not even ionloss!).
        // If we want to implement this, we have to solve the L where xp_MCS == t_VI. This is highly non-trivial
        // because t_VI depends on the exact x position, which changes during MCS.
        // Furthermore, the VRAM region implementation is not compatible with MCS to the VI point.
        //
//         // MCS to volume interaction
//         pc = amorphous_transport(everest, part, pc, length_VI);
// #ifdef XCOLL_REFINE_ENERGY
//         calculate_VI_parameters(everest, part, pc);
// #endif
        int64_t i_slot = CollimatorImpactsData_log(record, record_index, part, XC_DRIFT_TO_VI);
        Drift_single_particle_4d(part, length_VI);
        CollimatorImpactsData_log_child(record, i_slot, part, length_VI);

        // Are we reflecting or captured?
        if (RandomUniform_generate(part) > everest->Vcapt) {
            // Volume Reflection
            volume_reflection(everest, part, 0);
            // We call the main Amorphous function for the leftover
            pc = Amorphous(everest, part, pc, length - length_VI);

        } else {
            // Volume Capture
            CollimatorImpactsData_log(record, record_index, part, XC_VOLUME_CAPTURE);
            // We call the main Channel function for the leftover
            pc = Channel(everest, part, pc, length - length_VI);
        }
    }

    return pc;
}

#endif /* XCOLL_EVEREST_AMORPHOUS_H */
