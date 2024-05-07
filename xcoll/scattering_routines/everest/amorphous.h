// copyright ############################### #
// This file is part of the Xcoll Package.   #
// Copyright (c) CERN, 2023.                 #
// ######################################### #

#ifndef XCOLL_EVEREST_AMORPHOUS_H
#define XCOLL_EVEREST_AMORPHOUS_H
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

// #define XCOLL_TRANSITION

/*gpufun*/
void volume_reflection(EverestData restrict everest, LocalParticle* part, int8_t transition) {

    InteractionRecordData record = everest->coll->record;
    RecordIndex record_index     = everest->coll->record_index;
    int8_t sc = everest->coll->record_scatterings;
    int64_t i_slot;


    double Ang_avr = everest->Ang_avr;
    double Ang_rms = everest->Ang_rms;

    if (transition == XC_VOLUME_REFLECTION_TRANS_CH){
        // We are in transition from CH to VR
        double xp = LocalParticle_get_xp(part);
        // TODO: we believe the original 0.45 comes from the 0.9 saturation factor, so we changed it to 0.5
        Ang_avr *= 0.5*((xp - everest->t_I)/everest->t_c + 1);
        Ang_rms  = 0;   // TODO: why does transition to CH not use any spread?
        if (sc) i_slot = InteractionRecordData_log(record, record_index, part, XC_VOLUME_REFLECTION_TRANS_CH);

    } else if (transition == XC_VOLUME_REFLECTION_TRANS_MCS){
        // We are in transition from VR to MCS
//         double t_c = everest->t_c;
//         double t_P = everest->t_P;
//         double xp_rel = LocalParticle_get_xp(part) - everest->t_I;
//         // TODO: where does 3 come from
//         Ang_rms *= -3.*(xp_rel-t_P)/(2.*t_c); // TODO: why no random number?
        Ang_rms *= RandomNormal_generate(part);
        if (sc) i_slot = InteractionRecordData_log(record, record_index, part, XC_VOLUME_REFLECTION_TRANS_MCS);

    } else {
        Ang_rms *= RandomNormal_generate(part);
        if (sc) i_slot = InteractionRecordData_log(record, record_index, part, XC_VOLUME_REFLECTION);
    }

    double t_VR = Ang_avr + Ang_rms;
    LocalParticle_add_to_xp(part, t_VR);
    if (sc) InteractionRecordData_log_child(record, i_slot, part, 0);
}


// Amorphous transport is just Multiple Coulomb scattering
/*gpufun*/
double amorphous_transport(EverestData restrict everest, LocalParticle* part, double pc, double length, int8_t transition) {

    InteractionRecordData record = everest->coll->record;
    RecordIndex record_index     = everest->coll->record_index;
    int8_t sc = everest->coll->record_scatterings;
    int64_t i_slot;

    // Accumulated effect of mcs on the angles (with initial energy)
    // TODO: pc is energy
    // TODO: missing factor (1 + 0.038*log( L / dlri) )  ( *Z if not protons)
    double dya = (13.6/pc)*sqrt(length/everest->coll->dlri)*1.0e-3; // RMS of coloumb scattering MCS (rad)
    double kxmcs, kymcs;

    if (transition == XC_MULTIPLE_COULOMB_TRANS_VR){
        // Transition MCS
        if (sc) i_slot = InteractionRecordData_log(record, record_index, part, XC_MULTIPLE_COULOMB_TRANS_VR);
//         double xp_rel = LocalParticle_get_xp(part) - everest->t_I;
//         double t_P = everest->t_P;
//         double t_c = everest->t_c;
//         dya *= 1 - (xp_rel-t_P)/(2.*t_c);
    } else {
        // Normal MCS
        if (sc) i_slot = InteractionRecordData_log(record, record_index, part, XC_MULTIPLE_COULOMB_SCATTERING);
    }
    kxmcs = dya*RandomNormal_generate(part);
    kymcs = dya*RandomNormal_generate(part);

    // Transport of Multiple Coulomb Scattering
    Drift_single_particle_4d(part, length);

    // Energy lost because of ionisation process[GeV]
    double energy_loss = calcionloss(everest, part, length);
    pc  = pc - energy_loss*length;

    // Store new angles
    LocalParticle_add_to_xp_yp(part, kxmcs, kymcs);

    // Finally log particle at end of mcs
    if (sc) InteractionRecordData_log_child(record, i_slot, part, length);

    return pc;
}

double Channel(EverestData restrict everest, LocalParticle* part, double pc, double length);

/*gpufun*/
double Amorphous(EverestData restrict everest, LocalParticle* part, double pc, double length) {

    if (LocalParticle_get_state(part) < 1){
        // Do nothing if already absorbed
        return pc;
    }

    InteractionRecordData record = everest->coll->record;
    RecordIndex record_index     = everest->coll->record_index;
    int8_t sc = everest->coll->record_scatterings;

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
    // Calculate extra length to transition region between VR and AM
    double length_VR_trans = everest->r * sin(xp - everest->t_I - 2.*everest->t_c) * cos(xp - 2.*everest->t_c);
    if (length_VR_trans < 0) {length_VR_trans = 1e10;}

    // ------------------------------------------------------------------------
    // Compare the 3 lengths: the first one encountered is what will be applied
    // ------------------------------------------------------------------------
    if (length_VI <= fmin(length_nucl, length_exit)){
        // MCS to volume interaction
        pc = amorphous_transport(everest, part, pc, length_VI, 0);
#ifdef XCOLL_REFINE_ENERGY
        calculate_VI_parameters(everest, part, pc);
#endif
        // Are we reflecting or captured?
        if (RandomUniform_generate(part) > everest->Vcapt) {
            // Volume Reflection
            volume_reflection(everest, part, 0);
            // We call the main Amorphous function for the leftover
            pc = Amorphous(everest, part, pc, length - length_VI);

        } else {
            // Volume Capture
            if (sc) InteractionRecordData_log(record, record_index, part, XC_VOLUME_CAPTURE);
            // We call the main Channel function for the leftover
            pc = Channel(everest, part, pc, length - length_VI);
        }

#ifdef XCOLL_TRANSITION
    } else if (length_VR_trans <= fmin(length_nucl, length_exit)){
        // Transition region between VR and AM for t_P < xp - tI < t_P + 2t_c
#ifdef XCOLL_REFINE_ENERGY
        calculate_critical_angle(everest, part, pc);
#endif
        double xp_rel = xp - everest->t_I;
        double t_P = everest->t_P;
        double t_c = everest->t_c;
        double prob_MCS = (xp_rel - t_P) / (2*t_c);
        if (RandomUniform_generate(part) > prob_MCS){
            // We are on the VR side
            pc = amorphous_transport(everest, part, pc, length_VR_trans, 0);
            volume_reflection(everest, part, XC_VOLUME_REFLECTION_TRANS_MCS);
            pc = amorphous_transport(everest, part, pc, length - length_VR_trans, 0);
        } else {
            // We are on the AM side
            pc = amorphous_transport(everest, part, pc, length, XC_MULTIPLE_COULOMB_TRANS_VR);
        }
#endif

    } else if (length_nucl < length_exit) {
        // MCS to nuclear interaction
        pc = amorphous_transport(everest, part, pc, length_nucl, 0);
        // interact
        pc = nuclear_interaction(everest, part, pc);
        if (LocalParticle_get_state(part) == XC_LOST_ON_EVEREST_COLL){
            LocalParticle_set_state(part, XC_LOST_ON_EVEREST_CRYSTAL);
        } else {
            // We call the main Amorphous function for the leftover
            pc = Amorphous(everest, part, pc, length - length_nucl);
        }

    } else {
        // Exit crystal
        // MCS to exit point
        pc = amorphous_transport(everest, part, pc, length_exit, 0);
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
            InteractionRecordData_log(record, record_index, part, XC_EXIT_JAW);
            Drift_single_particle_4d(part, length - length_exit);
        }
    }

    return pc;
}

#endif /* XCOLL_EVEREST_AMORPHOUS_H */
