// copyright ############################### #
// This file is part of the Xcoll package.   #
// Copyright (c) CERN, 2023.                 #
// ######################################### #

#ifndef XCOLL_EVEREST_AMORPHOUS_H
#define XCOLL_EVEREST_AMORPHOUS_H
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef MIN
#define MIN(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a < _b ? _a : _b; })
#endif

/*gpufun*/
void volume_reflection(EverestData restrict everest, LocalParticle* part, int8_t transition) {

    InteractionRecordData record = everest->coll->record;
    RecordIndex record_index     = everest->coll->record_index;
    int8_t sc = everest->coll->record_scatterings;
    int64_t i_slot = -1;


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
    if (sc) InteractionRecordData_log_child(record, i_slot, part);
}


// Amorphous transport is just Multiple Coulomb scattering
/*gpufun*/
double amorphous_transport(EverestData restrict everest, LocalParticle* part, double pc, double length, int8_t transition) {

    InteractionRecordData record = everest->coll->record;
    RecordIndex record_index     = everest->coll->record_index;
    int8_t sc = everest->coll->record_scatterings;
    int64_t i_slot = -1;

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
    if (sc) InteractionRecordData_log_child(record, i_slot, part);

    return pc;
}

double Channel(EverestData restrict everest, LocalParticle* part, CrystalGeometry restrict cg, double pc, double length);

/*gpufun*/
double Amorphous(EverestData restrict everest, LocalParticle* part, CrystalGeometry restrict cg, double pc, double length, int8_t allow_VI) {

    if (LocalParticle_get_state(part) < 1){
        // Do nothing if already absorbed
        return pc;
    }

    InteractionRecordData record = everest->coll->record;
    RecordIndex record_index     = everest->coll->record_index;
    int8_t sc = everest->coll->record_scatterings;

    calculate_initial_angle(everest, part, cg);

    // -----------------------------------------------
    // Calculate longitudinal position where we go out
    // -----------------------------------------------
    // TODO: rewrite better readable, and use general output function
    // TODO: check for all R and side
    // double s  = LocalParticle_get_s(part);
    // double x  = LocalParticle_get_x(part);
    // double xp = LocalParticle_get_xp(part);
    // double s_exit = get_s_of_crossing_after_s_with_vlimit(double part_x, double part_tan_x, \
    //                             double part_y, double part_tan_y, Segment* segments, \
    //                             int8_t n_segments, double y_min, double y_max, double current_s
    // s = get_s_of_first_crossing_with_vlimit(part_x, part_tan_x, part_y, part_tan_y, cg->segments, 4, -cg->height/2, cg->height/2);
    double R = cg->bending_radius;
    double d = cg->width;
    double sB = cg->s_B;
    double xB = cg->x_B;
    double s  = LocalParticle_get_s(part);
    double x  = LocalParticle_get_x(part);
    double xp = LocalParticle_get_xp(part);
    double cot = cos(cg->bending_angle)/sin(cg->bending_angle);
    double xx = xB - x + s*xp;
    // We either exit the crystal at the exit face
    double s1 = (xx + cot*sB) / (xp + cot);
    if (s1 < s) {s1 = 1e10;}
    // or the larger bend
    double s2 = 1e10;
    double s3 = 1e10;
    double xpxp = 1 + xp*xp;
    double cc = (xp*xx + sB) / xpxp;
    double ee = (xx - sB*xp)*(xx - sB*xp);
    double dd = xpxp * R*R - ee;
    if (dd >= 0){
        s2 = cc + sqrt(dd) / xpxp;
        if (s2 < s || s2 < sB || (xB > 0 && xp*s2 > xx) || (xB < 0 && xp*s2 < xx)){
            s2 = 1e10; // wrong quadrant
        }
        if (dd > 0){
            s3 = cc - sqrt(dd) / xpxp;
            if (s3 < s || s3 < sB || (xB > 0 && xp*s3 > xx) || (xB < 0 && xp*s3 < xx)){
                s3 = 1e10; // wrong quadrant
            }
        }
    }
    // or the smaller bend
    double s4 = 1e10;
    double s5 = 1e10;
    dd = xpxp * (R-d)*(R-d) - ee;
    if (dd >= 0){
        s4 = cc + sqrt(dd) / xpxp;
        if (s4 < s || s4 < sB || (xB > 0 && xp*s4 > xx) || (xB < 0 && xp*s4 < xx)){
            s4 = 1e10; // wrong quadrant
        }
        if (dd > 0){
            s5 = cc - sqrt(dd) / xpxp;
            if (s5 < s || s5 < sB || (xB > 0 && xp*s5 > xx) || (xB < 0 && xp*s5 < xx)){
                s5 = 1e10; // wrong quadrant
            }
        }
    }
    // whichever comes first after the current position.
    // This gives us an "exit length" \with respect to the current position
    double exit_point  = fmin(fmin(fmin(fmin(s1, s2), s3), s4), s5);
    double length_exit = fmin(exit_point - s, length);
    printf("exit length: %f\n", length);

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
    if (length_VI <= fmin(length_nucl, length_exit) && allow_VI == 1){
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
            pc = Amorphous(everest, part, cg, pc, length - length_VI, 0);

        } else {
            // Volume Capture
            if (sc) InteractionRecordData_log(record, record_index, part, XC_VOLUME_CAPTURE);
            // We call the main Channel function for the leftover
            pc = Channel(everest, part, cg, pc, length - length_VI);
        }

#ifdef XCOLL_TRANSITION
    } else if (length_VR_trans <= fmin(length_nucl, length_exit) && allow_VI == 1){
        // Transition region between VR and AM for t_P < xp - tI < t_P + 2t_c
#ifdef XCOLL_REFINE_ENERGY
        calculate_critical_angle(everest, part, cg, pc);
#endif
        double xp_rel = xp - everest->t_I;
        double t_P = everest->t_P;
        double t_c = everest->t_c;
        double prob_MCS = (xp_rel - t_P) / (2*t_c);
        if (RandomUniform_generate(part) > prob_MCS){
            // We are on the VR side
            pc = amorphous_transport(everest, part, pc, length_VR_trans, 0);
            volume_reflection(everest, part, XC_VOLUME_REFLECTION_TRANS_MCS);
            pc = Amorphous(everest, part, cg,  pc, length - length_VR_trans, 0);
        } else {
            // We are on the AM side
            // if (sc) InteractionRecordData_log(record, record_index, part, XC_MULTIPLE_COULOMB_TRANS_VR);
            pc = amorphous_transport(everest, part, pc, length_VR_trans, XC_MULTIPLE_COULOMB_TRANS_VR);
            pc = Amorphous(everest, part, cg, pc, length - length_VR_trans, 0);
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
            pc = Amorphous(everest, part, cg, pc, length - length_nucl, 1);
        }

    } else {
        // Exit crystal
        // MCS to exit point
        pc = amorphous_transport(everest, part, pc, length_exit, 0);
        // However, if we have exited at s3, and we encounter s4 before s2, we reenter:
        double s4 = dd*xp + sqrt( (R-d)*(R-d) / (1 + xp*xp) * dd*dd);  // second solution for smaller bend
        if (s3 < fmin(s1, s2) && s4 < s2){
            // We drift until re-entry
            Drift_single_particle_4d(part, s4 - exit_point);
            // We call the main Amorphous function for the leftover
            pc = Amorphous(everest, part, cg, pc, length - length_exit - s4 + exit_point, 1);
        }
    }

    return pc;
}

#endif /* XCOLL_EVEREST_AMORPHOUS_H */
