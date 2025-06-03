// copyright ############################### #
// This file is part of the Xcoll package.   #
// Copyright (c) CERN, 2023.                 #
// ######################################### #

#ifndef XCOLL_EVEREST_AMORPHOUS_H
#define XCOLL_EVEREST_AMORPHOUS_H
#include <math.h>
#include <stdio.h>
#include <stdlib.h>


/*gpufun*/
void volume_reflection(EverestData restrict everest, LocalParticle* part, int8_t transition) {

    InteractionRecordData record = everest->coll->record;
    RecordIndex record_index     = everest->coll->record_index;
    int8_t sc = everest->coll->record_scatterings;
    int64_t i_slot = -1;

    double Ang_avr = everest->Ang_avr;
    double Ang_rms = everest->Ang_rms * RandomNormal_generate(part);
    // TODO: Use standard MCS
    // TODO: should we automatically channel after VC? Now it will still roll a dice, and do VRCH if failed. Many VR are missed this way

    if (transition == XC_VOLUME_REFLECTION_TRANS_CH){
        // We are in transition from CH to VR
        double xp = LocalParticle_get_xp(part);
        // TODO: we believe the original 0.45 comes from the 0.9 saturation factor, so we changed it to 0.5
        Ang_avr *= 0.5*((xp - everest->t_I)/everest->t_c + 1);
        Ang_rms  = 0;   // TODO: why does transition to CH not use any spread?
        if (sc) i_slot = InteractionRecordData_log(record, record_index, part, XC_VOLUME_REFLECTION_TRANS_CH);

    } else if (transition == XC_VOLUME_REFLECTION_TRANS_MCS){
        if (sc) i_slot = InteractionRecordData_log(record, record_index, part, XC_VOLUME_REFLECTION_TRANS_MCS);

    } else {
        if (sc) i_slot = InteractionRecordData_log(record, record_index, part, XC_VOLUME_REFLECTION);

    }
    LocalParticle_add_to_xp(part, Ang_avr + Ang_rms);
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
    } else {
        // Normal MCS
        if (sc) i_slot = InteractionRecordData_log(record, record_index, part, XC_MULTIPLE_COULOMB_SCATTERING);
    }
    kxmcs = dya*RandomNormal_generate(part);
    kymcs = dya*RandomNormal_generate(part);

    // Transport of Multiple Coulomb Scattering
    Drift_single_particle_4d(part, length);

    // Energy lost because of ionisation process[GeV]
    pc = calcionloss(everest, part, length, pc, 1);

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
    // Assumption: crystal bends away from beam (this does not work when crystal bends towards beam)
    // All s are absolute in the collimator frame
    double R = cg->bending_radius;
    double t = cg->bending_angle;
    double d = cg->width;
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
    // hence the longitudinal length is r sin(xp - t_I) cos xp
    double length_VI = 1e10;
    if (xp > everest->t_I){
    // xp has to be larger than t_I to be able to VR (no matter which transition we consider)
#ifdef XCOLL_TRANSITION_VRAM_OLD
        // Transition region between VR and AM for t_B < xp - tI < t_B + 2t_c
        length_VI = everest->r * sin(xp - everest->t_I) * cos(xp);
#else
#ifdef XCOLL_TRANSITION_VRAM
        // Transition region between VR and AM for t_B - t_c < xp - tI < t_B + t_c
        // We need to start transitioning earlier than t_B
        length_VI = everest->r * sin(xp - everest->t_I + everest->t_c) * cos(xp);
#else
        // Normal behaviour, no transition
        length_VI = everest->r * sin(xp - everest->t_I) * cos(xp);
#endif
#endif
    }
    // Calculate extra length to transition region between VR and AM
    double length_VR_trans = 1e10;
#ifdef XCOLL_TRANSITION_VRAM_OLD
    if (xp - 2.*everest->t_c > everest->t_I){
        length_VR_trans = everest->r * sin(xp - everest->t_I - 2.*everest->t_c) * cos(xp);
    }
#else
#ifdef XCOLL_TRANSITION_VRAM
    if (xp - everest->t_c > everest->t_I){
        length_VR_trans = everest->r * sin(xp - everest->t_I - everest->t_c) * cos(xp);
    }
#endif
#endif

    // ------------------------------------------------------------------------
    // Compare the 3 lengths: the first one encountered is what will be applied
    // ------------------------------------------------------------------------

    if (length_nucl < fmin(length_VI, length_exit)) {
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

    } else if (length_VI <= length_exit && allow_VI == 1){
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

    } else if (length_VR_trans <= length_exit && allow_VI == 1){
#ifdef XCOLL_REFINE_ENERGY
        calculate_critical_angle(everest, part, cg, pc);
#endif
        // We estimate where we are in the transition region (rather on the VR side or rather on the AM side)
        // by looking at how much correction we needed on the length. If length_VI is closest to length_exit,
        // we are on the VR side. On the other hand, if length_VR_trans is closest to length_exit, we are on the
        // AM side. Note that, by reaching this part we automatically have length_VI > length_exit > length_VR_trans.
        double prob_MCS = (length_VI - length_exit) / (length_VI - length_VR_trans);
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
