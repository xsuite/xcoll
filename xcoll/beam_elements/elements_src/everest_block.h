// copyright ############################### #
// This file is part of the Xcoll package.   #
// Copyright (c) CERN, 2024.                 #
// ######################################### #

#ifndef XCOLL_EVEREST_BLOCK_H
#define XCOLL_EVEREST_BLOCK_H

#ifdef XO_CONTEXT_CPU
#include <stdint.h>  // for int64_t etc
#endif /* XO_CONTEXT_CPU */

#include "xobjects/headers/common.h"
#include "xobjects/headers/atomicadd.h"
#include "xtrack/headers/checks.h"
#include "xtrack/headers/track.h"
#include "xtrack/beam_elements/elements_src/track_drift.h"
#include "xcoll/headers/checks.h"
#include "xcoll/scattering_routines/geometry/objects.h"
#include "xcoll/scattering_routines/geometry/collimator_geometry.h"
#include "xcoll/scattering_routines/everest/everest.h"
#include "xcoll/scattering_routines/everest/properties.h"
#include "xcoll/scattering_routines/everest/ionisation_loss.h"
#include "xcoll/scattering_routines/everest/jaw.h"


GPUFUN
int8_t EverestBlockData_get_record_impacts(EverestBlockData el) {
    return EverestBlockData_get__record_interactions(el) % 2;
}

GPUFUN
int8_t EverestBlockData_get_record_exits(EverestBlockData el) {
    return (EverestBlockData_get__record_interactions(el) >> 1) % 2;
}

GPUFUN
int8_t EverestBlockData_get_record_scatterings(EverestBlockData el) {
    return (EverestBlockData_get__record_interactions(el) >> 2) % 2;
}

GPUFUN
int8_t EverestBlockData_get_mark_scattered_particles(EverestBlockData el) {
    return (EverestBlockData_get__record_interactions(el) >> 3) % 2;
}


GPUKERN
void EverestBlock_set_material(EverestBlockData el) {
    RandomRutherford_set_by_xcoll_material(
        EverestBlockData_getp_rutherford_rng(el),
        EverestBlockData_getp__material(el)
    );
}


// TODO: it would be great if we could set EverestData as an xofield, because then we could
// run this function at creation of the collimator instead of every turn
// Hmmmm this should be called whenever we change an xofield
GPUFUN
void EverestBlock_init(EverestBlockData el, LocalParticle* part0, EverestCollData RESTRICT coll) {
    // Random generator
    coll->rng = EverestBlockData_getp_rutherford_rng(el);
    // Impact table:  need it here to record interactions
    coll->record = EverestBlockData_getp_internal_record(el, part0);
    coll->record_index = NULL;
    coll->record_scatterings = 0;
    if (coll->record) {
        coll->record_index = InteractionRecordData_getp__index(coll->record);
        coll->record_scatterings = EverestBlockData_get_record_scatterings(el);
    }
    coll->orient = 0;
}


GPUFUN
void EverestBlock_init_data(LocalParticle* part, MaterialData RESTRICT material,
                            EverestCollData RESTRICT coll, EverestData RESTRICT everest) {
    everest->coll = coll;
    everest->rescale_scattering = 1;
    // Preinitialise scattering parameters
    double energy = LocalParticle_get_energy(part) / 1e9; // energy in GeV
    calculate_scattering(everest, material, energy);
    calculate_ionisation_properties(everest, material, energy);
}


GPUFUN
void EverestBlock_track_local_particle(EverestBlockData el, LocalParticle* part0) {
    int8_t active = EverestBlockData_get__tracking(el);
    active       *= EverestBlockData_get_active(el);
    double const length = EverestBlockData_get_length(el);

    // Initialise collimator data
    // TODO: we want this to happen before tracking (instead of every turn), as a separate kernel
    EverestCollData_ coll_;
    EverestCollData coll = &coll_; // pointer
    MaterialData material;         // pointer
    if (active) {
        EverestBlock_init(el, part0, coll);
        material = EverestBlockData_getp__material(el);
    }

    START_PER_PARTICLE_BLOCK(part0, part);
        if (!active) {
            // Drift full length
            Drift_single_particle(part, length);

        } else {
            // Check collimator initialisation
            int8_t is_valid = xcoll_check_particle_init(coll->rng, part);

            if (is_valid) {
                // Store s-location of start of block
                double const s_block = LocalParticle_get_s(part);
                LocalParticle_set_s(part, 0);

                // Store initial coordinates for updating later
                double const rvv_in  = LocalParticle_get_rvv(part);
#ifdef XCOLL_USE_EXACT
                double const xp_in   = LocalParticle_get_exact_xp(part);
                double const yp_in   = LocalParticle_get_exact_yp(part);
#else
                double const xp_in   = LocalParticle_get_xp(part);
                double const yp_in   = LocalParticle_get_yp(part);
#endif
                double const zeta_in = LocalParticle_get_zeta(part);
                double const p0c     = LocalParticle_get_p0c(part);
                double const delta   = LocalParticle_get_delta(part);
                double const qq0     = LocalParticle_get_charge_ratio(part);
                double const chi     = LocalParticle_get_chi(part);
                double const pc_in   = (1 + delta)*p0c*qq0/chi;
                double const e_in    = LocalParticle_get_energy(part);
                double pc_out;

                // Scatter
                EverestData_ everest_;
                EverestData everest = &everest_;
                EverestBlock_init_data(part, material, coll, everest);
                pc_out = jaw(everest, material, part, pc_in, length, 0);
                LocalParticle_add_to_s(part, s_block);

                LocalParticle_set_zeta(part, zeta_in);

                // Survived particles need correcting:
                if (LocalParticle_get_state(part) > 0) {
                    double const rpp_old  = LocalParticle_get_rpp(part);
                    LocalParticle_update_delta(part, pc_out*chi/p0c/qq0 - 1);
                    // Keep angles constant (this is also correct for exact angles): px_new = px_old*(1 + δ_new)/(1 + δ_old)
                    double const scale = rpp_old / LocalParticle_get_rpp(part);
                    LocalParticle_scale_px(part, scale);
                    LocalParticle_scale_py(part, scale);

                    // Update zeta
#ifdef XCOLL_USE_EXACT
                    double xp  = LocalParticle_get_exact_xp(part);
                    double yp  = LocalParticle_get_exact_yp(part);
#else
                    double xp  = LocalParticle_get_xp(part);
                    double yp  = LocalParticle_get_yp(part);
#endif
                    double rvv = LocalParticle_get_rvv(part);
                    // First we drift half the length with the old angles:
                    LocalParticle_add_to_zeta(part, drift_zeta_single(rvv_in, xp_in, yp_in, length/2) );
                    // then half the length with the new angles:
                    LocalParticle_add_to_zeta(part, drift_zeta_single(rvv, xp, yp, length/2) );

                    // Store deposited energy in the block
                    double e_out = LocalParticle_get_energy(part);
                    if (LocalParticle_get_state(part) == XC_SECONDARY_PARTICLE) {
                        GPUGLMEM double *acc_loss = EverestBlockData_getp__acc_ionisation_loss_sec(el);
                        atomicAdd(acc_loss, e_in - e_out);
                    } else {
                        GPUGLMEM double *acc_loss = EverestBlockData_getp__acc_ionisation_loss(el);
                        atomicAdd(acc_loss, e_in - e_out);
                    }

                    // Mark scattered particles as secondaries (if desired)
                    if (EverestBlockData_get_mark_scattered_particles(el)) {
                        LocalParticle_set_state(part, XC_SECONDARY_PARTICLE);
                    }
                }
            }
        }
    END_PER_PARTICLE_BLOCK;
}


#endif /* XCOLL_EVEREST_BLOCK_H */
