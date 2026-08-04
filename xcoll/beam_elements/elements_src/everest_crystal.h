// copyright ############################### #
// This file is part of the Xcoll package.   #
// Copyright (c) CERN, 2024.                 #
// ######################################### #

#ifndef XCOLL_EVEREST_CRYSTAL_H
#define XCOLL_EVEREST_CRYSTAL_H

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
#include "xcoll/scattering_routines/geometry/crystal_geometry.h"
#include "xcoll/scattering_routines/everest/everest.h"
#include "xcoll/scattering_routines/everest/properties.h"
#include "xcoll/scattering_routines/everest/ionisation_loss.h"
#include "xcoll/scattering_routines/everest/crystal_parameters.h"
#include "xcoll/scattering_routines/everest/channelling.h"


GPUFUN
int8_t EverestCrystalData_get_record_impacts(EverestCrystalData el) {
    return EverestCrystalData_get__record_interactions(el) % 2;
}

GPUFUN
int8_t EverestCrystalData_get_record_exits(EverestCrystalData el) {
    return (EverestCrystalData_get__record_interactions(el) >> 1) % 2;
}

GPUFUN
int8_t EverestCrystalData_get_record_scatterings(EverestCrystalData el) {
    return (EverestCrystalData_get__record_interactions(el) >> 2) % 2;
}

GPUFUN
int8_t EverestCrystalData_get_mark_scattered_particles(EverestCrystalData el) {
    return (EverestCrystalData_get__record_interactions(el) >> 3) % 2;
}


GPUKERN
void EverestCrystal_set_material(EverestCrystalData el) {
    RandomRutherford_set_by_xcoll_material(
        EverestCrystalData_getp_rutherford_rng(el),
        EverestCrystalData_getp__material(el)
    );
}


GPUFUN
void EverestCrystal_init_geometry(EverestCrystalData el, LocalParticle* part0,
                                  CrystalGeometry RESTRICT cg) {
    // Dimensions
    cg->length = EverestCrystalData_get_length(el);
    cg->side   = EverestCrystalData_get__side(el);
    cg->bending_radius = EverestCrystalData_get__bending_radius(el);
    cg->bending_angle  = EverestCrystalData_get__bending_angle(el);
    cg->miscut_angle   = EverestCrystalData_get_miscut(el);
    cg->width  = EverestCrystalData_get__width(el);
    cg->height = EverestCrystalData_get__height(el);
    if (cg->width<=0) {
        xcoll_kill_with_message(part0, XC_ERR_INVALID_XOFIELD,
                                "Crystal width should be positive!");
        return;
    }
    if (cg->height<=0) {
        xcoll_kill_with_message(part0, XC_ERR_INVALID_XOFIELD,
                                "Crystal height should be positive!");
        return;
    }
    if (cg->bending_radius==0) {
        xcoll_kill_with_message(part0, XC_ERR_INVALID_XOFIELD,
                                "Crystal bending radius should be non-zero!");
        return;
    }
    if (cg->side != 1 && cg->side != -1) {
        xcoll_kill_with_message(part0, XC_ERR_INVALID_XOFIELD,
                                "Crystal side should be either 1 or -1!");
        return;
    }
    // Jaw corners (with tilts)
    cg->jaw_U = EverestCrystalData_get__jaw_U(el);
    // Jaw angles
    cg->sin_z = EverestCrystalData_get__sin_z(el);
    cg->cos_z = EverestCrystalData_get__cos_z(el);
    // Jaw tilts
    cg->sin_y = EverestCrystalData_get__sin_y(el);
    cg->cos_y = EverestCrystalData_get__cos_y(el);
    // Segments
    int8_t status = 0;
    double jaw;
    if (cg->side == 1) {
        jaw = cg->jaw_U;
    } else if (cg->side == -1) {
        jaw = cg->jaw_U - cg->width; // To ensure that jaw_U is the inner corner
    }
    status = create_crystal(cg->segments, cg->bending_radius, cg->width,
                            cg->length, jaw, cg->sin_y, cg->cos_y);
    if (status != 0) {
        xcoll_kill_with_message(part0, XC_ERR_FAILED_GEOMETRY,
                                "Failed creating crystal geometry!");
        return;
    }
    // Bend centre
    double R = cg->bending_radius;
    cg->s_B = 0;
    cg->x_B = R;
    // Miscut centre
    cg->s_P = -R*sin(cg->miscut_angle);
    cg->x_P = R*cos(cg->miscut_angle);
    if (cg->side == 1 && R < 0) {
        // If R<0, a left-sided crystal bends towards the beam
        cg->x_P = cg->x_P + cg->width;
        cg->x_B = cg->x_B + cg->width;
    } else if (cg->side == -1 && R > 0) {
        // If R>0, a right-sided crystal bends towards the beam
        cg->x_P = cg->x_P - cg->width;
        cg->x_B = cg->x_B - cg->width;
    }
    if (cg->side == -1) {
        // Mirror the crystal geometry
        cg->bending_radius = -cg->bending_radius;
        cg->bending_angle  = -cg->bending_angle;
        cg->miscut_angle   = -cg->miscut_angle;
        cg->x_P            = -cg->x_P;
        cg->x_B            = -cg->x_B;
    }
    // From here on, crystal geometry parameters can always be treated as left-sided.
    // Note that the segments are not mirrored, which is fine as get_s_of_first_crossing_with_vlimit
    // is absolute (not in the jaw reference frame). It is only after a hit is registered, that we
    // need to transform the particle to the jaw reference frame.
    double Rb;
    if (cg->miscut_angle > 0) {
        Rb = R - cg->width;
    } else {
        Rb = R;
    }
    double t_R = cg->bending_angle;
    cg->t_VImax = atan( (Rb*sin(t_R) - cg->s_P) / (R - Rb*cos(t_R) - cg->x_P) );
    // Impact table
    cg->record = EverestCrystalData_getp_internal_record(el, part0);
    cg->record_index = NULL;
    cg->record_impacts = 0;
    cg->record_exits = 0;
    if (cg->record) {
        cg->record_index = InteractionRecordData_getp__index(cg->record);
        cg->record_impacts = EverestCrystalData_get_record_impacts(el);
        cg->record_exits = EverestCrystalData_get_record_exits(el);
    }
}


// TODO: it would be great if we could set EverestData as an xofield, because then we could
// run this function at creation of the collimator instead of every turn
// Hmmmm this should be called whenever we change an xofield
GPUFUN
void EverestCrystal_init(EverestCrystalData el, LocalParticle* part0, EverestCollData RESTRICT coll) {
    // Random generator
    coll->rng = EverestCrystalData_getp_rutherford_rng(el);
    // Impact table:  need it here to record interactions
    coll->record = EverestCrystalData_getp_internal_record(el, part0);
    coll->record_index = NULL;
    if (coll->record) {
        coll->record_index = InteractionRecordData_getp__index(coll->record);
        coll->record_scatterings = EverestCrystalData_get_record_scatterings(el);
    }
    coll->orient = EverestCrystalData_get__orient(el);
}


GPUFUN
void EverestCrystal_init_data(LocalParticle* part, MaterialData RESTRICT material,
                              EverestCollData RESTRICT coll, CrystalGeometry RESTRICT cg,
                              EverestData RESTRICT everest) {
    everest->coll = coll;
    everest->rescale_scattering = 1;
    // Preinitialise scattering parameters
    double energy = LocalParticle_get_energy(part) / 1e9; // energy in GeV
    calculate_scattering(everest, material, energy);
    calculate_ionisation_properties(everest, material, energy);
    calculate_critical_angle(everest, material, part, cg, energy);
    calculate_VI_parameters(everest, part, energy);
}


GPUFUN
void EverestCrystal_track_local_particle(EverestCrystalData el, LocalParticle* part0) {
    int8_t active = EverestCrystalData_get_active(el);
    active       *= EverestCrystalData_get__tracking(el);
    double const length = EverestCrystalData_get_length(el);

    // Initialise collimator data
    // TODO: we want this to happen before tracking (instead of every turn), as a separate kernel
    EverestCollData_ coll_;
    CrystalGeometry_ cg_;
    EverestCollData coll = &coll_; // pointer
    CrystalGeometry cg = &cg_;     // pointer
    MaterialData material;         // pointer
    if (active) {
        EverestCrystal_init(el, part0, coll);
        EverestCrystal_init_geometry(el, part0, cg);
        material = EverestCrystalData_getp__material(el);

        // For info
        double const e0 = LocalParticle_get_energy0(part0);
        double t_c0  = _critical_angle0(material, e0);
        double Rcrit = _critical_radius(material, e0);
        double t_c = _critical_angle(coll, t_c0, Rcrit / fabs(cg->bending_radius));
        EverestCrystalData_set__critical_radius(el, Rcrit);
        EverestCrystalData_set__critical_angle(el, t_c);
    }

    START_PER_PARTICLE_BLOCK(part0, part);
        if (!active) {
            // Drift full length
            Drift_single_particle(part, length);

        } else {
            // Check collimator initialisation
            int8_t is_valid = xcoll_check_particle_init(coll->rng, part);

            if (is_valid) {
                // Store s-location of start of crystal
                double const s_coll = LocalParticle_get_s(part);
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

                // Check if hit on jaws
                int8_t is_hit = hit_crystal_check_and_transform(part, cg);

                if (is_hit != 0) {
                    // Hit one of the jaws, so scatter
                    double remaining_length = length - LocalParticle_get_s(part);
                    // Scatter
                    EverestData_ everest_;
                    EverestData everest = &everest_;
                    EverestCrystal_init_data(part, material, coll, cg, everest);
                    everest->shape_id = is_hit;
                    pc_out = do_crystal(everest, material, part, cg, pc_in/1.e9, remaining_length)*1.e9;
                }

                // Transform back to the lab frame
                hit_crystal_transform_back(is_hit, part, cg);
                LocalParticle_add_to_s(part, s_coll);

                LocalParticle_set_zeta(part, zeta_in);

                // Hit and survived particles need correcting:
                if (is_hit!=0 && LocalParticle_get_state(part) > 0) {
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

                    // Store deposited energy in the crystal
                    double e_out = LocalParticle_get_energy(part);
                    if (LocalParticle_get_state(part) == XC_SECONDARY_PARTICLE) {
                        GPUGLMEM double *acc_loss = EverestCrystalData_getp__acc_ionisation_loss_sec(el);
                        atomicAdd(acc_loss, e_in - e_out);
                    } else {
                        GPUGLMEM double *acc_loss = EverestCrystalData_getp__acc_ionisation_loss(el);
                        atomicAdd(acc_loss, e_in - e_out);
                    }

                    // Mark scattered particles as secondaries (if desired)
                    if (EverestCrystalData_get_mark_scattered_particles(el)) {
                        LocalParticle_set_state(part, XC_SECONDARY_PARTICLE);
                    }
                }
            }
        }
    END_PER_PARTICLE_BLOCK;
}


#endif /* XCOLL_EVEREST_CRYSTAL_H */
