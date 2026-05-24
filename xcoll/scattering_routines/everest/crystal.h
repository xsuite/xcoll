// copyright ############################### #
// This file is part of the Xcoll package.   #
// Copyright (c) CERN, 2026.                 #
// ######################################### #

#ifndef XCOLL_EVEREST_CRYSTAL_H
#define XCOLL_EVEREST_CRYSTAL_H

#ifdef XO_CONTEXT_CPU
#include <math.h>
#include <stdint.h>  // for int64_t etc
#include <stdlib.h>  // for malloc and free
#endif  // XO_CONTEXT_CPU


double do_crystal(EverestData restrict everest, MaterialData restrict material,
                  LocalParticle* part, CrystalGeometry restrict cg, double pc, double length) {
    calculate_initial_angle(everest, part, cg);
    calculate_critical_angle(everest, material, part, cg, pc);
#ifndef XCOLL_REFINE_ENERGY
    // Calculate once at start
    calculate_scattering(everest, material, pc);
    calculate_ionisation_properties(everest, material, pc);
    calculate_VI_parameters(everest, part, pc);
#endif

#ifdef XCOLL_USE_EXACT
    double const xp = LocalParticle_get_exact_xp(part);
#else
    double const xp = LocalParticle_get_xp(part);
#endif
    if (fabs(xp - everest->t_I) < everest->t_c) {
        double alpha = fabs(xp - everest->t_I) / everest->t_c;
        double ratio = everest->Rc_over_R;
        double eta = MaterialData_get__eta(material);
        double xi = RandomUniform_generate(part)/(1 - ratio)/sqrt(eta);
        if (xi > 1 || alpha > 2*sqrt(xi)*sqrt(1-xi)) {
#ifdef XCOLL_TRANSITION_VRCH
            volume_reflection(everest, part, XC_VOLUME_REFLECTION_TRANS_CH);
#endif
            pc = Amorphous(everest, material, part, cg, pc, length, 1);
        } else {
            pc = Channel(everest, material, part, cg, pc, length);
        }
    } else {
        pc = Amorphous(everest, material, part, cg, pc, length, 1);
    }
    return pc;
}


/*gpufun*/
double do_crystal_precise(EverestData restrict everest, MaterialData restrict material,
                          LocalParticle* part, CrystalGeometry restrict cg, double pc, double length) {
    calculate_initial_angle(everest, part, cg);
#ifndef XCOLL_REFINE_ENERGY
    // Calculate once at start
    calculate_scattering(everest, material, pc);
    calculate_ionisation_properties(everest, material, pc);
    calculate_critical_angle(everest, material, part, cg, pc);
    calculate_VI_parameters(everest, part, pc);
    calculate_potentials(everest, cg, pc);
#endif

    double* result = calculate_local_coords(everest, material, part, cg);
    double x_local = result[0];
    double xp_local = result[1];
    free(result);

    if (can_channel_precise(everest, part, cg, x_local, xp_local, pc)) {
        pc = Channel(everest, material, part, cg, pc, length);
    } else {
        pc = Amorphous(everest, material, part, cg, pc, length, 1);
    }
    return pc;
}

#endif /* XCOLL_EVEREST_CRYSTAL_H */
