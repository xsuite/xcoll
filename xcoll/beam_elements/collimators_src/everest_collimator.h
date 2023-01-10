#ifndef XCOLL_EVEREST_H
#define XCOLL_EVEREST_H

#include <math.h>

// TODO:
//    Do not split 4d and zeta in drifts
//    Use drift function from xtrack Drift element (call its C function)

/*gpufun*/
void drift_6d(LocalParticle* part0, double length) {
    //start_per_particle_block (part0->part)
        double const rpp    = LocalParticle_get_rpp(part);
        double const rv0v   = 1./LocalParticle_get_rvv(part);
        double const xp     = LocalParticle_get_px(part) * rpp;
        double const yp     = LocalParticle_get_py(part) * rpp;
        double const dzeta  = 1 - rv0v * ( 1. + ( xp*xp + yp*yp ) / 2. );

        LocalParticle_add_to_x(part, xp * length );
        LocalParticle_add_to_y(part, yp * length );
        LocalParticle_add_to_s(part, length);
        LocalParticle_add_to_zeta(part, length * dzeta );
    //end_per_particle_block
}

/*gpufun*/
void drift_4d(LocalParticle* part0, double length) {
    //start_per_particle_block (part0->part)
        double const rpp    = LocalParticle_get_rpp(part);
        double const xp     = LocalParticle_get_px(part) * rpp;
        double const yp     = LocalParticle_get_py(part) * rpp;
        LocalParticle_add_to_x(part, xp * length );
        LocalParticle_add_to_y(part, yp * length );
    //end_per_particle_block
}

/*gpufun*/
double drift_zeta_single(double rvv, double xp, double yp, double length){
    double const rv0v = 1./rvv;
    double const dzeta = 1 - rv0v * (1. + (pow(xp,2.) + pow(yp,2.))/2.);
    return length * dzeta;
}

/*gpufun*/
void shift_4d(LocalParticle* part0, double dx, double dpx, double dy, double dpy) {
    //start_per_particle_block (part0->part)
        LocalParticle_add_to_x(part, dx);
        LocalParticle_add_to_px(part, dpx);
        LocalParticle_add_to_y(part, dy);
        LocalParticle_add_to_py(part, dpy);
    //end_per_particle_block
}


// TODO: 
// Write impacts



/*gpufun*/
void track_collimator(CollimatorData el, LocalParticle* part0) {

    // Store initial coordinates for updating later
    double *x_in, *px_in, *y_in, *py_in, *zeta_in, *ptau_in, *rpp_in, *rvv_in;
    //start_per_particle_block (part0->part)
        long long pid = LocalParticle_get_particle_id(part);
        x_in[pid]    = LocalParticle_get_x(part);
        px_in[pid]   = LocalParticle_get_px(part);
        y_in[pid]    = LocalParticle_get_y(part);
        py_in[pid]   = LocalParticle_get_py(part);
        zeta_in[pid] = LocalParticle_get_zeta(part);
        ptau_in[pid] = LocalParticle_get_ptau(part);
        rpp_in[pid]  = LocalParticle_get_rpp(part);
        rvv_in[pid]  = LocalParticle_get_rvv(part);
    //end_per_particle_block
    double e0_ref    = LocalParticle_get_energy0(&part0[0]) / 1e9; // Reference energy in GeV
    double beta0_ref = LocalParticle_get_beta0(&part0[0]);

    double const length  = CollimatorData_get_active_length(el);

    // Drift to centre of collimator
    drift_4d(part0, length/2);
    // Move to closed orbit  (dpx = dxp because ref. particle has delta = 0)
    shift_4d(part0, -CollimatorData_get_dx(el), -CollimatorData_get_dpx(el), -CollimatorData_get_dy(el), -CollimatorData_get_dpy(el) );
    // Backtrack to start of collimator
    drift_4d(part0, -length/2);

    // Initialise status arrays
    double *part_hit, *part_abs, *part_impact, *part_indiv, *part_linteract, *nabs_type, *linside;
    //start_per_particle_block (part0->part)
        long long pid = LocalParticle_get_particle_id(part);
        part_hit[pid] = 0;
        part_abs[pid] = 0;
        part_impact[pid] = 0;
        part_indiv[pid] = 0;
        part_linteract[pid] = 0;
        nabs_type[pid] = 0;
        linside[pid] = 0;
    //end_per_particle_block

    // if collimator.jaw_F_L != collimator.jaw_B_L or collimator.jaw_F_R != collimator.jaw_B_R:
    //     raise NotImplementedError
    double opening = CollimatorData_get_jaw_F_L(el) - CollimatorData_get_jaw_F_R(el);
    double offset  = CollimatorData_get_offset(el) + ( CollimatorData_get_jaw_F_L(el) + CollimatorData_get_jaw_F_R(el) )/2;

    // Get material properties
    MaterialData material = CollimatorData_get_material(el);
    double zatom    = MaterialData_get_Z(material);
    double anuc     = MaterialData_get_A(material);
    double rho      = MaterialData_get_density(material);
    double exenergy = MaterialData_get_excitation_energy(material);
    double emr      = MaterialData_get_nuclear_radius(material);
    double bnref    = MaterialData_get_nuclear_elastic_slope(material);
    double csref0   = MaterialData_get_cross_section(material)[0];
    double csref1   = MaterialData_get_cross_section(material)[1];
    double csref5   = MaterialData_get_cross_section(material)[5];
    double hcut     = MaterialData_get_hcut(material);
    double radl     = MaterialData_get_radiation_length(material);

    // Calculate scattering parameters
    double* result = calculate_scattering(e0_ref,anuc,rho,zatom,emr,csref0,csref1,csref5,bnref);
    double cprob0 = result[0];
    double cprob1 = result[1];
    double cprob2 = result[2];
    double cprob3 = result[3];
    double cprob4 = result[4];
    double cprob5 = result[5];
    double xintl = result[6];
    double bn = result[7];
    double ecmsq = result[8];
    double xln15s = result[9];
    double bpp = result[10];

    set_rutherford_parameters(zatom, emr, hcut);

    // Initialise accumulated variables
    double p0 = e0_ref;
    double x0 = 0;
    double xp0 = 0;
    double nhit   = 0;
    double nabs   = 0;
    double fracab = 0;
    // Set energy and nucleon change variables as with the coupling
    // ien0, ien1: particle energy entering/leaving the collimator
    double nnuc0 = 0;
    double ien0  = 0;
    double nnuc1 = 0;
    double ien1  = 0;

    double *energy_out, *s_out;

    //start_per_particle_block (part0->part)
        long long pid = LocalParticle_get_particle_id(part);

        if (part_abs[pid] != 0){
            continue;
        }

        double rpp = LocalParticle_get_rpp(part);
        double angle = atan2(CollimatorData_get_sin_z(el), Collimator_get_cos_z(el) );

        double* result = scatter(
                LocalParticle_get_x(part),
                LocalParticle_get_px(part)*rpp,
                LocalParticle_get_y(part),
                LocalParticle_get_py(part)*rpp,
                0,   // s at start of collimator
                LocalParticle_get_energy(part) / 1e9,      // energy, not momentum, in GeV
                part_hit[pid],
                part_abs[pid],
                part_impact[pid],         // impact parameter
                part_indiv[pid],          // particle divergence
                part_linteract[pid],      // interaction length
                nabs_type[pid],
                linside[pid],
                exenergy,
                anuc,
                zatom,
                emr,
                rho,
                hcut,
                bnref,
                csref0,
                csref1,
                csref5,
                radl,
                0,   // dlri
                0,   // dlyi
                0,   // eUm
                0,   // ai
                0,   // collnt
                cprob0,
                cprob1,
                cprob2,
                cprob3,
                cprob4,
                cprob5,
                xintl,
                bn,
                ecmsq,
                xln15s,
                bpp,
                0,   // is_crystal
                length,
                angle,
                opening,
                offset,
                CollimatorData_get_tilt(el)[0],
                CollimatorData_get_tilt(el)[1],
                CollimatorData_get_onesided(el),
                0,   // cry_tilt
                0,   // cry_rcurv
                0,   // cry_bend
                0,   // cry_alayer
                0,   // cry_xmax
                0,   // cry_ymax
                0,   // cry_orient
                0,   // cry_miscut
                0,   // cry_cBend
                0,   // cry_sBend
                0,   // cry_cpTilt
                0,   // cry_spTilt
                0,   // cry_cnTilt
                0,   // cry_snTilt
                p0,
                x0,
                xp0,
                nhit,
                nabs,
                fracab,
                nnuc0,
                ien0,
                nnuc1,
                ien1,
                0,   // iProc
                0,   // n_chan
                0,   // n_VR
                0,   // n_amorphous
                0,   // s_imp
                )

        LocalParticle_set_x(part, result[0]);
        LocalParticle_set_px(part/rpp_in, result[1]);
        LocalParticle_set_y(part, result[2]);
        LocalParticle_set_py(part/rpp_in, result[3]);
        s_out[pid] = result[4];
        energy_out[pid] = result[5];       //  Cannot assign energy directly to LocalParticle as it would update dependent variables, but needs to be corrected first!
        part_hit[pid]=result[6];
        part_abs[pid]=result[7];
        part_impact[pid]=result[8];
        part_indiv[pid]=result[9];
        part_linteract[pid]=result[10];
        nabs_type[pid]=result[11];
        linside[pid]=result[12];
        p0=result[13];
        x0=result[14];
        xp0=result[15];
        nhit=result[16];
        nabs=result[17];
        fracab=result[18];
        nnuc0=result[19];
        ien0=result[20];
        nnuc1=result[21];
        ien1=result[22];

    //end_per_particle_block

    # Masks of hit and survived particles
    lost = part_abs > 0
    hit = part_hit > 0
    not_hit = ~hit
    not_lost = ~lost
    survived_hit = hit & (~lost)

    // Backtrack to centre of collimator
    drift_4d(part0, -length/2);
    // Return from closed orbit  (dpx = dxp because ref. particle has delta = 0)
    shift_4d(part0, CollimatorData_get_dx(el), CollimatorData_get_dpx(el), CollimatorData_get_dy(el), CollimatorData_get_dpy(el) );

    // Update energy    ---------------------------------------------------
    // Only particles that hit the jaw and survived need to be updated
    //start_per_particle_block (part0->part)
        long long pid = LocalParticle_get_particle_id(part);
        double ptau_out = ptau_in[pid];
        int survived_hit = part_hit[pid] + (1-part_abs[pid]);
        if (survived_hit){
            ptau_out = (energy_out[pid] * 1e9 - e0_ref) / e0_ref / beta0_ref;
        }
        LocalParticle_set_ptau(part, ptau_out);
        // rpp_out[pid] = LocalParticle_get_rpp(part);
    //end_per_particle_block

    // // Rescale angles, because K2 did not update angles with new energy!
    // // So need to do xp = xp * p_in / p_out = xp * rpp_out / rpp_in
    // // (see collimation.f90 line 1709 and mod_particles.f90 line 210)
    // xp_part *= rpp_out/rpp_in
    // yp_part *= rpp_out/rpp_in

    // Drift to end of collimator
    drift_4d(part0, length/2);

    // Update 4D coordinates    -------------------------------------------
    // Absorbed particles get their coordinates set to the entrance of collimator
    //start_per_particle_block (part0->part)
        long long pid = LocalParticle_get_particle_id(part);
        if (part_abs[pid]>0){
            LocalParticle_set_x(part,x_in[pid]);
            LocalParticle_set_px(part,px_in[pid]);
            LocalParticle_set_y(part,y_in[pid]);
            LocalParticle_set_py(part,py_in[pid]);
        }
    //end_per_particle_block

    // Update longitudinal coordinate zeta    -----------------------------
    // Absorbed particles keep coordinates at the entrance of collimator, others need correcting:
    //start_per_particle_block (part0->part)
        long long pid = LocalParticle_get_particle_id(part);
        // Non-hit particles are just drifting (zeta not yet drifted in K2, so do here)
        if (part_hit[pid]==0){
            LocalParticle_add_to_zeta(part, drift_zeta_single(rvv_in[pid], px_in[pid]*rpp_in[pid], py_in[pid]*rpp_in[pid], length) );
        }
        // Hit and survived particles need correcting:
        if (part_hit[pid]>0 && part_abs[pid]==0){
            // First we drift half the length with the old angles:
            LocalParticle_add_to_zeta(part, drift_zeta_single(rvv_in[pid], px_in[pid]*rpp_in[pid], py_in[pid]*rpp_in[pid], length/2) );
            // then half the length with the new angles:
            double px  = LocalParticle_get_px(part);
            double py  = LocalParticle_get_py(part);
            double rvv = LocalParticle_get_rvv(part);
            double rpp = LocalParticle_get_rpp(part);
            LocalParticle_add_to_zeta(part, drift_zeta_single(rvv, px*rpp, py*rpp, length/2) );
        }
    //end_per_particle_block

    // Update s    --------------------------------------------------------
    // TODO: move absorbed particles to last impact location
    //start_per_particle_block (part0->part)
        long long pid = LocalParticle_get_particle_id(part);
        if (part_abs[pid]==0){
            LocalParticle_add_to_s(part, length);
        }
    //end_per_particle_block
 
    // Update state    ----------------------------------------------------
    //start_per_particle_block (part0->part)
        long long pid = LocalParticle_get_particle_id(part);
        if (part_abs[pid] > 0){
            LocalParticle_set_state(part, -333);
        }
    //end_per_particle_block
}

/*gpufun*/
void Collimator_track_local_particle(CollimatorData el, LocalParticle* part0) {
    int8_t const is_active      = CollimatorData_get__active(el);
    double const inactive_front = CollimatorData_get_inactive_front(el);
    double const active_length  = CollimatorData_get_active_length(el);
    double const inactive_back  = CollimatorData_get_inactive_back(el);

    if (!is_active){
        // Drift full length
        drift_6d(part0, inactive_front+active_length+inactive_back);
    } else {
        drift_6d(part0, inactive_front);
        track_collimator(el, part0);
        drift_6d(part0, inactive_back);
    }
}

#endif
