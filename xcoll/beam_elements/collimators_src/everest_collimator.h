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
void drift_4d_single(LocalParticle* part, double length) {
    double const rpp    = LocalParticle_get_rpp(part);
    double const xp     = LocalParticle_get_px(part) * rpp;
    double const yp     = LocalParticle_get_py(part) * rpp;
    LocalParticle_add_to_x(part, xp * length );
    LocalParticle_add_to_y(part, yp * length );
}

/*gpufun*/
double drift_zeta_single(double rvv, double xp, double yp, double length){
    double const rv0v = 1./rvv;
    double const dzeta = 1 - rv0v * (1. + (pow(xp,2.) + pow(yp,2.))/2.);
    return length * dzeta;
}

/*gpufun*/
void shift_4d_single(LocalParticle* part, double dx, double dpx, double dy, double dpy) {
    LocalParticle_add_to_x(part, dx);
    LocalParticle_add_to_px(part, dpx);
    LocalParticle_add_to_y(part, dy);
    LocalParticle_add_to_py(part, dpy);
}


// TODO: 
// Write impacts



/*gpufun*/
void track_collimator(EverestCollimatorData el, LocalParticle* part0) {

    double const e0_ref    = LocalParticle_get_energy0(&part0[0]) / 1e9; // Reference energy in GeV
    double const p0c_ref   = LocalParticle_get_p0c(&part0[0]) / 1e9;
    double const beta0_ref = LocalParticle_get_beta0(&part0[0]);

    // Collimator properties
    double const length  = EverestCollimatorData_get_active_length(el);
    // if collimator.jaw_F_L != collimator.jaw_B_L or collimator.jaw_F_R != collimator.jaw_B_R:
    //     raise NotImplementedError
    double const opening  = EverestCollimatorData_get_jaw_F_L(el) - EverestCollimatorData_get_jaw_F_R(el);
    double const offset   = EverestCollimatorData_get_offset(el) + ( EverestCollimatorData_get_jaw_F_L(el) + EverestCollimatorData_get_jaw_F_R(el) )/2;
    double const tilt0    = EverestCollimatorData_get_tilt(el, 0);
    double const tilt1    = EverestCollimatorData_get_tilt(el, 1);
    double const onesided = EverestCollimatorData_get_onesided(el);
    double const angle    = atan2(EverestCollimatorData_get_sin_z(el), EverestCollimatorData_get_cos_z(el) );
    double const co_x     = EverestCollimatorData_get_dx(el);
    double const co_px    = EverestCollimatorData_get_dpx(el);
    double const co_y     = EverestCollimatorData_get_dy(el);
    double const co_py    = EverestCollimatorData_get_dpy(el);

    // Material properties
    MaterialData material = EverestCollimatorData_getp_material(el);
    double const zatom    = MaterialData_get_Z(material);
    double const anuc     = MaterialData_get_A(material);
    double const rho      = MaterialData_get_density(material);
    double const exenergy = MaterialData_get_excitation_energy(material);
    double const emr      = MaterialData_get_nuclear_radius(material);
    double const bnref    = MaterialData_get_nuclear_elastic_slope(material);
    double const csref0   = MaterialData_get_cross_section(material, 0);
    double const csref1   = MaterialData_get_cross_section(material, 1);
    double const csref5   = MaterialData_get_cross_section(material, 5);
    double const hcut     = MaterialData_get_hcut(material);
    double const radl     = MaterialData_get_radiation_length(material);

    // Calculate scattering parameters
    double* result_scat_init = calculate_scattering(e0_ref,anuc,rho,zatom,emr,csref0,csref1,csref5,bnref);
    double const cprob0 = result_scat_init[0];
    double const cprob1 = result_scat_init[1];
    double const cprob2 = result_scat_init[2];
    double const cprob3 = result_scat_init[3];
    double const cprob4 = result_scat_init[4];
    double const cprob5 = result_scat_init[5];
    double const xintl  = result_scat_init[6];
    double const bn     = result_scat_init[7];
    double const ecmsq  = result_scat_init[8];
    double const xln15s = result_scat_init[9];
    double const bpp    = result_scat_init[10];

    set_rutherford_parameters(zatom, emr, hcut);

    // Initialise accumulated variables   TODO: this is NOT GPU-proof...
    double p0     = e0_ref;
    double x0     = 0;
    double xp0    = 0;
    double nhit   = 0;
    double nabs   = 0;
    double fracab = 0;
    // Set energy and nucleon change variables as with the coupling
    // ien0, ien1: particle energy entering/leaving the collimator
    double nnuc0 = 0;
    double ien0  = 0;
    double nnuc1 = 0;
    double ien1  = 0;

    //start_per_particle_block (part0->part)

        // Store initial coordinates for updating later
        double const x_in    = LocalParticle_get_x(part);
        double const px_in   = LocalParticle_get_px(part);
        double const y_in    = LocalParticle_get_y(part);
        double const py_in   = LocalParticle_get_py(part);
//         double const zeta_in = LocalParticle_get_zeta(part);
        double const ptau_in = LocalParticle_get_ptau(part);
        double const rpp_in  = LocalParticle_get_rpp(part);
        double const rvv_in  = LocalParticle_get_rvv(part);

        // Drift to centre of collimator
        drift_4d_single(part, length/2);
        // Move to closed orbit  (dpx = dxp because ref. particle has delta = 0)
        shift_4d_single(part, -co_x, -co_px/rpp_in, -co_y, -co_py/rpp_in );
        // Backtrack to start of collimator
        drift_4d_single(part, -length/2);

        // Status variables
        int part_hit = 0;
        int part_abs = 0;
        int part_impact = 0;
        double part_indiv = 0;
        double part_linteract = 0;
        double nabs_type = 0;
        double linside = 0;

        // TODO: missing correction due to m/m0 (but also wrong in xpart...)
        double energy = p0c_ref*ptau_in + e0_ref;  // energy in GeV

        double* result_scat = scatter(
                LocalParticle_get_x(part),
                LocalParticle_get_px(part)*rpp_in,
                LocalParticle_get_y(part),
                LocalParticle_get_py(part)*rpp_in,
                0,   // s at start of collimator
                energy,                   // energy, not momentum
                part_hit,
                part_abs,
                part_impact,         // impact parameter
                part_indiv,          // particle divergence
                part_linteract,      // interaction length
                nabs_type,
                linside,
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
                tilt0,
                tilt1,
                onesided,
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
                0    // s_imp
        );

        LocalParticle_set_x(part, result_scat[0]);
        LocalParticle_set_px(part, result_scat[1]/rpp_in);
        LocalParticle_set_y(part, result_scat[2]);
        LocalParticle_set_py(part, result_scat[3]/rpp_in);
//         double s_out      = result_scat[4];       //  Unused
        double energy_out = result_scat[5];       //  Cannot assign energy directly to LocalParticle as it would update dependent variables, but needs to be corrected first!
        part_hit          = result_scat[6];
        part_abs          = result_scat[7];
        part_impact       = result_scat[8];
        part_indiv        = result_scat[9];
        part_linteract    = result_scat[10];
        nabs_type         = result_scat[11];
        linside           = result_scat[12];
        p0                = result_scat[13];
        x0                = result_scat[14];
        xp0               = result_scat[15];
        nhit              = result_scat[16];
        nabs              = result_scat[17];
        fracab            = result_scat[18];
        nnuc0             = result_scat[19];
        ien0              = result_scat[20];
        nnuc1             = result_scat[21];
        ien1              = result_scat[22];

        // Backtrack to centre of collimator
        drift_4d_single(part, -length/2);
        // Return from closed orbit  (dpx = dxp because ref. particle has delta = 0)
        shift_4d_single(part, co_x, co_px/rpp_in, co_y, co_py/rpp_in );

        // Update energy    ---------------------------------------------------
        // Only particles that hit the jaw and survived need to be updated
        if (part_hit>0 && part_abs==0){
            double ptau_out = (energy_out - e0_ref) / (e0_ref * beta0_ref);
            LocalParticle_update_ptau(part, ptau_out);
        }

        // Drift to end of collimator
        drift_4d_single(part, length/2);

        // Update 4D coordinates    -------------------------------------------
        // Absorbed particles get their coordinates set to the entrance of collimator
        if (part_abs>0){
            LocalParticle_set_x(part, x_in);
            LocalParticle_set_px(part, px_in);
            LocalParticle_set_y(part, y_in);
            LocalParticle_set_py(part, py_in);
        }

        // Update longitudinal coordinate zeta    -----------------------------
        // Absorbed particles keep coordinates at the entrance of collimator, others need correcting:
        // Non-hit particles are just drifting (zeta not yet drifted in K2, so do here)
        if (part_hit==0){
            LocalParticle_add_to_zeta(part, drift_zeta_single(rvv_in, px_in*rpp_in, py_in*rpp_in, length) );
        }
        // Hit and survived particles need correcting:
        if (part_hit>0 && part_abs==0){
            double px  = LocalParticle_get_px(part);
            double py  = LocalParticle_get_py(part);
            double rvv = LocalParticle_get_rvv(part);
            double rpp = LocalParticle_get_rpp(part);
            // First we drift half the length with the old angles:
            LocalParticle_add_to_zeta(part, drift_zeta_single(rvv_in, px_in*rpp_in, py_in*rpp_in, length/2) );
            // then half the length with the new angles:
            LocalParticle_add_to_zeta(part, drift_zeta_single(rvv, px*rpp, py*rpp, length/2) );
        }

        // Update s    --------------------------------------------------------
        // TODO: move absorbed particles to last impact location
        if (part_abs==0){
            LocalParticle_add_to_s(part, length);
        }
 
        // Update state    ----------------------------------------------------
        if (part_abs > 0){
            LocalParticle_set_state(part, -333);
        }
    //end_per_particle_block
}

/*gpufun*/
void EverestCollimator_track_local_particle(EverestCollimatorData el, LocalParticle* part0) {
    int8_t const is_active      = EverestCollimatorData_get__active(el);
    double const inactive_front = EverestCollimatorData_get_inactive_front(el);
    double const active_length  = EverestCollimatorData_get_active_length(el);
    double const inactive_back  = EverestCollimatorData_get_inactive_back(el);

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
