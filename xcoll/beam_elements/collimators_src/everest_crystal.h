#ifndef XCOLL_EVEREST_CRYSTAL_H
#define XCOLL_EVEREST_CRYSTAL_H

#include <math.h>

// TODO:
//    Do not split 4d and zeta in drifts
//    Use drift function from xtrack Drift element (call its C function)

/*gpufun*/
void cry_drift_6d(LocalParticle* part0, double length) {
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
void cry_drift_4d_single(LocalParticle* part, double length) {
    double const rpp    = LocalParticle_get_rpp(part);
    double const xp     = LocalParticle_get_px(part) * rpp;
    double const yp     = LocalParticle_get_py(part) * rpp;
    LocalParticle_add_to_x(part, xp * length );
    LocalParticle_add_to_y(part, yp * length );
}

/*gpufun*/
double cry_drift_zeta_single(double rvv, double xp, double yp, double length){
    double const rv0v = 1./rvv;
    double const dzeta = 1 - rv0v * (1. + (pow(xp,2.) + pow(yp,2.))/2.);
    return length * dzeta;
}

/*gpufun*/
void cry_shift_4d_single(LocalParticle* part, double dx, double dpx, double dy, double dpy) {
    LocalParticle_add_to_x(part, dx);
    LocalParticle_add_to_px(part, dpx);
    LocalParticle_add_to_y(part, dy);
    LocalParticle_add_to_py(part, dpy);
}


// TODO: 
// Write impacts



/*gpufun*/
void track_crystal(EverestCrystalData el, LocalParticle* part0) {

    double const e0_ref    = LocalParticle_get_energy0(&part0[0]) / 1e9; // Reference energy in GeV
    double const p0c_ref   = LocalParticle_get_p0c(&part0[0]) / 1e9;
    double const beta0_ref = LocalParticle_get_beta0(&part0[0]);

    // Crystal properties
    double length  = EverestCrystalData_get_active_length(el);
    // if collimator.jaw_F_L != collimator.jaw_B_L or collimator.jaw_F_R != collimator.jaw_B_R:
    //     raise NotImplementedError
    double const opening  = EverestCrystalData_get_jaw_F_L(el) - EverestCrystalData_get_jaw_F_R(el);
    double const offset   = EverestCrystalData_get_offset(el) + ( EverestCrystalData_get_jaw_F_L(el) + EverestCrystalData_get_jaw_F_R(el) )/2;
    double const tilt0    = EverestCrystalData_get_tilt(el, 0);
    double const tilt1    = EverestCrystalData_get_tilt(el, 1);
    double const onesided = EverestCrystalData_get_onesided(el);
    double const angle    = atan2(EverestCrystalData_get_sin_z(el), EverestCrystalData_get_cos_z(el) );
    double const co_x     = EverestCrystalData_get_dx(el);
    double const co_px    = EverestCrystalData_get_dpx(el);
    double const co_y     = EverestCrystalData_get_dy(el);
    double const co_py    = EverestCrystalData_get_dpy(el);
    double const bend     = EverestCrystalData_get_bend(el);
    double const cry_tilt = EverestCrystalData_get_align_angle(el) + EverestCrystalData_get_crytilt(el);
    double const bend_ang = length/bend;    // temporary value
    if (cry_tilt >= -bend_ang) {
        length = bend*(sin(bend_ang + cry_tilt) - sin(cry_tilt));
    } else {
        length = bend*(sin(bend_ang - cry_tilt) + sin(cry_tilt));
    }
    double const cry_rcurv  = bend;
    double const cry_bend   = length/cry_rcurv; //final value (with corrected length) 
    double const cry_alayer = EverestCrystalData_get_thick(el);
    double const cry_xmax   = EverestCrystalData_get_xdim(el);
    double const cry_ymax   = EverestCrystalData_get_ydim(el);
    double const cry_orient = EverestCrystalData_get_orient(el);
    double const cry_miscut = EverestCrystalData_get_miscut(el);
    double const cry_cBend  = cos(cry_bend); 
    double const cry_sBend  = sin(cry_bend); 
    double const cry_cpTilt = cos(cry_tilt);
    double const cry_spTilt = sin(cry_tilt);
    double const cry_cnTilt = cos(-cry_tilt);
    double const cry_snTilt = sin(-cry_tilt);

    // Material properties
    CrystalMaterialData material = EverestCrystalData_getp_material(el);
    double const zatom    = CrystalMaterialData_get_Z(material);
    double const anuc     = CrystalMaterialData_get_A(material);
    double const rho      = CrystalMaterialData_get_density(material);
    double const exenergy = CrystalMaterialData_get_excitation_energy(material);
    double const emr      = CrystalMaterialData_get_nuclear_radius(material);
    double const bnref    = CrystalMaterialData_get_nuclear_elastic_slope(material);
    double const csref0   = CrystalMaterialData_get_cross_section(material, 0);
    double const csref1   = CrystalMaterialData_get_cross_section(material, 1);
    double const csref5   = CrystalMaterialData_get_cross_section(material, 5);
    double const hcut     = CrystalMaterialData_get_hcut(material);
    double const dlri     = CrystalMaterialData_get_crystal_radiation_length(material);
    double const dlyi     = CrystalMaterialData_get_crystal_nuclear_length(material);
    double const eUm      = CrystalMaterialData_get_crystal_potential(material);
    double const ai       = CrystalMaterialData_get_crystal_plane_distance(material);
    double const collnt   = CrystalMaterialData_get_nuclear_collision_length(material);

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
    // specifically for crystals:
    double iProc  = 0;
    double n_chan = 0;
    double n_VR   = 0;
    double n_amorphous = 0;
    double s_imp  = 0;

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
        cry_drift_4d_single(part, length/2);
        // Move to closed orbit  (dpx = dxp because ref. particle has delta = 0)
        cry_shift_4d_single(part, -co_x, -co_px/rpp_in, -co_y, -co_py/rpp_in );
        // Backtrack to start of collimator
        cry_drift_4d_single(part, -length/2);

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
                0,   // radl not used
                dlri,
                dlyi,
                eUm,
                ai,
                collnt,
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
                1,   // is_crystal
                length,
                angle,
                opening,
                offset,
                tilt0,
                tilt1,
                onesided,
                cry_tilt, 
                cry_rcurv,  
                cry_bend,
                cry_alayer, 
                cry_xmax,   
                cry_ymax,   
                cry_orient, 
                cry_miscut, 
                cry_cBend,  
                cry_sBend,  
                cry_cpTilt, 
                cry_spTilt, 
                cry_cnTilt, 
                cry_snTilt, 
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
                iProc,
                n_chan,
                n_VR,
                n_amorphous,
                s_imp
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
        // p0                = result_scat[13];
        // x0                = result_scat[14];
        // xp0               = result_scat[15];
        // nhit              = result_scat[16];
        // nabs              = result_scat[17];
        // fracab            = result_scat[18];
        // nnuc0             = result_scat[19];
        // ien0              = result_scat[20];
        // nnuc1             = result_scat[21];
        // ien1              = result_scat[22];
        // iProc             = result_scat[23];
        // n_chan            = result_scat[24];
        // n_VR              = result_scat[25];
        // n_amorphous       = result_scat[26];
        // s_imp             = result_scat[27];

        // Backtrack to centre of collimator
        cry_drift_4d_single(part, -length/2);
        // Return from closed orbit  (dpx = dxp because ref. particle has delta = 0)
        cry_shift_4d_single(part, co_x, co_px/rpp_in, co_y, co_py/rpp_in );

        // Update energy    ---------------------------------------------------
        // Only particles that hit the jaw and survived need to be updated
        if (part_hit>0 && part_abs==0){
            double ptau_out = (energy_out - e0_ref) / (e0_ref * beta0_ref);
            LocalParticle_update_ptau(part, ptau_out);
        }

        // Drift to end of collimator
        cry_drift_4d_single(part, length/2);

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
            LocalParticle_add_to_zeta(part, cry_drift_zeta_single(rvv_in, px_in*rpp_in, py_in*rpp_in, length) );
        }
        // Hit and survived particles need correcting:
        if (part_hit>0 && part_abs==0){
            double px  = LocalParticle_get_px(part);
            double py  = LocalParticle_get_py(part);
            double rvv = LocalParticle_get_rvv(part);
            double rpp = LocalParticle_get_rpp(part);
            // First we drift half the length with the old angles:
            LocalParticle_add_to_zeta(part, cry_drift_zeta_single(rvv_in, px_in*rpp_in, py_in*rpp_in, length/2) );
            // then half the length with the new angles:
            LocalParticle_add_to_zeta(part, cry_drift_zeta_single(rvv, px*rpp, py*rpp, length/2) );
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
void EverestCrystal_track_local_particle(EverestCrystalData el, LocalParticle* part0) {
    int8_t const is_active      = EverestCrystalData_get__active(el);
    double const inactive_front = EverestCrystalData_get_inactive_front(el);
    double const active_length  = EverestCrystalData_get_active_length(el);
    double const inactive_back  = EverestCrystalData_get_inactive_back(el);

    if (!is_active){
        // Drift full length
        cry_drift_6d(part0, inactive_front+active_length+inactive_back);
    } else {
        cry_drift_6d(part0, inactive_front);
        track_crystal(el, part0);
        cry_drift_6d(part0, inactive_back);
    }
}

#endif
