// copyright ############################### #
// This file is part of the Xcoll package.   #
// Copyright (c) CERN, 2024.                 #
// ######################################### #

#ifndef XCOLL_EVEREST_ENGINE_H
#define XCOLL_EVEREST_ENGINE_H

#define XCOLL_TRANSITION

typedef struct EverestCollData_ {
    // Collimator properties
    RandomRutherfordData restrict rng;
    InteractionRecordData record;
    RecordIndex record_index;
    int8_t record_scatterings;
    // Material properties
    // TODO: can we use pointers for the MaterialData? It then gets a bit difficult to read them, ie *coll->exenergy
    double exenergy;
    double rho;
    double anuc;
    double zatom;
    double bnref;
    double csref[6];
    double radl;  // TODO: is this the same physically as collnt or dlri ?
    double dlri;
    double dlyi;
    double ai;
    double eum;
    double collnt;
    double eta;
    int8_t orient;
    int8_t only_mcs;
} EverestCollData_;
typedef EverestCollData_ *EverestCollData;

typedef struct EverestData_ {
    EverestCollData coll;
    // Dynamic parameters
    double cprob[6];
    double xintl;
    double bn;
    double ecmsq;
    double xln15s;
    double bpp;
    double prob_tail_c1;
    double prob_tail_c2;
    double prob_tail_c3;
    double prob_tail_c4;
    double energy_loss;
    double energy_loss_tail;
    // Crystal data
    double rescale_scattering;
    double t_c;
    double t_c0;
    double Rc_over_R;
    double Ang_rms;
    double Ang_avr;
    double Vcapt;
    double t_I;
    double t_P;
    double r;
} EverestData_;
typedef EverestData_ *EverestData;


/*gpufun*/
double LocalParticle_get_energy(LocalParticle* part){
    double mass_ratio = LocalParticle_get_charge_ratio(part) / LocalParticle_get_chi(part);
    return (LocalParticle_get_ptau(part)*LocalParticle_get_p0c(part) \
            + LocalParticle_get_energy0(part)) * mass_ratio;
}


/*gpufun*/
double drift_zeta_single(double rvv, double xp, double yp, double length){
    double const rv0v = 1./rvv;
    double const dzeta = 1 - rv0v * (1. + (pow(xp,2.) + pow(yp,2.))/2.);
    return length * dzeta;
}

/*gpufun*/
void Drift_single_particle_4d(LocalParticle* part, double length){
    double zeta = LocalParticle_get_zeta(part);
    Drift_single_particle(part, length);
    LocalParticle_set_zeta(part, zeta);
}

/*gpukern*/
void RandomRutherford_set_by_xcoll_material(RandomRutherfordData ran, GeneralMaterialData material){
    if (GeneralMaterialData_get__only_mcs(material)){
        RandomRutherford_set(ran, 1, 1, 0.0001, 0.01);
        return;
    }
    double const zatom    = GeneralMaterialData_get_Z(material);
    double const emr      = GeneralMaterialData_get_nuclear_radius(material);
    double const hcut     = GeneralMaterialData_get_hcut(material);
    double const lcut     = 0.0009982;
    double const c = 0.8561e3; // TODO: Where tha fuck does this come from??
    double A = pow(zatom,2);
    double B = c*pow(emr,2);
    RandomRutherford_set(ran, A, B, lcut, hcut);
}

#endif /* XCOLL_EVEREST_ENGINE_H */