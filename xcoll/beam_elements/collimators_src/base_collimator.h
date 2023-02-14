// copyright ############################### #
// This file is part of the Xcoll Package.   #
// Copyright (c) CERN, 2023.                 #
// ######################################### #

#ifndef XCOLL_BASE_H
#define XCOLL_BASE_H

/*gpukern*/
void RandomRutherfordData_set_by_xcoll_material(RandomRutherfordData ran, GeneralMaterialData material){
    double const zatom    = GeneralMaterialData_get_Z(material);
    double const emr      = GeneralMaterialData_get_nuclear_radius(material);
    double const hcut     = GeneralMaterialData_get_hcut(material);
    double const lcut     = 0.0009982;
    double const c = 0.8561e3; // TODO: Where tha fuck does this come from??
    double A = pow(zatom,2);
    double B = c*pow(emr,2);
    RandomRutherfordData_set(ran, A, B, lcut, hcut);
}

/*gpufun*/
void BaseCollimator_track_local_particle(BaseCollimatorData el, LocalParticle* part0) {
    kill_all_particles(part0, xcoll_state_invalid_tracking);
}

#endif /* XCOLL_BASE_H */
