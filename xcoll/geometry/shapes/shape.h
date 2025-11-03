// copyright ############################### #
// This file is part of the Xcoll package.   #
// Copyright (c) CERN, 2025.                 #
// ######################################### #

/*gpufun*/
void Shape2D_crossing_drift(Shape2D shape, FindRoot finder, double s, double x, double px){
    DriftTrajectory_  _dri;   // Assigns memory
    DriftTrajectory dri = (DriftTrajectory)&_dri; // safe: DriftTrajectory_ has same layout as DriftTrajectory
    DriftTrajectory_set_params(dri, s, x, px);
    int64_t n_segments = Shape2D_len_segments(shape);
    for (int8_t i=0; i<n_segments;i++) {
        LocalSegment seg = Shape2D_getp1_segments(shape, i);
        FindRoot_find_crossing(finder, (LocalSegment) seg, (LocalTrajectory) dri);
    }
}


/*gpufun*/
void Shape2D_crossing_mcs(Shape2D shape, FindRoot finder, double s, double x, double px, ...){
    MultipleCoulombTrajectory_  _mcs;   // Assigns memory
    MultipleCoulombTrajectory mcs = (MultipleCoulombTrajectory)&_mcs; // safe: MultipleCoulombTrajectory_ has same layout as MultipleCoulombTrajectory
    MultipleCoulombTrajectory_set_params(mcs, s, x, px, ...);  // TODO
    int64_t n_segments = Shape2D_len_segments(shape);
    for (int8_t i=0; i<n_segments;i++) {
        LocalSegment seg = Shape2D_getp1_segments(shape, i);
        FindRoot_find_crossing(finder, (LocalSegment) seg, (LocalTrajectory) mcs);
    }
}


EverestCollimator_track_local_particle(EverestCollimatorData el, LocalParticle part){
    double length = EverestCollimator_get_length(el);
    Shape2D shape = EverestCollimator_get_shape(el);
    Material material = EverestCollimator_get_material(el);

    // Temporary FindRoot object
    FindRoot_  _finder;   // Assigns memory
    FindRoot finder = (FindRoot)&_finder; // safe: FindRoot_ has same layout as FindRoot

    // Find hit on shape
    double s = LocalParticle_get_s(part);
    double x = LocalParticle_get_x(part);
    double px = LocalParticle_get_exact_xp(part);
    Shape2D_crossing_drift(shape, finder, s, x, px);
    double s_cross = FindRoot_get_first_solution_s(finder);
    if (s_cross > XC_GEOM_LARGE_VALUE){
        // No crossing found
        Drift_single_particle(part, length);
        return;
    } else {
        Drift_single_particle(part, s_cross - s);

        // Move in material
        // Find exit point
        // get MCS parameters + random numbers for this material
        double s = LocalParticle_get_s(part);
        double x = LocalParticle_get_x(part);
        double px = LocalParticle_get_exact_xp(part);
        Shape2D_crossing_mcs(shape, finder, s, x, px, ...);
        double s_cross = FindRoot_get_next_solution_s(finder);
        double l_cross = FindRoot_get_path_length(finder, solution_index);

        // Find distance to next nuclear interaction
        double l_nucl = find_nuclear_interaction_point(material, part);

        if (l_nucl < l_cross){
            // Nuclear interaction occurs before exit
            apply_mcs(l_nucl, part, ...); // drift + kicks + ionisation loss
            next step
        } else {
            // Exit occurs before nuclear interaction
            apply_mcs(l_cross, part, ...); // drift + kicks + ionisation loss
            drift to end position
        }

}
