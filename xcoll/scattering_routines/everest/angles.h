// copyright ############################### #
// This file is part of the Xcoll Package.   #
// Copyright (c) CERN, 2023.                 #
// ######################################### #

#ifndef XCOLL_EVEREST_ANGLES_H
#define XCOLL_EVEREST_ANGLES_H


/*gpufun*/
double LocalParticle_get_xp(LocalParticle* part){
    double const px = LocalParticle_get_px(part);
#ifndef XTRACK_USE_EXACT_DRIFTS
    double const rpp = LocalParticle_get_rpp(part);
#else
    double const py = LocalParticle_get_py(part);
    double const one_plus_delta = 1. + LocalParticle_get_delta(part);
    double const rpp = 1./sqrt(one_plus_delta*one_plus_delta - px*px - py*py);
#endif
    return px*rpp;    // TODO: this is not the angle, but sin(angle)
}

/*gpufun*/
double LocalParticle_get_yp(LocalParticle* part){
    double const py = LocalParticle_get_py(part);
#ifndef XTRACK_USE_EXACT_DRIFTS
    double const rpp = LocalParticle_get_rpp(part);
#else
    double const px = LocalParticle_get_px(part);
    double const one_plus_delta = 1. + LocalParticle_get_delta(part);
    double const rpp = 1./sqrt(one_plus_delta*one_plus_delta - px*px - py*py);
#endif
    return py*rpp;    // TODO: this is not the angle, but sin(angle)
}

/*gpufun*/
void LocalParticle_set_xp(LocalParticle* part, double xp){
    double const rpp = LocalParticle_get_rpp(part);
#ifdef XTRACK_USE_EXACT_DRIFTS
    // Careful! If yp also changes, use a different function!
    double const yp   = LocalParticle_get_yp(part);
    double const rpp *= sqrt(1 + xp*xp + yp*yp);
#endif
    LocalParticle_set_px(part, xp/rpp);    // TODO: xp is not the angle, but sin(angle)
}

/*gpufun*/
void LocalParticle_set_yp(LocalParticle* part, double yp){
    double const rpp = LocalParticle_get_rpp(part);
#ifdef XTRACK_USE_EXACT_DRIFTS
    // Careful! If xp also changes, use a different function!
    double const xp   = LocalParticle_get_xp(part);
    double const rpp *= sqrt(1 + xp*xp + yp*yp);
#endif
    LocalParticle_set_py(part, yp/rpp);    // TODO: yp is not the angle, but sin(angle)
}

/*gpufun*/
void LocalParticle_set_xp_yp(LocalParticle* part, double xp, double yp){
    double const rpp = LocalParticle_get_rpp(part);
#ifdef XTRACK_USE_EXACT_DRIFTS
    double const rpp *= sqrt(1 + xp*xp + yp*yp);
#endif
    LocalParticle_set_px(part, xp/rpp);    // TODO: xp is not the angle, but sin(angle)
    LocalParticle_set_py(part, yp/rpp);    // TODO: yp is not the angle, but sin(angle)
}

/*gpufun*/
void LocalParticle_add_to_xp(LocalParticle* part, double xp){
    LocalParticle_set_xp(part, LocalParticle_get_xp(part) + xp);
}

/*gpufun*/
void LocalParticle_add_to_yp(LocalParticle* part, double yp){
    LocalParticle_set_yp(part, LocalParticle_get_yp(part) + yp);
}

/*gpufun*/
void LocalParticle_add_to_xp_yp(LocalParticle* part, double xp, double yp){
    LocalParticle_set_xp_yp(part, LocalParticle_get_xp(part) + xp, LocalParticle_get_yp(part) + yp);
}

#endif /* XCOLL_EVEREST_ANGLES_H */
