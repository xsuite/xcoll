#ifndef XCOLL_EVEREST_RAN_H
#define XCOLL_EVEREST_RAN_H
#include <stdlib.h>
#include <math.h>
#include <time.h>

/*
=========================================
========  Main Random Generator  ========
=========================================
*/

// Generate a random uniform value in [0,1]
// Note that Ranlux (in old K2) generated randoms in ]0,1[
/*gpufun*/
double get_random(LocalParticle* part){
  double r = LocalParticle_generate_random_double(part);
  return r;
}

/*gpufun*/
void EverestRandomData_get_random(ParticlesData particles, ArrNFloat64 samples, int n){
    LocalParticle thispart;
    Particles_to_LocalParticle(particles, &thispart, 0);
    LocalParticle* part0 = &thispart;

    //start_per_particle_block (part0->part)
    int i;
    for (i=0; i<n; ++i){
        ArrNFloat64_set(samples, n*LocalParticle_get_particle_id(part)+ i, LocalParticle_generate_random_double(part));
    }
    //end_per_particle_block
}


/*
=========================================
======  Gaussian Random Generator  ======
=========================================
*/

// Generate a random value weighted within the normal (Gaussian) distribution
/*gpufun*/
double get_random_gauss(LocalParticle* part){
  double r = LocalParticle_generate_random_double_gauss(part);
  return r;
}

/*gpufun*/
void EverestRandomData_get_random_gauss(ParticlesData particles, ArrNFloat64 samples, int n){
    LocalParticle thispart;
    Particles_to_LocalParticle(particles, &thispart, 0);
    LocalParticle* part0 = &thispart;

    //start_per_particle_block (part0->part)
    int i;
    for (i=0; i<n; ++i){
        ArrNFloat64_set(samples, n*LocalParticle_get_particle_id(part)+ i, LocalParticle_generate_random_double_gauss(part));
    }
    //end_per_particle_block
}

/*gpufun*/
void EverestRandomData_get_random_exp(ParticlesData particles, ArrNFloat64 samples, int n){
    LocalParticle thispart;
    Particles_to_LocalParticle(particles, &thispart, 0);
    LocalParticle* part0 = &thispart;

    //start_per_particle_block (part0->part)
    int i;
    for (i=0; i<n; ++i){
        ArrNFloat64_set(samples, n*LocalParticle_get_particle_id(part)+ i, LocalParticle_generate_random_double_exp(part));
    }
    //end_per_particle_block
}

/*
=========================================
=====  Rutherford Random Generator  =====
=========================================
*/

// TODO: how to optimise Newton's method??

// PDF of Rutherford distribution
/*gpufun*/
double ruth_PDF(double t, double A, double B){
    return (A/pow(t,2))*(exp(-B*t));
}

// CDF of Rutherford distribution
/*gpufun*/
double ruth_CDF(double t, double A, double B, double t0){
    return A*B*Exponential_Integral_Ei(-B*t0) + t0*ruth_PDF(t0, A, B)
         - A*B*Exponential_Integral_Ei(-B*t)  - t*ruth_PDF(t, A, B);
        
}

/*gpufun*/
void EverestRandomData_set_rutherford(EverestRandomData ran, double z, double emr, double upper_val){
    double c = 0.8561e3; // TODO: Where tha fuck does this come from??
    double A = pow(z,2);
    double B = c*pow(emr,2);
    double lower_val = EverestRandomData_get_rutherford_lower_val(ran);

    // Normalise PDF
    double N = ruth_CDF(upper_val, A, B, lower_val);
    EverestRandomData_set_rutherford_A(ran, A/N);
    EverestRandomData_set_rutherford_B(ran, B);
    EverestRandomData_set_rutherford_upper_val(ran, upper_val);
}

// Generate a random value weighted with a Rutherford distribution
/*gpufun*/
double get_random_ruth(EverestRandomData ran, LocalParticle* part){

    // get the parameters
    double x0     = EverestRandomData_get_rutherford_lower_val(ran);
    int8_t n_iter = EverestRandomData_get_rutherford_iterations(ran);
    double A      = EverestRandomData_get_rutherford_A(ran);
    double B      = EverestRandomData_get_rutherford_B(ran);

    // sample a random uniform
    double t = get_random(part);

    // initial estimate is the lower border
    double x = x0;

    // HACK to let iterations depend on sample to improve speed
    // based on Berylium being worst performing and hcut as in materials table
    // DOES NOT WORK
//     if (n_iter==0){
//         if (t<0.1) {
//             n_iter = 3;
//         } else if (t<0.35) {
//             n_iter = 4;
//         } else if (t<0.63) {
//             n_iter = 5;
//         } else if (t<0.8) {
//             n_iter = 6;
//         } else if (t<0.92) {
//             n_iter = 7;
//         } else if (t<0.96) {
//             n_iter = 8;
//         } else if (t<0.98) {
//             n_iter = 9;
//         } else {
//             n_iter = 10;
//         }
//     }

    // solve CDF(x) == t for x
    int8_t i = 1;
    while(i <= n_iter) {
        x = x - (ruth_CDF(x, A, B, x0)-t)/ruth_PDF(x, A, B);
        i++;
    }

    return x;
}


/*gpufun*/
void EverestRandomData_get_random_ruth(EverestRandomData ran, ParticlesData particles, ArrNFloat64 samples, int n){
    LocalParticle thispart;
    Particles_to_LocalParticle(particles, &thispart, 0);
    LocalParticle* part0 = &thispart;

    //start_per_particle_block (part0->part)
    int i;
    for (i=0; i<n; ++i){
        ArrNFloat64_set(samples, n*LocalParticle_get_particle_id(part)+ i, get_random_ruth(ran, part));
    }
    //end_per_particle_block
}

#endif
