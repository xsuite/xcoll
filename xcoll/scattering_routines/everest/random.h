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
double get_random(LocalParticle* part){
  double r = LocalParticle_generate_random_double(part);
  return r;
}


/*
=========================================
======  Gaussian Random Generator  ======
=========================================
*/

// Generate a random value weighted within the normal (Gaussian) distribution
double get_random_gauss(LocalParticle* part){
  double r = LocalParticle_generate_random_double_gauss(part);
  return r;
}


/*
=========================================
=====  Rutherford Random Generator  =====
=========================================
*/

/*
    Iterations for Newton's method
*/
// unsigned int n_iter = 0;   // automatic
unsigned int n_iter = 7;   // good enough precision
// TODO: how to optimise Newton's method??

void set_rutherford_iterations(unsigned int n){
    n_iter = n;
}
unsigned int get_rutherford_iterations(void){
    return n_iter;
}


/*
    Rutherford parameters
*/
float ruth_lower_val = 0.0009982;
float ruth_upper_val;
float ruth_A;
float ruth_B;
double ruth_CDF(float t);
double ruth_PDF(float t);

void set_rutherford_parameters(float z, float emr, float upper_val){
    double c = 0.8561e3; // TODO: Where tha fuck does this come from??
    ruth_A = pow(z,2);
    ruth_B = c*pow(emr,2);
    ruth_upper_val = upper_val;

    // Normalise PDF
    double N = ruth_CDF(ruth_upper_val);
    ruth_A = pow(z,2) / N;
}


// PDF of Rutherford distribution
double ruth_PDF(float t){
    return (ruth_A/pow(t,2))*(exp(-ruth_B*t));
}

// CDF of Rutherford distribution
double ruth_CDF(float t){
    return - ruth_A*ruth_B*Exponential_Integral_Ei(-ruth_B*t) - t*ruth_PDF(t)
           + ruth_A*ruth_B*Exponential_Integral_Ei(-ruth_B*ruth_lower_val) + ruth_lower_val*ruth_PDF(ruth_lower_val);
        
}

// Generate a random value weighted with a Rutherford distribution
double get_random_ruth(LocalParticle* part){

    // sample a random uniform
    double t = get_random(part);
    
    // initial estimate the lower border
    double x = ruth_lower_val;

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
    unsigned short i = 1;
    while(i <= n_iter) {
        x = x - (ruth_CDF(x)-t)/ruth_PDF(x);
        i++;
    }

    return x;
}

