#include <stdlib.h>
#include <math.h>
#include <time.h>

unsigned int this_seed;

static void set_random_seed(unsigned int seed){
    /* Intializes random number generator */
    if (seed == 0) {
        time_t t;
        this_seed = (unsigned) time(&t);
    } else {
        this_seed = seed;
    }
    srand(this_seed);
}

// This wrong if not manually set first: it will return 0, but not sure what the default seed is..
static unsigned int get_random_seed(){
    return this_seed;
}

// Generate a random uniform value in [0,1]
// Note that Ranlux (in old K2) generated randoms in ]0,1[
static double get_random(void){
  double x = (double)rand() / RAND_MAX;
  return x;
}

// Generate a random value weighted within the normal (Gaussian) distribution
static double get_random_gauss(void){
  double x = get_random(),
         y = get_random(),
         z = sqrt(-2 * log(x)) * cos(2 * M_PI * y);
  return z;
}

// PDF of Rutherford distribution
static double PDF(float t){

    // Constants
    double z = 29;
    double c = 0.8561e3;
    double N = 2.607e-5;
    double e = 0.302;
    double A = N*pow(z,2);
    double B = c*pow(e,2);

    // CDF from Rutherford PDF
    double F = (A/pow(t,2))*(exp(-B*t));

    return F;
        
}

// CDF of Rutherford distribution
static double CDF(float t, float t0){

    // Constants
    double z = 29;
    double c = 0.8561e-3;
    double N = 2.607e-5;
    double e = 0.302;
    double A = N*pow(z,2);
    double B = c*pow(e,2);

    // CDF from Rutherford PDF
    double F = -A*B*Exponential_Integral_Ei(-B*t) - t*PDF(t) + A*B*Exponential_Integral_Ei(-B*t0) + t0*PDF(t0);

    return F;
        
}

// Generate a random value weighted with a Rutherford distribution
static double get_random_ruth(int n){

    // range of rutherford random
    double lower_CDF = 0.0009982;
    double upper_CDF = 0.01;

    // corresponding range in uniform
    // double lower_uni = 0;
    double upper_uni = CDF(upper_CDF, lower_CDF);

    // sample a random uniform
    double t = upper_uni * get_random() ;
    
    // initial estimate the lower border
    double x = lower_CDF;

    // solve CDF(x) == t for x
    unsigned short i = 1;
    while(i <= n) {
        x = x - (CDF(x, lower_CDF)-t)/PDF(x);
        i++;
    }

    return x;
}

