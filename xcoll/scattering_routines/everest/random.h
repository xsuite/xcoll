#include <stdlib.h>
#include <math.h>
#include <time.h>

unsigned int this_seed;
unsigned int n_sampled = 0;

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

// This is wrong if not manually set first: it will return 0, but not sure what the default seed is..
static unsigned int get_random_seed(){
    return this_seed;
}

static unsigned int get_n_sampled(){
    return n_sampled;
}

// Generate a random uniform value in [0,1]
// Note that Ranlux (in old K2) generated randoms in ]0,1[
static double get_random(void){
  n_sampled = n_sampled + 1;
  double x = (double)rand()/RAND_MAX;
  return x;
}







static double get_random_normal(void) {
    double u = ((double) rand() / (RAND_MAX)) * 2 - 1;
    double v = ((double) rand() / (RAND_MAX)) * 2 - 1;
    double r = u * u + v * v;
    if (r == 0 || r > 1) return get_random_normal();
    double c = sqrt(-2 * log(r) / r);
    return u * c;
}

// Generate a random value weighted within the normal (Gaussian) distribution
static double get_random_gauss(void){
  double x = get_random(),
         y = get_random(),
         z = sqrt(-2 * log(x)) * cos(2 * M_PI * y);
  return z;
}


// PDF of Rutherford distribution
static double ruth_PDF(float t, float z, float emr){

    // Constants
    double c = 0.8561e3;
    double N = 2.607e-5;
    double A = N*pow(z,2);
    double B = c*pow(emr,2);

    // CDF from Rutherford PDF
    double F = (A/pow(t,2))*(exp(-B*t));

    return F;
        
}

// CDF of Rutherford distribution
static double ruth_CDF(float t, float z, float emr, float t0){

    // Constants
    double c = 0.8561e-3;
    double N = 2.607e-5;
    double A = N*pow(z,2);
    double B = c*pow(emr,2);

    // CDF from Rutherford PDF
    double F = -A*B*Exponential_Integral_Ei(-B*t) - t*ruth_PDF(t,z,emr) + A*B*Exponential_Integral_Ei(-B*t0) +       t0*ruth_PDF(t0,z,emr);

    return F;
        
}

// Generate a random value weighted with a Rutherford distribution
static double get_random_ruth(float z, float emr, double lower_val, double upper_val, unsigned int n){

//     unsigned int n = 20; //iterations

    // corresponding range in uniform
    // double lower_uni = 0;
    double upper_uni = ruth_CDF(upper_val,z,emr,lower_val);

    // sample a random uniform
    double t = upper_uni * get_random() ;
    
    // initial estimate the lower border
    double x = lower_val;

    // solve CDF(x) == t for x
    unsigned short i = 1;
    while(i <= n) {
        x = x - (ruth_CDF(x,z,emr,lower_val)-t)/ruth_PDF(x,z,emr);
        i++;
    }

    return x;
}

