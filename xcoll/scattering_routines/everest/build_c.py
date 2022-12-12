from cffi import FFI
from pathlib import Path

ffibuilder = FFI()

# ffibuilder.cdef("double jaw(double run_exenergy, double run_anuc, double run_zatom, double run_rho, double run_radl, double run_cprob[], double run_xintl, double run_bn, double run_ecmsq, double run_xln15s, double run_bpp, double p0, double nabs, double s, double zlm, double x, double xp, double z, double zp, double dpop);")

# Random generator
ffibuilder.cdef("static void set_random_seed(unsigned int seed);")
ffibuilder.cdef("static unsigned int get_random_seed(void);")
ffibuilder.cdef("static unsigned int get_n_sampled(void);")
ffibuilder.cdef("static double get_random(void);")
ffibuilder.cdef("static double get_random_gauss(void);")
ffibuilder.cdef("static void set_rutherford_parameters(float z, float emr, float upper_val);")
ffibuilder.cdef("static void set_rutherford_iterations(unsigned int n);")
ffibuilder.cdef("static unsigned int get_rutherford_iterations(void);")
ffibuilder.cdef("static double get_random_ruth(void);")

# Jaw code
ffibuilder.cdef("double soln3(double a, double b, double dh, double smax, double s);")
ffibuilder.cdef("double* scamcs(double x0, double x0p, double s);")
ffibuilder.cdef("double* mcs(double s, double mc_radl, double mc_zlm1, double mc_p0, double mc_x, double mc_xp, double mc_z, double mc_zp, double mc_dpop);")
ffibuilder.cdef("double calcionloss(double p, double rlen, double il_exenergy, double il_anuc, double il_zatom, double il_rho, double enlo);")
ffibuilder.cdef("double* gettran(double inter, double p, double tt_bn, double tt_ecmsq, double tt_xln15s, double tt_bpp);")
ffibuilder.cdef("double* tetat(double t, double p, double tx, double tz);")
ffibuilder.cdef("int ichoix(double ich_cprob0, double ich_cprob1, double ich_cprob2, double ich_cprob3, double ich_cprob4, double ich_cprob5);")


ffibuilder.set_source("_everest",  # name of the output C extension
"""
    #include "exponential_integral_Ei.h"
    #include "random.h"
    #include "jaw.h"
""",
    libraries=['m'])    # on Unix, link with the math library

if __name__ == "__main__":
    ffibuilder.compile(verbose=True)

# clean up
file = Path.cwd() / "_everest.c"
file.unlink()
file = Path.cwd() / "_everest.o"
file.unlink()
