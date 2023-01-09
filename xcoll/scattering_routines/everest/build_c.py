from cffi import FFI
from pathlib import Path

ffibuilder = FFI()

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

# Jaw function
ffibuilder.cdef("double* jaw(double run_exenergy, double run_anuc, double run_zatom, double run_rho, double run_radl, double ich_cprob0, double ich_cprob1, double ich_cprob2, double ich_cprob3, double ich_cprob4, double ich_cprob5, double run_xintl, double run_bn, double run_ecmsq, double run_xln15s, double run_bpp, double p0, double nabs, double s, double zlm, double x, double xp, double z, double zp, double dpop);")

# Crystal function
ffibuilder.cdef("double* crystal(double x, double xp, double z, double zp, double s, double p, double x0, double xp0, double zlm, double s_imp, double isimp, double val_part_hit, double val_part_abs, double val_part_impact, double val_part_indiv, double c_length, double exenergy, double rho, double anuc, double zatom, double emr, double dlri, double dlyi, double ai, double eum, double collnt, double hcut, double bnref, double csref0, double csref1, double csref5, double nhit, double nabs, double cry_tilt, double  cry_rcurv, double  cry_bend, double  cry_alayer, double  cry_xmax, double  cry_ymax, double  cry_orient, double  cry_miscut, double  cry_cBend, double  cry_sBend, double cry_cpTilt, double cry_spTilt, double  cry_cnTilt, double  cry_snTilt, double  iProc, double  n_chan, double  n_VR, double  n_amorphous);")


ffibuilder.set_source("_everest",  # name of the output C extension
"""
    #include "exponential_integral_Ei.h"
    #include "random.h"
    #include "jaw.h"
    #include "crystal.h"
""",
    libraries=['m'])    # on Unix, link with the math library

if __name__ == "__main__":
    ffibuilder.compile(verbose=True)

# clean up
file = Path.cwd() / "_everest.c"
file.unlink()
file = Path.cwd() / "_everest.o"
file.unlink()
