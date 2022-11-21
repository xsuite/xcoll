from cffi import FFI
from pathlib import Path

ffibuilder = FFI()

# ffibuilder.cdef("float iterat(float a, float b, float dh, float s);")
ffibuilder.cdef("double soln3(double a, double b, double dh, double smax, double s);")
ffibuilder.cdef("static double get_random(void);")
ffibuilder.cdef("static double get_random_gauss(void);")
ffibuilder.cdef("static double get_random_ruth(int n);")
ffibuilder.cdef("static void set_random_seed(unsigned int seed);")
ffibuilder.cdef("static unsigned int get_random_seed();")

ffibuilder.set_source("_everest",  # name of the output C extension
"""
    #include "jaw.h"
    #include "exponential_integral_Ei.h"
    #include "random.h"
""",
    libraries=['m'])    # on Unix, link with the math library

if __name__ == "__main__":
    ffibuilder.compile(verbose=True)

# clean up
file = Path.cwd() / "_everest.c"
file.unlink()
file = Path.cwd() / "_everest.o"
file.unlink()
