from cffi import FFI
from pathlib import Path

ffibuilder = FFI()

# ffibuilder.cdef("float iterat(float a, float b, float dh, float s);")
ffibuilder.cdef("float soln3(float a, float b, float dh, float smax, float s);")

ffibuilder.set_source("_jaw",  # name of the output C extension
"""
    #include "jaw.h"
""",
    libraries=['m'])    # on Unix, link with the math library

if __name__ == "__main__":
    ffibuilder.compile(verbose=True)

# clean up
file = Path.cwd() / "_jaw.c"
file.unlink()
file = Path.cwd() / "_jaw.o"
file.unlink()
