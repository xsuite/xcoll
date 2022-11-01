from ctypes import *
so_file = "/Users/macbookpro/Documents/CERN/xsuite/xcoll/xcoll/scattering_routines/everest/my_functions.so"
my_functions = CDLL(so_file)

print(type(my_functions))

print(my_functions.square(10))