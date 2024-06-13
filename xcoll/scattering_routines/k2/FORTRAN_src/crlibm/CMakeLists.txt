#CMakeLists.txt for crlibm

SET(CRLIBM1_SOURCES csh_fast.c rem_pio2.c trigo_fast.c tan.c sine.c cosine.c exp.c exp_fast.c log.c log_fast.c log10.c atan.c atan_fast.c)
SET(CRLIBM2_SOURCES addition_scs.c double2scs.c multiplication_scs.c division_scs.c scs2double.c zero_scs.c)

add_library(crlibm1 OBJECT ${CRLIBM1_SOURCES})
add_library(crlibm2 OBJECT ${CRLIBM2_SOURCES})

set(lessWarnings  -Wno-sign-conversion -Wno-conversion -Wno-long-long)

target_compile_options(crlibm1 PRIVATE -mfpmath=sse -msse2 -fPIC -std=c99 -Wall -Wshadow -Wpointer-arith -Wcast-align -Wconversion -Waggregate-return -Wstrict-prototypes -Wnested-externs -Wlong-long -Winline -pedantic -fno-strict-aliasing ${lessWarnings})
target_compile_options(crlibm2 PRIVATE -mfpmath=sse -msse2 -fPIC -Wall -Wshadow -Wpointer-arith -Wcast-align -Wconversion -Waggregate-return -Wstrict-prototypes -Wnested-externs -Wlong-long -Winline ${lessWarnings})

#glue the parts together
add_library(crlibm $<TARGET_OBJECTS:crlibm1> $<TARGET_OBJECTS:crlibm2>)

