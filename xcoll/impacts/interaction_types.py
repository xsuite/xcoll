# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2023.                 #
# ######################################### #


source = r'''
#ifndef XCOLL_INTERACTIONS_H
#define XCOLL_INTERACTIONS_H

#define  XC_ABSORBED                      -1        // point (no children)
#define  XC_MULTIPLE_COULOMB_SCATTERING    1        // continous
#define  XC_PN_ELASTIC                     2        // point (no children)
#define  XC_PP_ELASTIC                     3        // point (no children)
#define  XC_SINGLE_DIFFRACTIVE             5        // point (no children)
#define  XC_COULOMB                        6        // point (no children)
#define  XC_RUTHERFORD                     7        // point (no children)
#define  XC_CRYSTAL_CHANNELING           100        // continous
#define  XC_CRYSTAL_VOLUME_REFLECTION    101        // point (no children)
#define  XC_CRYSTAL_VOLUME_CAPTURE       102        // point (no children)
#define  XC_CRYSTAL_DECHANNELING         103        // point (no children)
#define  XC_CRYSTAL_AMORPHOUS            104        // continous
#define  XC_CRYSTAL_DRIFT                105        // continous

#endif /* XCOLL_INTERACTIONS_H */
'''

#int const proc_TRVR        = 100;     // Volume reflection in VR-AM transition region
#int const proc_TRAM        = 101;     // Amorphous in VR-AM transition region


interactions = {
    int(line.split()[2]): line.split()[1][3:].replace('_',' ').title()
    for line in source.split('\n')
    if len(line.split()) > 1 and line.split()[1][:3] == 'XC_' # select the source lines with the definitions
}
