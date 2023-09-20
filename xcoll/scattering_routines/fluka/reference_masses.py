# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2023.                 #
# ######################################### #


source = r'''
#ifndef XCOLL_FLUKA_MASSES_H
#define XCOLL_FLUKA_MASSES_H

// 1u = 931494102.42(28)    2018 CODATA

//       name                              mass [eV]        PDG ID
#define  PROTON_MASS_EV_FLUKA          938272310.000          2212   // best known value:     938272088.15(29)
#define  Pb208_MASS_EV_FLUKA        193687690162.648    1000822080   // best known value: 193.7290249(13)e9   (207.9766521(13) u)
                                

#endif /* XCOLL_FLUKA_MASSES_H */
'''

masses = {
    int(line.split()[3]): [line.split()[1].split('_')[0], float(line.split()[2])]
    for line in source.split('\n')
    if len(line.split()) > 1 and len(line.split()[1]) > 6 and line.split()[1][-6:] == '_FLUKA'
}
