# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #


source = r'''
ifndef XCOLL_FLUKA_MASSES_H
define XCOLL_FLUKA_MASSES_H

// 1u = 931494102.42(28)    2018 CODATA

//       name                           mass [eV]          PDG ID
define  PROTON_MASS_EV_FLUKA       938272310.0000           2212   // best known value:     938272088.15(29)
define  Pb208_MASS_EV_FLUKA     193687690162.6481     1000822080   // best known value: 193.7290249(13)e9   (207.9766521(13) u)

define  DEUTERIUM_MASS_EV_FLUKA   1875613390.0000     1000010020   //
define  TRITIUM_MASS_EV_FLUKA     2808921774.1526     1000010030   //
define  He3_MASS_EV_FLUKA         2808392236.6048     1000020030   //
define  He4_MASS_EV_FLUKA         3727380256.9289     1000020040   //
define  He6_MASS_EV_FLUKA         5605538106.5358     1000020060   //
define  He8_MASS_EV_FLUKA         7482530706.1729     1000020080   //
define  Li6_MASS_EV_FLUKA         5601519470.4097     1000030060   //
define  Li7_MASS_EV_FLUKA         6533835190.3885     1000030070   //
define  Li8_MASS_EV_FLUKA         7471367020.2320     1000030080   //
define  Li11_MASS_EV_FLUKA       10285692859.7178     1000030110   //
define  Be7_MASS_EV_FLUKA         6534186187.0551     1000040070   //
define  Be10_MASS_EV_FLUKA        9325506337.1371     1000040100   //
define  B10_MASS_EV_FLUKA         9324439705.4834     1000050100   //
define  C12_MASS_EV_FLUKA        11174866883.2345     1000060120   //
define  N14_MASS_EV_FLUKA        13040208386.2960     1000070140   //
define  N15_MASS_EV_FLUKA        13968940806.3676     1000070150   //
define  O16_MASS_EV_FLUKA        14895086191.2700     1000080160   //
define  O17_MASS_EV_FLUKA        15830508491.1682     1000080170   //
define  Na24_MASS_EV_FLUKA       22341829488.0393     1000110240   //

endif /* XCOLL_FLUKA_MASSES_H */
'''

fluka_masses = {
    int(line.split()[3]): [line.split()[1].split('_')[0], float(line.split()[2])]
    for line in source.split('\n')
    if len(line.split()) > 1 and len(line.split()[1]) > 6 and line.split()[1][-6:] == '_FLUKA'
}
