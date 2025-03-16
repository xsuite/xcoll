# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #


source = r'''
#ifndef XCOLL_FLUKA_MASSES_H
#define XCOLL_FLUKA_MASSES_H

// 1u = 931494102.42(28)    2018 CODATA

//       name                            mass [eV]         PDG ID
#define  ELECTRON_MASS_EV_FLUKA        510998.9506918          11   // best known value:        510998.95069(16)
#define  MUON_MASS_EV_FLUKA         105658389.0                13   //
#define  PION_MASS_EV_FLUKA         139569950.0000            211   //
#define  KAON_MASS_EV_FLUKA         493646000.0000            321   //
#define  PROTON_MASS_EV_FLUKA       938272310.0000           2212   // best known value:     938272088.15(29)
#define  DEUTERIUM_MASS_EV_FLUKA   1875613390.0000     1000010020   //
#define  TRITIUM_MASS_EV_FLUKA     2808921774.1526     1000010030   //
#define  He3_MASS_EV_FLUKA         2808392236.6048     1000020030   //
#define  He4_MASS_EV_FLUKA         3727380256.9289     1000020040   //
#define  He6_MASS_EV_FLUKA         5605538106.5358     1000020060   //
#define  He8_MASS_EV_FLUKA         7482530706.1729     1000020080   //
#define  Li6_MASS_EV_FLUKA         5601519470.4097     1000030060   //
#define  Li7_MASS_EV_FLUKA         6533835190.3885     1000030070   //
#define  Li8_MASS_EV_FLUKA         7471367020.2320     1000030080   //
#define  Li9_MASS_EV_FLUKA         8406870100.1281     1000030090   //
#define  Li11_MASS_EV_FLUKA       10285692859.7178     1000030110   //
#define  Be7_MASS_EV_FLUKA         6534186187.0551     1000040070   //
#define  Be9_MASS_EV_FLUKA         8392753017.1697     1000040090   //
#define  Be10_MASS_EV_FLUKA        9325506337.1371     1000040100   //
#define  Be11_MASS_EV_FLUKA       10264568076.9410     1000040110   //
#define  B8_MASS_EV_FLUKA          7472321245.2018     1000050080   //
#define  B10_MASS_EV_FLUKA         9324439705.4834     1000050100   //
#define  B11_MASS_EV_FLUKA        10252551235.5711     1000050110   //
#define  B12_MASS_EV_FLUKA        11188746425.4493     1000050120   //
#define  C11_MASS_EV_FLUKA        10254022852.9586     1000060110   //
#define  C12_MASS_EV_FLUKA        11174866883.2345     1000060120   //
#define  C13_MASS_EV_FLUKA        12109486203.1536     1000060130   //
#define  C14_MASS_EV_FLUKA        13040875423.1563     1000060140   //
#define  N14_MASS_EV_FLUKA        13040208386.2960     1000070140   //
#define  N15_MASS_EV_FLUKA        13968940806.3676     1000070150   //
#define  O16_MASS_EV_FLUKA        14895086191.2700     1000080160   //
#define  O17_MASS_EV_FLUKA        15830508491.1682     1000080170   //
#define  O18_MASS_EV_FLUKA        16762029651.1675     1000080180   //
#define  O19_MASS_EV_FLUKA        17697638381.0609     1000080190   //
#define  Ne20_MASS_EV_FLUKA       18617738027.1804     1000100200   //
#define  Ne21_MASS_EV_FLUKA       19550542547.1464     1000100210   //
#define  Ne22_MASS_EV_FLUKA       20479744147.2058     1000100220   //
#define  Na23_MASS_EV_FLUKA       21409223268.0681     1000110230   //
#define  Na24_MASS_EV_FLUKA       22341829488.0393     1000110240   //
#define  Mg26_MASS_EV_FLUKA       24196511250.3284     1000120260   //
#define  Si28_MASS_EV_FLUKA       26053202014.9776     1000140280   //
#define  Si29_MASS_EV_FLUKA       26984294034.9880     1000140290   //
#define  P31_MASS_EV_FLUKA        28844227175.2577     1000150310   //
#define  S34_MASS_EV_FLUKA        31632709902.4699     1000160340   //
#define  S35_MASS_EV_FLUKA        32565289722.4418     1000160350   //
#define  Cl36_MASS_EV_FLUKA       33495599129.8856     1000170360   //
#define  Cl37_MASS_EV_FLUKA       34424853859.9436     1000170370   //
#define  Ar36_MASS_EV_FLUKA       33494381341.0459     1000180360   //
#define  Ar40_MASS_EV_FLUKA       37215549251.1705     1000180400   //
#define  K40_MASS_EV_FLUKA        37216545099.2299     1000190400   //
#define  Ca44_MASS_EV_FLUKA       40934079347.6229     1000200440   //
#define  Ti44_MASS_EV_FLUKA       40936982876.2058     1000220440   //
#define  Cr50_MASS_EV_FLUKA       46512225920.8323     1000240500   //
#define  Mn51_MASS_EV_FLUKA       47445229663.1891     1000250510   //
#define  Mn54_MASS_EV_FLUKA       50232398113.3786     1000250540   //
#define  Fe56_MASS_EV_FLUKA       52089829152.9667     1000260560   //
#define  Co57_MASS_EV_FLUKA       53022077022.0146     1000270570   //
#define  Cu65_MASS_EV_FLUKA       60465097019.2569     1000290650   //
#define  Ga65_MASS_EV_FLUKA       60468689064.1419     1000310650   //
#define  Ga68_MASS_EV_FLUKA       63258742124.2567     1000310680   //
#define  Ge72_MASS_EV_FLUKA       66978710348.2507     1000320720   //
#define  As67_MASS_EV_FLUKA       62336673285.6217     1000330670   //
#define  Y87_MASS_EV_FLUKA        80937153217.2207     1000390870   //
#define  Te118_MASS_EV_FLUKA      109802221414.320     1000521180   //
#define  Cs123_MASS_EV_FLUKA     114464850781.1856     1000551230   //
#define  Hf164_MASS_EV_FLUKA     152676900144.6193     1000721640   //
#define  Tl193_MASS_EV_FLUKA     179710102492.4562     1000811930   //
#define  Pb208_MASS_EV_FLUKA     193687690162.6481     1000822080   // best known value: 193.7290249(13)e9   (207.9766521(13) u)

#endif /* XCOLL_FLUKA_MASSES_H */
'''

fluka_masses = {
    int(line.split()[3]): [line.split()[1].split('_')[0], float(line.split()[2])]
    for line in source.split('\n')
    if len(line.split()) > 1 and len(line.split()[1]) > 6 and line.split()[1][-6:] == '_FLUKA'
}
