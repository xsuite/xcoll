# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

import xtrack.particles.pdg as pdg


source = r'''
#ifndef XCOLL_FLUKA_MASSES_H
#define XCOLL_FLUKA_MASSES_H

// 1u = 931494102.42(28)    2018 CODATA

//       name                            mass [eV]                      PDG ID
#define  ELECTRON_MASS_EV_FLUKA        510999.0600000001      //         11   // best known value:     510998.95069(16) eV
#define  MUON_MASS_EV_FLUKA         105658389.0000            //         13   // best known value:     105658375.5(23) eV
#define  PION_MASS_EV_FLUKA         139569950.0000            //        211   //
#define  KAON_MASS_EV_FLUKA         493646000.0000            //        321   //
#define  PROTON_MASS_EV_FLUKA       938272310.0000            //       2212   // best known value:     938272088.15(29) eV
#define  DEUTERON_MASS_EV_FLUKA    1875613390.0000            // 1000010020   //
#define  TRITON_MASS_EV_FLUKA      2808921774.1526            // 1000010030   //

#define  He3_MASS_EV_FLUKA         2808392236.6048            // 1000020030   //
#define  He4_MASS_EV_FLUKA         3727380256.9289            // 1000020040   //
#define  He5_MASS_EV_FLUKA         4667835926.6967            // 1000020050   //
#define  He6_MASS_EV_FLUKA         5605538106.5358            // 1000020060   //
#define  He7_MASS_EV_FLUKA         6545548646.3151            // 1000020070   //
#define  He8_MASS_EV_FLUKA         7482530706.1729            // 1000020080   //

#define  Li6_MASS_EV_FLUKA         5601519470.4097            // 1000030060   //
#define  Li7_MASS_EV_FLUKA         6533835190.3885            // 1000030070   //
#define  Li8_MASS_EV_FLUKA         7471367020.2320            // 1000030080   //
#define  Li9_MASS_EV_FLUKA         8406870100.1281            // 1000030090   //
#define  Li10_MASS_EV_FLUKA        9346855049.9081            // 1000030100   //
#define  Li11_MASS_EV_FLUKA       10285692859.7178            // 1000030110   //

#define  Be7_MASS_EV_FLUKA         6534186187.0551            // 1000040070   //
#define  Be8_MASS_EV_FLUKA         7454852657.3357            // 1000040080   //
#define  Be9_MASS_EV_FLUKA         8392753017.1697            // 1000040090   //
#define  Be10_MASS_EV_FLUKA        9325506337.1371            // 1000040100   //
#define  Be11_MASS_EV_FLUKA       10264568076.9410            // 1000040110   //
#define  Be12_MASS_EV_FLUKA       11200964756.8139            // 1000040120   //
#define  Be13_MASS_EV_FLUKA       12142541006.5527            // 1000040130   //
#define  Be14_MASS_EV_FLUKA       13078759286.4303            // 1000040140   //

#define  B8_MASS_EV_FLUKA          7472321245.2018            // 1000050080   //
#define  B9_MASS_EV_FLUKA          8393310365.4740            // 1000050090   //
#define  B10_MASS_EV_FLUKA         9324439705.4834            // 1000050100   //
#define  B11_MASS_EV_FLUKA        10252551235.5711            // 1000050110   //
#define  B12_MASS_EV_FLUKA        11188746425.4493            // 1000050120   //
#define  B13_MASS_EV_FLUKA        12123434145.3666            // 1000050130   //
#define  B14_MASS_EV_FLUKA        13062029875.1825            // 1000050140   //
#define  B15_MASS_EV_FLUKA        13998827395.0451            // 1000050150   //
#define  B16_MASS_EV_FLUKA        14938494194.8333            // 1000050160   //
#define  B17_MASS_EV_FLUKA        15876565424.6629            // 1000050170   //
#define  B18_MASS_EV_FLUKA        16816666234.4399            // 1000050180   //
#define  B19_MASS_EV_FLUKA        17753646308.4133            // 1000050190   //
#define  B20_MASS_EV_FLUKA        18695610342.2641            // 1000050200   //
#define  B21_MASS_EV_FLUKA        19627796224.8218            // 1000050210   //
#define  B22_MASS_EV_FLUKA        20561795084.8209            // 1000050220   //
#define  B23_MASS_EV_FLUKA        21491344416.3474            // 1000050230   //
#define  B24_MASS_EV_FLUKA        22425823568.7792            // 1000050240   //
#define  B25_MASS_EV_FLUKA        23361517466.3034            // 1000050250   //

#define  C9_MASS_EV_FLUKA          8409297942.4853            // 1000060090   //
#define  C10_MASS_EV_FLUKA         9327576872.8277            // 1000060100   //
#define  C11_MASS_EV_FLUKA        10254022852.9586            // 1000060110   //
#define  C12_MASS_EV_FLUKA        11174866883.2345            // 1000060120   //
#define  C13_MASS_EV_FLUKA        12109486203.1536            // 1000060130   //
#define  C14_MASS_EV_FLUKA        13040875423.1563            // 1000060140   //
#define  C15_MASS_EV_FLUKA        13979223022.9787            // 1000060150   //
#define  C16_MASS_EV_FLUKA        14914538302.8797            // 1000060160   //
#define  C17_MASS_EV_FLUKA        15853375102.6894            // 1000060170   //
#define  C18_MASS_EV_FLUKA        16788757102.5887            // 1000060180   //
#define  C30_MASS_EV_FLUKA        28026329089.7061            // 1000060300   //

#define  N12_MASS_EV_FLUKA        11191694455.9209            // 1000070120   //
#define  N13_MASS_EV_FLUKA        12111196196.2317            // 1000070130   //
#define  N14_MASS_EV_FLUKA        13040208386.2960            // 1000070140   //
#define  N15_MASS_EV_FLUKA        13968940806.3676            // 1000070150   //
#define  N16_MASS_EV_FLUKA        14906015646.2229            // 1000070160   //
#define  N17_MASS_EV_FLUKA        15839698756.1662            // 1000070170   //
#define  N18_MASS_EV_FLUKA        16776439386.0303            // 1000070180   //
#define  N19_MASS_EV_FLUKA        17710677005.9592            // 1000070190   //
#define  N20_MASS_EV_FLUKA        18648077455.8061            // 1000070200   //
#define  N21_MASS_EV_FLUKA        19583037215.7163            // 1000070210   //

#define  O13_MASS_EV_FLUKA        12128451020.5484            // 1000080130   //
#define  O14_MASS_EV_FLUKA        13044841910.9398            // 1000080140   //
#define  O15_MASS_EV_FLUKA        13971184421.0733            // 1000080150   //
#define  O16_MASS_EV_FLUKA        14895086191.2700            // 1000080160   //
#define  O17_MASS_EV_FLUKA        15830508491.1682            // 1000080170   //
#define  O18_MASS_EV_FLUKA        16762029651.1675            // 1000080180   //
#define  O19_MASS_EV_FLUKA        17697638381.0609            // 1000080190   //
#define  O20_MASS_EV_FLUKA        18629597381.0489            // 1000080200   //
#define  O21_MASS_EV_FLUKA        19565356500.9384            // 1000080210   //
#define  O22_MASS_EV_FLUKA        20498073470.9067            // 1000080220   //

#define  F17_MASS_EV_FLUKA        15832758902.9589            // 1000090170   //
#define  F18_MASS_EV_FLUKA        16763174962.9869            // 1000090180   //
#define  F19_MASS_EV_FLUKA        17692308413.0480            // 1000090190   //
#define  F20_MASS_EV_FLUKA        18625272733.0099            // 1000090200   //
#define  F21_MASS_EV_FLUKA        19556736843.0107            // 1000090210   //
#define  F22_MASS_EV_FLUKA        20491072592.9371            // 1000090220   //
#define  F23_MASS_EV_FLUKA        21423102682.9232            // 1000090230   //
#define  F24_MASS_EV_FLUKA        22358811992.8140            // 1000090240   //
#define  F25_MASS_EV_FLUKA        23294028222.7175            // 1000090250   //

#define  Ne18_MASS_EV_FLUKA       16767110436.8600            // 1000100180   //
#define  Ne19_MASS_EV_FLUKA       17695036666.9525            // 1000100190   //
#define  Ne20_MASS_EV_FLUKA       18617738027.1804            // 1000100200   //
#define  Ne21_MASS_EV_FLUKA       19550542547.1464            // 1000100210   //
#define  Ne22_MASS_EV_FLUKA       20479744147.2058            // 1000100220   //
#define  Ne23_MASS_EV_FLUKA       21414109157.1314            // 1000100230   //
#define  Ne24_MASS_EV_FLUKA       22344809697.1520            // 1000100240   //
#define  Ne25_MASS_EV_FLUKA       23280192787.0512            // 1000100250   //
#define  Ne26_MASS_EV_FLUKA       24214175696.9867            // 1000100260   //
#define  Ne27_MASS_EV_FLUKA       25152333756.8140            // 1000100270   //

#define  Na19_MASS_EV_FLUKA       17705704117.4861            // 1000110190   //
#define  Na20_MASS_EV_FLUKA       18631114797.6438            // 1000110200   //
#define  Na21_MASS_EV_FLUKA       19553579777.8778            // 1000110210   //
#define  Na22_MASS_EV_FLUKA       20482076227.9555            // 1000110220   //
#define  Na23_MASS_EV_FLUKA       21409223268.0681            // 1000110230   //
#define  Na24_MASS_EV_FLUKA       22341829488.0393            // 1000110240   //
#define  Na25_MASS_EV_FLUKA       23272383918.0637            // 1000110250   //
#define  Na26_MASS_EV_FLUKA       24206333238.0001            // 1000110260   //
#define  Na27_MASS_EV_FLUKA       25139149127.9658            // 1000110270   //
#define  Na28_MASS_EV_FLUKA       26075190707.8480            // 1000110280   //
#define  Na29_MASS_EV_FLUKA       27010337597.7533            // 1000110290   //
#define  Na30_MASS_EV_FLUKA       27947807517.5985            // 1000110300   //
#define  Na31_MASS_EV_FLUKA       28883371167.4930            // 1000110310   //

#define  Mg20_MASS_EV_FLUKA       18641330379.4529            // 1000120200   //
#define  Mg21_MASS_EV_FLUKA       19566165889.6254            // 1000120210   //
#define  Mg22_MASS_EV_FLUKA       20486351679.9185            // 1000120220   //
#define  Mg23_MASS_EV_FLUKA       21412770100.0500            // 1000120230   //
#define  Mg24_MASS_EV_FLUKA       22335803710.2693            // 1000120240   //
#define  Mg25_MASS_EV_FLUKA       23268038730.2501            // 1000120250   //
#define  Mg26_MASS_EV_FLUKA       24196511250.3284            // 1000120260   //
#define  Mg27_MASS_EV_FLUKA       25129633570.2862            // 1000120270   //
#define  Mg28_MASS_EV_FLUKA       26060695590.2974            // 1000120280   //
#define  Mg29_MASS_EV_FLUKA       26996547500.1845            // 1000120290   //
#define  Mg30_MASS_EV_FLUKA       27929820780.1384            // 1000120300   //
#define  Mg31_MASS_EV_FLUKA       28866982189.9915            // 1000120310   //
#define  Mg32_MASS_EV_FLUKA       29800895989.9288            // 1000120320   //

#define  Al23_MASS_EV_FLUKA       21424500177.2724            // 1000130230   //
#define  Al24_MASS_EV_FLUKA       22349172137.4492            // 1000130240   //
#define  Al25_MASS_EV_FLUKA       23271805727.6788            // 1000130250   //
#define  Al26_MASS_EV_FLUKA       24200005587.7642            // 1000130260   //
#define  Al27_MASS_EV_FLUKA       25126513317.8934            // 1000130270   //
#define  Al28_MASS_EV_FLUKA       26058353937.8845            // 1000130280   //
#define  Al29_MASS_EV_FLUKA       26988483367.9198            // 1000130290   //
#define  Al30_MASS_EV_FLUKA       27922320787.8591            // 1000130300   //
#define  Al31_MASS_EV_FLUKA       28854733297.8353            // 1000130310   //
#define  Al32_MASS_EV_FLUKA       29790119647.7345            // 1000130320   //

#define  Si25_MASS_EV_FLUKA       23284037184.3215            // 1000140250   //
#define  Si26_MASS_EV_FLUKA       24204561464.6058            // 1000140260   //
#define  Si27_MASS_EV_FLUKA       25130815494.7416            // 1000140270   //
#define  Si28_MASS_EV_FLUKA       26053202014.9776            // 1000140280   //
#define  Si29_MASS_EV_FLUKA       26984294034.9880            // 1000140290   //
#define  Si30_MASS_EV_FLUKA       27913250555.0538            // 1000140300   //
#define  Si31_MASS_EV_FLUKA       28846228775.0153            // 1000140310   //
#define  Si32_MASS_EV_FLUKA       29776591195.0447            // 1000140320   //
#define  Si33_MASS_EV_FLUKA       30711674014.9517            // 1000140330   //
#define  Si34_MASS_EV_FLUKA       31643704144.9378            // 1000140340   //
#define  Si35_MASS_EV_FLUKA       32580795234.7928            // 1000140350   //

#define  P27_MASS_EV_FLUKA        25141937994.6439            // 1000150270   //
#define  P28_MASS_EV_FLUKA        26067024104.8099            // 1000150280   //
#define  P29_MASS_EV_FLUKA        26988727585.0637            // 1000150290   //
#define  P30_MASS_EV_FLUKA        27916973245.1478            // 1000150300   //
#define  P31_MASS_EV_FLUKA        28844227175.2577            // 1000150310   //
#define  P32_MASS_EV_FLUKA        29775857195.2542            // 1000150320   //
#define  P33_MASS_EV_FLUKA        30705319135.3069            // 1000150330   //
#define  P34_MASS_EV_FLUKA        31638593515.2608            // 1000150340   //
#define  P35_MASS_EV_FLUKA        32569787865.2685            // 1000150350   //
#define  P36_MASS_EV_FLUKA        33505888995.1491            // 1000150360   //
#define  P37_MASS_EV_FLUKA        34438639385.1166            // 1000150370   //
#define  P38_MASS_EV_FLUKA        35374662334.9992            // 1000150380   //
#define  P39_MASS_EV_FLUKA        36307973054.9522            // 1000150390   //

#define  S29_MASS_EV_FLUKA        27002011261.7761            // 1000160290   //
#define  S30_MASS_EV_FLUKA        27922601602.0587            // 1000160300   //
#define  S31_MASS_EV_FLUKA        28849113842.1878            // 1000160310   //
#define  S32_MASS_EV_FLUKA        29773637162.3685            // 1000160320   //
#define  S33_MASS_EV_FLUKA        30704561182.3832            // 1000160330   //
#define  S34_MASS_EV_FLUKA        31632709902.4699            // 1000160340   //
#define  S35_MASS_EV_FLUKA        32565289722.4418            // 1000160350   //
#define  S36_MASS_EV_FLUKA        33494966332.4889            // 1000160360   //
#define  S37_MASS_EV_FLUKA        34430228442.3913            // 1000160370   //
#define  S38_MASS_EV_FLUKA        35361757822.3904            // 1000160380   //
#define  S39_MASS_EV_FLUKA        36296951962.2945            // 1000160390   //
#define  S40_MASS_EV_FLUKA        37228758112.2864            // 1000160400   //
#define  S41_MASS_EV_FLUKA        38164500032.1763            // 1000160410   //

#define  Cl32_MASS_EV_FLUKA       29785812989.4660            // 1000170320   //
#define  Cl33_MASS_EV_FLUKA       30709634529.6648            // 1000170330   //
#define  Cl34_MASS_EV_FLUKA       31637692779.7539            // 1000170340   //
#define  Cl35_MASS_EV_FLUKA       32564613219.8724            // 1000170350   //
#define  Cl36_MASS_EV_FLUKA       33495599129.8856            // 1000170360   //
#define  Cl37_MASS_EV_FLUKA       34424853859.9436            // 1000170370   //
#define  Cl38_MASS_EV_FLUKA       35358311669.8927            // 1000170380   //
#define  Cl39_MASS_EV_FLUKA       36289804179.8928            // 1000170390   //
#define  Cl40_MASS_EV_FLUKA       37223540589.8347            // 1000170400   //
#define  Cl41_MASS_EV_FLUKA       38155253479.8290            // 1000170410   //
#define  Cl42_MASS_EV_FLUKA       39089099649.7681            // 1000170420   //
#define  Cl43_MASS_EV_FLUKA       40021551859.7432            // 1000170430   //

#define  Ar34_MASS_EV_FLUKA       31643244900.7387            // 1000180340   //
#define  Ar35_MASS_EV_FLUKA       32570069270.8597            // 1000180350   //
#define  Ar36_MASS_EV_FLUKA       33494381341.0459            // 1000180360   //
#define  Ar37_MASS_EV_FLUKA       34425158161.0645            // 1000180370   //
#define  Ar38_MASS_EV_FLUKA       35352885661.1621            // 1000180380   //
#define  Ar39_MASS_EV_FLUKA       36285852981.1239            // 1000180390   //
#define  Ar40_MASS_EV_FLUKA       37215549251.1705            // 1000180400   //
#define  Ar41_MASS_EV_FLUKA       38149016101.1194            // 1000180410   //
#define  Ar42_MASS_EV_FLUKA       39079155651.1545            // 1000180420   //
#define  Ar43_MASS_EV_FLUKA       40013094441.0911            // 1000180430   //
#define  Ar44_MASS_EV_FLUKA       40944304311.0985            // 1000180440   //
#define  Ar45_MASS_EV_FLUKA       41878341291.0326            // 1000180450   //
#define  Ar46_MASS_EV_FLUKA       42809834231.0327            // 1000180460   //

#define  K37_MASS_EV_FLUKA        34430797839.0035            // 1000190370   //
#define  K38_MASS_EV_FLUKA        35358289719.1072            // 1000190380   //
#define  K39_MASS_EV_FLUKA        36284778979.2369            // 1000190390   //
#define  K40_MASS_EV_FLUKA        37216545099.2299            // 1000190400   //
#define  K41_MASS_EV_FLUKA        38146015519.2823            // 1000190410   //
#define  K42_MASS_EV_FLUKA        39078047439.2684            // 1000190420   //
#define  K43_MASS_EV_FLUKA        40007969189.3091            // 1000190430   //
#define  K44_MASS_EV_FLUKA        40940247179.2888            // 1000190440   //
#define  K45_MASS_EV_FLUKA        41870943619.3095            // 1000190450   //
#define  K46_MASS_EV_FLUKA        42803627039.2787            // 1000190460   //
#define  K47_MASS_EV_FLUKA        43734842569.2859            // 1000190470   //

#define  Ca39_MASS_EV_FLUKA       36290800617.2551            // 1000200390   //
#define  Ca40_MASS_EV_FLUKA       37214725127.4512            // 1000200400   //
#define  Ca41_MASS_EV_FLUKA       38145928037.4588            // 1000200410   //
#define  Ca42_MASS_EV_FLUKA       39074013057.5471            // 1000200420   //
#define  Ca43_MASS_EV_FLUKA       40005645767.5435            // 1000200430   //
#define  Ca44_MASS_EV_FLUKA       40934079347.6229            // 1000200440   //
#define  Ca45_MASS_EV_FLUKA       41866230267.6058            // 1000200450   //
#define  Ca46_MASS_EV_FLUKA       42795402147.6660            // 1000200460   //
#define  Ca47_MASS_EV_FLUKA       43727691777.6454            // 1000200470   //
#define  Ca48_MASS_EV_FLUKA       44657310977.6940            // 1000200480   //

#define  Sc41_MASS_EV_FLUKA       38151914628.6004            // 1000210410   //
#define  Sc42_MASS_EV_FLUKA       39079930238.6905            // 1000210420   //
#define  Sc43_MASS_EV_FLUKA       40007357878.7959            // 1000210430   //
#define  Sc44_MASS_EV_FLUKA       40937223998.8381            // 1000210440   //
#define  Sc45_MASS_EV_FLUKA       41865464728.9224            // 1000210450   //
#define  Sc46_MASS_EV_FLUKA       42796269748.9403            // 1000210460   //
#define  Sc47_MASS_EV_FLUKA       43725191059.0069            // 1000210470   //
#define  Sc48_MASS_EV_FLUKA       44656524149.0111            // 1000210480   //
#define  Sc49_MASS_EV_FLUKA       45585959079.0645            // 1000210490   //
#define  Sc50_MASS_EV_FLUKA       46519468119.0123            // 1000210500   //
#define  Sc51_MASS_EV_FLUKA       47452281138.9781            // 1000210510   //

#define  Ti42_MASS_EV_FLUKA       39086421665.8838            // 1000220420   //
#define  Ti43_MASS_EV_FLUKA       40013716565.9926            // 1000220430   //
#define  Ti44_MASS_EV_FLUKA       40936982876.2058            // 1000220440   //
#define  Ti45_MASS_EV_FLUKA       41867018666.2436            // 1000220450   //
#define  Ti46_MASS_EV_FLUKA       42793394486.3763            // 1000220460   //
#define  Ti47_MASS_EV_FLUKA       43724082406.3972            // 1000220470   //
#define  Ti48_MASS_EV_FLUKA       44652021436.4893            // 1000220480   //
#define  Ti49_MASS_EV_FLUKA       45583444746.4911            // 1000220490   //
#define  Ti50_MASS_EV_FLUKA       46512071266.5655            // 1000220500   //
#define  Ti51_MASS_EV_FLUKA       47445264586.5214            // 1000220510   //
#define  Ti52_MASS_EV_FLUKA       48377021646.5146            // 1000220520   //
#define  Ti53_MASS_EV_FLUKA       49311155226.4462            // 1000220530   //

#define  V46_MASS_EV_FLUKA        42799937472.4909            // 1000230460   //
#define  V47_MASS_EV_FLUKA        43726501802.6187            // 1000230470   //
#define  V48_MASS_EV_FLUKA        44655525402.6827            // 1000230480   //
#define  V49_MASS_EV_FLUKA        45583538242.7729            // 1000230490   //
#define  V50_MASS_EV_FLUKA        46513771062.8056            // 1000230500   //
#define  V51_MASS_EV_FLUKA        47442285482.8828            // 1000230510   //
#define  V52_MASS_EV_FLUKA        48374539902.8631            // 1000230520   //
#define  V53_MASS_EV_FLUKA        49305626802.8737            // 1000230530   //
#define  V54_MASS_EV_FLUKA        50239078932.8230            // 1000230540   //
#define  V55_MASS_EV_FLUKA        51171312762.8038            // 1000230550   //
#define  V56_MASS_EV_FLUKA        52105797572.7263            // 1000230560   //

#define  Cr48_MASS_EV_FLUKA       44656676410.6395            // 1000240480   //
#define  Cr49_MASS_EV_FLUKA       45585660580.7045            // 1000240490   //
#define  Cr50_MASS_EV_FLUKA       46512225920.8323            // 1000240500   //
#define  Cr51_MASS_EV_FLUKA       47442529940.8631            // 1000240510   //
#define  Cr52_MASS_EV_FLUKA       48370056060.9659            // 1000240520   //
#define  Cr53_MASS_EV_FLUKA       49301682480.9625            // 1000240530   //
#define  Cr54_MASS_EV_FLUKA       50231529101.0052            // 1000240540   //
#define  Cr55_MASS_EV_FLUKA       51164848520.9579            // 1000240550   //
#define  Cr56_MASS_EV_FLUKA       52096157450.9627            // 1000240560   //
#define  Cr60_MASS_EV_FLUKA       55830597500.7434            // 1000240600   //

#define  Mn50_MASS_EV_FLUKA       46519350843.0435            // 1000250500   //
#define  Mn51_MASS_EV_FLUKA       47445229663.1891            // 1000250510   //
#define  Mn52_MASS_EV_FLUKA       48374259773.2529            // 1000250520   //
#define  Mn53_MASS_EV_FLUKA       49301771403.3561            // 1000250530   //
#define  Mn54_MASS_EV_FLUKA       50232398113.3786            // 1000250540   //
#define  Mn55_MASS_EV_FLUKA       51161737343.4344            // 1000250550   //
#define  Mn56_MASS_EV_FLUKA       52094032463.4137            // 1000250560   //
#define  Mn57_MASS_EV_FLUKA       53024947463.4287            // 1000250570   //
#define  Mn58_MASS_EV_FLUKA       53958024583.3877            // 1000250580   //
#define  Mn59_MASS_EV_FLUKA       54889947903.3766            // 1000250590   //
#define  Mn60_MASS_EV_FLUKA       55824140833.3066            // 1000250600   //

#define  Fe50_MASS_EV_FLUKA       46526992882.2895            // 1000260500   //
#define  Fe52_MASS_EV_FLUKA       48376123872.6487            // 1000260520   //
#define  Fe53_MASS_EV_FLUKA       49305006082.7163            // 1000260530   //
#define  Fe54_MASS_EV_FLUKA       50231193312.8539            // 1000260540   //
#define  Fe55_MASS_EV_FLUKA       51161461032.8857            // 1000260550   //
#define  Fe56_MASS_EV_FLUKA       52089829152.9667            // 1000260560   //
#define  Fe57_MASS_EV_FLUKA       53021748772.9557            // 1000260570   //
#define  Fe58_MASS_EV_FLUKA       53951269883.0068            // 1000260580   //
#define  Fe59_MASS_EV_FLUKA       54884254612.9681            // 1000260590   //
#define  Fe60_MASS_EV_FLUKA       55815000402.9875            // 1000260600   //
#define  Fe61_MASS_EV_FLUKA       56748984242.9230            // 1000260610   //
#define  Fe62_MASS_EV_FLUKA       57680498162.9225            // 1000260620   //

#define  Co53_MASS_EV_FLUKA       49312800341.5818            // 1000270530   //
#define  Co54_MASS_EV_FLUKA       50238928571.7209            // 1000270540   //
#define  Co55_MASS_EV_FLUKA       51164404491.8768            // 1000270550   //
#define  Co56_MASS_EV_FLUKA       52093887401.9290            // 1000270560   //
#define  Co57_MASS_EV_FLUKA       53022077022.0146            // 1000270570   //
#define  Co58_MASS_EV_FLUKA       53953069642.0276            // 1000270580   //
#define  Co59_MASS_EV_FLUKA       54882181762.0893            // 1000270590   //
#define  Co60_MASS_EV_FLUKA       55814255492.0743            // 1000270600   //
#define  Co61_MASS_EV_FLUKA       56744498902.1068            // 1000270610   //
#define  Co62_MASS_EV_FLUKA       57677460222.0687            // 1000270620   //
#define  Co63_MASS_EV_FLUKA       58608545642.0793            // 1000270630   //
#define  Co64_MASS_EV_FLUKA       59542087662.0263            // 1000270640   //
#define  Co65_MASS_EV_FLUKA       60474207092.0101            // 1000270650   //

#define  Ni56_MASS_EV_FLUKA       52095515343.0931            // 1000280560   //
#define  Ni57_MASS_EV_FLUKA       53024833743.1495            // 1000280570   //
#define  Ni58_MASS_EV_FLUKA       53952180483.2570            // 1000280580   //
#define  Ni59_MASS_EV_FLUKA       54882746703.2811            // 1000280590   //
#define  Ni60_MASS_EV_FLUKA       55810924023.3670            // 1000280600   //
#define  Ni61_MASS_EV_FLUKA       56742669633.3605            // 1000280610   //
#define  Ni62_MASS_EV_FLUKA       57671638063.4260            // 1000280620   //
#define  Ni63_MASS_EV_FLUKA       58604365883.3940            // 1000280630   //
#define  Ni64_MASS_EV_FLUKA       59534273503.4351            // 1000280640   //
#define  Ni65_MASS_EV_FLUKA       60467741113.3840            // 1000280650   //
#define  Ni66_MASS_EV_FLUKA       61398329433.4074            // 1000280660   //

#define  Cu58_MASS_EV_FLUKA       53960236178.8527            // 1000290580   //
#define  Cu59_MASS_EV_FLUKA       54887038898.9743            // 1000290590   //
#define  Cu60_MASS_EV_FLUKA       55816543509.0258            // 1000290600   //
#define  Cu61_MASS_EV_FLUKA       56744399339.1201            // 1000290610   //
#define  Cu62_MASS_EV_FLUKA       57675078839.1412            // 1000290620   //
#define  Cu63_MASS_EV_FLUKA       58603791479.2133            // 1000290630   //
#define  Cu64_MASS_EV_FLUKA       59535441209.2093            // 1000290640   //
#define  Cu65_MASS_EV_FLUKA       60465097019.2569            // 1000290650   //
#define  Cu66_MASS_EV_FLUKA       61397596739.2309            // 1000290660   //
#define  Cu67_MASS_EV_FLUKA       62328044799.2580            // 1000290670   //
#define  Cu68_MASS_EV_FLUKA       63261297359.2124            // 1000290680   //
#define  Cu70_MASS_EV_FLUKA       65126867629.1455            // 1000290700   //
#define  Cu71_MASS_EV_FLUKA       66058558029.1404            // 1000290710   //

#define  Zn61_MASS_EV_FLUKA       56749529231.7960            // 1000300610   //
#define  Zn62_MASS_EV_FLUKA       57676198651.9210            // 1000300620   //
#define  Zn63_MASS_EV_FLUKA       58606651061.9480            // 1000300630   //
#define  Zn64_MASS_EV_FLUKA       59534355082.0462            // 1000300640   //
#define  Zn65_MASS_EV_FLUKA       60465941202.0438            // 1000300650   //
#define  Zn66_MASS_EV_FLUKA       61394447132.1213            // 1000300660   //
#define  Zn67_MASS_EV_FLUKA       62326960642.0949            // 1000300670   //
#define  Zn68_MASS_EV_FLUKA       63256328062.1500            // 1000300680   //
#define  Zn69_MASS_EV_FLUKA       64189411492.1088            // 1000300690   //
#define  Zn70_MASS_EV_FLUKA       65119761292.1385            // 1000300700   //
#define  Zn71_MASS_EV_FLUKA       66053493432.0805            // 1000300710   //
#define  Zn72_MASS_EV_FLUKA       66984183402.1013            // 1000300720   //
#define  Zn73_MASS_EV_FLUKA       67918394042.0309            // 1000300730   //

#define  Ga62_MASS_EV_FLUKA       57684862593.8657            // 1000310620   //
#define  Ga63_MASS_EV_FLUKA       58611664033.9873            // 1000310630   //
#define  Ga64_MASS_EV_FLUKA       59541012824.0429            // 1000310640   //
#define  Ga65_MASS_EV_FLUKA       60468689064.1419            // 1000310650   //
#define  Ga66_MASS_EV_FLUKA       61399115074.1695            // 1000310660   //
#define  Ga67_MASS_EV_FLUKA       62327454004.2513            // 1000310670   //
#define  Ga68_MASS_EV_FLUKA       63258742124.2567            // 1000310680   //
#define  Ga69_MASS_EV_FLUKA       64187998434.3146            // 1000310690   //
#define  Ga70_MASS_EV_FLUKA       65119909054.3038            // 1000310700   //
#define  Ga71_MASS_EV_FLUKA       66050173684.3357            // 1000310710   //
#define  Ga72_MASS_EV_FLUKA       66983218304.2955            // 1000310720   //
#define  Ga73_MASS_EV_FLUKA       67913593174.3246            // 1000310730   //
#define  Ga74_MASS_EV_FLUKA       68846737294.2818            // 1000310740   //
#define  Ga75_MASS_EV_FLUKA       69777821424.2924            // 1000310750   //

#define  Ge65_MASS_EV_FLUKA       60474424427.8315            // 1000320650   //
#define  Ge66_MASS_EV_FLUKA       61400708217.9665            // 1000320660   //
#define  Ge67_MASS_EV_FLUKA       62331169917.9933            // 1000320670   //
#define  Ge68_MASS_EV_FLUKA       63258341128.1053            // 1000320680   //
#define  Ge69_MASS_EV_FLUKA       64189718878.1083            // 1000320690   //
#define  Ge70_MASS_EV_FLUKA       65117746418.1982            // 1000320700   //
#define  Ge71_MASS_EV_FLUKA       66049896228.1812            // 1000320710   //
#define  Ge72_MASS_EV_FLUKA       66978710348.2507            // 1000320720   //
#define  Ge73_MASS_EV_FLUKA       67911493078.2173            // 1000320730   //
#define  Ge74_MASS_EV_FLUKA       68840862498.2723            // 1000320740   //
#define  Ge75_MASS_EV_FLUKA       69773922918.2317            // 1000320750   //
#define  Ge76_MASS_EV_FLUKA       70704060328.2669            // 1000320760   //
#define  Ge77_MASS_EV_FLUKA       71637553348.2151            // 1000320770   //

#define  As66_MASS_EV_FLUKA       61410000975.4967            // 1000330660   //
#define  As67_MASS_EV_FLUKA       62336673285.6217            // 1000330670   //
#define  As68_MASS_EV_FLUKA       63265934505.6795            // 1000330680   //
#define  As69_MASS_EV_FLUKA       64193225205.7885            // 1000330690   //
#define  As70_MASS_EV_FLUKA       65123459705.8211            // 1000330700   //
#define  As71_MASS_EV_FLUKA       66051402225.9132            // 1000330710   //
#define  As72_MASS_EV_FLUKA       66982559755.9219            // 1000330720   //
#define  As73_MASS_EV_FLUKA       67911327275.9925            // 1000330730   //
#define  As74_MASS_EV_FLUKA       68842918215.9900            // 1000330740   //
#define  As75_MASS_EV_FLUKA       69772239736.0464            // 1000330750   //
#define  As76_MASS_EV_FLUKA       70704476966.0271            // 1000330760   //
#define  As77_MASS_EV_FLUKA       71634344676.0693            // 1000330770   //
#define  As78_MASS_EV_FLUKA       72566938916.0408            // 1000330780   //
#define  As78_MASS_EV_FLUKA       72566938916.0408            // 1000330780   //
#define  As81_MASS_EV_FLUKA       75361705426.0334            // 1000330810   //
#define  As82_MASS_EV_FLUKA       76295493635.9740            // 1000330820   //

#define  Se69_MASS_EV_FLUKA       64199501809.5504            // 1000340690   //
#define  Se70_MASS_EV_FLUKA       65125352759.6967            // 1000340700   //
#define  Se71_MASS_EV_FLUKA       66055695089.7265            // 1000340710   //
#define  Se72_MASS_EV_FLUKA       66982388299.8509            // 1000340720   //
#define  Se73_MASS_EV_FLUKA       67913560819.8593            // 1000340730   //
#define  Se74_MASS_EV_FLUKA       68841058729.9628            // 1000340740   //
#define  Se75_MASS_EV_FLUKA       69772596849.9617            // 1000340750   //
#define  Se76_MASS_EV_FLUKA       70701008480.0416            // 1000340760   //
#define  Se77_MASS_EV_FLUKA       71633155290.0247            // 1000340770   //
#define  Se78_MASS_EV_FLUKA       72562223020.0876            // 1000340780   //
#define  Se79_MASS_EV_FLUKA       73494825740.0588            // 1000340790   //
#define  Se80_MASS_EV_FLUKA       74424477950.1066            // 1000340800   //
#define  Se81_MASS_EV_FLUKA       75357342570.0710            // 1000340810   //

#define  Br71_MASS_EV_FLUKA       66061688691.8296            // 1000350710   //
#define  Br72_MASS_EV_FLUKA       66990596491.8966            // 1000350720   //
#define  Br73_MASS_EV_FLUKA       67917710922.0101            // 1000350730   //
#define  Br74_MASS_EV_FLUKA       68847459132.0553            // 1000350740   //
#define  Br75_MASS_EV_FLUKA       69775120562.1546            // 1000350750   //
#define  Br76_MASS_EV_FLUKA       70705464992.1845            // 1000350760   //
#define  Br77_MASS_EV_FLUKA       71634014082.2608            // 1000350770   //
#define  Br78_MASS_EV_FLUKA       72565290392.2664            // 1000350780   //
#define  Br79_MASS_EV_FLUKA       73494168732.3342            // 1000350790   //
#define  Br80_MASS_EV_FLUKA       74425842252.3296            // 1000350800   //
#define  Br81_MASS_EV_FLUKA       75355251162.3836            // 1000350810   //
#define  Br82_MASS_EV_FLUKA       76287223882.3712            // 1000350820   //
#define  Br83_MASS_EV_FLUKA       77217205482.4104            // 1000350830   //
#define  Br88_MASS_EV_FLUKA       81882957692.1959            // 1000350880   //

#define  Kr73_MASS_EV_FLUKA       67923879414.5838            // 1000360730   //
#define  Kr74_MASS_EV_FLUKA       68850090684.7207            // 1000360740   //
#define  Kr75_MASS_EV_FLUKA       69779513054.7744            // 1000360750   //
#define  Kr76_MASS_EV_FLUKA       70706270274.8972            // 1000360760   //
#define  Kr77_MASS_EV_FLUKA       71636571614.9281            // 1000360770   //
#define  Kr78_MASS_EV_FLUKA       72564077865.0314            // 1000360780   //
#define  Kr79_MASS_EV_FLUKA       73495288405.0388            // 1000360790   //
#define  Kr80_MASS_EV_FLUKA       74423331825.1282            // 1000360800   //
#define  Kr81_MASS_EV_FLUKA       75355025755.1230            // 1000360810   //
#define  Kr82_MASS_EV_FLUKA       76283625175.1981            // 1000360820   //
#define  Kr83_MASS_EV_FLUKA       77215727085.1823            // 1000360830   //
#define  Kr84_MASS_EV_FLUKA       78144772315.2458            // 1000360840   //
#define  Kr85_MASS_EV_FLUKA       79077218535.2211            // 1000360850   //
#define  Kr86_MASS_EV_FLUKA       80006929135.2673            // 1000360860   //

#define  Rb76_MASS_EV_FLUKA       70714264690.0032            // 1000370760   //
#define  Rb77_MASS_EV_FLUKA       71641410210.1158            // 1000370770   //
#define  Rb78_MASS_EV_FLUKA       72570796330.1705            // 1000370780   //
#define  Rb79_MASS_EV_FLUKA       73498431770.2705            // 1000370790   //
#define  Rb80_MASS_EV_FLUKA       74428548780.3062            // 1000370800   //
#define  Rb81_MASS_EV_FLUKA       75356757810.3913            // 1000370810   //
#define  Rb82_MASS_EV_FLUKA       76287520620.4102            // 1000370820   //
#define  Rb83_MASS_EV_FLUKA       77216130210.4850            // 1000370830   //
#define  Rb84_MASS_EV_FLUKA       78146947700.5025            // 1000370840   //
#define  Rb85_MASS_EV_FLUKA       79076025620.5652            // 1000370850   //
#define  Rb86_MASS_EV_FLUKA       80006940040.5802            // 1000370860   //
#define  Rb87_MASS_EV_FLUKA       80936586160.6281            // 1000370870   //
#define  Rb90_MASS_EV_FLUKA       83736310670.4923            // 1000370900   //

#define  Sr79_MASS_EV_FLUKA       73503244220.1075            // 1000380790   //
#define  Sr80_MASS_EV_FLUKA       74429911050.2326            // 1000380800   //
#define  Sr81_MASS_EV_FLUKA       75360183670.2643            // 1000380810   //
#define  Sr82_MASS_EV_FLUKA       76287195210.3804            // 1000380820   //
#define  Sr83_MASS_EV_FLUKA       77217900700.4009            // 1000380830   //
#define  Sr84_MASS_EV_FLUKA       78145547680.5006            // 1000380840   //
#define  Sr85_MASS_EV_FLUKA       79076584690.5124            // 1000380850   //
#define  Sr86_MASS_EV_FLUKA       80004659630.6011            // 1000380860   //
#define  Sr87_MASS_EV_FLUKA       80935797150.6103            // 1000380870   //
#define  Sr88_MASS_EV_FLUKA       81864250160.6891            // 1000380880   //
#define  Sr89_MASS_EV_FLUKA       82797457080.6447            // 1000380890   //
#define  Sr90_MASS_EV_FLUKA       83729215200.6379            // 1000380900   //
#define  Sr91_MASS_EV_FLUKA       84663001380.5785            // 1000380910   //
#define  Sr95_MASS_EV_FLUKA       88397469020.3585            // 1000380950   //

#define  Y79_MASS_EV_FLUKA        73509858736.5816            // 1000390790   //
#define  Y80_MASS_EV_FLUKA        74438545206.6544            // 1000390800   //
#define  Y81_MASS_EV_FLUKA        75365188566.7801            // 1000390810   //
#define  Y82_MASS_EV_FLUKA        76294505736.8365            // 1000390820   //
#define  Y83_MASS_EV_FLUKA        77221861226.9438            // 1000390830   //
#define  Y84_MASS_EV_FLUKA        78151452166.9931            // 1000390840   //
#define  Y85_MASS_EV_FLUKA        79079333777.0867            // 1000390850   //
#define  Y86_MASS_EV_FLUKA        80009394117.1239            // 1000390860   //
#define  Y87_MASS_EV_FLUKA        80937153217.2207            // 1000390870   //
#define  Y88_MASS_EV_FLUKA        81867367237.2538            // 1000390880   //
#define  Y89_MASS_EV_FLUKA        82795454967.3421            // 1000390890   //
#define  Y90_MASS_EV_FLUKA        83728163477.3107            // 1000390900   //
#define  Y91_MASS_EV_FLUKA        84659796597.3071            // 1000390910   //
#define  Y93_MASS_EV_FLUKA        86524888857.2526            // 1000390930   //

#define  Zr82_MASS_EV_FLUKA       76298000372.0780            // 1000400820   //
#define  Zr84_MASS_EV_FLUKA       78153687272.2672            // 1000400840   //
#define  Zr85_MASS_EV_FLUKA       79083521472.3102            // 1000400850   //
#define  Zr86_MASS_EV_FLUKA       80010361772.4308            // 1000400860   //
#define  Zr87_MASS_EV_FLUKA       80940313242.4708            // 1000400870   //
#define  Zr88_MASS_EV_FLUKA       81867531732.5817            // 1000400880   //
#define  Zr89_MASS_EV_FLUKA       82797781932.6139            // 1000400890   //
#define  Zr90_MASS_EV_FLUKA       83725376152.7150            // 1000400900   //
#define  Zr91_MASS_EV_FLUKA       84657747182.6922            // 1000400910   //
#define  Zr92_MASS_EV_FLUKA       85588678102.7068            // 1000400920   //
#define  Zr93_MASS_EV_FLUKA       86521509622.6722            // 1000400930   //
#define  Zr94_MASS_EV_FLUKA       87452855142.6760            // 1000400940   //
#define  Zr95_MASS_EV_FLUKA       88385958162.6343            // 1000400950   //
#define  Zr96_MASS_EV_FLUKA       89317669872.6287            // 1000400960   //

#define  Nb84_MASS_EV_FLUKA       78162795038.0215            // 1000410840   //
#define  Nb86_MASS_EV_FLUKA       80017834598.2275            // 1000410860   //
#define  Nb87_MASS_EV_FLUKA       80944977138.3402            // 1000410870   //
#define  Nb88_MASS_EV_FLUKA       81874226518.3984            // 1000410880   //
#define  Nb89_MASS_EV_FLUKA       82801566698.5061            // 1000410890   //
#define  Nb90_MASS_EV_FLUKA       83730982008.5600            // 1000410900   //
#define  Nb91_MASS_EV_FLUKA       84658495448.6631            // 1000410910   //
#define  Nb92_MASS_EV_FLUKA       85590178668.6582            // 1000410920   //
#define  Nb93_MASS_EV_FLUKA       86520913298.6779            // 1000410930   //
#define  Nb94_MASS_EV_FLUKA       87453251518.6560            // 1000410940   //
#define  Nb95_MASS_EV_FLUKA       88384328638.6668            // 1000410950   //
#define  Nb96_MASS_EV_FLUKA       89317001138.6363            // 1000410960   //
#define  Nb97_MASS_EV_FLUKA       90248492768.6364            // 1000410970   //
#define  Nb99_MASS_EV_FLUKA       92114760918.5514            // 1000410990   //

#define  Mo87_MASS_EV_FLUKA       80950958556.7755            // 1000420870   //
#define  Mo88_MASS_EV_FLUKA       81877445376.9053            // 1000420880   //
#define  Mo89_MASS_EV_FLUKA       82806636896.9650            // 1000420890   //
#define  Mo90_MASS_EV_FLUKA       83732966077.0988            // 1000420900   //
#define  Mo91_MASS_EV_FLUKA       84662424347.1516            // 1000420910   //
#define  Mo92_MASS_EV_FLUKA       85589317647.2708            // 1000420920   //
#define  Mo93_MASS_EV_FLUKA       86520813667.2708            // 1000420930   //
#define  Mo94_MASS_EV_FLUKA       87450701497.3124            // 1000420940   //
#define  Mo95_MASS_EV_FLUKA       88382898117.2942            // 1000420950   //
#define  Mo96_MASS_EV_FLUKA       89313309437.3223            // 1000420960   //
#define  Mo97_MASS_EV_FLUKA       90246053957.2898            // 1000420970   //
#define  Mo99_MASS_EV_FLUKA       92110617397.2491            // 1000420990   //
#define  Mo100_MASS_EV_FLUKA      93041893477.2547            // 1000421000   //
#define  Mo101_MASS_EV_FLUKA      93976060597.1855            // 1000421010   //
#define  Mo102_MASS_EV_FLUKA      94907508957.1867            // 1000421020   //

#define  Tc90_MASS_EV_FLUKA       83741601419.9785            // 1000430900   //
#define  Tc91_MASS_EV_FLUKA       84668139620.1069            // 1000430910   //
#define  Tc92_MASS_EV_FLUKA       85596683100.1834            // 1000430920   //
#define  Tc93_MASS_EV_FLUKA       86523509740.3043            // 1000430930   //
#define  Tc94_MASS_EV_FLUKA       87454452560.3186            // 1000430940   //
#define  Tc95_MASS_EV_FLUKA       88384084760.3669            // 1000430950   //
#define  Tc96_MASS_EV_FLUKA       89315777980.3617            // 1000430960   //
#define  Tc97_MASS_EV_FLUKA       90245869520.3981            // 1000430970   //
#define  Tc98_MASS_EV_FLUKA       91178156440.3775            // 1000430980   //
#define  Tc99_MASS_EV_FLUKA       92108755380.4007            // 1000430990   //
#define  Tc100_MASS_EV_FLUKA      93041556700.3669            // 1000431000   //
#define  Tc101_MASS_EV_FLUKA      93972731520.3752            // 1000431010   //
#define  Tc102_MASS_EV_FLUKA      94905994260.3293            // 1000431020   //

#define  Ru91_MASS_EV_FLUKA       84675040179.4302            // 1000440910   //
#define  Ru93_MASS_EV_FLUKA       86529342189.6553            // 1000440930   //
#define  Ru94_MASS_EV_FLUKA       87455540969.7926            // 1000440940   //
#define  Ru95_MASS_EV_FLUKA       88386152699.8154            // 1000440950   //
#define  Ru96_MASS_EV_FLUKA       89315025349.8834            // 1000440960   //
#define  Ru97_MASS_EV_FLUKA       90246479459.8844            // 1000440970   //
#define  Ru98_MASS_EV_FLUKA       91175856109.9393            // 1000440980   //
#define  Ru99_MASS_EV_FLUKA       92107957369.9236            // 1000440990   //
#define  Ru100_MASS_EV_FLUKA      93037849889.9651            // 1000441000   //
#define  Ru101_MASS_EV_FLUKA      93970613409.9322            // 1000441010   //
#define  Ru102_MASS_EV_FLUKA      94900959529.9620            // 1000441020   //
#define  Ru103_MASS_EV_FLUKA      95834292749.9143            // 1000441030   //
#define  Ru104_MASS_EV_FLUKA      96764954849.9358            // 1000441040   //
#define  Ru106_MASS_EV_FLUKA      98629711249.8900            // 1000441060   //

#define  Rh94_MASS_EV_FLUKA       87464665737.3159            // 1000450940   //
#define  Rh95_MASS_EV_FLUKA       88390758397.4559            // 1000450950   //
#define  Rh96_MASS_EV_FLUKA       89320967507.4892            // 1000450960   //
#define  Rh97_MASS_EV_FLUKA       90249498207.5660            // 1000450970   //
#define  Rh98_MASS_EV_FLUKA       91180409257.5811            // 1000450980   //
#define  Rh99_MASS_EV_FLUKA       92109556077.6419            // 1000450990   //
#define  Rh100_MASS_EV_FLUKA      93040975587.6439            // 1000451000   //
#define  Rh101_MASS_EV_FLUKA      93970650717.6910            // 1000451010   //
#define  Rh102_MASS_EV_FLUKA      94902777797.6746            // 1000451020   //
#define  Rh103_MASS_EV_FLUKA      95833025137.7069            // 1000451030   //
#define  Rh104_MASS_EV_FLUKA      96765591757.6791            // 1000451040   //
#define  Rh105_MASS_EV_FLUKA      97696189057.7024            // 1000451050   //
#define  Rh106_MASS_EV_FLUKA      98629167547.6639            // 1000451060   //
#define  Rh107_MASS_EV_FLUKA      99560164227.6768            // 1000451070   //
#define  Rh108_MASS_EV_FLUKA     100493503257.6290            // 1000451080   //

#define  Pd96_MASS_EV_FLUKA       89323913415.2637            // 1000460960   //
#define  Pd97_MASS_EV_FLUKA       90253784125.3057            // 1000460970   //
#define  Pd98_MASS_EV_FLUKA       91181777755.3965            // 1000460980   //
#define  Pd99_MASS_EV_FLUKA       92112417285.4186            // 1000460990   //
#define  Pd100_MASS_EV_FLUKA      93040834405.4984            // 1000461000   //
#define  Pd101_MASS_EV_FLUKA      93972126615.5036            // 1000461010   //
#define  Pd102_MASS_EV_FLUKA      94901123225.5683            // 1000461020   //
#define  Pd103_MASS_EV_FLUKA      95833064145.5567            // 1000461030   //
#define  Pd104_MASS_EV_FLUKA      96762646645.6063            // 1000461040   //
#define  Pd105_MASS_EV_FLUKA      97695118275.5810            // 1000461050   //
#define  Pd106_MASS_EV_FLUKA      98625122385.6196            // 1000461060   //

#define  Ag98_MASS_EV_FLUKA       91189693874.9419            // 1000470980   //
#define  Ag99_MASS_EV_FLUKA       92117343405.0415            // 1000470990   //
#define  Ag100_MASS_EV_FLUKA      93047403955.0787            // 1000471000   //
#define  Ag101_MASS_EV_FLUKA      93975826545.1583            // 1000471010   //
#define  Ag102_MASS_EV_FLUKA      94906542515.1785            // 1000471020   //
#define  Ag103_MASS_EV_FLUKA      95835248275.2507            // 1000471030   //
#define  Ag104_MASS_EV_FLUKA      96766421455.2590            // 1000471040   //
#define  Ag105_MASS_EV_FLUKA      97695960025.3097            // 1000471050   //
#define  Ag106_MASS_EV_FLUKA      98627583805.3064            // 1000471060   //
#define  Ag107_MASS_EV_FLUKA      99557612115.3444            // 1000471070   //
#define  Ag108_MASS_EV_FLUKA     100489908245.3236            // 1000471080   //
#define  Ag109_MASS_EV_FLUKA     101420286685.3525            // 1000471090   //
#define  Ag110_MASS_EV_FLUKA     102353043105.3198            // 1000471100   //
#define  Ag111_MASS_EV_FLUKA     103283777625.3395            // 1000471110   //
#define  Ag114_MASS_EV_FLUKA     106081533085.2547            // 1000471140   //

#define  Cd101_MASS_EV_FLUKA      93980799478.4646            // 1000481010   //
#define  Cd102_MASS_EV_FLUKA      94908625958.5596            // 1000481020   //
#define  Cd103_MASS_EV_FLUKA      95838886508.5916            // 1000481030   //
#define  Cd104_MASS_EV_FLUKA      96767054648.6778            // 1000481040   //
#define  Cd105_MASS_EV_FLUKA      97698195058.6870            // 1000481050   //
#define  Cd106_MASS_EV_FLUKA      98626885528.7596            // 1000481060   //
#define  Cd107_MASS_EV_FLUKA      99558525438.7558            // 1000481070   //
#define  Cd108_MASS_EV_FLUKA     100487755068.8145            // 1000481080   //
#define  Cd109_MASS_EV_FLUKA     101419996808.7951            // 1000481090   //
#define  Cd110_MASS_EV_FLUKA     102349647338.8429            // 1000481100   //
#define  Cd111_MASS_EV_FLUKA     103282237158.8145            // 1000481110   //
#define  Cd112_MASS_EV_FLUKA     104212404678.8489            // 1000481120   //
#define  Cd113_MASS_EV_FLUKA     105145430098.8092            // 1000481130   //
#define  Cd114_MASS_EV_FLUKA     106075953018.8344            // 1000481140   //
#define  Cd116_MASS_EV_FLUKA     107940243158.8007            // 1000481160   //
#define  Cd119_MASS_EV_FLUKA     110739540488.6759            // 1000481190   //

#define  In103_MASS_EV_FLUKA      95844433057.3294            // 1000491030   //
#define  In104_MASS_EV_FLUKA      96774459987.3674            // 1000491040   //
#define  In105_MASS_EV_FLUKA      97702540597.4559            // 1000491050   //
#define  In106_MASS_EV_FLUKA      98632903427.4852            // 1000491060   //
#define  In107_MASS_EV_FLUKA      99561448047.5616            // 1000491070   //
#define  In108_MASS_EV_FLUKA     100492399267.5757            // 1000491080   //
#define  In109_MASS_EV_FLUKA     101421513637.6374            // 1000491090   //
#define  In110_MASS_EV_FLUKA     102353021807.6370            // 1000491100   //
#define  In111_MASS_EV_FLUKA     103282599587.6867            // 1000491110   //
#define  In112_MASS_EV_FLUKA     104214487407.6765            // 1000491120   //
#define  In113_MASS_EV_FLUKA     105144610747.7120            // 1000491130   //
#define  In114_MASS_EV_FLUKA     106076902057.6913            // 1000491140   //
#define  In115_MASS_EV_FLUKA     107007428677.7164            // 1000491150   //
#define  In116_MASS_EV_FLUKA     107940209987.6831            // 1000491160   //
#define  In118_MASS_EV_FLUKA     109804219297.6566            // 1000491180   //
#define  In120_MASS_EV_FLUKA     111668704877.6178            // 1000491200   //

#define  Sn106_MASS_EV_FLUKA      98635584553.4832            // 1000501060   //
#define  Sn107_MASS_EV_FLUKA      99565944443.5126            // 1000501070   //
#define  Sn108_MASS_EV_FLUKA     100493988113.6020            // 1000501080   //
#define  Sn109_MASS_EV_FLUKA     101424860763.6181            // 1000501090   //
#define  Sn110_MASS_EV_FLUKA     102353156263.7010            // 1000501100   //
#define  Sn111_MASS_EV_FLUKA     103284541233.7039            // 1000501110   //
#define  Sn112_MASS_EV_FLUKA     104213320683.7742            // 1000501120   //
#define  Sn113_MASS_EV_FLUKA     105145143403.7657            // 1000501130   //
#define  Sn114_MASS_EV_FLUKA     106074410133.8234            // 1000501140   //
#define  Sn115_MASS_EV_FLUKA     107006430153.8098            // 1000501150   //
#define  Sn117_MASS_EV_FLUKA     108869053493.8193            // 1000501170   //

#define  Sb109_MASS_EV_FLUKA     101430737828.4374            // 1000511090   //
#define  Sb110_MASS_EV_FLUKA     102360952768.4705            // 1000511100   //
#define  Sb111_MASS_EV_FLUKA     103289138088.5563            // 1000511110   //
#define  Sb112_MASS_EV_FLUKA     104219872688.5760            // 1000511120   //
#define  Sb113_MASS_EV_FLUKA     105148545898.6491            // 1000511130   //
#define  Sb114_MASS_EV_FLUKA     106079788628.6556            // 1000511140   //
#define  Sb115_MASS_EV_FLUKA     107008957148.7159            // 1000511150   //
#define  Sb116_MASS_EV_FLUKA     107940636728.7111            // 1000511160   //
#define  Sb118_MASS_EV_FLUKA     109802446488.7416            // 1000511180   //
#define  Sb120_MASS_EV_FLUKA     111665009198.7527            // 1000511200   //

#define  Te112_MASS_EV_FLUKA     104223713694.0493            // 1000521120   //
#define  Te113_MASS_EV_FLUKA     105154142934.0768            // 1000521130   //
#define  Te114_MASS_EV_FLUKA     106082028254.1704            // 1000521140   //
#define  Te115_MASS_EV_FLUKA     107013093064.1815            // 1000521150   //
#define  Te116_MASS_EV_FLUKA     107941633804.2580            // 1000521160   //
#define  Te117_MASS_EV_FLUKA     108873339094.2526            // 1000521170   //
#define  Te118_MASS_EV_FLUKA     109802221414.3200            // 1000521180   //
#define  Te119_MASS_EV_FLUKA     110734254164.3063            // 1000521190   //
#define  Te120_MASS_EV_FLUKA     111663528464.3638            // 1000521200   //
#define  Te121_MASS_EV_FLUKA     112595868964.3419            // 1000521210   //
#define  Te122_MASS_EV_FLUKA     113525613084.3872            // 1000521220   //

#define  I113_MASS_EV_FLUKA      105160840511.7544            // 1000531130   //
#define  I114_MASS_EV_FLUKA      106090662581.7977            // 1000531140   //
#define  I115_MASS_EV_FLUKA      107018549701.8912            // 1000531150   //
#define  I116_MASS_EV_FLUKA      107948876511.9214            // 1000531160   //
#define  I117_MASS_EV_FLUKA      108877490171.9961            // 1000531170   //
#define  I118_MASS_EV_FLUKA      109808763282.0018            // 1000531180   //
#define  I119_MASS_EV_FLUKA      110737265222.0793            // 1000531190   //
#define  I120_MASS_EV_FLUKA      111668640882.0824            // 1000531200   //
#define  I121_MASS_EV_FLUKA      112597637212.1472            // 1000531210   //
#define  I123_MASS_EV_FLUKA      114458978922.1898            // 1000531230   //

#define  Xe116_MASS_EV_FLUKA     107953033923.6005            // 1000541160   //
#define  Xe117_MASS_EV_FLUKA     108883433973.6289            // 1000541170   //
#define  Xe118_MASS_EV_FLUKA     109811205703.7253            // 1000541180   //
#define  Xe119_MASS_EV_FLUKA     110741764523.7496            // 1000541190   //
#define  Xe120_MASS_EV_FLUKA     111670098513.8315            // 1000541200   //
#define  Xe121_MASS_EV_FLUKA     112600867153.8503            // 1000541210   //
#define  Xe122_MASS_EV_FLUKA     113529736613.9183            // 1000541220   //
#define  Xe123_MASS_EV_FLUKA     114461152593.9203            // 1000541230   //
#define  Xe124_MASS_EV_FLUKA     115390241913.9827            // 1000541240   //
#define  Xe125_MASS_EV_FLUKA     116322204233.9706            // 1000541250   //
#define  Xe126_MASS_EV_FLUKA     117251714714.0219            // 1000541260   //
#define  Xe131_MASS_EV_FLUKA     121909945154.0023            // 1000541310   //

#define  Cs118_MASS_EV_FLUKA     109820003570.8584            // 1000551180   //
#define  Cs119_MASS_EV_FLUKA     110747590290.9597            // 1000551190   //
#define  Cs120_MASS_EV_FLUKA     111677518311.0003            // 1000551200   //
#define  Cs121_MASS_EV_FLUKA     112605765041.0844            // 1000551210   //
#define  Cs122_MASS_EV_FLUKA     113536289261.1095            // 1000551220   //
#define  Cs123_MASS_EV_FLUKA     114464850781.1856            // 1000551230   //
#define  Cs124_MASS_EV_FLUKA     115395656801.2034            // 1000551240   //
#define  Cs125_MASS_EV_FLUKA     116324794241.2645            // 1000551250   //
#define  Cs127_MASS_EV_FLUKA     118185635791.3201            // 1000551270   //
#define  Cs128_MASS_EV_FLUKA     119117444231.3120            // 1000551280   //
#define  Cs129_MASS_EV_FLUKA     120047367661.3527            // 1000551290   //
#define  Cs130_MASS_EV_FLUKA     120979465441.3371            // 1000551300   //

#define  Ba120_MASS_EV_FLUKA     111682016405.4391            // 1000561200   //
#define  Ba121_MASS_EV_FLUKA     112612077825.4762            // 1000561210   //
#define  Ba122_MASS_EV_FLUKA     113539629775.5784            // 1000561220   //
#define  Ba123_MASS_EV_FLUKA     114469810095.6125            // 1000561230   //
#define  Ba124_MASS_EV_FLUKA     115397801405.7032            // 1000561240   //
#define  Ba125_MASS_EV_FLUKA     116328852305.7147            // 1000561250   //
#define  Ba126_MASS_EV_FLUKA     117257209035.7960            // 1000561260   //
#define  Ba127_MASS_EV_FLUKA     118188583965.7991            // 1000561270   //
#define  Ba128_MASS_EV_FLUKA     119117463585.8669            // 1000561280   //
#define  Ba129_MASS_EV_FLUKA     120049299105.8580            // 1000561290   //
#define  Ba131_MASS_EV_FLUKA     121910662885.9002            // 1000561310   //

#define  La121_MASS_EV_FLUKA     112619509388.6356            // 1000571210   //
#define  La122_MASS_EV_FLUKA     113548861818.6911            // 1000571220   //
#define  La123_MASS_EV_FLUKA     114476192228.7990            // 1000571230   //
#define  La124_MASS_EV_FLUKA     115406093648.8403            // 1000571240   //
#define  La125_MASS_EV_FLUKA     116333993078.9335            // 1000571250   //
#define  La126_MASS_EV_FLUKA     117264276388.9648            // 1000571260   //
#define  La127_MASS_EV_FLUKA     118192780789.0423            // 1000571270   //
#define  La128_MASS_EV_FLUKA     119123611889.0595            // 1000571280   //
#define  La129_MASS_EV_FLUKA     120052517399.1266            // 1000571290   //
#define  La130_MASS_EV_FLUKA     120983686769.1350            // 1000571300   //
#define  La131_MASS_EV_FLUKA     121913121279.1884            // 1000571310   //
#define  La132_MASS_EV_FLUKA     122844617259.1884            // 1000571320   //
#define  La133_MASS_EV_FLUKA     123774515009.2297            // 1000571330   //
#define  La134_MASS_EV_FLUKA     124706096119.2275            // 1000571340   //

#define  Ce124_MASS_EV_FLUKA     115411171912.4421            // 1000581240   //
#define  Ce126_MASS_EV_FLUKA     117268180752.5970            // 1000581260   //
#define  Ce127_MASS_EV_FLUKA     118198417172.6297            // 1000581270   //
#define  Ce128_MASS_EV_FLUKA     119126297592.7233            // 1000581280   //
#define  Ce129_MASS_EV_FLUKA     120057066012.7421            // 1000581290   //
#define  Ce130_MASS_EV_FLUKA     120985392332.8242            // 1000581300   //
#define  Ce131_MASS_EV_FLUKA     121916639812.8306            // 1000581310   //
#define  Ce132_MASS_EV_FLUKA     122845399972.9014            // 1000581320   //
#define  Ce133_MASS_EV_FLUKA     123776950292.9000            // 1000581330   //
#define  Ce134_MASS_EV_FLUKA     124706094692.9609            // 1000581340   //
#define  Ce135_MASS_EV_FLUKA     125637700822.9580            // 1000581350   //
#define  Ce137_MASS_EV_FLUKA     127499414822.9910            // 1000581370   //
#define  Ce139_MASS_EV_FLUKA     129361349943.0183            // 1000581390   //

#define  Pr125_MASS_EV_FLUKA     116348465726.4607            // 1000591250   //
#define  Pr127_MASS_EV_FLUKA     118205444368.1303            // 1000591270   //
#define  Pr128_MASS_EV_FLUKA     119135046188.1794            // 1000591280   //
#define  Pr129_MASS_EV_FLUKA     120062870608.2745            // 1000591290   //
#define  Pr130_MASS_EV_FLUKA     120992986028.3102            // 1000591300   //
#define  Pr131_MASS_EV_FLUKA     121921388608.3903            // 1000591310   //
#define  Pr132_MASS_EV_FLUKA     122852006768.4130            // 1000591320   //
#define  Pr133_MASS_EV_FLUKA     123780781088.4835            // 1000591330   //
#define  Pr134_MASS_EV_FLUKA     124711800418.4958            // 1000591340   //
#define  Pr135_MASS_EV_FLUKA     125640919508.5574            // 1000591350   //
#define  Pr136_MASS_EV_FLUKA     126571954998.5693            // 1000591360   //
#define  Pr137_MASS_EV_FLUKA     127501615618.6168            // 1000591370   //
#define  Pr138_MASS_EV_FLUKA     128433175668.6151            // 1000591380   //

#define  Nd129_MASS_EV_FLUKA     120070193237.2685            // 1000601290   //
#define  Nd130_MASS_EV_FLUKA     120997514657.3766            // 1000601300   //
#define  Nd131_MASS_EV_FLUKA     121927447627.4171            // 1000601310   //
#define  Nd132_MASS_EV_FLUKA     122855231497.5132            // 1000601320   //
#define  Nd133_MASS_EV_FLUKA     123785877917.5352            // 1000601330   //
#define  Nd134_MASS_EV_FLUKA     124714069397.6208            // 1000601340   //
#define  Nd135_MASS_EV_FLUKA     125645168757.6310            // 1000601350   //
#define  Nd136_MASS_EV_FLUKA     126573665017.7087            // 1000601360   //
#define  Nd137_MASS_EV_FLUKA     127504804627.7179            // 1000601370   //
#define  Nd138_MASS_EV_FLUKA     128433774717.7833            // 1000601380   //
#define  Nd139_MASS_EV_FLUKA     129365263297.7835            // 1000601390   //
#define  Nd140_MASS_EV_FLUKA     130294322737.8466            // 1000601400   //

#define  Pm132_MASS_EV_FLUKA     122864632561.4900            // 1000611320   //
#define  Pm133_MASS_EV_FLUKA     123792372981.5872            // 1000611330   //
#define  Pm134_MASS_EV_FLUKA     124722451401.6240            // 1000611340   //
#define  Pm135_MASS_EV_FLUKA     125650685821.7085            // 1000611350   //
#define  Pm136_MASS_EV_FLUKA     126581014321.7386            // 1000611360   //
#define  Pm137_MASS_EV_FLUKA     127509883691.8067            // 1000611370   //
#define  Pm138_MASS_EV_FLUKA     128440173741.8379            // 1000611380   //
#define  Pm139_MASS_EV_FLUKA     129369284551.8996            // 1000611390   //
#define  Pm140_MASS_EV_FLUKA     130299910671.9222            // 1000611400   //
#define  Pm141_MASS_EV_FLUKA     131229306121.9765            // 1000611410   //
#define  Pm145_MASS_EV_FLUKA     134954492581.9971            // 1000611450   //

#define  Sm134_MASS_EV_FLUKA     124727371702.2735            // 1000621340   //
#define  Sm136_MASS_EV_FLUKA     126585032532.4116            // 1000621360   //
#define  Sm137_MASS_EV_FLUKA     127515436962.4398            // 1000621370   //
#define  Sm138_MASS_EV_FLUKA     128443587382.5265            // 1000621380   //
#define  Sm139_MASS_EV_FLUKA     129374244092.5482            // 1000621390   //
#define  Sm140_MASS_EV_FLUKA     130302430122.6339            // 1000621400   //
#define  Sm141_MASS_EV_FLUKA     131233348732.6488            // 1000621410   //
#define  Sm142_MASS_EV_FLUKA     132161799752.7277            // 1000621420   //
#define  Sm143_MASS_EV_FLUKA     133092754442.7417            // 1000621430   //
#define  Sm144_MASS_EV_FLUKA     134021800272.8051            // 1000621440   //
#define  Sm145_MASS_EV_FLUKA     134954608792.7711            // 1000621450   //

#define  Eu137_MASS_EV_FLUKA     127522463501.0964            // 1000631370   //
#define  Eu138_MASS_EV_FLUKA     128452317921.1389            // 1000631380   //
#define  Eu139_MASS_EV_FLUKA     129380421341.2268            // 1000631390   //
#define  Eu140_MASS_EV_FLUKA     130310329611.2679            // 1000631400   //
#define  Eu141_MASS_EV_FLUKA     131238398471.3566            // 1000631410   //
#define  Eu142_MASS_EV_FLUKA     132168659401.3886            // 1000631420   //
#define  Eu143_MASS_EV_FLUKA     133097422581.4594            // 1000631430   //
#define  Eu144_MASS_EV_FLUKA     134027628821.4928            // 1000631440   //
#define  Eu145_MASS_EV_FLUKA     134956768711.5538            // 1000631450   //
#define  Eu146_MASS_EV_FLUKA     135889136711.5311            // 1000631460   //
#define  Eu147_MASS_EV_FLUKA     136820203771.5423            // 1000631470   //
#define  Eu148_MASS_EV_FLUKA     137753013801.5081            // 1000631480   //

#define  Gd139_MASS_EV_FLUKA     129387625219.4307            // 1000641390   //
#define  Gd140_MASS_EV_FLUKA     130315289439.5299            // 1000641400   //
#define  Gd141_MASS_EV_FLUKA     131245146059.5724            // 1000641410   //
#define  Gd142_MASS_EV_FLUKA     132172659369.6755            // 1000641420   //
#define  Gd143_MASS_EV_FLUKA     133102929979.7073            // 1000641430   //
#define  Gd144_MASS_EV_FLUKA     134030868119.7994            // 1000641440   //
#define  Gd145_MASS_EV_FLUKA     134961323109.8263            // 1000641450   //
#define  Gd146_MASS_EV_FLUKA     135889666909.9080            // 1000641460   //
#define  Gd147_MASS_EV_FLUKA     136821891539.8890            // 1000641470   //
#define  Gd148_MASS_EV_FLUKA     137752473359.9127            // 1000641480   //
#define  Gd149_MASS_EV_FLUKA     138685112069.8830            // 1000641490   //
#define  Gd150_MASS_EV_FLUKA     139615969869.8995            // 1000641500   //
#define  Gd155_MASS_EV_FLUKA     144277132539.8039            // 1000641550   //

#define  Tb143_MASS_EV_FLUKA     133109823818.9469            // 1000651430   //
#define  Tb144_MASS_EV_FLUKA     134039287238.9995            // 1000651440   //
#define  Tb145_MASS_EV_FLUKA     134967331349.0889            // 1000651450   //
#define  Tb146_MASS_EV_FLUKA     135897242169.1300            // 1000651460   //
#define  Tb147_MASS_EV_FLUKA     136826003189.2008            // 1000651470   //
#define  Tb148_MASS_EV_FLUKA     137757666589.1964            // 1000651480   //
#define  Tb149_MASS_EV_FLUKA     138688248389.2201            // 1000651490   //
#define  Tb150_MASS_EV_FLUKA     139620126479.2101            // 1000651500   //
#define  Tb155_MASS_EV_FLUKA     144277454249.2138            // 1000651550   //

#define  Dy143_MASS_EV_FLUKA     133118089180.6405            // 1000661430   //
#define  Dy146_MASS_EV_FLUKA     135901902630.9170            // 1000661460   //
#define  Dy147_MASS_EV_FLUKA     136831876210.9564            // 1000661470   //
#define  Dy148_MASS_EV_FLUKA     137759845151.0478            // 1000661480   //
#define  Dy149_MASS_EV_FLUKA     138691560891.0420            // 1000661490   //
#define  Dy150_MASS_EV_FLUKA     139621420871.0844            // 1000661500   //
#define  Dy151_MASS_EV_FLUKA     140553473901.0699            // 1000661510   //
#define  Dy152_MASS_EV_FLUKA     141483603001.1053            // 1000661520   //
#define  Dy153_MASS_EV_FLUKA     142416073431.0800            // 1000661530   //
#define  Dy154_MASS_EV_FLUKA     143346319311.1123            // 1000661540   //
#define  Dy156_MASS_EV_FLUKA     145209173471.1158            // 1000661560   //

#define  Ho147_MASS_EV_FLUKA     136839525266.6038            // 1000671470   //
#define  Ho148_MASS_EV_FLUKA     137768745716.6627            // 1000671480   //
#define  Ho149_MASS_EV_FLUKA     138697075286.7448            // 1000671490   //
#define  Ho150_MASS_EV_FLUKA     139628161526.7553            // 1000671500   //
#define  Ho151_MASS_EV_FLUKA     140558102336.7956            // 1000671510   //
#define  Ho152_MASS_EV_FLUKA     141489577436.7961            // 1000671520   //
#define  Ho153_MASS_EV_FLUKA     142419702816.8316            // 1000671530   //
#define  Ho154_MASS_EV_FLUKA     143351571416.8219            // 1000671540   //
#define  Ho155_MASS_EV_FLUKA     144281651906.8585            // 1000671550   //
#define  Ho157_MASS_EV_FLUKA     146143810516.8800            // 1000671570   //
#define  Ho159_MASS_EV_FLUKA     148006352366.8915            // 1000671590   //

#define  Er149_MASS_EV_FLUKA     138704309787.7756            // 1000681490   //
#define  Er150_MASS_EV_FLUKA     139631770477.8801            // 1000681500   //
#define  Er151_MASS_EV_FLUKA     140562823597.8915            // 1000681510   //
#define  Er152_MASS_EV_FLUKA     141492183387.9469            // 1000681520   //
#define  Er153_MASS_EV_FLUKA     142423767427.9445            // 1000681530   //
#define  Er154_MASS_EV_FLUKA     143353104298.0004            // 1000681540   //
#define  Er155_MASS_EV_FLUKA     144284996327.9901            // 1000681550   //
#define  Er156_MASS_EV_FLUKA     145214605248.0390            // 1000681560   //
#define  Er157_MASS_EV_FLUKA     146146781428.0213            // 1000681570   //
#define  Er159_MASS_EV_FLUKA     148008621808.0511            // 1000681590   //
#define  Er162_MASS_EV_FLUKA     150801329878.0970            // 1000681620   //

#define  Tm151_MASS_EV_FLUKA     140569854735.7222            // 1000691510   //
#define  Tm152_MASS_EV_FLUKA     141500349115.7481            // 1000691520   //
#define  Tm153_MASS_EV_FLUKA     142429727715.8029            // 1000691530   //
#define  Tm154_MASS_EV_FLUKA     143360658945.8175            // 1000691540   //
#define  Tm155_MASS_EV_FLUKA     144290075955.8713            // 1000691550   //
#define  Tm156_MASS_EV_FLUKA     145221325335.8777            // 1000691560   //
#define  Tm157_MASS_EV_FLUKA     146150762705.9309            // 1000691570   //
#define  Tm158_MASS_EV_FLUKA     147082449195.9260            // 1000691580   //
#define  Tm159_MASS_EV_FLUKA     148011973045.9770            // 1000691590   //
#define  Tm160_MASS_EV_FLUKA     148944015565.9629            // 1000691600   //
#define  Tm161_MASS_EV_FLUKA     149873643106.0112            // 1000691610   //
#define  Tm162_MASS_EV_FLUKA     150805640475.9982            // 1000691620   //

#define  Yb154_MASS_EV_FLUKA     143364649401.9305            // 1000701540   //
#define  Yb155_MASS_EV_FLUKA     144295565521.9455            // 1000701550   //
#define  Yb156_MASS_EV_FLUKA     145224399712.0145            // 1000701560   //
#define  Yb157_MASS_EV_FLUKA     146155794512.0170            // 1000701570   //
#define  Yb158_MASS_EV_FLUKA     147084679682.0847            // 1000701580   //
#define  Yb159_MASS_EV_FLUKA     148016524482.0755            // 1000701590   //
#define  Yb160_MASS_EV_FLUKA     148945527222.1401            // 1000701600   //
#define  Yb161_MASS_EV_FLUKA     149877294522.1331            // 1000701610   //
#define  Yb162_MASS_EV_FLUKA     150806829852.1838            // 1000701620   //
#define  Yb163_MASS_EV_FLUKA     151738804172.1714            // 1000701630   //
#define  Yb164_MASS_EV_FLUKA     152668672602.2135            // 1000701640   //
#define  Yb165_MASS_EV_FLUKA     153600984502.1924            // 1000701650   //
#define  Yb166_MASS_EV_FLUKA     154531064762.2289            // 1000701660   //
#define  Yb173_MASS_EV_FLUKA     161055555752.1245            // 1000701730   //

#define  Lu156_MASS_EV_FLUKA     145233347537.5995            // 1000711560   //
#define  Lu159_MASS_EV_FLUKA     148022014047.7503            // 1000711590   //
#define  Lu160_MASS_EV_FLUKA     148952908887.7658            // 1000711600   //
#define  Lu161_MASS_EV_FLUKA     149882096197.8256            // 1000711610   //
#define  Lu162_MASS_EV_FLUKA     150813551527.8266            // 1000711620   //
#define  Lu163_MASS_EV_FLUKA     151742905857.8821            // 1000711630   //
#define  Lu164_MASS_EV_FLUKA     152674424257.8814            // 1000711640   //
#define  Lu165_MASS_EV_FLUKA     153604406137.9206            // 1000711650   //

#define  Hf161_MASS_EV_FLUKA     149887922014.4767            // 1000721610   //
#define  Hf162_MASS_EV_FLUKA     150816502594.5522            // 1000721620   //
#define  Hf163_MASS_EV_FLUKA     151747857704.5558            // 1000721630   //
#define  Hf164_MASS_EV_FLUKA     152676900144.6193            // 1000721640   //
#define  Hf165_MASS_EV_FLUKA     153608503294.6165            // 1000721650   //
#define  Hf166_MASS_EV_FLUKA     154537864684.6718            // 1000721660   //
#define  Hf167_MASS_EV_FLUKA     155469685104.6634            // 1000721670   //
#define  Hf168_MASS_EV_FLUKA     156399344524.7109            // 1000721680   //
#define  Hf169_MASS_EV_FLUKA     157331331664.6982            // 1000721690   //
#define  Hf170_MASS_EV_FLUKA     158261420074.7346            // 1000721700   //
#define  Hf171_MASS_EV_FLUKA     159193697384.7143            // 1000721710   //

#define  Ta159_MASS_EV_FLUKA     148036183583.3451            // 1000731590   //
#define  Ta162_MASS_EV_FLUKA     150825267773.4850            // 1000731620   //
#define  Ta164_MASS_EV_FLUKA     152684923143.5714            // 1000731640   //
#define  Ta165_MASS_EV_FLUKA     153613853643.6379            // 1000731650   //
#define  Ta166_MASS_EV_FLUKA     154545023883.6462            // 1000731660   //
#define  Ta167_MASS_EV_FLUKA     155474191063.7066            // 1000731670   //
#define  Ta168_MASS_EV_FLUKA     156405516453.7109            // 1000731680   //
#define  Ta169_MASS_EV_FLUKA     157335268933.7561            // 1000731690   //
#define  Ta170_MASS_EV_FLUKA     158266921263.7520            // 1000731700   //
#define  Ta171_MASS_EV_FLUKA     159196897573.7913            // 1000731710   //
#define  Ta172_MASS_EV_FLUKA     160128652783.7845            // 1000731720   //
#define  Ta174_MASS_EV_FLUKA     161991109663.7983            // 1000731740   //

#define  W163_MASS_EV_FLUKA      151761121366.2512            // 1000741630   //
#define  W165_MASS_EV_FLUKA      153620360626.3483            // 1000741650   //
#define  W166_MASS_EV_FLUKA      154548765926.4284            // 1000741660   //
#define  W167_MASS_EV_FLUKA      155479933636.4369            // 1000741670   //
#define  W168_MASS_EV_FLUKA      156408813076.5046            // 1000741680   //
#define  W169_MASS_EV_FLUKA      157340210286.5071            // 1000741690   //
#define  W170_MASS_EV_FLUKA      158269402446.5667            // 1000741700   //
#define  W171_MASS_EV_FLUKA      159200972956.5648            // 1000741710   //
#define  W172_MASS_EV_FLUKA      160130655286.6117            // 1000741720   //
#define  W173_MASS_EV_FLUKA      161062629506.5993            // 1000741730   //
#define  W174_MASS_EV_FLUKA      161992465896.6423            // 1000741740   //
#define  W175_MASS_EV_FLUKA      162924529326.6276            // 1000741750   //
#define  W176_MASS_EV_FLUKA      163854923636.6560            // 1000741760   //
#define  W177_MASS_EV_FLUKA      164787377856.6312            // 1000741770   //
#define  W178_MASS_EV_FLUKA      165718153866.6498            // 1000741780   //
#define  W182_MASS_EV_FLUKA      169446326326.5929            // 1000741820   //

#define  Re169_MASS_EV_FLUKA     157346299104.2974            // 1000751690   //
#define  Re170_MASS_EV_FLUKA     158277172214.3135            // 1000751700   //
#define  Re171_MASS_EV_FLUKA     159206145604.3789            // 1000751710   //
#define  Re172_MASS_EV_FLUKA     160137483954.3829            // 1000751720   //
#define  Re173_MASS_EV_FLUKA     161066904134.4366            // 1000751730   //
#define  Re174_MASS_EV_FLUKA     161998445494.4354            // 1000751740   //
#define  Re176_MASS_EV_FLUKA     163859997344.4727            // 1000751760   //
#define  Re177_MASS_EV_FLUKA     164790280664.5040            // 1000751770   //
#define  Re178_MASS_EV_FLUKA     165722316564.4900            // 1000751780   //
#define  Re179_MASS_EV_FLUKA     166652999954.5110            // 1000751790   //
#define  Re180_MASS_EV_FLUKA     167585245594.4916            // 1000751800   //
#define  Re181_MASS_EV_FLUKA     168516066334.5090            // 1000751810   //
#define  Re182_MASS_EV_FLUKA     169448629054.4813            // 1000751820   //
#define  Re183_MASS_EV_FLUKA     170379759904.4908            // 1000751830   //
#define  Re187_MASS_EV_FLUKA     174110328834.3718            // 1000751870   //
#define  Re192_MASS_EV_FLUKA     178777310464.1253            // 1000751920   //

#define  Os172_MASS_EV_FLUKA     160141448038.6361            // 1000761720   //
#define  Os173_MASS_EV_FLUKA     161072675238.6430            // 1000761730   //
#define  Os174_MASS_EV_FLUKA     162001682398.7075            // 1000761740   //
#define  Os175_MASS_EV_FLUKA     162933093888.7097            // 1000761750   //
#define  Os176_MASS_EV_FLUKA     163862664308.7595            // 1000761760   //
#define  Os177_MASS_EV_FLUKA     164794256548.7569            // 1000761770   //
#define  Os178_MASS_EV_FLUKA     165724154328.7983            // 1000761780   //
#define  Os179_MASS_EV_FLUKA     166656181228.7845            // 1000761790   //
#define  Os180_MASS_EV_FLUKA     167586214608.8224            // 1000761800   //
#define  Os181_MASS_EV_FLUKA     168518498888.8019            // 1000761810   //
#define  Os182_MASS_EV_FLUKA     169449039998.8266            // 1000761820   //
#define  Os183_MASS_EV_FLUKA     170381394648.8043            // 1000761830   //
#define  Os184_MASS_EV_FLUKA     171312312238.8192            // 1000761840   //

#define  Ir174_MASS_EV_FLUKA     162010204690.5783            // 1000771740   //
#define  Ir176_MASS_EV_FLUKA     163870129420.6577            // 1000771760   //
#define  Ir177_MASS_EV_FLUKA     164799439600.7144            // 1000771770   //
#define  Ir178_MASS_EV_FLUKA     165730854010.7164            // 1000771780   //
#define  Ir179_MASS_EV_FLUKA     166660546290.7631            // 1000771790   //
#define  Ir180_MASS_EV_FLUKA     167592134710.7607            // 1000771800   //
#define  Ir181_MASS_EV_FLUKA     168522071090.8010            // 1000771810   //
#define  Ir182_MASS_EV_FLUKA     169454152630.7858            // 1000771820   //
#define  Ir183_MASS_EV_FLUKA     170384347830.8195            // 1000771830   //
#define  Ir184_MASS_EV_FLUKA     171316377460.8056            // 1000771840   //
#define  Ir185_MASS_EV_FLUKA     172247128410.8249            // 1000771850   //
#define  Ir186_MASS_EV_FLUKA     173179890410.7921            // 1000771860   //
#define  Ir187_MASS_EV_FLUKA     174110834780.8063            // 1000771870   //
#define  Ir188_MASS_EV_FLUKA     175043718090.7703            // 1000771880   //
#define  Ir194_MASS_EV_FLUKA     180638481160.6201            // 1000771940   //
#define  Ir195_MASS_EV_FLUKA     181570814980.5983            // 1000771950   //

#define  Pt176_MASS_EV_FLUKA     163874743011.6830            // 1000781760   //
#define  Pt177_MASS_EV_FLUKA     164805727221.6962            // 1000781770   //
#define  Pt178_MASS_EV_FLUKA     165734665381.7625            // 1000781780   //
#define  Pt179_MASS_EV_FLUKA     166665783871.7722            // 1000781790   //
#define  Pt180_MASS_EV_FLUKA     167595330291.8227            // 1000781800   //
#define  Pt181_MASS_EV_FLUKA     168526798541.8233            // 1000781810   //
#define  Pt182_MASS_EV_FLUKA     169456505811.8696            // 1000781820   //
#define  Pt183_MASS_EV_FLUKA     170388429221.8585            // 1000781830   //
#define  Pt184_MASS_EV_FLUKA     171318213591.9028            // 1000781840   //
#define  Pt185_MASS_EV_FLUKA     172250449861.8836            // 1000781850   //
#define  Pt186_MASS_EV_FLUKA     173180772881.9140            // 1000781860   //
#define  Pt187_MASS_EV_FLUKA     174113446571.8834            // 1000781870   //
#define  Pt188_MASS_EV_FLUKA     175043728091.9148            // 1000781880   //
#define  Pt189_MASS_EV_FLUKA     175976560161.8802            // 1000781890   //
#define  Pt190_MASS_EV_FLUKA     176907214431.9019            // 1000781900   //
#define  Pt196_MASS_EV_FLUKA     182500842081.7811            // 1000781960   //
#define  Pt197_MASS_EV_FLUKA     183434561301.7235            // 1000781970   //
#define  Pt198_MASS_EV_FLUKA     184366570411.7101            // 1000781980   //

#define  Au181_MASS_EV_FLUKA     168532601103.1779            // 1000791810   //
#define  Au182_MASS_EV_FLUKA     169463789513.1858            // 1000791820   //
#define  Au183_MASS_EV_FLUKA     170393421793.2340            // 1000791830   //
#define  Au184_MASS_EV_FLUKA     171324844213.2360            // 1000791840   //
#define  Au185_MASS_EV_FLUKA     172254660603.2794            // 1000791850   //
#define  Au186_MASS_EV_FLUKA     173186317243.2752            // 1000791860   //
#define  Au187_MASS_EV_FLUKA     174116550323.3079            // 1000791870   //
#define  Au188_MASS_EV_FLUKA     175048531693.2953            // 1000791880   //
#define  Au189_MASS_EV_FLUKA     175978913913.3241            // 1000791890   //
#define  Au190_MASS_EV_FLUKA     176911160213.3046            // 1000791900   //
#define  Au191_MASS_EV_FLUKA     177841676903.3299            // 1000791910   //
#define  Au192_MASS_EV_FLUKA     178774252453.3019            // 1000791920   //
#define  Au193_MASS_EV_FLUKA     179705114703.3183            // 1000791930   //
#define  Au194_MASS_EV_FLUKA     180637733603.2891            // 1000791940   //
#define  Au195_MASS_EV_FLUKA     181568928903.2969            // 1000791950   //
#define  Au197_MASS_EV_FLUKA     183433346143.2599            // 1000791970   //
#define  Au198_MASS_EV_FLUKA     184366399563.2195            // 1000791980   //

#define  Hg180_MASS_EV_FLUKA     167608411036.1504            // 1000801800   //
#define  Hg185_MASS_EV_FLUKA     172259986556.3031            // 1000801850   //
#define  Hg186_MASS_EV_FLUKA     173189121536.3642            // 1000801860   //
#define  Hg187_MASS_EV_FLUKA     174120919236.3564            // 1000801870   //
#define  Hg188_MASS_EV_FLUKA     175050333616.4103            // 1000801880   //
#define  Hg189_MASS_EV_FLUKA     175982367836.3963            // 1000801890   //
#define  Hg190_MASS_EV_FLUKA     176912138286.4409            // 1000801900   //
#define  Hg191_MASS_EV_FLUKA     177844360866.4221            // 1000801910   //
#define  Hg192_MASS_EV_FLUKA     178774467356.4581            // 1000801920   //
#define  Hg193_MASS_EV_FLUKA     179706959076.4323            // 1000801930   //
#define  Hg194_MASS_EV_FLUKA     180637277596.4627            // 1000801940   //
#define  Hg195_MASS_EV_FLUKA     181569942886.4324            // 1000801950   //
#define  Hg196_MASS_EV_FLUKA     182500669716.4522            // 1000801960   //
#define  Hg197_MASS_EV_FLUKA     183433449936.4189            // 1000801970   //
#define  Hg198_MASS_EV_FLUKA     184364531166.4296            // 1000801980   //
#define  Hg199_MASS_EV_FLUKA     185297432786.3932            // 1000801990   //
#define  Hg200_MASS_EV_FLUKA     186228970206.3921            // 1000802000   //
#define  Hg202_MASS_EV_FLUKA     188094116946.3361            // 1000802020   //
#define  Hg204_MASS_EV_FLUKA     189959760386.2673            // 1000802040   //
#define  Hg206_MASS_EV_FLUKA     191826497036.1702            // 1000802060   //

#define  Tl188_MASS_EV_FLUKA     175057632842.3265            // 1000811880   //
#define  Tl189_MASS_EV_FLUKA     175987048032.3804            // 1000811890   //
#define  Tl190_MASS_EV_FLUKA     176918642272.3778            // 1000811900   //
#define  Tl191_MASS_EV_FLUKA     177848357802.4239            // 1000811910   //
#define  Tl192_MASS_EV_FLUKA     178780091652.4177            // 1000811920   //
#define  Tl193_MASS_EV_FLUKA     179710102492.4562            // 1000811930   //
#define  Tl195_MASS_EV_FLUKA     181572249202.4779            // 1000811950   //
#define  Tl196_MASS_EV_FLUKA     182504550492.4570            // 1000811960   //
#define  Tl197_MASS_EV_FLUKA     183435137012.4806            // 1000811970   //
#define  Tl198_MASS_EV_FLUKA     184367495382.4582            // 1000811980   //
#define  Tl199_MASS_EV_FLUKA     185298381682.4739            // 1000811990   //
#define  Tl200_MASS_EV_FLUKA     186230930442.4466            // 1000812000   //
#define  Tl201_MASS_EV_FLUKA     187162292702.4500            // 1000812010   //
#define  Tl202_MASS_EV_FLUKA     188094985622.4189            // 1000812020   //
#define  Tl203_MASS_EV_FLUKA     189026702132.4132            // 1000812030   //
#define  Tl204_MASS_EV_FLUKA     189959611952.3765            // 1000812040   //
#define  Tl205_MASS_EV_FLUKA     190891631172.3629            // 1000812050   //
#define  Tl206_MASS_EV_FLUKA     191824693192.3223            // 1000812060   //
#define  Tl207_MASS_EV_FLUKA     192757410582.2906            // 1000812070   //

#define  Pb192_MASS_EV_FLUKA     178782965692.6692            // 1000821920   //
#define  Pb193_MASS_EV_FLUKA     179714758002.6615            // 1000821930   //
#define  Pb194_MASS_EV_FLUKA     180644286362.7124            // 1000821940   //
#define  Pb196_MASS_EV_FLUKA     182506100752.7429            // 1000821960   //
#define  Pb197_MASS_EV_FLUKA     183438220372.7266            // 1000821970   //
#define  Pb198_MASS_EV_FLUKA     184368409862.7605            // 1000821980   //
#define  Pb199_MASS_EV_FLUKA     185300769542.7380            // 1000821990   //
#define  Pb200_MASS_EV_FLUKA     186231245422.7644            // 1000822000   //
#define  Pb201_MASS_EV_FLUKA     187163700122.7396            // 1000822010   //
#define  Pb202_MASS_EV_FLUKA     188094539772.7565            // 1000822020   //
#define  Pb203_MASS_EV_FLUKA     189027181622.7268            // 1000822030   //
#define  Pb204_MASS_EV_FLUKA     189958352782.7352            // 1000822040   //
#define  Pb205_MASS_EV_FLUKA     190891186902.7004            // 1000822050   //
#define  Pb206_MASS_EV_FLUKA     191822664522.7009            // 1000822060   //
#define  Pb207_MASS_EV_FLUKA     192755492342.6663            // 1000822070   //
#define  Pb208_MASS_EV_FLUKA     193687690162.6481            // 1000822080   //

#define  Bi197_MASS_EV_FLUKA     183442901148.4195            // 1000831970   //

#endif /* XCOLL_FLUKA_MASSES_H */
'''


fluka_masses = {
    int(line.split()[4]): [line.split()[1].split('_')[0], float(line.split()[2])]
    for line in source.split('\n')
    if len(line.split()) > 1 and len(line.split()[1]) > 6 and line.split()[1][-6:] == '_FLUKA'
}

for pdg_id, vals in fluka_masses.items():
    q, A, Z, name = pdg.get_properties_from_pdg_id(pdg_id, subscripts=False)
    if Z > 2:
        if vals[0] != name:
            raise ValueError(f"Name mismatch in FLUKA masses database for {pdg_id}: "
                           + f"{vals[0]} != {name}")
    if not pdg._mass_consistent(pdg_id, vals[1]):
        raise ValueError(f"Mass mismatch in FLUKA masses database for {pdg_id}: "
                       + f"{vals[1]} != {pdg.get_mass_from_pdg_id(pdg_id)}")
