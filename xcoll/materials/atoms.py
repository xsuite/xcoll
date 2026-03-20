# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

from xtrack.particles import pdg

from .material import Material, IsotopeData
from .database import db, _manually_add_material_to_db


# Atomic mass from https://iupac.qmul.ac.uk/AtWt/ rounded to 5 digits
# Densities (in g/cm3) from:
#     https://en.wikipedia.org/wiki/List_of_chemical_elements
#     https://en.wikipedia.org/wiki/Densities_of_the_elements_(data_page)
Hydrogen       = Material(Z=1   , A=1.008   , state='gas'    , density=89.88e-6,
                 isotopes=(IsotopeData(mass_number=  1, atomic_mass=1.00782503223, abundance=0.999885),
                           IsotopeData(mass_number=  2, atomic_mass=2.01410177812, abundance=0.000115),
                           IsotopeData(mass_number=  3, atomic_mass=3.01604927791, abundance=None),
                 ))
Helium         = Material(Z=2   , A=4.0026  , state='gas'    , density=178.5e-6,
                 isotopes=(IsotopeData(mass_number=  3, atomic_mass=3.01602932008, abundance=1.34e-06),
                           IsotopeData(mass_number=  4, atomic_mass=4.00260325413, abundance=0.99999866),
                 ))
Lithium        = Material(Z=3   , A=6.94    , state='solid'  , density=0.534,
                 isotopes=(IsotopeData(mass_number=  6, atomic_mass=6.01512228741, abundance=0.0759),
                           IsotopeData(mass_number=  7, atomic_mass=7.0160034366, abundance=0.9241),
                 ))
Beryllium      = Material(Z=4   , A=9.0122  , state='solid'  , density=1.848,
                 isotopes=(IsotopeData(mass_number=  9, atomic_mass=9.012183065, abundance=1.0),
                 ))
Boron          = Material(Z=5   , A=10.81   , state='solid'  , density=2.34,
                 isotopes=(IsotopeData(mass_number= 10, atomic_mass=10.01293695, abundance=0.199),
                           IsotopeData(mass_number= 11, atomic_mass=11.00930536, abundance=0.801),
                 ))
Carbon         = Material(Z=6   , A=12.011  , state='solid'  , density=2.265,
                 isotopes=(IsotopeData(mass_number= 12, atomic_mass=12.0, abundance=0.9893),
                           IsotopeData(mass_number= 13, atomic_mass=13.00335483507, abundance=0.0107),
                           IsotopeData(mass_number= 14, atomic_mass=14.0032419884, abundance=None),
                 ))   # Amorphous
Nitrogen       = Material(Z=7   , A=14.007  , state='gas'    , density=0.0012506,
                 isotopes=(IsotopeData(mass_number= 14, atomic_mass=14.00307400443, abundance=0.99636),
                           IsotopeData(mass_number= 15, atomic_mass=15.00010889888, abundance=0.00364),
                 ))
Oxygen         = Material(Z=8   , A=15.999  , state='gas'    , density=0.001429,
                 isotopes=(IsotopeData(mass_number= 16, atomic_mass=15.99491461957, abundance=0.99757),
                           IsotopeData(mass_number= 17, atomic_mass=16.9991317565, abundance=0.00038),
                           IsotopeData(mass_number= 18, atomic_mass=17.99915961286, abundance=0.00205),
                 ))
Fluorine       = Material(Z=9   , A=18.998  , state='gas'    , density=0.001696,
                 isotopes=(IsotopeData(mass_number= 19, atomic_mass=18.99840316273, abundance=1.0),
                 ))
Neon           = Material(Z=10  , A=20.18   , state='gas'    , density=0.0009002,
                 isotopes=(IsotopeData(mass_number= 20, atomic_mass=19.9924401762, abundance=0.9048),
                           IsotopeData(mass_number= 21, atomic_mass=20.9938853, abundance=0.0027),
                           IsotopeData(mass_number= 22, atomic_mass=21.991385114, abundance=0.0925),
                 ))
Sodium         = Material(Z=11  , A=22.99   , state='solid'  , density=0.968,
                 isotopes=(IsotopeData(mass_number= 23, atomic_mass=22.98976928, abundance=1.0),
                 ))
Magnesium      = Material(Z=12  , A=24.305  , state='solid'  , density=1.738,
                 isotopes=(IsotopeData(mass_number= 24, atomic_mass=23.985041697, abundance=0.7899),
                           IsotopeData(mass_number= 25, atomic_mass=24.985836976, abundance=0.1),
                           IsotopeData(mass_number= 26, atomic_mass=25.982592968, abundance=0.1101),
                 ))
Aluminium      = Material(Z=13  , A=26.982  , state='solid'  , density=2.70,
                 isotopes=(IsotopeData(mass_number= 27, atomic_mass=26.98153853, abundance=1.0),
                 ))
Silicon        = Material(Z=14  , A=28.085  , state='solid'  , density=2.3296,
                 isotopes=(IsotopeData(mass_number= 28, atomic_mass=27.97692653465, abundance=0.92223),
                           IsotopeData(mass_number= 29, atomic_mass=28.9764946649, abundance=0.04685),
                           IsotopeData(mass_number= 30, atomic_mass=29.973770136, abundance=0.03092),
                 ))
Phosphorus     = Material(Z=15  , A=30.974  , state='solid'  , density=2.2,
                 isotopes=(IsotopeData(mass_number= 31, atomic_mass=30.97376199842, abundance=1.0),
                 ))   # Red phosphorus
Sulfur         = Material(Z=16  , A=32.06   , state='solid'  , density=2.067,
                 isotopes=(IsotopeData(mass_number= 32, atomic_mass=31.9720711744, abundance=0.9499),
                           IsotopeData(mass_number= 33, atomic_mass=32.9714589098, abundance=0.0075),
                           IsotopeData(mass_number= 34, atomic_mass=33.967867004, abundance=0.0425),
                           IsotopeData(mass_number= 36, atomic_mass=35.96708071, abundance=0.0001),
                 ))
Chlorine       = Material(Z=17  , A=35.45   , state='gas'    , density=0.003214,
                 isotopes=(IsotopeData(mass_number= 35, atomic_mass=34.968852682, abundance=0.7576),
                           IsotopeData(mass_number= 37, atomic_mass=36.965902602, abundance=0.2424),
                 ))
Argon          = Material(Z=18  , A=39.948  , state='gas'    , density=0.001784,
                 isotopes=(IsotopeData(mass_number= 36, atomic_mass=35.967545105, abundance=0.003336),
                           IsotopeData(mass_number= 38, atomic_mass=37.96273211, abundance=0.000629),
                           IsotopeData(mass_number= 40, atomic_mass=39.9623831237, abundance=0.996035),
                 ))
Potassium      = Material(Z=19  , A=39.098  , state='solid'  , density=0.862,
                 isotopes=(IsotopeData(mass_number= 39, atomic_mass=38.9637064864, abundance=0.932581),
                           IsotopeData(mass_number= 40, atomic_mass=39.963998166, abundance=0.000117),
                           IsotopeData(mass_number= 41, atomic_mass=40.9618252579, abundance=0.067302),
                 ))
Calcium        = Material(Z=20  , A=40.078  , state='solid'  , density=1.54,
                 isotopes=(IsotopeData(mass_number= 40, atomic_mass=39.962590863, abundance=0.96941),
                           IsotopeData(mass_number= 42, atomic_mass=41.95861783, abundance=0.00647),
                           IsotopeData(mass_number= 43, atomic_mass=42.95876644, abundance=0.00135),
                           IsotopeData(mass_number= 44, atomic_mass=43.95548156, abundance=0.02086),
                           IsotopeData(mass_number= 46, atomic_mass=45.953689, abundance=4e-05),
                           IsotopeData(mass_number= 48, atomic_mass=47.95252276, abundance=0.00187),
                 ))
Scandium       = Material(Z=21  , A=44.956  , state='solid'  , density=2.99,
                 isotopes=(IsotopeData(mass_number= 45, atomic_mass=44.95590828, abundance=1.0),
                 ))
Titanium       = Material(Z=22  , A=47.867  , state='solid'  , density=4.506,
                 isotopes=(IsotopeData(mass_number= 46, atomic_mass=45.95262772, abundance=0.0825),
                           IsotopeData(mass_number= 47, atomic_mass=46.95175879, abundance=0.0744),
                           IsotopeData(mass_number= 48, atomic_mass=47.94794198, abundance=0.7372),
                           IsotopeData(mass_number= 49, atomic_mass=48.94786568, abundance=0.0541),
                           IsotopeData(mass_number= 50, atomic_mass=49.94478689, abundance=0.0518),
                 ))
Vanadium       = Material(Z=23  , A=50.942  , state='solid'  , density=6.11,
                 isotopes=(IsotopeData(mass_number= 50, atomic_mass=49.94715601, abundance=0.0025),
                           IsotopeData(mass_number= 51, atomic_mass=50.94395704, abundance=0.9975),
                 ))
Chromium       = Material(Z=24  , A=51.996  , state='solid'  , density=7.19,
                 isotopes=(IsotopeData(mass_number= 50, atomic_mass=49.94604183, abundance=0.04345),
                           IsotopeData(mass_number= 52, atomic_mass=51.94050623, abundance=0.83789),
                           IsotopeData(mass_number= 53, atomic_mass=52.94064815, abundance=0.09501),
                           IsotopeData(mass_number= 54, atomic_mass=53.93887916, abundance=0.02365),
                 ))
Manganese      = Material(Z=25  , A=54.938  , state='solid'  , density=7.26,
                 isotopes=(IsotopeData(mass_number= 55, atomic_mass=54.93804391, abundance=1.0),
                 ))
Iron           = Material(Z=26  , A=55.845  , state='solid'  , density=7.874,
                 isotopes=(IsotopeData(mass_number= 54, atomic_mass=53.93960899, abundance=0.05845),
                           IsotopeData(mass_number= 56, atomic_mass=55.93493633, abundance=0.91754),
                           IsotopeData(mass_number= 57, atomic_mass=56.93539284, abundance=0.02119),
                           IsotopeData(mass_number= 58, atomic_mass=57.93327443, abundance=0.00282),
                 ))
Cobalt         = Material(Z=27  , A=58.933  , state='solid'  , density=8.90,
                 isotopes=(IsotopeData(mass_number= 59, atomic_mass=58.93319429, abundance=1.0),
                 ))
Nickel         = Material(Z=28  , A=58.693  , state='solid'  , density=8.908,
                 isotopes=(IsotopeData(mass_number= 58, atomic_mass=57.93534241, abundance=0.68077),
                           IsotopeData(mass_number= 60, atomic_mass=59.93078588, abundance=0.26223),
                           IsotopeData(mass_number= 61, atomic_mass=60.93105557, abundance=0.011399),
                           IsotopeData(mass_number= 62, atomic_mass=61.92834537, abundance=0.036346),
                           IsotopeData(mass_number= 64, atomic_mass=63.92796682, abundance=0.009255),
                 ))
Copper         = Material(Z=29  , A=63.546  , state='solid'  , density=8.96,
                 isotopes=(IsotopeData(mass_number= 63, atomic_mass=62.92959772, abundance=0.6915),
                           IsotopeData(mass_number= 65, atomic_mass=64.9277897, abundance=0.3085),
                 ))
Zinc           = Material(Z=30  , A=65.38   , state='solid'  , density=7.134,
                 isotopes=(IsotopeData(mass_number= 64, atomic_mass=63.92914201, abundance=0.4917),
                           IsotopeData(mass_number= 66, atomic_mass=65.92603381, abundance=0.2773),
                           IsotopeData(mass_number= 67, atomic_mass=66.92712775, abundance=0.0404),
                           IsotopeData(mass_number= 68, atomic_mass=67.92484455, abundance=0.1845),
                           IsotopeData(mass_number= 70, atomic_mass=69.9253192, abundance=0.0061),
                 ))
Gallium        = Material(Z=31  , A=69.723  , state='solid'  , density=5.91,
                 isotopes=(IsotopeData(mass_number= 69, atomic_mass=68.9255735, abundance=0.60108),
                           IsotopeData(mass_number= 71, atomic_mass=70.92470258, abundance=0.39892),
                 ))
Germanium      = Material(Z=32  , A=72.63   , state='solid'  , density=5.323,
                 isotopes=(IsotopeData(mass_number= 70, atomic_mass=69.92424875, abundance=0.2057),
                           IsotopeData(mass_number= 72, atomic_mass=71.922075826, abundance=0.2745),
                           IsotopeData(mass_number= 73, atomic_mass=72.923458956, abundance=0.0775),
                           IsotopeData(mass_number= 74, atomic_mass=73.921177761, abundance=0.365),
                           IsotopeData(mass_number= 76, atomic_mass=75.921402726, abundance=0.0773),
                 ))
Arsenic        = Material(Z=33  , A=74.922  , state='solid'  , density=5.727,
                 isotopes=(IsotopeData(mass_number= 75, atomic_mass=74.92159457, abundance=1.0),
                 ))
Selenium       = Material(Z=34  , A=78.971  , state='solid'  , density=4.809,
                 isotopes=(IsotopeData(mass_number= 74, atomic_mass=73.922475934, abundance=0.0089),
                           IsotopeData(mass_number= 76, atomic_mass=75.919213704, abundance=0.0937),
                           IsotopeData(mass_number= 77, atomic_mass=76.919914154, abundance=0.0763),
                           IsotopeData(mass_number= 78, atomic_mass=77.91730928, abundance=0.2377),
                           IsotopeData(mass_number= 80, atomic_mass=79.9165218, abundance=0.4961),
                           IsotopeData(mass_number= 82, atomic_mass=81.9166995, abundance=0.0873),
                 ))
Bromine        = Material(Z=35  , A=79.904  , state='liquid' , density=3.1028,
                 isotopes=(IsotopeData(mass_number= 79, atomic_mass=78.9183376, abundance=0.5069),
                           IsotopeData(mass_number= 81, atomic_mass=80.9162897, abundance=0.4931),
                 ))
Krypton        = Material(Z=36  , A=83.798  , state='gas'    , density=0.003749,
                 isotopes=(IsotopeData(mass_number= 78, atomic_mass=77.92036494, abundance=0.000355),
                           IsotopeData(mass_number= 80, atomic_mass=79.91637808, abundance=0.002286),
                           IsotopeData(mass_number= 82, atomic_mass=81.91348273, abundance=0.11593),
                           IsotopeData(mass_number= 83, atomic_mass=82.91412716, abundance=0.115),
                           IsotopeData(mass_number= 84, atomic_mass=83.9115107282, abundance=0.56987),
                           IsotopeData(mass_number= 86, atomic_mass=85.9106106269, abundance=0.17279),
                 ))
Rubidium       = Material(Z=37  , A=85.468  , state='solid'  , density=1.532,
                 isotopes=(IsotopeData(mass_number= 85, atomic_mass=84.9117897379, abundance=0.7217),
                           IsotopeData(mass_number= 87, atomic_mass=86.909180531, abundance=0.2783),
                 ))
Strontium      = Material(Z=38  , A=87.62   , state='solid'  , density=2.64,
                 isotopes=(IsotopeData(mass_number= 84, atomic_mass=83.9134191, abundance=0.0056),
                           IsotopeData(mass_number= 86, atomic_mass=85.9092606, abundance=0.0986),
                           IsotopeData(mass_number= 87, atomic_mass=86.9088775, abundance=0.07),
                           IsotopeData(mass_number= 88, atomic_mass=87.9056125, abundance=0.8258),
                 ))
Yttrium        = Material(Z=39  , A=88.906  , state='solid'  , density=4.469,
                 isotopes=(IsotopeData(mass_number= 89, atomic_mass=88.9058403, abundance=1.0),
                 ))
Zirconium      = Material(Z=40  , A=91.224  , state='solid'  , density=6.506,
                 isotopes=(IsotopeData(mass_number= 90, atomic_mass=89.9046977, abundance=0.5145),
                           IsotopeData(mass_number= 91, atomic_mass=90.9056396, abundance=0.1122),
                           IsotopeData(mass_number= 92, atomic_mass=91.9050347, abundance=0.1715),
                           IsotopeData(mass_number= 94, atomic_mass=93.9063108, abundance=0.1738),
                           IsotopeData(mass_number= 96, atomic_mass=95.9082714, abundance=0.028),
                 ))
Niobium        = Material(Z=41  , A=92.906  , state='solid'  , density=8.57,
                 isotopes=(IsotopeData(mass_number= 93, atomic_mass=92.906373, abundance=1.0),
                 ))
Molybdenum     = Material(Z=42  , A=95.95   , state='solid'  , density=10.223,
                 isotopes=(IsotopeData(mass_number= 92, atomic_mass=91.90680796, abundance=0.1453),
                           IsotopeData(mass_number= 94, atomic_mass=93.9050849, abundance=0.0915),
                           IsotopeData(mass_number= 95, atomic_mass=94.90583877, abundance=0.1584),
                           IsotopeData(mass_number= 96, atomic_mass=95.90467612, abundance=0.1667),
                           IsotopeData(mass_number= 97, atomic_mass=96.90601812, abundance=0.096),
                           IsotopeData(mass_number= 98, atomic_mass=97.90540482, abundance=0.2439),
                           IsotopeData(mass_number=100, atomic_mass=99.9074718, abundance=0.0982),
                 ))
Technetium     = Material(Z=43  , A=97      , state='solid'  , density=11.5,
                 isotopes=(IsotopeData(mass_number= 97, atomic_mass=96.9063667, abundance=None),
                           IsotopeData(mass_number= 98, atomic_mass=97.9072124, abundance=None),
                           IsotopeData(mass_number= 99, atomic_mass=98.9062508, abundance=None),
                 ))   # Only unstable isotopes
Ruthenium      = Material(Z=44  , A=101.07  , state='solid'  , density=12.37,
                 isotopes=(IsotopeData(mass_number= 96, atomic_mass=95.90759025, abundance=0.0554),
                           IsotopeData(mass_number= 98, atomic_mass=97.9052868, abundance=0.0187),
                           IsotopeData(mass_number= 99, atomic_mass=98.9059341, abundance=0.1276),
                           IsotopeData(mass_number=100, atomic_mass=99.9042143, abundance=0.126),
                           IsotopeData(mass_number=101, atomic_mass=100.9055769, abundance=0.1706),
                           IsotopeData(mass_number=102, atomic_mass=101.9043441, abundance=0.3155),
                           IsotopeData(mass_number=104, atomic_mass=103.9054275, abundance=0.1862),
                 ))
Rhodium        = Material(Z=45  , A=102.91  , state='solid'  , density=12.41,
                 isotopes=(IsotopeData(mass_number=103, atomic_mass=102.905498, abundance=1.0),
                 ))
Palladium      = Material(Z=46  , A=106.42  , state='solid'  , density=12.02,
                 isotopes=(IsotopeData(mass_number=102, atomic_mass=101.9056022, abundance=0.0102),
                           IsotopeData(mass_number=104, atomic_mass=103.9040305, abundance=0.1114),
                           IsotopeData(mass_number=105, atomic_mass=104.9050796, abundance=0.2233),
                           IsotopeData(mass_number=106, atomic_mass=105.9034804, abundance=0.2733),
                           IsotopeData(mass_number=108, atomic_mass=107.9038916, abundance=0.2646),
                           IsotopeData(mass_number=110, atomic_mass=109.9051722, abundance=0.1172),
                 ))
Silver         = Material(Z=47  , A=107.87  , state='solid'  , density=10.49,
                 isotopes=(IsotopeData(mass_number=107, atomic_mass=106.9050916, abundance=0.51839),
                           IsotopeData(mass_number=109, atomic_mass=108.9047553, abundance=0.48161),
                 ))
Cadmium        = Material(Z=48  , A=112.41  , state='solid'  , density=8.65,
                 isotopes=(IsotopeData(mass_number=106, atomic_mass=105.9064599, abundance=0.0125),
                           IsotopeData(mass_number=108, atomic_mass=107.9041834, abundance=0.0089),
                           IsotopeData(mass_number=110, atomic_mass=109.90300661, abundance=0.1249),
                           IsotopeData(mass_number=111, atomic_mass=110.90418287, abundance=0.128),
                           IsotopeData(mass_number=112, atomic_mass=111.90276287, abundance=0.2413),
                           IsotopeData(mass_number=113, atomic_mass=112.90440813, abundance=0.1222),
                           IsotopeData(mass_number=114, atomic_mass=113.90336509, abundance=0.2873),
                           IsotopeData(mass_number=116, atomic_mass=115.90476315, abundance=0.0749),
                 ))
Indium         = Material(Z=49  , A=114.82  , state='solid'  , density=7.31,
                 isotopes=(IsotopeData(mass_number=113, atomic_mass=112.90406184, abundance=0.0429),
                           IsotopeData(mass_number=115, atomic_mass=114.903878776, abundance=0.9571),
                 ))
Tin            = Material(Z=50  , A=118.71  , state='solid'  , density=7.365,
                 isotopes=(IsotopeData(mass_number=112, atomic_mass=111.90482387, abundance=0.0097),
                           IsotopeData(mass_number=114, atomic_mass=113.9027827, abundance=0.0066),
                           IsotopeData(mass_number=115, atomic_mass=114.903344699, abundance=0.0034),
                           IsotopeData(mass_number=116, atomic_mass=115.9017428, abundance=0.1454),
                           IsotopeData(mass_number=117, atomic_mass=116.90295398, abundance=0.0768),
                           IsotopeData(mass_number=118, atomic_mass=117.90160657, abundance=0.2422),
                           IsotopeData(mass_number=119, atomic_mass=118.90311117, abundance=0.0859),
                           IsotopeData(mass_number=120, atomic_mass=119.90220163, abundance=0.3258),
                           IsotopeData(mass_number=122, atomic_mass=121.9034438, abundance=0.0463),
                           IsotopeData(mass_number=124, atomic_mass=123.9052766, abundance=0.0579),
                 ))
Antimony       = Material(Z=51  , A=121.76  , state='solid'  , density=6.697,
                 isotopes=(IsotopeData(mass_number=121, atomic_mass=120.903812, abundance=0.5721),
                           IsotopeData(mass_number=123, atomic_mass=122.9042132, abundance=0.4279),
                 ))
Tellurium      = Material(Z=52  , A=127.6   , state='solid'  , density=6.24,
                 isotopes=(IsotopeData(mass_number=120, atomic_mass=119.9040593, abundance=0.0009),
                           IsotopeData(mass_number=122, atomic_mass=121.9030435, abundance=0.0255),
                           IsotopeData(mass_number=123, atomic_mass=122.9042698, abundance=0.0089),
                           IsotopeData(mass_number=124, atomic_mass=123.9028171, abundance=0.0474),
                           IsotopeData(mass_number=125, atomic_mass=124.9044299, abundance=0.0707),
                           IsotopeData(mass_number=126, atomic_mass=125.9033109, abundance=0.1884),
                           IsotopeData(mass_number=128, atomic_mass=127.90446128, abundance=0.3174),
                           IsotopeData(mass_number=130, atomic_mass=129.906222748, abundance=0.3408),
                 ))
Iodine         = Material(Z=53  , A=126.9   , state='solid'  , density=4.933,
                 isotopes=(IsotopeData(mass_number=127, atomic_mass=126.9044719, abundance=1.0),
                 ))
Xenon          = Material(Z=54  , A=131.29  , state='gas'    , density=0.005894,
                 isotopes=(IsotopeData(mass_number=124, atomic_mass=123.905892, abundance=0.000952),
                           IsotopeData(mass_number=126, atomic_mass=125.9042983, abundance=0.00089),
                           IsotopeData(mass_number=128, atomic_mass=127.903531, abundance=0.019102),
                           IsotopeData(mass_number=129, atomic_mass=128.9047808611, abundance=0.264006),
                           IsotopeData(mass_number=130, atomic_mass=129.903509349, abundance=0.04071),
                           IsotopeData(mass_number=131, atomic_mass=130.90508406, abundance=0.212324),
                           IsotopeData(mass_number=132, atomic_mass=131.9041550856, abundance=0.269086),
                           IsotopeData(mass_number=134, atomic_mass=133.90539466, abundance=0.104357),
                           IsotopeData(mass_number=136, atomic_mass=135.907214484, abundance=0.088573),
                 ))
Cesium         = Material(Z=55  , A=132.91  , state='solid'  , density=1.93,
                 isotopes=(IsotopeData(mass_number=133, atomic_mass=132.905451961, abundance=1.0),
                 ))
Barium         = Material(Z=56  , A=137.33  , state='solid'  , density=3.62,
                 isotopes=(IsotopeData(mass_number=130, atomic_mass=129.9063207, abundance=0.00106),
                           IsotopeData(mass_number=132, atomic_mass=131.9050611, abundance=0.00101),
                           IsotopeData(mass_number=134, atomic_mass=133.90450818, abundance=0.02417),
                           IsotopeData(mass_number=135, atomic_mass=134.90568838, abundance=0.06592),
                           IsotopeData(mass_number=136, atomic_mass=135.90457573, abundance=0.07854),
                           IsotopeData(mass_number=137, atomic_mass=136.90582714, abundance=0.11232),
                           IsotopeData(mass_number=138, atomic_mass=137.905247, abundance=0.71698),
                 ))
Lanthanum      = Material(Z=57  , A=138.91  , state='solid'  , density=6.15,
                 isotopes=(IsotopeData(mass_number=138, atomic_mass=137.9071149, abundance=0.0008881),
                           IsotopeData(mass_number=139, atomic_mass=138.9063563, abundance=0.9991119),
                 ))
Cerium         = Material(Z=58  , A=140.12  , state='solid'  , density=6.77,
                 isotopes=(IsotopeData(mass_number=136, atomic_mass=135.90712921, abundance=0.00185),
                           IsotopeData(mass_number=138, atomic_mass=137.905991, abundance=0.00251),
                           IsotopeData(mass_number=140, atomic_mass=139.9054431, abundance=0.8845),
                           IsotopeData(mass_number=142, atomic_mass=141.9092504, abundance=0.11114),
                 ))
Praseodymium   = Material(Z=59  , A=140.91  , state='solid'  , density=6.77,
                 isotopes=(IsotopeData(mass_number=141, atomic_mass=140.9076576, abundance=1.0),
                 ))
Neodymium      = Material(Z=60  , A=144.24  , state='solid'  , density=7.01,
                 isotopes=(IsotopeData(mass_number=142, atomic_mass=141.907729, abundance=0.27152),
                           IsotopeData(mass_number=143, atomic_mass=142.90982, abundance=0.12174),
                           IsotopeData(mass_number=144, atomic_mass=143.910093, abundance=0.23798),
                           IsotopeData(mass_number=145, atomic_mass=144.9125793, abundance=0.08293),
                           IsotopeData(mass_number=146, atomic_mass=145.9131226, abundance=0.17189),
                           IsotopeData(mass_number=148, atomic_mass=147.9168993, abundance=0.05756),
                           IsotopeData(mass_number=150, atomic_mass=149.9209022, abundance=0.05638),
                 ))
Promethium     = Material(Z=61  , A=145     , state='solid'  , density=7.26,
                 isotopes=(IsotopeData(mass_number=145, atomic_mass=144.9127559, abundance=None),
                           IsotopeData(mass_number=147, atomic_mass=146.915145, abundance=None),
                 ))   # Only unstable isotopes
Samarium       = Material(Z=62  , A=150.36  , state='solid'  , density=7.52,
                 isotopes=(IsotopeData(mass_number=144, atomic_mass=143.9120065, abundance=0.0307),
                           IsotopeData(mass_number=147, atomic_mass=146.9149044, abundance=0.1499),
                           IsotopeData(mass_number=148, atomic_mass=147.9148292, abundance=0.1124),
                           IsotopeData(mass_number=149, atomic_mass=148.9171921, abundance=0.1382),
                           IsotopeData(mass_number=150, atomic_mass=149.9172829, abundance=0.0738),
                           IsotopeData(mass_number=152, atomic_mass=151.9197397, abundance=0.2675),
                           IsotopeData(mass_number=154, atomic_mass=153.9222169, abundance=0.2275),
                 ))
Europium       = Material(Z=63  , A=151.96  , state='solid'  , density=5.24,
                 isotopes=(IsotopeData(mass_number=151, atomic_mass=150.9198578, abundance=0.4781),
                           IsotopeData(mass_number=153, atomic_mass=152.921238, abundance=0.5219),
                 ))
Gadolinium     = Material(Z=64  , A=157.25  , state='solid'  , density=7.90,
                 isotopes=(IsotopeData(mass_number=152, atomic_mass=151.9197995, abundance=0.002),
                           IsotopeData(mass_number=154, atomic_mass=153.9208741, abundance=0.0218),
                           IsotopeData(mass_number=155, atomic_mass=154.9226305, abundance=0.148),
                           IsotopeData(mass_number=156, atomic_mass=155.9221312, abundance=0.2047),
                           IsotopeData(mass_number=157, atomic_mass=156.9239686, abundance=0.1565),
                           IsotopeData(mass_number=158, atomic_mass=157.9241123, abundance=0.2484),
                           IsotopeData(mass_number=160, atomic_mass=159.9270624, abundance=0.2186),
                 ))
Terbium        = Material(Z=65  , A=158.93  , state='solid'  , density=8.23,
                 isotopes=(IsotopeData(mass_number=159, atomic_mass=158.9253547, abundance=1.0),
                 ))
Dysprosium     = Material(Z=66  , A=162.5   , state='solid'  , density=8.55,
                 isotopes=(IsotopeData(mass_number=156, atomic_mass=155.9242847, abundance=5.6e-05),
                           IsotopeData(mass_number=158, atomic_mass=157.9244159, abundance=9.5e-05),
                           IsotopeData(mass_number=160, atomic_mass=159.9252046, abundance=0.02329),
                           IsotopeData(mass_number=161, atomic_mass=160.9269405, abundance=0.18889),
                           IsotopeData(mass_number=162, atomic_mass=161.9268056, abundance=0.25475),
                           IsotopeData(mass_number=163, atomic_mass=162.9287383, abundance=0.24896),
                           IsotopeData(mass_number=164, atomic_mass=163.9291819, abundance=0.2826),
                 ))
Holmium        = Material(Z=67  , A=164.93  , state='solid'  , density=8.80,
                 isotopes=(IsotopeData(mass_number=165, atomic_mass=164.9303288, abundance=1.0),
                 ))
Erbium         = Material(Z=68  , A=167.26  , state='solid'  , density=9.07,
                 isotopes=(IsotopeData(mass_number=162, atomic_mass=161.9287884, abundance=0.00139),
                           IsotopeData(mass_number=164, atomic_mass=163.9292088, abundance=0.01601),
                           IsotopeData(mass_number=166, atomic_mass=165.9302995, abundance=0.33503),
                           IsotopeData(mass_number=167, atomic_mass=166.9320546, abundance=0.22869),
                           IsotopeData(mass_number=168, atomic_mass=167.9323767, abundance=0.26978),
                           IsotopeData(mass_number=170, atomic_mass=169.9354702, abundance=0.1491),
                 ))
Thulium        = Material(Z=69  , A=168.93  , state='solid'  , density=9.32,
                 isotopes=(IsotopeData(mass_number=169, atomic_mass=168.9342179, abundance=1.0),
                 ))
Ytterbium      = Material(Z=70  , A=173.05  , state='solid'  , density=6.97,
                 isotopes=(IsotopeData(mass_number=168, atomic_mass=167.9338896, abundance=0.00123),
                           IsotopeData(mass_number=170, atomic_mass=169.9347664, abundance=0.02982),
                           IsotopeData(mass_number=171, atomic_mass=170.9363302, abundance=0.1409),
                           IsotopeData(mass_number=172, atomic_mass=171.9363859, abundance=0.2168),
                           IsotopeData(mass_number=173, atomic_mass=172.9382151, abundance=0.16103),
                           IsotopeData(mass_number=174, atomic_mass=173.9388664, abundance=0.32026),
                           IsotopeData(mass_number=176, atomic_mass=175.9425764, abundance=0.12996),
                 ))
Lutetium       = Material(Z=71  , A=174.97  , state='solid'  , density=9.84,
                 isotopes=(IsotopeData(mass_number=175, atomic_mass=174.9407752, abundance=0.97401),
                           IsotopeData(mass_number=176, atomic_mass=175.9426897, abundance=0.02599),
                 ))
Hafnium        = Material(Z=72  , A=178.49  , state='solid'  , density=13.31,
                 isotopes=(IsotopeData(mass_number=174, atomic_mass=173.9400461, abundance=0.0016),
                           IsotopeData(mass_number=176, atomic_mass=175.9414076, abundance=0.0526),
                           IsotopeData(mass_number=177, atomic_mass=176.9432277, abundance=0.186),
                           IsotopeData(mass_number=178, atomic_mass=177.9437058, abundance=0.2728),
                           IsotopeData(mass_number=179, atomic_mass=178.9458232, abundance=0.1362),
                           IsotopeData(mass_number=180, atomic_mass=179.946557, abundance=0.3508),
                 ))
Tantalum       = Material(Z=73  , A=180.95  , state='solid'  , density=16.65,
                 isotopes=(IsotopeData(mass_number=180, atomic_mass=179.9474648, abundance=0.0001201),
                           IsotopeData(mass_number=181, atomic_mass=180.9479958, abundance=0.9998799),
                 ))
Tungsten       = Material(Z=74  , A=183.84  , state='solid'  , density=19.3,
                 isotopes=(IsotopeData(mass_number=180, atomic_mass=179.9467108, abundance=0.0012),
                           IsotopeData(mass_number=182, atomic_mass=181.94820394, abundance=0.265),
                           IsotopeData(mass_number=183, atomic_mass=182.95022275, abundance=0.1431),
                           IsotopeData(mass_number=184, atomic_mass=183.95093092, abundance=0.3064),
                           IsotopeData(mass_number=186, atomic_mass=185.9543628, abundance=0.2843),
                 ))
Rhenium        = Material(Z=75  , A=186.21  , state='solid'  , density=21.02,
                 isotopes=(IsotopeData(mass_number=185, atomic_mass=184.9529545, abundance=0.374),
                           IsotopeData(mass_number=187, atomic_mass=186.9557501, abundance=0.626),
                 ))
Osmium         = Material(Z=76  , A=190.23  , state='solid'  , density=22.59,
                 isotopes=(IsotopeData(mass_number=184, atomic_mass=183.9524885, abundance=0.0002),
                           IsotopeData(mass_number=186, atomic_mass=185.953835, abundance=0.0159),
                           IsotopeData(mass_number=187, atomic_mass=186.9557474, abundance=0.0196),
                           IsotopeData(mass_number=188, atomic_mass=187.9558352, abundance=0.1324),
                           IsotopeData(mass_number=189, atomic_mass=188.9581442, abundance=0.1615),
                           IsotopeData(mass_number=190, atomic_mass=189.9584437, abundance=0.2626),
                           IsotopeData(mass_number=192, atomic_mass=191.961477, abundance=0.4078),
                 ))
Iridium        = Material(Z=77  , A=192.22  , state='solid'  , density=22.56,
                 isotopes=(IsotopeData(mass_number=191, atomic_mass=190.9605893, abundance=0.373),
                           IsotopeData(mass_number=193, atomic_mass=192.9629216, abundance=0.627),
                 ))
Platinum       = Material(Z=78  , A=195.08  , state='solid'  , density=21.45,
                 isotopes=(IsotopeData(mass_number=190, atomic_mass=189.9599297, abundance=1.2e-05),
                           IsotopeData(mass_number=192, atomic_mass=191.9610387, abundance=0.000782),
                           IsotopeData(mass_number=194, atomic_mass=193.9626809, abundance=0.3286),
                           IsotopeData(mass_number=195, atomic_mass=194.9647917, abundance=0.3378),
                           IsotopeData(mass_number=196, atomic_mass=195.96495209, abundance=0.2521),
                           IsotopeData(mass_number=198, atomic_mass=197.9678949, abundance=0.07356),
                 ))
Gold           = Material(Z=79  , A=196.97  , state='solid'  , density=19.32,
                 isotopes=(IsotopeData(mass_number=197, atomic_mass=196.96656879, abundance=1.0),
                 ))
Mercury        = Material(Z=80  , A=200.59  , state='liquid' , density=13.534,
                 isotopes=(IsotopeData(mass_number=196, atomic_mass=195.9658326, abundance=0.0015),
                           IsotopeData(mass_number=198, atomic_mass=197.9667686, abundance=0.0997),
                           IsotopeData(mass_number=199, atomic_mass=198.96828064, abundance=0.1687),
                           IsotopeData(mass_number=200, atomic_mass=199.96832659, abundance=0.231),
                           IsotopeData(mass_number=201, atomic_mass=200.97030284, abundance=0.1318),
                           IsotopeData(mass_number=202, atomic_mass=201.9706434, abundance=0.2986),
                           IsotopeData(mass_number=204, atomic_mass=203.97349398, abundance=0.0687),
                 ))
Thallium       = Material(Z=81  , A=204.38  , state='solid'  , density=11.85,
                 isotopes=(IsotopeData(mass_number=203, atomic_mass=202.9723446, abundance=0.2952),
                           IsotopeData(mass_number=205, atomic_mass=204.9744278, abundance=0.7048),
                 ))
Lead           = Material(Z=82  , A=207.19  , state='solid'  , density=11.348,
                 isotopes=(IsotopeData(mass_number=204, atomic_mass=203.973044, abundance=0.014),
                           IsotopeData(mass_number=206, atomic_mass=205.9744657, abundance=0.241),
                           IsotopeData(mass_number=207, atomic_mass=206.9758973, abundance=0.221),
                           IsotopeData(mass_number=208, atomic_mass=207.9766525, abundance=0.524),
                 ))
Bismuth        = Material(Z=83  , A=208.98  , state='solid'  , density=9.78,
                 isotopes=(IsotopeData(mass_number=209, atomic_mass=208.9803991, abundance=1.0),
                 ))
Polonium       = Material(Z=84  , A=209     , state='solid'  , density=9.32,
                 isotopes=(IsotopeData(mass_number=209, atomic_mass=208.9824308, abundance=None),
                           IsotopeData(mass_number=210, atomic_mass=209.9828741, abundance=None),
                 ))   # Only unstable isotopes
Astatine       = Material(Z=85  , A=210     , state='solid'  , density=7.0,
                 isotopes=(IsotopeData(mass_number=210, atomic_mass=209.9871479, abundance=None),
                           IsotopeData(mass_number=211, atomic_mass=210.9874966, abundance=None),
                 ))   # Only unstable isotopes
Radon          = Material(Z=86  , A=222     , state='gas'    , density=0.00973,
                 isotopes=(IsotopeData(mass_number=211, atomic_mass=210.9906011, abundance=None),
                           IsotopeData(mass_number=220, atomic_mass=220.0113941, abundance=None),
                           IsotopeData(mass_number=222, atomic_mass=222.0175782, abundance=None),
                 ))   # Only unstable isotopes
Francium       = Material(Z=87  , A=223     , state='solid'  , density=2.48,
                 isotopes=(IsotopeData(mass_number=223, atomic_mass=223.019736, abundance=None),
                 ))   # Only unstable isotopes
Radium         = Material(Z=88  , A=226     , state='solid'  , density=5.50,
                 isotopes=(IsotopeData(mass_number=223, atomic_mass=223.0185023, abundance=None),
                           IsotopeData(mass_number=224, atomic_mass=224.020212, abundance=None),
                           IsotopeData(mass_number=226, atomic_mass=226.0254103, abundance=None),
                           IsotopeData(mass_number=228, atomic_mass=228.0310707, abundance=None),
                 ))   # Only unstable isotopes
Actinium       = Material(Z=89  , A=227     , state='solid'  , density=10.07,
                 isotopes=(IsotopeData(mass_number=227, atomic_mass=227.0277523, abundance=None),
                 ))   # Only unstable isotopes
Thorium        = Material(Z=90  , A=232.04  , state='solid'  , density=11.72,
                 isotopes=(IsotopeData(mass_number=230, atomic_mass=230.0331341, abundance=None),
                           IsotopeData(mass_number=232, atomic_mass=232.0380558, abundance=1.0),
                 ))
Protactinium   = Material(Z=91  , A=231.04  , state='solid'  , density=15.37,
                 isotopes=(IsotopeData(mass_number=231, atomic_mass=231.0358842, abundance=1.0),
                 ))
Uranium        = Material(Z=92  , A=238.03  , state='solid'  , density=18.95,
                 isotopes=(IsotopeData(mass_number=233, atomic_mass=233.0396355, abundance=None),
                           IsotopeData(mass_number=234, atomic_mass=234.0409523, abundance=5.4e-05),
                           IsotopeData(mass_number=235, atomic_mass=235.0439301, abundance=0.007204),
                           IsotopeData(mass_number=236, atomic_mass=236.0455682, abundance=None),
                           IsotopeData(mass_number=238, atomic_mass=238.0507884, abundance=0.992742),
                 ))
Neptunium      = Material(Z=93  , A=237     , state='solid'  , density=20.45,
                 isotopes=(IsotopeData(mass_number=236, atomic_mass=236.04657, abundance=None),
                           IsotopeData(mass_number=237, atomic_mass=237.0481736, abundance=None),
                 ))   # Only unstable isotopes
Plutonium      = Material(Z=94  , A=244     , state='solid'  , density=19.84,
                 isotopes=(IsotopeData(mass_number=238, atomic_mass=238.0495601, abundance=None),
                           IsotopeData(mass_number=239, atomic_mass=239.0521636, abundance=None),
                           IsotopeData(mass_number=240, atomic_mass=240.0538138, abundance=None),
                           IsotopeData(mass_number=241, atomic_mass=241.0568517, abundance=None),
                           IsotopeData(mass_number=242, atomic_mass=242.0587428, abundance=None),
                           IsotopeData(mass_number=244, atomic_mass=244.0642053, abundance=None),
                 ))   # Only unstable isotopes
Americium      = Material(Z=95  , A=243     , state='solid'  , density=13.69)   # Only unstable isotopes
Curium         = Material(Z=96  , A=247     , state='solid'  , density=13.51)   # Only unstable isotopes
Berkelium      = Material(Z=97  , A=247     , state='solid'  , density=14.78)   # Only unstable isotopes
Californium    = Material(Z=98  , A=251     , state='solid'  , density=15.1)   # Only unstable isotopes
Einsteinium    = Material(Z=99  , A=252     , state='solid'  , density=8.84)   # Only unstable isotopes
Fermium        = Material(Z=100 , A=257     , state='solid'  , density=9.7)   # Only unstable isotopes
Mendelevium    = Material(Z=101 , A=258     , state='solid'  , density=10.3)   # Only unstable isotopes
Nobelium       = Material(Z=102 , A=259     , state='solid'  , density=9.9)   # Only unstable isotopes
Lawrencium     = Material(Z=103 , A=266     , state='solid'  , density=15.6)   # Only unstable isotopes
Rutherfordium  = Material(Z=104 , A=267     , state='solid'  , density=23.2)   # Only unstable isotopes
Dubnium        = Material(Z=105 , A=268     , state='solid'  , density=29.3)   # Only unstable isotopes
Seaborgium     = Material(Z=106 , A=269     , state='solid'  , density=35.0)   # Only unstable isotopes
Bohrium        = Material(Z=107 , A=270     , state='solid'  , density=37.1)   # Only unstable isotopes
Hassium        = Material(Z=108 , A=277     , state='solid'  , density=41.0)   # Only unstable isotopes
Meitnerium     = Material(Z=109 , A=278     , state='solid'  , density=37.4)   # Only unstable isotopes
Darmstadtium   = Material(Z=110 , A=281     , state='solid'  , density=34.8)   # Only unstable isotopes
Roentgenium    = Material(Z=111 , A=282     , state='solid'  , density=28.7)   # Only unstable isotopes
Copernicium    = Material(Z=112 , A=285     , state='solid'  , density=14.0)   # Only unstable isotopes
Nihonium       = Material(Z=113 , A=286     , state='solid'  , density=16.0)   # Only unstable isotopes
Flerovium      = Material(Z=114 , A=289     , state='solid'  , density=14.0)   # Only unstable isotopes
Moscovium      = Material(Z=115 , A=290     , state='solid'  , density=13.5)   # Only unstable isotopes
Livermorium    = Material(Z=116 , A=293     , state='solid'  , density=12.9)   # Only unstable isotopes
Tennessine     = Material(Z=117 , A=294     , state='solid'  , density=7.2)   # Only unstable isotopes
Oganesson      = Material(Z=118 , A=294     , state='solid'  , density=5.0)   # Only unstable isotopes


# Extra parameters for Everest
# ============================

Beryllium.adapt(inplace=True,  nuclear_radius=0.22, nuclear_elastic_slope=74.7,
                               cross_section=[0.271, 0.192, 0, 0, 0, 0.0035e-2], hcut=0.02)
Aluminium.adapt(inplace=True,  nuclear_radius=0.302, nuclear_elastic_slope=120.3,
                               cross_section=[0.643, 0.418, 0, 0, 0, 0.0340e-2], hcut=0.02)
Silicon.adapt(inplace=True,    nuclear_radius=0.441, nuclear_elastic_slope=120.14,
                               cross_section=[0.664, 0.430, 0, 0, 0, 0.0390e-2], hcut=0.02,
                               crystal_plane_distance=0.96e-7, crystal_potential=21.34,
                               nuclear_collision_length=0.3016, eta=0.9)
Copper.adapt(inplace=True,     nuclear_radius=0.366, nuclear_elastic_slope=217.8,
                               cross_section=[1.253, 0.769, 0, 0, 0, 0.1530e-2], hcut=0.01)
Germanium.adapt(inplace=True,  nuclear_radius=0.605, nuclear_elastic_slope=226.35,
                               cross_section=[1.388, 0.844, 0, 0, 0, 0.1860e-2], hcut=0.02,
                               crystal_plane_distance=1.0e-7, crystal_potential=40.0,
                               nuclear_collision_length=0.1632, eta=0.9)
Molybdenum.adapt(inplace=True, nuclear_radius=0.481, nuclear_elastic_slope=273.9,
                               cross_section=[1.713, 1.023, 0, 0, 0, 0.2650e-2], hcut=0.02)
Tungsten.adapt(inplace=True,   nuclear_radius=0.520, nuclear_elastic_slope=440.3,
                               cross_section=[2.765, 1.591, 0, 0, 0, 0.7680e-2], hcut=0.01,
                               crystal_plane_distance=0.56e-7, crystal_potential=21.0,
                               nuclear_collision_length=1.e-12, eta=0.9)
Lead.adapt(inplace=True,       nuclear_radius=0.542, nuclear_elastic_slope=455.3,
                               cross_section=[3.016, 1.724, 0, 0, 0, 0.9070e-2], hcut=0.01)


# Metadata for database
# =====================

for name, obj in list(globals().items()):  # Have to wrap in list to take a snapshot (avoid updating in-place)
    if isinstance(obj, Material) and obj.name is None:
        # Check consistency between name and Z
        if name == 'Aluminium':
            assert obj.Z == 13
        else:
            assert pdg.get_element_full_name_from_Z(obj.Z).lower() == name.lower(), \
            f'Inconsistency between material name {name} and Z={obj.Z} ({pdg.get_element_full_name_from_Z(obj.Z)})'

        short_name  = pdg.get_element_name_from_Z(obj.Z)
        fluka_name  = None
        geant4_name = short_name
        # The following elements are pre-defined in FLUKA
        if obj.Z in [1, 2, 4, 5, 6, 7, 8, 11, 12, 14, 15, 16, 18, 20, 22, 23, 24, 25, 26, 27, 28, 29, 30, 40,
                     41, 42, 47, 50, 72, 73, 74, 79, 80, 82]:
            fluka_name = name[:8].upper()
        if name == 'Aluminium': fluka_name = 'ALUMINUM'
        if name == 'Phosphorus': fluka_name = 'PHOSPHO'
        _manually_add_material_to_db(obj, name, short_name=short_name, fluka_name=fluka_name, geant4_name=geant4_name)
        db.geant4._names[f'G4_{obj._geant4_name}'] = db._strip(name)  # Additional Geant4 aliases

        obj.info = f'Elemental {obj.name.lower()}'
        if obj.state == 'gas':
            obj.temperature = 273.15
            obj.pressure = 1
        elif obj.state == 'liquid' or obj.state == 'solid':
            obj.temperature = 293.15

db._set_alias('Aluminum', 'Aluminium') # Allow American spelling


# Clean up namespace
del name, obj
del pdg
del Material
del db
