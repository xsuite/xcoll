#
#                                                                      #
# ==================================================================== #
# IMPORT built-in modules                                              #
# ==================================================================== #
#                                                                      #
import os
import os.path
import sys
from copy import deepcopy
#                                                                      #
# ==================================================================== #
# IMPORT modules made on purpose                                       #
# ==================================================================== #
#                                                                      #
#
# get path to generic_frame.fluka
from find_paths import GENERIC_FRAME
print(GENERIC_FRAME)
#
loc_path = os.path.realpath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0,loc_path)
from structure import present_FEDB
#
from FEDB_ELEMENTs import parse_prototypes_lbp, synch_against_BASIC_ELEMENTs, \
                          checkFlukaNameUniqueness, REPLICA
#
from COMPLEX_STRING import error_message, secure_open_file_4_w, red, grn, \
                           secure_open_file_4_r, secure_open_file_4_a
#
from FLUKA_lib import INCLUDE_STATEMENT_FMT
#
from LB_sty import HEADERs

#                                                                      #
# ==================================================================== #
#                               FUNCTIONs                              #
# ==================================================================== #
#                                                                      #

def start():

    print("")
    print(" " + sys.argv[0].split( "/" )[-1] + " <fedb_series> <fedb_tag> <ass_name>")
    print(" by A. Mereghetti")
    print(" ")
    print(" for creating the .inp file to test an assembly")
    print(" NB: if <ass_name> is not given, then <fedb_tag> will be used")
    print("")
    if ( len( sys.argv ) < 3 ):
        msg = "Please specify:\n"
        msg = msg + "a  <fedb_series> (mandatory); \n"
        msg = msg + "a  <fedb_tag>    (mandatory); \n"
        msg = msg + "an <ass_name>    (optional)."
        error_message( msg, True )

    # get input information:
    fedb_series = sys.argv[1]
    fedb_tag    = sys.argv[2]
    if ( len( sys.argv ) == 3 ):
        ass_name = fedb_tag[:]
    else:
        ass_name = sys.argv[3]

    ass_file_name = present_FEDB.inflect( "name", fedb_series, fedb_tag )
    ass_full_path = present_FEDB.inflect( "smbl", fedb_series, fedb_tag )

    # output files:
    out_file_name  = present_FEDB.inflect( "name", fedb_series, ass_name ) + ".inp"
    log_file_name  = present_FEDB.inflect( "name", fedb_series, ass_name ) + ".log"

    # initialise log file:
    logfile = secure_open_file_4_w( log_file_name )
    logfile.write( " called from:\n %s \n" % ( sys.argv[0] ) )
    logfile.close()

    return out_file_name, log_file_name, ass_full_path, ass_name

#                                                                      #
# ==================================================================== #
#   MAIN MAIN MAIN MAIN MAIN MAIN MAIN MAIN MAIN MAIN MAIN MAIN MAIN   #
# ==================================================================== #
#                                                                      #

# some initialisation
out_file_name, log_file_name, ass_full_path, ass_name = start()
position = [ 0.0, 0.0, 0.0 ]
eulers   = [ 0.0, 0.0, 0.0 ] # [deg]
ADDITIONAL_TRANSF = []
#
l_Bl = 2 # Max number of beams
# A.Mereghetti, 09/08/2012
# with devel, this dictonary should be replaced by a more complex
#      dictionary in customize_geo_builder.py
flukaIndices = { "ilatt": 0, "irot": 0 }
#
# exceptions in the mapping of prototypes:
# these are arrays of Dictionaries, with a Dictionary per beam;
# each dictionary stores the name of an element as key,
#   and the name of the corresponding ELEMENT_MODEL as value;
# eg:
#   mappingExceptions["map_entries"][0]={
#                "MCBY.A4R2":"MCBYrob",
#                "MCBY.B4R2":"MCBYang",
#                "MCBY.C4R2":"MCBYtri",
#                }
mappingExceptions={ \
    "map_entries" : [ {} for col in range( l_Bl ) ] ,\
    "map_vetos"   : [ {} for col in range( l_Bl ) ]  \
        }

# log FEDB
logfile = secure_open_file_4_a( log_file_name )
logfile.write( present_FEDB.log_structure() )
logfile.close()

# parse the assembly file, and return the BASIC_ELEMENTs and the
#     ELEMENT_MODELs
BASIC_ELEMENTs, ELEMENT_MODELs, mappingExceptions = \
    parse_prototypes_lbp( ass_full_path, present_FEDB, \
                              mappingExceptions, log_file_name )
if ( ass_name not in ELEMENT_MODELs ):
    msg = "No '%s' ASSEMBLY declared in file: \n" % ( red( ass_name ) )
    msg = msg + "'%s'" % ( red( ass_full_path ) )
    error_message( msg, True )

# log BASIC_ELEMENTs:
logfile = secure_open_file_4_a( log_file_name )
for element in BASIC_ELEMENTs.values():
    logfile.write( element.log() )
logfile.close()

# synchronise the ELEMENT_MODELs against the BASIC_ELEMENTs
ELEMENT_MODELs = synch_against_BASIC_ELEMENTs( ELEMENT_MODELs, BASIC_ELEMENTs, \
                                                   present_FEDB, log_file_name )

# build internal rotdefis, in order to set up the assembly
for ele_name, model in ELEMENT_MODELs.items():
    model.build_rotdefis()
    continue

# check uniqueness of FLUKA names:
checkFlukaNameUniqueness( ELEMENT_MODELs, log_file_name )

# log ELEMENT_MODELs:
logfile = secure_open_file_4_a( log_file_name )
for element in ELEMENT_MODELs.values():
    logfile.write( element.log() )
logfile.close()

# take the only ELEMENT_MODEL present in the dictionary:
current_replica = REPLICA()
current_replica.inherit_from( ELEMENT_MODELs[ ass_name ] )
current_replica.names["full"] = ass_name
current_replica.names["prto"] = ass_name
current_replica.names["FLKA"] = ass_name
current_replica.build_rotdefis( position, eulers, "", ADDITIONAL_TRANSF )
# A.Mereghetti, 28/06/2012
# rename only the transformations, in order to avoid
#     names longer than 10 characters (geoedit complains);
# bodies and regions will keep the names found in the .lbp assembly file,
#     so that geometry debugging becomes easier;
# A.Mereghetti, 09/08/2012
# starting from devel, these two calls should be changed into:
#    current_replica.update_names_and_indeces( key="TRANSFs" )
current_replica.change_names( key="TRANSFs" )
current_replica.change_indices( key="TRANSFs", flukaIndices=flukaIndices )

BODIEs  = current_replica.ret( "BODIEs"  )
REGIONs = current_replica.ret( "REGIONs" )
TRANSFs = current_replica.ret( "TRANSFs" )

# vacuum around assembly:
EXTERNAL_REG = deepcopy( current_replica.ret( "EXTERNAL" ) )
EXTERNAL_REG.name     = "VAROUND"
EXTERNAL_REG.material = "VACUUM"
EXTERNAL_REG.comment = "* vacuum around assembly"
zones = EXTERNAL_REG.defin.split( "|" )
EXTERNAL_REG.defin = ""
if ( len( zones ) > 1 ):
    for zone in zones:
        EXTERNAL_REG.defin = EXTERNAL_REG.defin + " | +%-8s %s" % ( "cont", zone.strip() )
        continue
else:
    EXTERNAL_REG.defin = " +%-8s %s" % ( "cont", zones[0].strip() )
REGIONs.append( EXTERNAL_REG )

# write everything
print("")
print(" dumping the geometry in the '%s' file..." % \
    ( grn( out_file_name ) ))
ofile = secure_open_file_4_w( out_file_name )
ifile = secure_open_file_4_r( GENERIC_FRAME )
flag = 0
fastforward_def = 'the_big_bang_theory_is_a_silly_sit_com'
fastforward = fastforward_def
#
for line in ifile.readlines():

    #
    # understand what to do:
    # NB: multi-line actions
    #

    if ( line.find("$START:build_line:BODIEs$") > -1 ):
        ofile.write( line )
        flag = 1
    elif ( line.find("$START:build_line:REGIONs$") > -1 ):
        ofile.write( line )
        flag = 2
    elif ( line.find("$START:build_line:LATTICEs$") > -1 ):
        ofile.write( line )
        flag = 3
    elif ( line.find("$START:build_line:ROT-DEFIs$") > -1 ):
        ofile.write( line )
        flag = 4
    elif ( line.find("$START:build_line:ASSIGNMAs$") > -1 ):
        ofile.write( line )
        flag = 6
    elif ( line.find("$START:build_line:PARKING_region$") > -1 ):
        ofile.write( line )
        flag = 8
    elif ( line.find( "MGNFIELD" ) > -1 ):
        # don't copy MGNFIELD card
        continue
    elif ( line.find( fastforward ) > -1 ):
        ofile.write( line )
        fastforward = fastforward_def
        flag = 0
        continue

    #
    # ...and do it!
    #

    # flag == -1: skip line, and do not check the other
    #   values of flag
    if ( flag == -1 ):
        continue

    # flag == 0: copy line
    elif ( flag == 0 ):
        ofile.write( line )

    # flag == 1: bodies declarations
    elif ( flag == 1 ):

        print(" ...dumping FLUKA BODIEs;")
        ofile.write( "* \n" )

        # BASIC ELEMENTs:
        for basic_element in BASIC_ELEMENTs.values():
            line = basic_element.ret_FLUKA_include_statement( "body" )
            ofile.write( "%s \n" % ( line ) )
            ofile.write( "* \n" )

        # a header:
        ofile.write( HEADERs["assembly"] % ( "bodies for assembly" ) + "\n" )

        # BODIEs:
        for body in BODIEs:
            body.convert2infinite()
            ofile.write( body.returnSTRING() + "\n" )
            ofile.write( "* \n" )

        ofile.write( "* \n" )
        fastforward = "$END:build_line:BODIEs$"
        flag = -1

    # flag == 2: regions declarations
    elif ( flag == 2 ):

        print(" ...dumping FLUKA REGIONs;")
        ofile.write( "* \n" )

        # BASIC elements:
        for basic_element in BASIC_ELEMENTs.values():
            line = basic_element.ret_FLUKA_include_statement( "regn" )
            ofile.write( "%s \n" % ( line ) )
            ofile.write( "* \n" )

        # a header:
        ofile.write( HEADERs["assembly"] % ( "regions for assembly" ) + "\n" )

        # REGIONs:
        for region in REGIONs:
            ofile.write( "%s \n" % ( region.returnSTRING() ) )
            ofile.write( "* \n" )

        ofile.write( "* \n" )
        fastforward = "$END:build_line:REGIONs$"
        flag = -1

    # flag == 3: LATTICE declarations
    elif ( flag == 3 ):

        print(" ...dumping FLUKA LATTICEs;")
        ofile.write( "* \n" )

        for region in REGIONs:
            if ( region.ilatt > 0 ):
                ofile.write( "%s \n" % ( region.returnLATTICE() ) )

        ofile.write( "* \n" )
        fastforward = "$END:build_line:LATTICEs$"
        flag = -1

    # flag == 4: ROT-DEFI declarations
    elif ( flag == 4 ):

        print(" ...dumping FLUKA ROT-DEFIs;")
        ofile.write( "* \n" )

        for transformation in TRANSFs.values():
            for rotdefi in transformation:
                ofile.write( "%s \n" % ( rotdefi.returnSTRING( insert_comment=True ) ) )

        ofile.write( "* \n" )
        fastforward = "$END:build_line:ROT-DEFIs$"
        flag = -1

    # flag == 6: material assignment
    elif ( flag == 6 ):

        print(" ...dumping FLUKA ASSIGNMA cards;")
        ofile.write( "* \n" )

        # FEDB elements:
        # - materials.inp file:
        for tmpMatFile in present_FEDB.ret_material_files():
            line = INCLUDE_STATEMENT_FMT % ( tmpMatFile )
            ofile.write( "%s \n" % ( line ) )
            ofile.write( "* \n" )
        # - element declarations:
        for basic_element in BASIC_ELEMENTs.values():
            line = basic_element.ret_FLUKA_include_statement( "assm" )
            ofile.write( "%s \n" % ( line ) )
            ofile.write( "* \n" )

        # a header:
        ofile.write( HEADERs["assembly"] % ( "material assignment for assembly" ) + "\n" )

        # REGIONs:
        for region in REGIONs:
            ofile.write( "%s \n" % ( region.returnASSIGNMA() ) )

        ofile.write( "* \n" )
        fastforward = "$END:build_line:ASSIGNMAs$"
        flag = -1

    # flag == 8: new prototypes to be subtracted from the parking:
    elif ( flag == 8 ):

        # FEDB lines:

        if ( len( BASIC_ELEMENTs ) > 0 ):

            print(" ...subtracting new prototypes from parking region;")
            ofile.write( "* \n" )

            for basic_element in BASIC_ELEMENTs.values():
                line = basic_element.ret_FLUKA_subtract_bounding_box_from_parking()
                ofile.write( "%s \n" % ( line ) )

            ofile.write( "* \n" )

        fastforward = "$END:build_line:PARKING_region$"
        flag = -1

    continue

ofile.close()
ifile.close()

