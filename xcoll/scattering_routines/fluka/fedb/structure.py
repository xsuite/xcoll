import os

def create_FEDB():
    """
    definition of fedb structure
    """
    # global path ( i.e. to present directory):
    global_path = os.path.realpath(os.path.dirname( __file__ ))

    # ...in order to import FEDB class:
    from SUPPORT import FEDB

    # instance of FEDB:
    present_FEDB = FEDB()

    # set global path:
    present_FEDB.set_global_path( global_path )

    # - naming convention:
    present_FEDB.set_naming_convention( "%s_%s" )

    # - materials file:
    present_FEDB.set_material_file_subpath( "materials/materials.inp" )

    # - SUBPATHs:
    present_FEDB.set_SUBPATHs( "body", "bodies/"     + present_FEDB.ret_naming_convention() + ".bodies" )
    present_FEDB.set_SUBPATHs( "regn", "regions/"    + present_FEDB.ret_naming_convention() + ".regions" )
    present_FEDB.set_SUBPATHs( "assm", "materials/"  + present_FEDB.ret_naming_convention() + ".assignmat" )
    present_FEDB.set_SUBPATHs( "stpz", "stepsizes/"  + present_FEDB.ret_naming_convention() + ".stepsizes" )
    present_FEDB.set_SUBPATHs( "smbl", "assemblies/" + present_FEDB.ret_naming_convention() + ".lbp" )

    # - inflect material files:
    present_FEDB.inflect_materials()

    # check integrity of declared structure:
    present_FEDB.check()

    return present_FEDB

# ------------------------------------------------------------------------------
present_FEDB = create_FEDB()
# ------------------------------------------------------------------------------
