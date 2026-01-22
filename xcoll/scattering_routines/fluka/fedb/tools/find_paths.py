import os
import sys
#
def set_LB_path():
    """
    1. check that LB_PATH is set as enviroment variable
    2. generate the full path+name of the 'generic_frame.fluka'
    Output: LB_PATH, GENERIC_FRAME
    """
    # is LB_PATH declared?
    assert os.environ["LB_PATH"], "Please set the environment variable LB_PATH"

    # path to LineBuilder + generic_frame.fluka
    LB_PATH = os.environ["LB_PATH"]

    # check it
    msg = ( " \n"
            " the following path to the LineBuilder doesn't exist:\n"
            "    '%s'\n"
            " please set LB_PATH correctly\n") % ( LB_PATH )
    assert os.path.isdir( LB_PATH ), msg

    # build paths to subdirs and append them to sys.path
    for subdir in [ "src", "lib", "additionals" ]:
        loc_path = LB_PATH[:] + "/" + subdir[:]
        msg = ( " \n"
                " the following path doesn't exist:"
                "    '%s'"
                " please set LB_PATH correctly") % ( loc_path )
        assert os.path.isdir( loc_path ), msg
        sys.path.append( loc_path )
        if ( subdir == "additionals" ):
            GENERIC_FRAME = loc_path + "/generic_frame.fluka"

    return LB_PATH, GENERIC_FRAME

# ------------------------------------------------------------------------------
LB_PATH, GENERIC_FRAME = set_LB_path()
# ------------------------------------------------------------------------------
