# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2023.                 #
# ######################################### #


source = r'''
#ifndef XCOLL_INTERACTIONS_H
#define XCOLL_INTERACTIONS_H

#define  XC_ABSORBED                      -1
    
#endif /* XCOLL_INTERACTIONS_H */
'''

interactions = {
    int(line.split()[2]): line.split()[1][3:].replace('_',' ').title()
    for line in source.split('\n')
    if len(line.split()) > 1 and line.split()[1][:3] == 'XC_' # select the source lines with the definitions
}
