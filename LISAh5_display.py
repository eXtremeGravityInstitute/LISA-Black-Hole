#!/usr/bin/env python3

##################################################
#                                                #
# Simple script to display the content LISA hdf5 #
#                  version 1.0                   #
#                                                #
#      A. Petiteau, ...                          #
#      for the LISA Data Challenge team          #
#                                                #
##################################################

import os
from optparse import OptionParser

from LISAhdf5 import LISAhdf5

parser = OptionParser(usage="usage: %prog  [options] YYY.hdf5\n\
                      YYY.hdf5 : name of the hdf5 file where everything is recorded [required]",
                      version="A. Petiteau - 23/12/2017 - $Id$")

parser.add_option("-p", "--disppath",
                  action="store_true", dest="DispPath", default=False,
                  help="display the path [off by default]")

parser.add_option("-v", "--verbose",
                  action="store_true", dest="verbose", default=False,
                  help="display parameter values [off by default]")

(options, args) = parser.parse_args()

if len(args) < 1:
    parser.error("You must specify the hdf5 file!")

if not os.path.isfile(args[0]):
    print("ERROR in LISA5_display.py : hdf5 file", args[0], "not found")

LH = LISAhdf5(args[0])
LH.display(ReturnStr=False, Print=True, DispPath=options.DispPath)
