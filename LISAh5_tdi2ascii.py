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

import os, sys
import numpy as np
from optparse import OptionParser

from LISAhdf5 import LISAhdf5

parser = OptionParser(usage="usage: %prog  [options] YYY.hdf5 ZZZ.txt\n\
                      YYY.hdf5 : name of the hdf5 file where everything is recorded [required]\n\
                      ZZZ.hdf5 : name of the output ascii file that will contain time the TDI generators",
                      version="A. Petiteau - 23/12/2017 - $Id$")

parser.add_option("-v", "--verbose",
                  action="store_true", dest="verbose", default=False,
                  help="display parameter values [off by default]")

(options, args) = parser.parse_args()

if len(args) < 2:
    parser.error("You must specify the hdf5 file and the ascii output file!")

if not os.path.isfile(args[0]):
    print("ERROR in LISA5_display.py : hdf5 file", args[0], "not found")

LH = LISAhdf5(args[0])
tdi = LH.getPreProcessTDI()
#tdi = (LH.get('/H5LISA/PreProcess/TDIGenerator')).split(',')

print (np.shape(tdi))

fOut = open(args[1],'w')

### Header
r = "#time(s)   X    Y   Z"
# for x in tdi:
#     r = r + " " + x
fOut.write(r+"\n")

### Data
for iT in range(tdi.shape[0]):
    r = "%.15e"%(tdi[iT,0])
    for iR in range(1, tdi.shape[1]):
        r = r + " %.15e"%(tdi[iT,iR])
    fOut.write(r+"\n")
fOut.close()
