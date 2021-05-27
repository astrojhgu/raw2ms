#!/usr/bin/env python3

from casacore.tables import table
import astropy.units as u
import astropy.coordinates as ac

import sys

obslist=table(sys.argv[1], readonly=False)
if '21CMA' in obslist.col('Name'):
    print('21CMA has alread patched')
    sys.exit(0)

ulastai_site={'MJD': 54466,
 'Name': '21CMA',
 'Type': 'ITRF',
 'Long': 86.71737156,
 'Lat': 42.92424534,
 'Height': 2598.13130234,
 'X': 267959.64783333,
 'Y': 4671913.08883333,
 'Z': 4323112.49433333,
 'Source': '',
 'Comment': '',
 'AntennaResponses': ''}

obslist.addrows()
obslist.row().put(obslist.nrows()-1, ulastai_site)
