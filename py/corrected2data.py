#!/usr/bin/env python2

import sys
from pyrap.tables import table
import pyrap.tables as tbl

t=table(sys.argv[1],readonly=False)

if not ('DATA1' in t.colnames()):
    t.addcols(tbl.makecoldesc("DATA1", t.getcoldesc('DATA')))

t.putcol('DATA1',t.getcol('CORRECTED_DATA'))

