#!/usr/bin/env python2
import scipy
import pyfits
import sys


f=pyfits.open(sys.argv[1])
f[0].data=scipy.sum(f[0].data,axis=0)
f.writeto(sys.argv[1],clobber=True)
