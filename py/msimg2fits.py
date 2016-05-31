#!/usr/bin/env python2

from pyrap.tables import table
import sys
import numpy as np
import pyfits

t=table(sys.argv[1])

imgcube=t.row()[0]['map']
img=np.zeros([imgcube.shape[2],imgcube.shape[3]])

for i in range(0,imgcube.shape[0]):
    img+=imgcube[i,0]

pyfits.PrimaryHDU(img).writeto('img.fits',clobber=True)
