#!/usr/bin/env python
import scipy.fftpack

import pyfits
import sys

mxr=pyfits.open(sys.argv[1])[0].data
mxi=pyfits.open(sys.argv[2])[0].data

mx=mxr+mxi*(1j)

img=scipy.fftpack.fftshift(scipy.fftpack.fft2(scipy.fftpack.fftshift(mx))).real

pyfits.PrimaryHDU(img).writeto('img.fits',clobber=True)

