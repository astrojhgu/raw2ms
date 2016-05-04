#!/usr/bin/env python2

import astropy.io.fits as pyfits
from pyrap.tables import table
import sys
import numpy as np

img_size=1024
max_freq=200E6
c=2.99792458E8
max_bl=2640.00

max_uv=max_bl/(c/max_freq)

if len(sys.argv)!=3:
    print("Usage:{0} <ms name> <outfile>".format(sys.argv[0]))
    sys.exit(-1)

msname=sys.argv[1]
outname=sys.argv[2]


tab=table(msname,readonly=True)
tdesc=table(tab.getkeyword('DATA_DESCRIPTION'),readonly=True)
tspw=table(tab.getkeyword('SPECTRAL_WINDOW'),readonly=True)

nrows=tab.nrows()

tdesc_vec=[i['SPECTRAL_WINDOW_ID'] for i in tdesc.row()]
ch_freq=[i['CHAN_FREQ'] for i in tspw]

mx=np.zeros([img_size,img_size])
wgt=np.zeros([img_size,img_size])

cnt=0
for i in tab.row():
    if cnt%1000==0:
        print("{0}% finished".format(100*cnt/float(nrows)))
    data=i['DATA']
    desc_id=i['DATA_DESC_ID']
    cf=ch_freq[tdesc_vec[desc_id]]
    uvw=i['UVW']
    lbd=c/cf
    u1=uvw[0]/lbd
    v1=uvw[1]/lbd
    w1=uvw[2]/lbd
    #print(u1,v1,w1)
    x=(u1*(img_size/2/max_uv)+img_size/2).astype('int')
    y=(v1*(img_size/2/max_uv)+img_size/2).astype('int')

    for j in range(data.shape[0]):
        mx[x[j],y[j]]+=data[j].real
        wgt[x[j],y[j]]+=1
    cnt+=1

for i in range(0,img_size):
    for j in range(0,img_size):
        if wgt[i,j]>0:
            mx[i,j]/=wgt[i,j]


pyfits.PrimaryHDU(mx).writeto(outname,clobber=True)

