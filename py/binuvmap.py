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

if len(sys.argv)!=8:
    print("Usage:{0} <ant table> <path> <date> <ant1> <ant2> <delay> <outfile>".format(sys.argv[0]))
    sys.exit(-1)

prefix=sys.argv[2]
date=sys.argv[3]
ant1_name=sys.argv[4]
ant2_name=sys.argv[5]
binfile=prefix+"/"+ant1_name+ant2_name+"-0-"+date+".bin"
timefile=prefix+"/"+"time-0-"+date+".txt"
delay=float(sys.argv[6])
outname=sys.argv[7]

ant_tab=table(sys.argv[1],readonly=True)

nrows=tab.nrows()

ant1_row=-1
ant2_row=-1

for i in range(nrows):
    if tab.row()[i]['NANME']==ant1_name:
        ant1_row=i
        break

if ant1_row==-1:
    raise RuntimeException("antenna {0} not found".format(ant1_name))

for i in range(nrows):
    if tab.row()[i]['NANME']==ant2_name:
        ant2_row=i
        break

if ant2_row==-1:
    raise RuntimeException("antenna {0} not found".format(ant2_name))



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

