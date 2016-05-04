#!/usr/bin/env python
from scipy.optimize import leastsq
from scipy.fftpack import fft
import numpy as np
import scipy as sp
import sys
import cmath
from cmath import pi

channels=[]
phases=[]
c=2.99792458E8
nch=8192
max_freq=200E6

def func(lp):
    result=0
    for i in range(0,len(channels)):
        nu=channels[i]/float(nch)*200E6
        result+=abs(cmath.exp(1j*phases[i])-cmath.exp(1j*lp[0]/c*nu*2*cmath.pi))**2
        #result+=(phases[i]-lp[0]/(c/nu)*2*cmath.pi)**2
    return result

if len(sys.argv)!=2:
    print("Usage:{0} infile".format(sys.argv[0]))
    sys.exit(-10)

for l in open(sys.argv[1]):
    ch,ph=l.split()
    channels.append(int(ch))
    phases.append(float(ph))

    
fphase=fft(phases)

max_value=0
idx=0
for i in range(1,int(len(fphase)/2)):
    v=abs(fphase[i])**2
    if max_value<v:
        max_value=v
        idx=i

lp0=(idx*2*pi/len(fphase)/(2*pi/c/nch*max_freq))
#print(lp0/1.47)
#print(func([lp0]))
#print(func([-lp0]))

pos_slopes=[]
neg_slopes=[]


for i in range(1,len(phases)):
    d=phases[i]-phases[i-1]
    if d>0:
        pos_slopes.append(d)
    else:
        neg_slopes.append(d)

if len(pos_slopes)<len(neg_slopes):
    lp0=-lp0

lp1=leastsq(func,[lp0])[0][0]
print(-lp1/1.47,-lp0/1.47)
