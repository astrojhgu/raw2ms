#!/usr/bin/env python

import numpy
import struct
import sys
import cmath

if len(sys.argv)!=3:
    print("Usage:{0} <node list> <gain_prefix>".format(sys.argv[0]))
    sys.exit(-1)


ch_limits=(2048,8192)

gain_prefix=sys.argv[2]

size_of_float=4

nodes=[i.strip() for i in open(sys.argv[1])]
print(nodes)

nantennas=len(nodes)

phases=[]#phase for each channel

for i in range(ch_limits[0],ch_limits[1]):
    fname=gain_prefix+"{0}.dat".format(i)
    fgain=open(fname,'rb')
    print(fname)
    phase_sum=[0.0 for i in range(0,nantennas)]
    cnt=[0.0 for i in range(0,nantennas)]
    while True:
        buffer=fgain.read(nantennas*2*size_of_float)
        if buffer:
            data=struct.unpack('<{0}f'.format(2*nantennas),buffer)
            for j in range(0,nantennas):
                if cmath.isfinite(data[j]):
                    pf=cmath.exp(1j*data[j])
                    phase_sum[j]+=pf
                    cnt[j]+=1
            #print(data[0:nantennas])
        else:
            break

    for j in range(0,nantennas):
        if cnt[j]>0:
            phase_sum[j]/=cnt[j]
    
    phases.append([cmath.phase(j) for j in phase_sum])
    

for i in range(0,nantennas):
    fout=open(nodes[i]+"_gain.txt",'w')
    for j in range(ch_limits[0],ch_limits[1]):
        fout.write("{0} {1}\n".format(j,phases[j-ch_limits[0]][i]))
