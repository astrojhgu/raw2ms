#vis_prefix="6144:6400"
"""
casa -c selfcal.py <prefix> [solint] [refant]
"""

import sys

vis_prefix=sys.argv[3]
if len(sys.argv)>4:
    solint=sys.argv[4]
else:
    solint="5min"

if len(sys.argv)>5:
    refant=sys.argv[5]
else:
    refant="E01"

vis=vis_prefix+".MS"
flagversion="init.flag"

flagmanager(vis=vis,mode="restore",versionname=flagversion,oldname="",comment="",merge="replace")

img_prefix=vis_prefix+"_img"
model=img_prefix+".model"

#ft(vis=vis,field="",spw="",model=model,nterms=1,reffreq="",complist="",incremental=False,usescratch=True)


caltable=vis_prefix+".cal"

bandpass(vis=vis,caltable=caltable,field="",spw="",intent="",selectdata=False,timerange="",uvrange="",antenna="",scan="",observation="",msselect="",solint=solint,combine="scan",refant=refant,minblperant=4,minsnr=2,solnorm=False,bandtype="B",smodel=[],append=False,fillgaps=0,degamp=3,degphase=3,visnorm=False,maskcenter=0,maskedge=5,docallib=False,callib="",gaintable=[],gainfield=[],interp=[],spwmap=[],parang=False)

applycal(vis=vis,field="",spw="",intent="",selectdata=False,timerange="",uvrange="",antenna="",scan="",observation="",msselect="",docallib=False,callib="",gaintable=[caltable],gainfield=[],interp=[],spwmap=[],calwt=[True],parang=False,applymode="",flagbackup=True)

imagename=img_prefix
tclean(vis=vis,selectdata=False,field="",spw="",timerange="",uvrange="",antenna="",scan="",observation="",intent="",datacolumn="corrected",imagename=imagename,imsize=[4096, 4096],cell="0.25arcmin",phasecenter="J2000 00h00m00s +90d00m00s",stokes="I",projection="ARC",startmodel="",specmode="mfs",reffreq="",nchan=-1,start="",width=16,outframe="",veltype="radio",restfreq=[],interpolation="linear",gridder="standard",facets=1,wprojplanes=1,aterm=None,psterm=False,wbawp=True,conjbeams=True,cfcache="",computepastep=360.0,rotatepastep=360.0,pblimit=0.05,normtype="flatnoise",deconvolver="hogbom",scales=[],nterms=2,restoringbeam=[],outlierfile="",weighting="uniform",robust=0.5,npixels=0,uvtaper=[''],niter=1000,gain=0.1,threshold=0.0,cycleniter=100,cyclefactor=0.5,minpsffraction=0.05,maxpsffraction=0.8,interactive=False,mask="",savemodel="modelcolumn",calcres=True,calcpsf=True,parallel=False)


exportfits(imagename=img_prefix+".image",fitsimage="img.fits",velocity=False,optical=False,bitpix=-32,minpix=0,maxpix=-1,overwrite=True,dropstokes=False,stokeslast=True,history=True,dropdeg=True)
