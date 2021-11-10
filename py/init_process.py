#vis_prefix="6144:6400"
complist="components.cl"
model=""

vis=vis_prefix+".MS"
flagversion="init.flag"

flagdata(vis=vis,mode="tfcrop",autocorr=False,inpfile="",reason="any",tbuff=0.0,spw="",field="",antenna="",uvrange="",timerange="",correlation="",scan="",intent="",array="",observation="",feed="",clipminmax=[],datacolumn="DATA",clipoutside=True,channelavg=False,timeavg=False,timebin="0s",clipzeros=False,quackinterval=1.0,quackmode="beg",quackincrement=False,tolerance=0.0,addantenna="",lowerlimit=0.0,upperlimit=90.0,ntime="scan",combinescans=False,timecutoff=4.0,freqcutoff=3.0,timefit="line",freqfit="poly",maxnpieces=7,flagdimension="freqtime",usewindowstats="none",halfwin=1,extendflags=True,winsize=3,timedev="",freqdev="",timedevscale=5.0,freqdevscale=5.0,spectralmax=1000000.0,spectralmin=0.0,extendpols=True,growtime=50.0,growfreq=50.0,growaround=False,flagneartime=False,flagnearfreq=False,minrel=0.0,maxrel=1.0,minabs=0,maxabs=-1,spwchan=False,spwcorr=False,basecnt=False,fieldcnt=None,name="Summary",action="apply",display="",flagbackup=True,savepars=False,cmdreason="",outfile="",writeflags=None)

flagmanager(vis=vis,mode="save",versionname="init.flag",oldname="",comment="",merge="replace")

img_prefix=vis_prefix+"_img"
ft(vis=vis,field="",spw="",model=model,nterms=1,reffreq="",complist=complist,incremental=False,usescratch=True)


caltable=vis_prefix+".cal"

bandpass(vis=vis,caltable=caltable,field="",spw="",intent="",selectdata=False,timerange="",uvrange="",antenna="",scan="",observation="",msselect="",solint="5min",combine="scan",refant="E01",minblperant=4,minsnr=2,solnorm=False,bandtype="B",smodel=[],append=False,fillgaps=0,degamp=3,degphase=3,visnorm=False,maskcenter=0,maskedge=5,docallib=False,callib="",gaintable=[],gainfield=[],interp=[],spwmap=[],parang=False)

applycal(vis=vis,field="",spw="",intent="",selectdata=False,timerange="",uvrange="",antenna="",scan="",observation="",msselect="",docallib=False,callib="",gaintable=[caltable],gainfield=[],interp=[],spwmap=[],calwt=[True],parang=False,applymode="",flagbackup=True)

imagename=img_prefix
tclean(vis=vis,selectdata=False,field="",spw="",timerange="",uvrange="",antenna="",scan="",observation="",intent="",datacolumn="corrected",imagename=imagename,imsize=[4096, 4096],cell="0.25arcmin",phasecenter="J2000 00h00m00s +90d00m00s",stokes="I",projection="ARC",startmodel="",specmode="cube",reffreq="",nchan=-1,start="",width=16,outframe="",veltype="radio",restfreq=[],interpolation="linear",gridder="standard",facets=1,wprojplanes=1,aterm=None,psterm=False,wbawp=True,conjbeams=True,cfcache="",computepastep=360.0,rotatepastep=360.0,pblimit=0.2,normtype="flatnoise",deconvolver="hogbom",scales=[],nterms=2,restoringbeam=[],outlierfile="",weighting="uniform",robust=0.5,npixels=0,uvtaper=[''],niter=1000,gain=0.1,threshold=0.0,cycleniter=100,cyclefactor=0.5,minpsffraction=0.05,maxpsffraction=0.8,interactive=False,mask="",savemodel="none",calcres=True,calcpsf=True,parallel=False)


exportfits(imagename=img_prefix+".image",fitsimage="img.fits",velocity=False,optical=False,bitpix=-32,minpix=0,maxpix=-1,overwrite=True,dropstokes=False,stokeslast=True,history=True,dropdeg=True)
