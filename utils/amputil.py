#!/usr/bin/env python

# $Id: amputil.py,v 1.3 2006/07/21 20:17:22 anderson Exp $
# ---------------------------------------------------------------------
# functions for measuring the amplifier discontinuity

__version__      = '$Revision: 1.3 $ '[11:-3]
__version_date__ = '$Date: 2006/07/21 20:17:22 $ '[7:-3]
__author__       = 'J. Blakeslee <jpb@pha.jhu.edu>'

import os,math,fUtil
from numpy import where,resize,less_equal,arange

def getstep(ff, sky, sig, A1maxcol=-1, ext=0, nsig=3.5, verb=1, dowrite=0):
    """find size of step across amplifier boundary.
    Returns a tuple:  (stepsize,error,Nr),
    where error is uncertainty and Nr tells how many values went
    into median estimate.  This one uses 3 cols on each side.
    """
    retnull = (0,-1,0)
    dat = ff[ext].data
    if len(dat.shape) != 2:
        print "amp.getstep: data array not 2-d!"
        return retnull

    if A1maxcol <= 0:
        A1maxcol = dat.shape[1]/2 - 1
        if verb:
            print "Taking amp boundary at %d/%d." %(A1maxcol,A1maxcol+1)

    ny = dat.shape[0]
    thresh = nsig*sig

    steplist=[]
    for j in range(ny):
        a1q = dat[j,A1maxcol-5]   # 2042
        a1p = dat[j,A1maxcol-4]   # 2043
        a1o = dat[j,A1maxcol-3]   # 2044
        # a1n = dat[j,A1maxcol-2]  # 2045  is this one ok to use?  better not
        a1m = dat[j,A1maxcol-1]   # 2046
        #a1m = dat[j,A1maxcol]     2047
        #a2m = dat[j,A1maxcol+1]   2048
        a2m = dat[j,A1maxcol+2]   # 2049
        a2n = dat[j,A1maxcol+3]   # 2050
        a2o = dat[j,A1maxcol+4]   # 2051
        a2p = dat[j,A1maxcol+5]   # 2052
        a2q = dat[j,A1maxcol+6]   # 2053

        okLeft  = okcheck(a1m,sky,thresh) + \
                  okcheck(a1o,sky,thresh) + okcheck(a1p,sky,thresh)
        okRight = okcheck(a2m,sky,thresh) + okcheck(a2n,sky,thresh) + \
                  okcheck(a2o,sky,thresh) + okcheck(a2p,sky,thresh)
        if okLeft < 2 or okRight < 2:
            continue

        leftamp  = getsideMed([a1m,a1o,a1p,a1q],sky,thresh)
        rightamp = getsideMed([a2m,a2n,a2o,a2p,a2q],sky,thresh)
        
        step = rightamp - leftamp
        steplist.append(step)
        del step,rightamp,leftamp,a1m,a1o,a1p,a1q,a2m,a2n,a2o,a2p,a2q

    Ns = len(steplist)
    if verb:  print "Ns = "+str(Ns)+" used in median step estimate."

    if Ns < min(500,ny/4):
        return retnull

    steplist.sort()
    medstep = (steplist[(Ns-1)/2] + steplist[Ns/2])/2.0

    err=0.0
    for i in range(Ns):
        err += (steplist[i]-medstep)**2
    err = math.sqrt(err)/(Ns-1.0)

    if dowrite:
        if verb: print "writing new fits file"
        _writenew(ff[ext],medstep,A1maxcol,verb=verb)

    return (medstep,err,Ns)

def getstep3col(ff, sky, sig, A1maxcol=-1, ext=0, nsig=3.5, verb=1, dowrite=0):
    """find size of step across amplifier boundary.
    Returns a tuple:  (stepsize,error,Nr),
    where error is uncertainty and Nr tells how many values went
    into median estimate.  This one uses 3 cols on each side.
    """
    retnull = (0,-1,0)
    dat = ff[ext].data
    if len(dat.shape) != 2:
        print "amp.getstep: data array not 2-d!"
        return retnull

    if A1maxcol <= 0:
        A1maxcol = dat.shape[1]/2 - 1
        if verb:
            print "Taking amp boundary at %d/%d." %(A1maxcol,A1maxcol+1)

    ny = dat.shape[0]
    thresh = nsig*sig

    steplist=[]
    for j in range(ny):
        a1p = dat[j,A1maxcol-2]
        a1n = dat[j,A1maxcol-1]
        a1m = dat[j,A1maxcol-3]
        #a1m = dat[j,A1maxcol]
        #a2m = dat[j,A1maxcol+1]
        a2m = dat[j,A1maxcol+4]
        a2n = dat[j,A1maxcol+2]
        a2p = dat[j,A1maxcol+3]

        if not (okcheck(a1m,sky,thresh) or okcheck(a1n,sky,thresh)):
            continue
        if not (okcheck(a2m,sky,thresh) or okcheck(a2n,sky,thresh)):
            continue

        leftamp,nn = getside([a1p,a1n,a1m],sky,thresh)
        rightamp,nn = getside([a2p,a2n,a2m],sky,thresh)
        
        step = rightamp - leftamp
        steplist.append(step)
        del step,rightamp,leftamp,nn,a1p,a1n,a1m,a2m,a2n,a2p

    Ns = len(steplist)
    if verb:  print "Ns = "+str(Ns)+" used in median step estimate."

    if Ns < min(500,ny/4):
        return retnull

    steplist.sort()
    medstep = steplist[(Ns-1)/2]

    err=0.0
    for i in range(Ns):
        err += (steplist[i]-medstep)**2
    err = math.sqrt(err)/(Ns-1.0)

    if dowrite:
        if verb: print "writing new fits file"
        _writenew(ff[ext],medstep,A1maxcol,verb=verb)

    return (medstep,err,Ns)


def getstep2col(ff, sky, sig, A1maxcol=-1, ext=0, nsig=2.7, verb=1, dowrite=0):
    """find size of step across amplifier boundary.
    Returns a tuple:  (stepsize,error,Nr),  where error is uncertainty
    and Nr tells how many values went into median estimate.  This version
    uses 2 cols on each side of amplifier boundary.
    """
    retnull = (0,-1,0)
    dat = ff[ext].data
    if len(dat.shape) != 2:
        print "amp.getstep: data array not 2-d!"
        return retnull

    if A1maxcol <= 0:
        A1maxcol = dat.shape[1]/2 - 1
        if verb:
            print "Taking amp boundary at %d/%d." %(A1maxcol,A1maxcol+1)

    ny = dat.shape[0]
    thresh = nsig*sig

    steplist=[]
    for j in range(ny):
        a1n = dat[j,A1maxcol-1]
        a1m = dat[j,A1maxcol]
        a2m = dat[j,A1maxcol+1]
        a2n = dat[j,A1maxcol+2]
        ok=1

        if (not okcheck(a1m,sky,thresh)) or (not okcheck(a2m,sky,thresh)):
            ok=0
        elif (not okcheck(a1n,sky,thresh)) and (not okcheck(a2n,sky,thresh)):
            step = a2m - a1m

        elif (okcheck(a1n,sky,thresh)) and (okcheck(a2n,sky,thresh)):
            step = 0.5*((a2m+a2n) - (a1m+a1n))

        elif okcheck(a1n,sky,thresh):
            step = a2m - 0.5*(a1m+a1n)
        else:
            step = 0.5*(a2m+a2n) - a1m

        if ok:
            steplist.append(step)

    Ns = len(steplist)
    if verb:
        print "Ns = "+str(Ns)+" used in median step estimate."

    if Ns < min(1000,ny/4):
        return retnull

    steplist.sort()
    medstep = steplist[(Ns-1)/2]

    err=0.0
    for i in range(Ns):
        err += (steplist[i]-medstep)**2
    err = math.sqrt(err)/(Ns-1.0)

    if dowrite:
        if verb: print "writing new fits file"
        _writenew(ff[ext],medstep,A1maxcol,verb=verb)

    return (medstep,err,Ns)

def getsideMed(colist,sky,thresh):
    " return clipped median of colist"
    N = len(colist)
    newList = []
    for i in range(N):
        if okcheck(colist[i],sky,thresh):
            newList.append(colist[i])

    nn = len(newList)
    if nn < 1:
        print "Bug in amputil: amputil.getsideMed found zero-length list!"
        raise Exception,"amputil: getsideMed found zero-length list!"
    
    newList.sort()
    med = (newList[nn/2] + newList[(nn-1)/2])/2.0
    return med
        
def getside(colist,sky,thresh):
    " return clipped ave of colist"
    N = len(colist)
    if N < 1:
        return (0,0)

    sum = 0.0
    ngood = 0
    for i in range(N):
        if okcheck(colist[i],sky,thresh):
            sum += colist[i]
            ngood += 1

    ave = sum/ngood
    return (ave,ngood)
        

def _writenew(hdu,step,A1maxcol,newname="new.fits",verb=0):
    "write the new, de-amp'ed fitsfile"

    NX = hdu.header['NAXIS1']
    NY = hdu.header['NAXIS2']
    if (NY,NX) != hdu.data.shape:
        raise Exception,"writenew: input header/data inconsistency!"

    if os.path.isfile(newname): os.remove(newname)
    nf = fUtil.createfits(newname, NX, NY, hdu.header['BITPIX'])
    nf[0].header = hdu.header

    colvals = resize(arange(NX),(NY,NX))
    nf[0].data = where(less_equal(colvals,A1maxcol), hdu.data + 0.5*step, hdu.data - 0.5*step)

    nf.close()
    if verb:
        print "wrote "+newname
    del nf,colvals,NX,NY
    return


def okcheck(x,mean,thresh):
    if math.fabs(x-mean) < thresh:
        return 1
    else:
        return 0

