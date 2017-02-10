#!/usr/bin/env python

# $Id: pyblot.py,v 1.20 2004/03/17 22:00:05 anderson Exp $
# ---------------------------------------------------------------------

__version__      = '$Revision: 1.20 $ '[11:-3]
__version_date__ = '$Date: 2004/03/17 22:00:05 $ '[7:-3]
__author__       = "J Blakeslee, jpb@pha.jhu.edu"

import os,popen2
import string,glob,numpy,math
import pyfits,fUtil
from   pyraf import iraf
import pdb

class blotter:
    """ blot back a pydrizzled image:
           pb = pyblot.blotter(drob.parlists, image_dir, [verbose=<0|1> ...])
           pb.runblots()
    """
    def __init__(self, shortparlists, parlists, FitsDir, logfile, verbose=1, clean_up=1, skyKey='ALIGNSKY', hdrGain=0,
                 crlower=None, imNsci=1):
        self.modName   = string.split(string.split(str(self))[0],'.')[0][1:]
        self.shortparlists  = shortparlists
        self.parlists  = parlists
        self.Fits      = FitsDir
        self.verbose   = verbose
        self.crmasks   = {}           # cosmic ray masks names
        self.removeList = []
        self.clean_up   = clean_up
        self.skyKey     = skyKey
        self.hdrGain    = hdrGain
        self.crlower    = crlower
        if imNsci < 1:
            raise ValueError,'Error: pyblot got imNsci = '+imNsci
        self.imNsci     = imNsci
        self.logfile    = logfile
        print self.modName,'version',__version__
        self.logfile.write('Instantiating '+self.modName+' version '+__version__)
        
        # make sure these packages are loaded
        iraf.stsdas()
        iraf.toolbox()
        iraf.imgtool()
        iraf.fourier()
        iraf.fitting()
        iraf.ttools()
        iraf.analysis()
        iraf.dither()
        # flush the cash! twice!
        iraf.flpr()
        iraf.flpr()
        iraf.reset(imtype='fits')  # seems to make deriv task a bit happier
        iraf.set(tmp='./')

    def runblots(self, doblots=1, doderivs=1, doCRrej=1, domedian=1, useMedMask=1):
        """Blot drizzled images from mulitple assocations back to
        the positions of the input images from the ASN table, looping
        over the parlists attribute of drizzleImage object.
        Input parameters are:
          doblots  (0/1): do the blotting? [def: 1(yes)]
          doderivs (0/1): produce deriv images? [def:1]
          doCRrej  (0/1): do CR rejection to produce masks? [def:1]
          domedian (0/1): median stack separate drizzled images for CR comp [def:1]
          useMedMask(0/1): use context images as masks in median stack [def:1]
        """
        # clean out the crmasks dictionary if it has junk in it
        while self.crmasks:
            a = self.crmasks.keys()[-1]
            del self.crmasks[a]
        curdir = os.getcwd()
        os.chdir(self.Fits)
        
        iraf.flpr('imcombine')
        iraf.flpr('blot')
        iraf.flpr('deriv')
        iraf.flpr('driz_cr')

        self.runNum = 0
        i = 0
        for parList in self.parlists:
            shortparList = self.shortparlists[i]
            if self.clean_up:
                self.removeList = []
            self.runNum += 1
            self._blot_asn(shortparList, parList, doblots=doblots, doderivs=doderivs, \
                          doCRrej=doCRrej,domedian=domedian,useMedMask=useMedMask)
            if self.clean_up:
                for file in self.removeList:
                    try: os.remove(file)
                    except: pass
            i = i + 1

        iraf.flpr('imcombine')
        iraf.flpr('blot')
        iraf.flpr('deriv')
        iraf.flpr('driz_cr')
        iraf.flpr()

        os.chdir(curdir)
        return self.crmasks

    def _blot_asn(self,shortparList, parList, doblots=1, doderivs=1, doCRrej=1, domedian=1, \
                 useMedMask=1, wdir=None):
        """Blot a drizzled image back to the positions of the input
        images from the ASN table.  Input parameters are:
          parList:  PyDrizzle parlist attribute
          doblots  (0/1): do the blotting?
          doderivs (0/1): produce deriv images?
          doCRrej  (0/1): do CR rejection to produce masks?
          domedian (0/1): median stack separate drizzled images for CR comp [def:1]
          useMedMask(0/1): use context images as masks in median stack [def:1]
          wdir:  optional directory to change to
        """
        if wdir:
            curdir = os.getcwd()
            os.chdir(wdir)

        # first, check if we want to median-stack drizzle outputs of images
        #pdb.set_trace()
        if domedian:
            medDrizzleImage = self._medDriz2(shortparList,usemask=useMedMask)

        iraf.unlearn(iraf.blot)
        iraf.unlearn(iraf.deriv)
        for parDict in parList:
            # Get name of orig MEF fits file and extname,extver of extension
            # origfile:  name of raw MEF file
            # basename:  name of simple-fits file, minus the .fits part
            # simplefits:  full name of simplefits file (with .fits part)
            origfile  = string.split(parDict['data'],'[')[0] 
            tmpstring = string.split(string.split(parDict['data'],'[')[1],']')[0]
            Xname = (string.split(tmpstring,',')[0]).upper()
            Xver  = int(string.split(tmpstring,',')[1])
            basename   = string.split(origfile,'.')[0]+"_"+Xname+"_"+str(Xver)
            simplefits = basename + '.fits' 
            #pdb.set_trace()
            # open the orig file and get needed info for blot from header
            mef = pyfits.open(origfile)
            ###> set blot parameters
            iraf.blot.expout = mef[0].header['EXPTIME']          # Exp time for output image
            iraf.blot.outnx  = mef[Xname,Xver].header['NAXIS1']   # X size for created images
            iraf.blot.outny  = mef[Xname,Xver].header['NAXIS2']   # Y size for created images
            
            if domedian:
                iraf.blot.data  = medDrizzleImage    # Input image
            else:
                iraf.blot.data  = parDict['outdata'] # Input image
            iraf.blot.scale   = parDict['scale']  # size of input pix in terms of output pix
            iraf.blot.coeffs  = parDict['coeffs'] # Coefficients file name
            iraf.blot.xsh     = parDict['xsh']    # X shift to be applied to input image
            iraf.blot.ysh     = parDict['ysh']    # Y shift to be applied to input image
            iraf.blot.rot     = parDict['rot']    # Rotn (deg CCW) of output im w.r.t. input
            # iraf.blot.shft_un = 'output'          # Units of shifts (input or output pixels)
            # iraf.blot.shft_fr = 'output'          # Frame in which shifts are applied
            ##TEST: 15/June/2002: Warren has these set as 'input'.... try it out
            # yes! this fixes the blot rotation problem!
            iraf.blot.shft_fr = 'input'
            iraf.blot.shft_un = 'input'
            iraf.blot.in_un   = parDict['units']  # Units of input image (counts or cps)
            iraf.blot.out_un  = 'counts'          # Units for output image (counts or cps)
            iraf.drizzle.align = 'center'

            iraf.blot.mode    = 'h'               # don't ask for verification!
            # Output data image name construction:
            ### I think "bl_" needs to be appended, not prepended, so other stuff will work
            iraf.blot.outdata = basename + "_bl.fits"
            # print parDict['data'], origfile, Xname, Xver, iraf.blot.outdata

            if doblots:
                if os.path.isfile(iraf.blot.outdata):
                    os.remove(iraf.blot.outdata)
                iraf.blot()
                self.removeList.append(basename + "_bl.fits")
            #pdb.set_trace()
            ###> set deriv parameters
            iraf.deriv.inlist = iraf.blot.outdata
            # following is the name of the output deriv image; can't be changed
            derivname = string.split(iraf.blot.outdata,'.')[0]+'_deriv.fits'
            iraf.deriv.mode    = 'h'               # don't ask for verification!
                
            if doderivs:
                try:
                    if os.path.exists(derivname):
                        os.remove(derivname)
                except: pass
                iraf.deriv()
                self.removeList.append(derivname)

            # make CR masks - get number of input raw images for noise estimate
            Ncombed = mef[Xname,Xver].header.get('NCOMBINE')   # number of input raw images
            if Ncombed == None or Ncombed < 1:
                Ncombed = 1
                print 'Warning: NCOMBINE not in',origfile,'sci,'+str(Xver)+' header. Setting to',Ncombed

            if doCRrej:
                self._do_driz_cr(simplefits, mef, Ncombed)

            mef.close()
            del mef

        # tidy-up: drz* & deriv* file names are arbitrarily chosen by IRAF,
        # so cannot be removed as part of self.removeList
        if doderivs: junkem('drz*')
        if doCRrej:  junkem('deriv*')
        if wdir:     os.chdir(curdir)
        return


    def _run_medianfilter(self,medianfilter_cmd):
        # procedure to run superalign program and check for errors
        self.logfile.write(medianfilter_cmd)
        sproc  = popen2.Popen3(medianfilter_cmd,1)
        output = sproc.fromchild.readlines()
        errs   = sproc.childerr.readlines()
        if errs:
            self.logfile.write('medianfilter choked on '+medianfilter_cmd)
            for f in errs:
                self.logfile.write(string.strip(f))
                self.errorList.append(('medianfilter',string.strip(f)))
        sproc.fromchild.close()
        sproc.childerr.close()
        return 0

    def _medDriz2(self,parList,usemask=1):
        # rewrite of _medDriz designed to call the _C program medianfilter.
        #   Medianfilter was written to get around a problem with imcombine
        #   where it crashed on a large number of files.

        self.logfile.write("Entered _medDriz2 for asn number %d." %(self.runNum))
        NumIm = len(parList)
        if NumIm == 1:                  #1,2
            nlow = 0
            nhigh = 0
        elif NumIm == 2:                 # 2
            nlow = 0
            nhigh = 1
        elif NumIm == 3:                 # 3    
            nlow = 0
            nhigh = 2
        elif NumIm == 4:                 # 4 HRC;   2 WFC
            nlow = 0
            nhigh = 3
        elif NumIm == 5:                 # 5 HRC
            nlow = 1
            nhigh = 3
        elif NumIm == 6:                 # 6 HRC;   3 WFC
            nlow = 1
            nhigh = 4
        elif NumIm == 7:                 # 7 HRC
            nlow = 1
            nhigh = 5
        elif NumIm < 10:                 # 8,9 HRC;  4 WFC
            nlow = 1
            nhigh = 6
        elif NumIm < 12:                 # 10,11 HRC;   5 WFC
            nlow = 2
            nhigh = 7
        # next added/changed
        else:
            nlow = (NumIm+2)/4 - 1
            nhigh = 3*NumIm/4

        output = open('medianfilter_input','w')
        S = '%i %i %i 0\n' % (len(parList),nlow,nhigh)
        output.write(S)        
        for ii in range(len(parList)):
            for jj in range(ii+1,len(parList)):
                if parList[ii]['outdata'] == parList[jj]['outdata']:
                    raise KeyError,"requested to median stack same images!"
            S = "%s %s\n" % (parList[ii]['outdata'],parList[ii]['outcontext'])
            output.write(S)  
        output.close()

        outfile = 'medriz_'+str(self.runNum)+'.fits'
        medianfilter_cmd = "medianfilter medianfilter_input %s" % (outfile)
        self.logfile.write(medianfilter_cmd)
        #os.system(medianfilter_cmd)
        self._run_medianfilter(medianfilter_cmd)

        if self.verbose:
            print 'NumIm = %d  nlow,high:  %d %d' %(NumIm,nlow,nhigh)
            print 'median image',outfile,'created'
        self.logfile.write('median image '+outfile+' created.')
        return outfile
        # pdb.set_trace()

    def _medDriz(self,parList,usemask=1):
        """median-combined separate drizzle images.  there is an iraf dependant kludge
       in here involving the inputString which the a string of input files for imcombine.
       iraf cannot apparently handle a string over a certain size (guessing 512 char).
       So, we now write out a temp file, imcombine_input, which is just a list of the
       input files imcombine is use.  We use the iraf idiom input = "@file" to get this
       task to run.  So far it seems to work.
       """        
        drizList = []
        maskList = []
        iraf.flpr('imcombine')
        iraf.unlearn(iraf.imcalc)
        iraf.imcalc.pixtype = 'short'
        self.logfile.write("Entered _medDriz for asn number %d." %(self.runNum))
        for ii in range(len(parList)):
            for jj in range(ii+1,len(parList)):
                if parList[ii]['outdata'] == parList[jj]['outdata']:
                    raise KeyError,"requested to median stack same images!"
            drizList.append(parList[ii]['outdata'])
            if usemask:
                plmask = string.split(parList[ii]['outcontext'],'.')[0] + '.pl'
                try:   os.remove(plmask)
                except: pass
                if self.verbose:
                    print 'making',plmask,'from',parList[ii]['outcontext'],'...'
                iraf.imcalc(parList[ii]['outcontext'], plmask, "if im1 .eq. 0 then 0 else 1")
                maskList.append(plmask)
                self.removeList.append(plmask)
                
        # construct input list and add masks info to the headers
        inputString = drizList[0]

        # ***> NOTE!  If list too big, imcombine crashes!
        # ***> NOTE!  define 76 as a safe maximum, or 80...
        MAX_IM = 80
        NumTot = len(drizList)
        NumIm = min(NumTot,MAX_IM)
        for ii in range(NumIm):
            if ii>0:
                inputString = inputString +','+drizList[ii]
            if usemask:
                fUtil.fixHeader(drizList[ii],[('BPM',maskList[ii])])
        if self.verbose:
            print 'median stacking: ',inputString
            if usemask: print '  with masks: ', maskList

        #if that all checks out, go ahead and median
        iraf.unlearn(iraf.imcombine)

        # 15/Apr/2002, jpb: want to keep all the medriz's around for CR-rej debugging
        outfile = 'medriz_'+str(self.runNum)+'.fits'
        # self.removeList.append(outfile)
        try:
            os.remove(outfile)
        except:
            pass
        
        # temp file for iraf input because the list might be too big.
        filekludge   = open("imcombine_input","w")
        newinputList = inputString.split(',')

        if(NumIm  != len(newinputList)):
            errtxt = "ERROR! Error:  NumIm != len(newinputList) in _medDriz ?!"
            print errtxt
            self.logfile.write(errtxt)
            
        for item in newinputList:
            filekludge.write(item+"\n")
        filekludge.close()

        #iraf.imcombine.input = inputString         # this is what we used to do.
        iraf.imcombine.input = "@imcombine_input"
        iraf.imcombine.output = outfile
        # iraf.imcombine.plfile = ''
        iraf.imcombine.sigma = ''
        iraf.imcombine.combine = 'median'
        iraf.imcombine.reject = 'minmax'
        iraf.imcombine.outtype = 'real'
        iraf.imcombine.offsets = 'none'
        if usemask:
            iraf.imcombine.masktype = 'badvalue'
        else:
            iraf.imcombine.masktype = 'none'
        iraf.imcombine.maskvalue = 0.
        iraf.imcombine.scale = 'exposure'
        iraf.imcombine.expname = 'EXPTIME'
        iraf.imcombine.nkeep = 1

        # paradoxically, this is not what we want
        # NumIm = len(drizList)/self.imNsci
        # imcombine considers the total number of images being
        # stacked, not the number at any given point!
        
        # NOTE: nhigh must be >= NumIm/2 if cr rejection to be done everywhere
        if NumIm == 1:                  #1,2
            iraf.imcombine.nlow = 0
            iraf.imcombine.nhigh = 0
        elif NumIm == 2:                 # 2
            iraf.imcombine.nlow = 0
            iraf.imcombine.nhigh = 1
        elif NumIm == 3:                 # 3    
            iraf.imcombine.nlow = 0
            iraf.imcombine.nhigh = 2
        elif NumIm == 4:                 # 4 HRC;   2 WFC
            iraf.imcombine.nlow = 0
            iraf.imcombine.nhigh = 3
        elif NumIm == 5:                 # 5 HRC   
            iraf.imcombine.nlow = 1
            iraf.imcombine.nhigh = 3
        elif NumIm == 6:                 # 6 HRC;   3 WFC
            iraf.imcombine.nlow = 1
            iraf.imcombine.nhigh = 4
        elif NumIm == 7:                 # 7 HRC
            iraf.imcombine.nlow = 1
            iraf.imcombine.nhigh = 5
        elif NumIm < 10:                 # 8,9 HRC;  4 WFC
            iraf.imcombine.nlow = 1
            iraf.imcombine.nhigh = 6
        elif NumIm < 12:                 # 10,11 HRC;   5 WFC
            iraf.imcombine.nlow = 2
            iraf.imcombine.nhigh = 7
        # next added/changed
        else:
            iraf.imcombine.nlow = (NumIm+2)/4 - 1
            iraf.imcombine.nhigh = 3*NumIm/4
            
        iraf.imcombine.mode = 'h'
        self.logfile.write(self.modName+' calling imcombine. come in imcombine.  NumIm/nlo/nhi: '+\
                           str(NumIm)+' '+str(iraf.imcombine.nlow)+' '+str(iraf.imcombine.nhigh)+\
                           (' [NumTot=%d]'%(NumTot)))
        iraf.imcombine()
        if self.verbose:
            print 'NumIm = %d  nlow,high:  %d %d' %(NumIm,iraf.imcombine.nlow,iraf.imcombine.nhigh)
            print 'median image',outfile,'created'
        self.logfile.write('median image '+outfile+' created.  Removing imcombine_input temp file.')
        print "removing imcombine_input temp file."
        try:    os.remove("imcombine_input")
        except: pass
        return outfile

    def _do_driz_cr(self, simplefits, mef, Ncombed):
        """ find cosmic rays using iraf.dither.driz_cr """
        iraf.unlearn(iraf.driz_cr)
        inputfits = '_dz_' + simplefits # WZ
        #os.rename(simplefits,'temp.fits') #WZ
        #os.rename(inputfits,simplefits)
        iraf.driz_cr.inlist = simplefits
        iraf.driz_cr.group  = 0

        #pdb.set_trace()
        if self.crlower:
            iraf.driz_cr.SNR    = "4.4 2.0"
            iraf.driz_cr.scale  = "0.6 0.4"
        else:
            # iraf.driz_cr.SNR    = "5.0 3.0"    
            # iraf.driz_cr.SNR    = "4.8 3.0"    #<-- changed from 5.0 to 4.8
            # iraf.driz_cr.SNR    = "4.7 2.7"
            
            #iraf.driz_cr.SNR    = "4.7 2.6"
            #iraf.driz_cr.SNR    = "4.7 2.4"
            
            # iraf.driz_cr.scale  = "0.9 0.5"
            # 15/Apr/2002, jpb: changing yet again - does this help?
            # iraf.driz_cr.scale  = "1.0 0.6"
            # iraf.driz_cr.scale  = "1.2 0.8"
            # 18/Jun/2002 now that the blot rotation problem is fixed...
            # iraf.driz_cr.scale  = "1.0 0.6"
            # iraf.driz_cr.scale  = "1.2 0.6"
            #iraf.driz_cr.scale  = "1.3 0.6"

            # ARGH, still problems with centers of some stars
            # has to do with the smooth of the derivative in driz/blotting
            # going with this 14/Nov/2002:
            # iraf.driz_cr.SNR    = "4.5 2.4"   tweak again 18/March/2003
            iraf.driz_cr.SNR    = "4.5 2.1"
            iraf.driz_cr.scale  = "1.7 0.7"

        if self.verbose:
            print "Using driz_cr S/N rejection thresholds: ", iraf.driz_cr.SNR
            print "and derivative scales: ", iraf.driz_cr.scale
        self.logfile.write(' pyblot driz_cr rej thresh: '+iraf.driz_cr.SNR+'   deriv scales: '+iraf.driz_cr.scale)
            
        iraf.driz_cr.backg  = self.skyKey

        # get gain, readnoise from image header
        _tmp_hdr_gain_,_tmp_single_rn = _gain_rn(mef,self.logfile,ext=0)
        # increase readnoise by sqrt of # of constituent raw images
        iraf.driz_cr.rn = math.sqrt(1.0*Ncombed) * _tmp_single_rn

        #  only use the header gain if hdrGain has been set
        if self.hdrGain:
            iraf.driz_cr.gain = _tmp_hdr_gain_
        else:
            iraf.driz_cr.gain = 1.0
        del _tmp_hdr_gain_,_tmp_single_rn
        iraf.driz_cr.mode   = 'h'
        print '\n ***> driz_cr using gain,rn_tot = ', iraf.driz_cr.gain, iraf.driz_cr.rn
        print ' ***> (NCOMBINE = '+str(Ncombed)+')'

        plmask   = string.split(simplefits,'.')[0] + '_cr.pl'
        fitsmask = string.split(simplefits,'.')[0] + '_cr.fits'
        try:    os.remove(plmask)
        except: pass
        try:    os.remove(fitsmask)
        except: pass

        if self.verbose:
            # print ' running driz_cr() on ',simplefits
            print ' running driz_cr() on ',inputfits #WZ
            print '   output: ', plmask
        #pdb.set_trace()
        iraf.driz_cr()

        #os.rename(simplefits,inputfile) # WZ
        #os.rename('temp.fits',simplefits) 

        # convert mask from pixel list to short integer fits
        iraf.unlearn(iraf.imcalc)
        iraf.imcalc.pixtype = 'short'
        
        # following line is not right -- drizzle uses 1=good, 0=bad
        #iraf.imcalc(plmask, fitsmask, "if im1 .gt. 0 then 0 else 1")
        if self.verbose:
            print ' converting',plmask,'to',fitsmask
        iraf.imcalc(plmask, fitsmask, "if im1 .eq. 0 then 0 else 1")
        # the pixel-list version can be removed immediately
        if self.clean_up:
            os.remove(plmask)

        # finally, save the mask name in a dictionary
        nzero = numMaskZeros(fitsmask)
        self.crmasks[simplefits] = (fitsmask,nzero)
        if self.verbose:
            print 'Masked '+str(nzero)+' CR pixels in '+simplefits
        return fitsmask


# the following are a couple little utility functions not specific
# to blotting, etc, so I leave them out of the class defintion

def _gain_rn(ff, logfile, ext=0, xname=None, xver=None):
    """ takes a fits file and either extenstion number or
    extname,extver pair and returns gain,readnoise """

    if xname and xver:
        header = ff[xname,xver].header
    else:
        header = ff[ext].header
        
    try: gain = header['CCDGAIN']
    except KeyError:
        gain = 1.0
        logfile.write("WARNING: No CCDGAIN keyword found. Using gain = 1.0")
        
    try: amps   = header['CCDAMP']
    except KeyError:
        amps = []
        logfile.write("WARNING: No CCDAMP keyword found.")
        
    rnvals = []
    try:    rnvals.append(header['READNSE'])
    except: pass

    for amp in amps:
        try:    rnvals.append(header['READNSE'+amp])
        except: pass

    if len(rnvals) < 1:
        logfile.write("WARNING: No readnoise found in header!")

    # take the biggest one
    try: rn = max(rnvals)
    except ValueError: rn = 0.0
    
    return gain,rn

def junkem(junkstring):
    removelist = glob.glob(junkstring)
    for removefile in removelist:
        os.remove(removefile)

def numMaskZeros(fitsmask):
    ff = pyfits.open(fitsmask)
    hh = ff[0].data
    # this is a good deal slower in numpy than Numeric,
    # but have to use numpy for pyfits arrays.
    nzero = numpy.multiply.reduce(numpy.empty(hh.shape)) - int(numpy.sum(1.*numpy.sum(hh)))

    ff.close()
    del ff,hh
    return nzero
