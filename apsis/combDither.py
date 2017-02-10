#!/usr/bin/env python

# $Id: combDither.py,v 3.10 2012/06/16 06:50:33 wz Exp $
# ---------------------------------------------------------------------

__version__      = '$Revision: 3.10 $ '[11:-3]
__version_date__ = '$Date: 2012/06/16 06:50:33 $ '[7:-3]
__author__       = 'W Zheng <zheng@pha.jhu.edu>'
# change to multidrizzle

import os,sys,string,glob,math
import path
import numpy 
import re #WZ
import pyfits,matutil
import xydrizzle
import xmlUtil,fUtil,pyblot,augmask,tableio #WZ
import astrometer
import pdb, popen2 #WZ
from   xydrizzle import wcsutil,drutil #WZ
from   pUtil import ptime
from   msg   import pMessage
from   sys   import version
pyversion = version
from   pyraf import iraf
import shutil
import subprocess

# get the version of drizzle by pretending to run it
tmplist = iraf.drizzle("None","None",Stdout=1)
drversion = "unknown"
for line in tmplist:
    wordlist = line.split()
    for __ii in range(len(wordlist)):
        word = wordlist[__ii]
        if word.lower() == 'version':
            drversion = wordlist[__ii+1]
            break
        del word
    del wordlist
del line,tmplist

try:
    pydriz_version = xydrizzle.__version__
except:
    pydriz_version = xydrizzle.version

# Good bits marked by a flag value in the DQ array are
#
#    0  good pixel
#    2  data replaced by fill value
#   64  pre-existing hot pixel
#  128  bias level pixel
#  256  saturation (full well or a-to-d)
# 2048  a-to-d saturation*
# 4096  reserved (possibly for Multidrizzle CR rejected pixels)
# 8192  cosmic ray rejected during cr-split image combination

# Recent changes to DQ bit flags values has indicated that a DQ pixel value of
# 128 should not now be considered "good" as that value, a bias level pixel, now
# indicates bad columns as well.
# 128 bias level pixel (bad columns); was in BPIXTAB, will be in superbias DQ

#_goodBits_ = 2+64+128+256+2048+4096+8192
#_goodBits_ = 2+64+256+512+1024+2048+4096+8192   # good-pixel summed value WFC3
#_goodBits_ = 2+64+256+2048+4096+8192   # good-pixel summed value ACS
_goodBits_ = 32+64+1024   # WZ Use this

#_goodBits_ = 96                     # CALACS good-pixel summed value (see run_drizzle docstring)

_rnVarSuppFac_ = 1.38                  # boost N*rn^2 by this factor for RMS images
_exVarPerSec_  = 0.02222               # add this much extra var per sec for RMS ims

class drizzleImage:
    """Some notes on the combDither module and its implementation of multidrizzle.

    combDither currently contains one class, drizzleImage, which defines
    methods for drizzling dithered image associations to a final images,
    including cosmic ray rejection, etc.  The original intention was to
    create other classes with methods for combining the images using other
    routines (e.g., SWARP, or something based on python numerical arrays).

    The basic procedure for processing a dataset through drizzle:

    drob = drizzleImage(obs,al.MatchDict) - create instance of this class
    drob.run_all()          -   run the drizzle/blot/drizzle sequence
    drob.makeFlagImage()    -   create Sextractor flag image for each asn
    drob.makeRmsImage()     -   create Sextractor rms image for each asn
    drob.fixAstrometry(obs,noGSC) - fix the astrometry of the images if possible.

    [The fixAstrometry method must receive a DataSet object (obs) and will generally
    be called by the apsis pipeline.  This behaviour is not internally consistent
    but the method is a bit of hack to removed the need for the pipeline to
    deal with a new object (an instance of astrometer.gscMatchup) and all the
    attendant messaging for that.  N.B. This method should not be called internally.]

    Alternatively, the user can replace run_all() with run_drizzle() if no
    cosmic ray rejection is desired.  This module no longer makes its own
    association tables, but just uses the original ones that are given as
    inputs to the pipeline.

    The run_all method does three things:

    1. run_drizzle() drizzles images onto separate output images
    2. runs pyblot (median stacks, blots, derivs, cr_rej) which
        returns the list of cr masks for 2nd drizzling
    3. run_drizzle() to drizzles images all together using CR masks
        from pyblot

    run_drizzle() is the heart of the drizzleImage class.  It runs pydrizzle
    on all of the association tables in the asnTableList and takes a number
    of optional parameters (outshape=None, deltmp=1, units='counts',
    bits=_goodBits_, separDriz=0, crmasks=None, exptimePyDrKludge=0).
    Here's what the params are for

    outshape:  x,y sizes of output drizzled image (def: determined 'otf')
    deltmp:    delete temporary (masks, coeffs) files? (def: 1 [yes])
    units:     units of output pixel values (def: counts)
    bits:      sum of all possible good-pixel values in dq array
          ('None' means don't use masks).
          The default value bits=8450 is the sum of the pix values
          that we accept as being good, bits = 0+2+256+8192 = 8450,
          which includes 'good', 'replaced', 'saturated',and
          'CR-repaired', but NOTE THAT THESE VALUES ARE PARTICULAR 
          TO CALACS! [ISR ACS-99-03]

    The remainder of the params relate to pyblot setup and results:
    crmasks:            List of (mask,nzap) tuples from pyblot
    separDriz:          Drizzle onto separate output images?

    run_drizzle loops over all associations and runs pydrizzle on them.  
    It measures the shift 'zeropoint' based on the mean of all the measured
    shifts and then sets the relative shifts accordingly, ensuring that all
    associations will have the same absolute shift zeropoint, as well as the
    same image size.  This *seems* to be ok, although the output WCS may then
    disagree with pydrizzle's idea of the WCS of the output product, but drizzle
    itself constructs the WCS image in the output product from the input and
    whatever image headers it gets.

    For each association, run_drizzle creates the pydrizzle object, then
    goes through the parlist (list of dictionaries holding input drizzle
    parameters).  It modifies these parlists in various way, e.g., it makes
    the output images different if 'separDriz' has been set.  It also sets
    the x,y shifts and rotations to whatever is in the MatchDict (thus,
    making delta shifts in the asn tables no longer necessary).  The list of
    modified parlists is saved as a class attribute.  It then runs drizzle.

    After the 2nd drizzle (separdriz = 0), the code checks if:
    (1) the number of science extensions in the multi-extension FITS (mef)
    file is greater than 1 and units='counts' (the default); if so, it
    divides the count levels by the number of extensions.
    (2) if the pydrizzle version is less than 1.4; if so and Nsci>1, then it
    divides the output exposure times by the number of science extensions.

    If it's the final drizzle, then it also goes through all the input images
    and sums the alignSky values in the headers and writes the sum to the
    science image header.  That's it.  The makeFlagImage and makeRmsImage
    methods are pretty straightforward.

    The drizzle kernel, for the *final drizzle only* can now be set when
    run_drizzle or run_all is called. The drizzle-default linear ('square')
    kernel is the fastest, and tends to smooth over bad pixels better,
    but the psf is not as tight and the noise correlation is much worse.

    The damped sinc interpolant ('lanczos3') is sometimes called the
    optimal kernel, though beware that you may get holes (negative
    pixels) in the output.  These holes can happen anywhere, though more
    commonly near stars, etc.  Because of the masking of bad pixels and
    columns (based on the DQ array), occasionally a pixel in the output
    image will have very little data from the input image, and by a
    quirk of the interpolation kernel with its negative sidelobes, will
    end up with a negative value.  Richard Hook said to check the
    drizzle weight image to make sure that those negative pixels are
    given very little weight; if they have low weight, just don't use
    them, but if their weight is high, then there might be a problem.
    I checked several of these pixels, and found they were given low
    weight, about 1/10 of the median, so it seems ok.  Of course, if
    there are multiple dithered images, these negative pixels get filled
    in quickly by other pixels with much higher weight.  For the linear
    kernel you don't get this effect because there are no negative
    sidelobes.

    One caveat: run_drizzle overwrites shifts, rotations, and output product
    names in the pydrizzle parlist rather than setting them through methods.
    This was historically necessary for the shifts (and *still* the most
    straightforward way), and it's still necessary for drizzling onto
    separate outputs for blotting.  So, I wouldn't recommend changing it,
    but we have to be careful in case later pydrizzles change the parlist
    elements or their meanings.

    """

    def __init__(self, obs, MatchDict, alignSky='ALIGNSKY',hdrGain=0, suppInstVar=1, crlower=None,
                 smartstack=0, notrim=0, padfac=None, outsize=None, outshift=None,origscale=None,noContext=None,
                 maskFile=None,dfilts=(None,None)):
        
        self.modName     = string.split(string.split(str(self))[0],'.')[0][1:]
        self.root        = obs.newobspath     # root path of the observation dir
        self.obsName     = obs.newobs
        self.obsAsnDict  = obs.asnDict
        self.obsAlign    = obs.newalign
        self.obsPars     = obs.newpar
        self.obsFits     = obs.newfits        # path to fits files
        self.obsFitsList = obs.fitslist
        self.messagedir  = obs.messagedir     # where the module message will go.
        self.logfile     = obs.logfile
        self.allShifts   = obs.allShifts      # list of tuples; 2nd element is dictionary of shifts
        self.alignSky    = alignSky           # header keyword for subtracted Sky's
        self.detector    = obs.detector       # detector matters for drizzle shifts
        self.refcdmat    = obs.refcdmat       # dictionary with final cd matrix
        self.refotherkeys = obs.refotherkeys  # dictionary with other keys relating to the wcs solution
        self.edgeMaskName= obs.edgeMaskName = 'Edgemask.fits'
        self.reflogfile_add = obs.reflogfile_add
        self.errorList   = []
        self.inputList   = []
        self.outputList  = {}
        self.hdrGain     = hdrGain
        self.crlower     = crlower
        self.smartstack  = smartstack
        self.notrim      = notrim
        self.padfac      = padfac
        self.outsize     = outsize
        self.outshift    = outshift
        self.maskFile    = maskFile		
        self.origscale   = origscale
        self.noContext   = noContext
        self.dfilts      = dfilts
        self.suppInstVar = suppInstVar
        self.MatchDict   = MatchDict          # MatchDict, from the alignImage object is
                                              # a dictionary of dict's for each image
        self.pardir      = obs.pardir #WZ
        # lists of the various image types . . .
        self.sciImageList     = obs.sciImageList     = []   
        self.contextImageList = obs.contextImageList = []
        self.weightImageList  = obs.weightImageList  = []
        self.flagImageList    = obs.flagImageList    = []
        self.rmsImageList     = obs.rmsImageList     = []
        self.removeList = []
        self.logfile.write('Instantiating drizzleImage object for observation: '+self.obsName)

        if not self.origscale:   # by xingxing    the self.origscale must be defined, or there will be error!! defined from alimage.outscale!!
            if obs.detector =='WFC':
                self.origscale = 0.05
            elif obs.detector == 'HRC':
                self.origscale = 0.025
            elif obs.detector == 'SBC':
                self.origscale = 0.025
            else:
                errtxt = "Detector type "+ obs.detector + " origscale unknown.\n"
                self.errorList.append((self.modName,errtxt))
                self.logfile.write(errtxt)
                raise Exception,errtxt
                
        # set up the inputList for mkMsg() method
        for i in self.obsFitsList:
            self.inputList.append(i)
            self.NumSci = 1
            _firstN = 0
            for im in self.MatchDict.keys():
                if not _firstN:
                    _firstN = self.MatchDict[im]['NumSci']
                    firstImage = im #WZ dec 2010
                else:
                    if self.MatchDict[im]['NumSci'] != _firstN:
                        self.NumSci = 0
                        errtxt = 'WARNING MatchDict purports that images to be combined'+\
                                 ' have different numbers of sci extensions. \n'+\
                                 ' CANNOT fix drizzle scale factor problem in this case!!!!'
                        self.logfile.write(errtxt)
                        self.errorList.append((self.modName,errtxt))
                        break
            if self.NumSci: self.NumSci = _firstN
            #self.logfile.write('setting iraf min_lenuserarea=640000')
            iraf.set(min_lenuserarea=640000)
        #self.base       = path.Env()
        #reffile=os.path.join(self.base.getenv('INGEST'), 'ref.fits')
        # reffile = string.replace(self.root,'datasets','ingest')+'/ref.fits'
        #substring=string.split(self.root,'datasets/')
        #reffile = substring[0]+'ingest/'+string.split(substring[1],'/')[0]+'/ref.fits'
	reffile = os.environ['INGEST']+self.obsName+'/ref.fits'  # by xingxing
        # pdb.set_trace()
        tem=string.split(obs.root,'/')[-1]
	tmp='tweak_'+str.lower(tem)+'.par'
	# tmp='tweak_'+str.lower(obs.filter)+'.par'
	tmp1=os.path.join(obs.ingestdir,'par')
        tmp2=os.path.join(tmp1,tmp)
	tweakfile=os.path.join(obs.ingestdir,'par','tweak_'+str.lower(tem)+'.par')
	# tweakfile=os.path.join(obs.ingestdir,'par','tweak_'+str.lower(obs.filter)+'.par')
	#tweakfile= string.replace(self.root,'datasets','ingest')+'/tweak.par'
        #firstImage = self.inputList[1]
        self.ref=self.setRef()
        self.ref.tweak = tweakfile # WZ piggyback for later use
        if (os.path.exists(reffile)):
            self.ref.prod=self.setgrid(reffile) #WZ: The output frame
            self.ref.model=self.ref.prod.copy()
            self.ref.model.cd11=-self.ref.pscale/3600. # "model" holds the native pixel scale
            self.ref.model.cd12=0.
            self.ref.model.cd21=0.
            self.ref.model.cd22=self.ref.pscale/3600.
            self.ref.model.pscale=self.ref.pscale
        #temname = string.upper(string.split(self.inputList[0],'_asn')[0])
        #ref.band=string.split(temname[-8:],'_')[1]
        # ref.band is the only band subject to astrometric correction if noGSC
        self.wcs=self.ref.copy()
        obs.ref = self.ref.copy()

    def run_all(self, clean_up=0, userKeepBits=None, units='counts',
                asecpix=None, pixfrac=None, kernel=None,
                noRej=None):
        
        """ Run_all()
        This is it, the whole enchilada: drizzle, blot, drizzle again.
        """
        loop=0 #WZ
        self.chi2=9999.
        chi2=9999.
        self.outname=''
        while (loop < 8):
            print '>>>>>>>>  Loop  >>>>>>>>>> ', loop
            #pdb.set_trace()  
            # Remove any files w/the output names so drizzle isn't confused
            if (self.outname and loop>0):
                bck = 'bck'+str(loop-1)
                if os.path.isfile(self.outname): #WZ                                    
                    back=string.replace(self.outname,'drz',bck)
                    os.rename(self.outname,back)
                if not self.noContext:
                    if os.path.isfile(self.outcontext):
                        back=string.replace(self.outcontext,'drz',bck)
                        os.rename(self.outcontext,back)
                if os.path.isfile(self.outweight):
                    back=string.replace(self.outweight,'drz',bck)
                    os.rename(self.outweight,back)
            if noRej: 
                # then just drizzle once
                self.logfile.write('Skipping cosmic ray rejection in combDither!')
                if not kernel:
                    # # kernel = 'square' # default to lanczos3, Feb 25 02, anderson
                    kernel = "lanczos3"
                if userKeepBits:
                    print "Using user supplied good bits value of",userKeepBits
                    self.run_drizzle(deltmp=clean_up, bits=userKeepBits, separDriz=0,
                                     crmasks=None, units=units, asecpix=asecpix, pixfrac=pixfrac, \
                                     kernel=kernel, delwght=0)
                else:
                    print "Using default good bits value of",_goodBits_
                    self.run_drizzle(deltmp=clean_up, bits=_goodBits_, separDriz=0, \
                                     crmasks=None, units=units, asecpix=asecpix, pixfrac=pixfrac, \
                                     kernel=kernel, delwght=0)
                self.logfile.write('Drizzling complete.')

            elif loop == 0: 
	    ##else:
                self.logfile.write('Starting initial drizzle process...')  
		print "loop= ",loop,": run_drizzle first!!! ___________by xingxing"
                if userKeepBits:
                    print "using user supplied good bits value of",userKeepBits
                    self.run_drizzle(deltmp=0, bits=userKeepBits, separDriz=1, \
                                    crmasks=None, units='counts', delwght=clean_up) # xingxing
                else:
                    print "using default good bits value of",_goodBits_
                    fname=os.path.join(self.obsAlign,'alignTweak.txt') #WZ              
                    if (not os.path.exists(fname)): #WZ Tweak alignment                 
                        self.run_drizzle(deltmp=0, bits=_goodBits_, separDriz=1, \
                                    crmasks=None, units='counts', delwght=clean_up)
                    self.run_drizzle(deltmp=0, bits=_goodBits_, separDriz=1, \
                                    crmasks=None, units='counts', delwght=clean_up)

		print "loop=",loop," : pyblot.blotter and .runblots  ___________by xingxing"                
		print "Creating pyblot object."
                self.logfile.write('getting blotter object.')        
                if self.crlower:
                    self.logfile.write('crlower flag set; will use tighter driz_cr SNR thresholds.')
                PyBlOb = pyblot.blotter(self.shortparlists,self.parlists, \
                        self.obsFits,self.logfile, skyKey=self.alignSky, \
                        hdrGain=self.hdrGain,crlower=self.crlower, \
                        clean_up=clean_up, imNsci=self.NumSci)

                self.logfile.write('Blotting...')
                try:
                    crmasks = PyBlOb.runblots()
                    # elements of crmasks dict are now (maskname,nzero) tuples!
                except Exception,err:
                    errtxt = "ERROR: runblots method failed.\n\t"+str(err)
                    self.errorList.append((self.modName,errtxt))
                    self.logfile.write(errtxt)
                    raise Exception,err
                self.logfile.write('blotting complete.')
 
                # get the crmasks listed in the outputList
                for key in crmasks.keys():
                    self.outputList[crmasks[key][0]] = [key]   # can't get at any predecessor info here.

                if not self.notrim:
                    self.logfile.write('Assembling Edge mask...')
                    auger = augmask.augmask(self.augimdict,self.edgeMaskName,logfile=self.logfile,clean_up=clean_up)
                    auger.make()
                    # test to see if edge mask has been made.
                    if os.path.isfile(os.path.join(self.obsFits,self.edgeMaskName)):
                        self.outputList[self.edgeMaskName] = []
                else:
                    self.logfile.write('Skipping Edge mask construction.')
  		
		#if loop > 0 :		
		#    del PyBlOb
                # otherwise, do the full cycle
		self.logfile.write('loop = '+str(loop)+' Skipping the first run_drizzle and pyblot.blotter to save time')
                # self.makeEdgemask() #WZ
                # Remove any files w/the output names so drizzle isn't confused
                if (self.outname):
                    if os.path.isfile(self.outname): #WZ
                        back=string.replace(self.outname,'drz','bck')
                        os.rename(self.outname,back)
                    if not self.noContext:
                        if os.path.isfile(self.outcontext):
                            back=string.replace(self.outcontext,'drz','bck')
                            os.rename(self.outcontext,back)
                    if os.path.isfile(self.outweight):
                        back=string.replace(self.outweight,'drz','bck')
                        os.rename(self.outweight,back)
                # pdb.set_trace()
                # if string.find(self.outname,'drz') > -1:
                #    os.remove(self.parlists[0][0]['outdata'])
                #    os.remove(self.parlists[0][0]['outweight'])
                #    os.remove(self.parlists[0][0]['outcontext'])
            if (not noRej):
                self.logfile.write('Starting second drizzle process...')
                if asecpix:
                    self.logfile.write('  using requested output scale: '+str(asecpix)+' asec/pix.')
                if pixfrac:
                    self.logfile.write('  using requested drizzle pixfrac: '+str(pixfrac))
                if not kernel:
                    kernel = "lanczos3"
                # pdb.set_trace()
		print "loop=",loop," : second run_drizzle !!!!  ___________by xingxing"                
                if userKeepBits:
                    self.logfile.write('  using bits parameter = '+str(userKeepBits))
                    self.run_drizzle(deltmp=clean_up, bits=userKeepBits,separDriz=0,\
                                 crmasks=crmasks, units=units, asecpix=asecpix, pixfrac=pixfrac, \
                                 kernel=kernel, delwght=0)
                else:
                    self.logfile.write('  using bits parameter = '+str(_goodBits_))
                    self.run_drizzle(deltmp=clean_up, bits=_goodBits_,separDriz=0,\
                                 crmasks=crmasks,units=units, asecpix=asecpix, pixfrac=pixfrac, \
                                 kernel=kernel, delwght=0)
                                                                
                if clean_up:
                    self.logfile.write('Removing 1st pass output files.')
                    self._clean_driz()
                else:
                    self.logfile.write('Keeping 1st pass drizzle output files.')

            # alignment calls #WZ
            #substring=string.split(self.root,'datasets/')
            #reffile = substring[0]+'ingest/'+string.split(substring[1],'/')[0]+'/ref.fits'
            # reffile= string.replace(self.root,'datasets','ingest')+'/ref.fits'
	    reffile = os.environ['INGEST']+self.obsName+'/ref.fits'  # by xingxing	
	    #print "Check the reffile here!"
	    #pdb.set_trace() ####  
            self.prepMatch(reffile)
            chi2=self.run_Match(reffile,loop)
            ratio = numpy.abs(self.chi2 - chi2)/self.chi2
            if (chi2 < self.chi2):
                self.chi2 = chi2
		print 'Match result: chi2 = ',chi2,', loop = ',loop   #xingxing
            else:
                print "WARNING: chi2 > self.chi2  !!! ________by xingxing"
		print "chi2 = ",chi2
		print "self.chi2 = ",self.chi2
                break
            # if ((chi2 < 0.1 and self.chi2 < 0.1) or ratio < 0.01)  and loop > 2: #WZ Mar 2013
	    #	print "loop = ", loop     # by xingxing
            # if (chi2 < 0.1 and self.chi2 < 0.1): 
	    #	   print "chi2 = ",chi2,". self.chi2 < 0.1"
	    #	if (ratio < 0.01):
            #       print "ratio < 0.01: ", ratio
            #    break
            # if (match < 0):
            #     print "Something wrong with tweaking. Check tweak.log"
        #pdb.set_trace()
        #print "Second Pass. Tweak =1"
            self.logfile.write('Starting drizzle process with a tweak..')
        #self.logfile.write('Starting initial drizzle process...')
        # if userKeepBits:
        #    print "using user supplied good bits value of",userKeepBits
        #    self.run_drizzle(deltmp=0, bits=userKeepBits, separDriz=1,
        #                 crmasks=None,units='counts',
        #                 delwght=clean_up)
        # else:
        #    print "using default good bits value of",_goodBits_
        #    self.run_drizzle(deltmp=0, bits=_goodBits_, separDriz=1, \
        #     crmasks=None,units='counts', delwght=clean_up)
        
        # Check the final alignment of sci images.

        # self.prepCheck(reffile)
        # match=self.checkFinalMatch(reffile,loop)
        # if (match<0):
        #     print ("Problem with tweaking. Check results")
        # else:
        #    self.logfile.write('Drizzling complete.')
        # f = open('tweak.res','r')
        # for line in f.readlines():
        #    print line[:-1]
        # f.close()
            # print "Happy with the tweaking results? Y or N"
            # read = str.lower(sys.stdin.readline())
            # if (string.find(read,'y')> -1):
            #     loop=-1
            #     break
            # else:
            #    self.updateTweak()
            loop = loop + 1                    
            # pdb.set_trace()
        # 5 > loop  > -1 #WZ
        # pdb.set_trace()
        self.updateTweak()
        self.logfile.write('A full drizzle cycle complete.')
        return

    def _clean_driz(self):
        curdir = os.getcwd()
        os.chdir(self.obsFits)
        for file in self.removeList:
            try:
                os.remove(file)
                self.logfile.write("removed "+file)
            except Exception,err:
                self.logfile.write(str(err))
                print err
        os.chdir(curdir)
        return

    def _get_refAsn(self):
        "choose the reference association for the output image size"
        refAsn = None
        maxshift = -9999.
        
        for tab in self.obsAsnDict.keys():
            xshifts = []
            yshifts = []
            for im in self.obsAsnDict[tab]:
                xshifts.append(self.MatchDict[im]['xpix_shift'])
                yshifts.append(self.MatchDict[im]['ypix_shift'])
            xshifts.sort()
            yshifts.sort()
            _maxi = max((xshifts[-1] - xshifts[0]),(yshifts[-1] - yshifts[0]))
            if _maxi  > maxshift:
                maxshift = _maxi
                refAsn = tab
            del _maxi,xshifts,yshifts

        del tab
        if not refAsn:
            self.logfile.write('Error: ref Asn could not be determined.')
            raise RuntimeError,'Bug in _get_refAsn method.'
        self.logfile.write(refAsn+(' chosen as reference asn table.  maxshift: %.4f'%maxshift))
        print refAsn+' chosen as reference asn table.  maxshift: %.4f'%maxshift
        return refAsn

    def run_drizzle(self,outshape=None, deltmp=1, units='counts',
                    bits=_goodBits_, separDriz=0, crmasks=None,
                    asecpix=None, pixfrac=None, kernel='square', delwght=0): 
        
        """ run_drizzle():
        Run drizzle on all asn tables in the obsAsnDict.keys()
           Input parameters (all optional) include:
             outshape: x,y sizes of output drizzled image (def: determined 'otf')
             deltmp:   delete temporary (masks, coeffs) files (def: 1 [yes])
             units:    units of output pixel values (def: counts)
             bits:     sum of possible good-pixel values in dq array (None => don't use masks).
                       The default bits=8578 is the sum of the pix values that we accept
                       as being good, bits = 0+2+256+8192+128 = 8578, which includes
                       'good', 'replaced', 'saturated', 'CR-repaired' and 'bias-level pixel' but
                       NOTE THAT THESE VALUES ARE PARTICULAR TO CALACS! [ISR ACS-99-03]
             The remainder of the params relate to pyblot setup and results:
             crmasks:      List of masks from pyblot
             separDriz:    Drizzle onto separate output images
             asecpix:      Scale in arsec/pix of final drizzled images.
             pixfrac:      Drizzle 'pixfrac' or 'dropsize' parameter.
             kernel:       interpolation kernel used by drizzle.
        """
        
        curdir = os.getcwd()
        self.logfile.write('Moving to directory '+self.obsFits+' for drizzling.')
        os.chdir(self.obsFits)       # move into the FITS dir.
        self.logfile.write('Using '+kernel+' drizzle kernel.')
        # will pick the output shape on-the-fly later if none specified
        self.outshape = outshape
        del outshape
        if bits==None:
            self.logfile.write('will drizzle without masks')
            if crmasks:
                warntxt='Warning: sent crmask dictionary but specified bits=None'
                self.errorList.append((self.modName,warntxt))
                self.logfile.write(warntxt)
        else:
            self.logfile.write('Will drizzle with masks; summed goodpix vals = '+str(bits))
            print '  Will drizzle with masks; summed goodpix vals =',bits
        # use drizzle pixfrac of 1 if not specified
        if not pixfrac:
            self.pixfrac = 1
            self.logfile.write('setting pixfrac = '+str(pixfrac)+' [DEFAULT VALUE]')  # by xingxing
        else:
            self.pixfrac = pixfrac
            self.logfile.write('setting pixfrac = '+str(pixfrac))
        self.asecpix = asecpix
        self.kernel  = kernel
        del pixfrac,asecpix,kernel

        # This bunch of thrashing is done to prevent multiple calls of this
        # method appending redundant or incorrect file names to the image
        # lists.  It is done this way so that the obs object passed to the
        # constructor also gets updated.  Just the setting lists to zero,
        # i.e.  self.sciImageList = [], and then appending items does not
        # result in the obs.sciImageList being updated as well.  
        while self.sciImageList:
            del self.sciImageList[0]
        while self.contextImageList:
            del self.contextImageList[0]
        while self.weightImageList:
            del self.weightImageList[0]

        # important to nuke any lingering inmasks!
        pyblot.junkem("*inmask*")
        self.parlists = []
        self.shortparlists = []
        self.wcslist  = []
        iraf.flpr()
        iraf.flpr()

        ### all this aug stuff is for making the mask for trimming the edges
        #pdb.set_trace()
        if separDriz:
            self.augimdict = {}
            for tab in self.obsAsnDict.keys():
                self.augimdict[tab] = {}
                self.augimdict[tab]['maskname'] = tab.split('.')[0][:-3]+'augmask.fits'
                self.augimdict[tab]['ctxlist'] = []

        #pdb.set_trace()
        for tab in self.obsAsnDict.keys():
                
            if self.outshape == None:
                # choose the reference Asn!
                self.refAsn = self._get_refAsn()
                REFtab = self.refAsn
                print "*** First step"
                # pdb.set_trace()  ##
                if self.asecpix: 
                    _skyf = xydrizzle.SkyField()
                    _skyf.set(psize=self.asecpix)
                    # if string.find(self.dfilts[0],'CLEAR' > -1 or string.find(self.dfilts[0],'CLEAR' > -1):
                    PyDrOb = xydrizzle.PyDrizzle(REFtab, bits=None, field=_skyf, pixfrac=self.pixfrac, kernel=self.kernel,
                            filter1=self.dfilts[0], filter2=self.dfilts[1],ref=self.ref) #WZ
                    del _skyf
                else:
                    PyDrOb = xydrizzle.PyDrizzle(REFtab, bits=None, pixfrac=self.pixfrac, kernel=self.kernel,
                            filter1=self.dfilts[0], filter2=self.dfilts[1],ref=self.ref) #WZ

                # OK, we've just created a PyDrizzle object to serve
                # as a prototype for later ones, but not be used itself...
                # pad the axis sizes a little 
                _maxAng = 0
                for imkey in self.MatchDict.keys():
                    delta45ang = 45 - abs(45 - abs(divmod(self.MatchDict[imkey]['angle'],90)[1]))
                    if delta45ang > _maxAng:
                        _maxAng = delta45ang
                if(_maxAng > 5):
                    _padx,_pady = self._rotateRect(_maxAng,PyDrOb.parlist[0]['outnx'],PyDrOb.parlist[0]['outny'])
                    if self.padfac:
                        _padx *= self.padfac
                        _pady *= self.padfac
                    self.logfile.write("padding x,y by factors "+str(_padx)+" "+str(_pady)+" including rotation.")
                else:
                    if self.padfac:
                        _padx = _pady = self.padfac
                    else:
                        _padx = _pady = 1.025
                # if self.outsize:
                #     _outNx,_outNy = self.outsize
                # else:
                #     _outNx = int(_padx * PyDrOb.parlist[0]['outnx'])
                #     _outNy = int(_pady * PyDrOb.parlist[0]['outny'])
                # self.logfile.write("Output image size will be %d x %d"%(_outNx,_outNy))
                # self.outshape  = (_outNx,_outNy)
                # del _maxAng,delta45ang,_padx,_pady
                #######################################
                self.protoWCS = PyDrOb.observation.product.geometry.wcs
                
                # check to see if scale of protoWCS is same as specified asecpix
                if self.asecpix:
                    if abs(self.protoWCS.pscale - self.asecpix) > 1e-9:
                        errtxt = 'Error: self.protoWCS.pscale: '+str(self.protoWCS.pscale)+\
                                 '  !=  self.asecpix: '+str(self.asecpix)
                        self.logfile.write(errtxt)
                        raise Exception,errtxt
                else:
                    self.asecpix = self.protoWCS.pscale
                #if self.outsize:
                #    _outNx = int(self.outsize[0] / self.asecpix + 0.5 + 2*separDriz)
                #    _outNy = int(self.outsize[1] / self.asecpix + 0.5 + 2*separDriz)
                #else:
                #    _outNx = int(_padx * PyDrOb.parlist[0]['outnx'])
                #    _outNy = int(_pady * PyDrOb.parlist[0]['outny'])
                _outNx=self.ref.prod.naxis1
                _outNy=self.ref.prod.naxis2

                self.logfile.write("Output image size will be %d x %d"%(_outNx,_outNy))
                self.outshape  = (_outNx,_outNy)
                del _maxAng,delta45ang,_padx,_pady
                # self.ra_coord  = PyDrOb.observation.product.geometry.wcs.crval1
                # self.dec_coord = PyDrOb.observation.product.geometry.wcs.crval2
                # self.crpix1  = PyDrOb.observation.product.geometry.wcs.crpix1
                # self.crpix2  = PyDrOb.observation.product.geometry.wcs.crpix2
                ## Set what we can using the SkyField object.
                print "*** Second step"
                # pdb.set_trace()
                _skyf = xydrizzle.SkyField()
                _skyf.set(shape=self.outshape, psize=self.asecpix, \
                          ra=self.protoWCS.crval1, dec=self.protoWCS.crval2) 

                # now, remake the PyDrOb object, just to be sure of the right
                # zeropoint shifts for an output image of this size (this SkyField)
                del PyDrOb, self.protoWCS
                PyDrOb = xydrizzle.PyDrizzle(REFtab, bits=None, field=_skyf, pixfrac=self.pixfrac, kernel=self.kernel,
                         filter1=self.dfilts[0], filter2=self.dfilts[1],ref=self.ref) #WZ
                self.protoWCS = PyDrOb.observation.product.geometry.wcs
                # paranoid check
                if self.asecpix:
                    if abs(self.protoWCS.pscale - self.asecpix) > 1e-9:
                        errtxt = 'Error: self.protoWCS.pscale: '+str(self.protoWCS.pscale)+\
                             '  !=  self.asecpix: '+str(self.asecpix)
                        self.logfile.write(errtxt)
                        raise Exception,errtxt
            
                # calculate zeropoints of the x,y shifts, using *all* images in MatchDict
                _dX_EMP = []
                _dY_EMP = []
                ###for imdict in PyDrOb.parlist:
                #pdb.set_trace() ##
                for i in range(len(PyDrOb.parlist)): #WZ reference rotation Dec 2010
                    if (string.find(PyDrOb.parlist[i]['data'],self.ref.rootname) > -1):
                        self.ref.rot = PyDrOb.parlist[i]['rot'] #WZ reference rotation Dec 2010
                        print "Reference rotation angle: ", PyDrOb.parlist[i]['rot']
                for imkey in self.MatchDict.keys():
                    # pdb.set_trace()
                    # deg = (self.ref.rot - self.MatchDict[imkey]['angle']) #WZ >>>
                    deg = self.ref.rot #WZ >>>
                    theta = -deg * numpy.pi /180. # Note the sign
                    x= self.MatchDict[imkey]['xarc_shift'] / self.asecpix
                    y= self.MatchDict[imkey]['yarc_shift'] / self.asecpix     
                    x1 =  x * numpy.cos (theta) + y * numpy.sin(theta)
                    y1 = -x * numpy.sin(theta)  + y * numpy.cos(theta)
                    x1,y1 = self.xyrot([x,y],self.ref.rot)
                    x_MsdShift = x1 
                    y_MsdShift = y1 
                    #x_MsdShift = self.MatchDict[imkey]['xarc_shift'] / self.asecpix
                    #y_MsdShift = self.MatchDict[imkey]['yarc_shift'] / self.asecpix
                    _dX_EMP.append(x_MsdShift)
                    _dY_EMP.append(y_MsdShift)
                del x_MsdShift, y_MsdShift, deg,theta,x,y,x1,y1
                self.logfile.write('dx_EMP '+ str(_dX_EMP))
                self.logfile.write('dy_EMP '+ str(_dY_EMP))

                ## this centers the refimage in the output image:
                #substring=string.split(self.root,'datasets/')
                #refimg = substring[0]+'ingest/'+string.split(substring[1],'/')[0]+'/ref.fits'
                # refimg= string.replace(self.root,'datasets','ingest')+'/ref.fits'
                refimg = os.environ['INGEST']+self.obsName+'/ref.fits'  # by xingxing
		im = refimg + '[0]'
                ref_wcs = wcsutil.WCSObject(im)
                print "Tweak values"
                #pdb.set_trace()
                #reffile= string.replace(self.root,'datasets','ingest')+'/ref.txt'<
                twf= 'tweak.log'
                if os.path.exists(twf): 
                    rdata=tableio.get_str(twf,(0,1,2,3))
                    j=len(rdata[0])-1
                    # im=rdata[0][i]+'.fits[0]'
                    # ref_wcs = wcsutil.WCSObject(im)
                    ref_wcs.dx=float(rdata[1][j])                
                    ref_wcs.dy=float(rdata[2][j])
                    ref_wcs.dr=float(rdata[3][j])/numpy.pi*180. #WZ degree 
                    del rdata,j
                else:   # this is the first pass. Use given ref image
                    rfile= string.replace(self.root,'datasets','ingest')+'/tweak.log'
                    if os.path.isfile(rfile): # Initial guess WZ
                        rdata=tableio.get_str(rfile,(0,1,2,3)) # Only the first line
                        ref_wcs.dx=float(rdata[1][0])
                        ref_wcs.dy=float(rdata[2][0])
                        ref_wcs.dr=float(rdata[3][0])/numpy.pi*180. #WZ degree 
                        del rdata,rfile
                    else:
                        ref_wcs.dx= 0.
                        ref_wcs.dy= 0.
                        ref_wcs.dr= 0.
                self.align=ref_wcs.copy()
                del ref_wcs, twf
                #self.align.xSHT = 0. # -26.22 #WZ
                #self.align.ySHT = 0. # -56.09
                if self.outshift:
                    self.xZPT = 1.0*int(0.5 + (self.outshift[0] / self.asecpix))
                    self.yZPT = 1.0*int(0.5 + (self.outshift[1] / self.asecpix))
                    self.logfile.write('using specified arcsec shift zeropoints '+\
                                       str(self.outshift[0])+' '+str(self.outshift[1]))
                else:
                    # this centers the avg image position:
                    self.xZPT = self.align.dx -1.0*int((max(_dX_EMP) + min(_dX_EMP))/2.0)   # preserves fractional
                    self.yZPT = self.align.dy -1.0*int((max(_dY_EMP) + min(_dY_EMP))/2.0)   # part of image shifts
                self.logfile.write('selected pixshift zeropoints '+str(self.xZPT)+' '+str(self.yZPT))
                del PyDrOb, _dX_EMP, _dY_EMP, REFtab

            _skyf.exptime = None
            print "*** Third step"
            # pdb.set_trace()
            PyDrOb = xydrizzle.PyDrizzle(tab, field=_skyf, clean=deltmp, bits=bits, units=units,\
                      pixfrac=self.pixfrac, kernel=self.kernel,
                      filter1=self.dfilts[0], filter2=self.dfilts[1], ref=self.ref) #WZ
            # verify all subsequent images get same pscale
            if abs(PyDrOb.observation.product.geometry.wcs.pscale - self.asecpix) > 1e-9:
                errtxt='Error: '+tab+' PyDrOb wcs.pscale: '+str(PyDrOb.observation.product.geometry.wcs.pscale)\
                        + '  !=  self.asecpix: '+str(self.asecpix)
                self.logfile.write(errtxt)
                raise Exception,errtxt

            # Now mask out specific regions if specified from the command line.
            # Usually, this is to explicitly remove figure 8's or other such
            # artifacts from the F850LP exposures...  

            if self.maskFile:
                maskreg_cmd = 'maskreg %s %s' % (self.maskFile,self.obsFits)
                self.logfile.write(maskreg_cmd)
                os.system(maskreg_cmd)

            # OK, PyDrOb should now have the desired shape and CRVAL's, and I'd have
            # thought we have to make sure the Drizzle Product has a consistent CRPIX
            # for this modified shifts, but drizzle seems to take care of it.
            # Now have to update the individual shifts accordingly...
            # Also, we need to change output product names if separDriz=1
            # and include crmasks if they've been sent.
            if separDriz and (self.MatchDict[self.MatchDict.keys()[0]]['NumSci'] == 2):
                if self.smartstack:
                    output = open('medianimages_input','w')
                    S = '%i\n' % (len(PyDrOb.parlist))
                    output.write(S)        
                    for imdict in PyDrOb.parlist:
                        imkey = string.split(imdict['data'],'[')[0]
                        x_MsdShift = self.MatchDict[imkey]['xarc_shift'] / self.asecpix
                        y_MsdShift = self.MatchDict[imkey]['yarc_shift'] / self.asecpix
                        xsh = (x_MsdShift + self.xZPT + _outNx/2)*self.asecpix
                        ysh = (y_MsdShift + self.yZPT + _outNy/2)*self.asecpix
                        S = '%s %g %g %g\n' % (imdict['data'],xsh,ysh,-(self.MatchDict[imkey]['angle']-self.align.dr))
                        output.write(S)
                    output.close()
                    medianfilter_cmd = "whmedian medianimages_input medianimages %g %g 10" % (int(float(_outNx)*self.asecpix+0.9),int(float(_outNy)*self.asecpix+0.9))
                    self.logfile.write(medianfilter_cmd)
                    os.system(medianfilter_cmd)
                else:
                    miout = open('medianimages', 'w')
                    for imdict in PyDrOb.parlist:
                        miout.write('%s\n' % imdict['data'])
                    miout.close()

            print "*** Fourth step"
            # pdb.set_trace()
            for i in range(len(PyDrOb.parlist)): #WZ reference rotation Mar 2011
                if (string.find(PyDrOb.parlist[i]['data'],self.ref.rootname) > -1):
                    self.ref.rot = PyDrOb.parlist[i]['rot'] #WZ reference rotation Dec 2010
                    print "Reference rotation angle: ", PyDrOb.parlist[i]['rot']
            # self.ref.rot = PyDrOb.parlist[0]['rot'] #WZ reference rotation Dec 2010   Key step
            for imdict in PyDrOb.parlist:
                imkey = string.split(imdict['data'],'[')[0]
                x_MsdShift = self.MatchDict[imkey]['xarc_shift'] / self.asecpix #WZ Feb 2011
                y_MsdShift = self.MatchDict[imkey]['yarc_shift'] / self.asecpix
                # x_MsdShift = self.MatchDict[imkey]['xarc_shift'] / self.ref.model.pscale
                # y_MsdShift = self.MatchDict[imkey]['yarc_shift'] / self.ref.model.pscale

                x1,y1=self.xyrot([x_MsdShift,y_MsdShift],self.ref.rot)
                imdict['xsh'] = x1 + self.xZPT # x_MsdShift + self.xZPT
                imdict['ysh'] = y1 + self.yZPT # y_MsdShift + self.yZPT
                ### changing to negative sign for drizzle convention (12/Apr/2002)
                # imdict['rot'] = -(self.MatchDict[imkey]['angle'])
                # pdb.set_trace()
                imdict['rot'] = self.ref.rot - (self.MatchDict[imkey]['angle'])-self.align.dr #WZ
                tmpstring = string.split(string.split(imdict['data'],'[')[1],']')[0]
                Xname = (string.split(tmpstring,',')[0]).upper()
                Xver  = int(string.split(tmpstring,',')[1])
                basefits   = string.split(imkey,'.')[0]+"_"+Xname+"_"+str(Xver)+'.fits'

                if separDriz:
                    imdict['outdata']    = '_dz_'+basefits
                    imdict['outcontext'] = '_cx_'+basefits
                    imdict['outweight']  = '_wt_'+basefits
                    # self.augimdict[tab]['ctxlist'].append(imdict['outcontext']) #WZ

                    # have to make sure exptime is right for a single image
                    _oldtexp = str(imdict['texptime'])
                    imdict['texptime']   = fUtil.getKeyVal(imdict['data'].split('[')[0],'EXPTIME')
                    self.logfile.write('changing '+imdict['outdata']+' parlist texptime from '+\
                                       _oldtexp+' to '+str(imdict['texptime']))
                    del _oldtexp
                    self.removeList.append(imdict['outdata'])
                    self.removeList.append(imdict['outcontext'])
                    # 'outweight' now done separately...
                    # self.removeList.append(imdict['outweight']) 
                    # Again, be sure to remove the files if they're already there!
                    if os.path.isfile(imdict['outdata']):
                        os.remove(imdict['outdata'])
                    if os.path.isfile(imdict['outcontext']):
                        os.remove(imdict['outcontext'])
                    if os.path.isfile(imdict['outweight']):
                        os.remove(imdict['outweight'])

                if crmasks:
                    if bits:
                        # then multiply the pixel in_masks by crmasks
                        iraf.unlearn(iraf.imcalc)
                        iraf.imcalc.pixtype = 'short'
                        imcalclist = imdict['in_mask']+','+crmasks[basefits][0]
                        print imdict['in_mask'],' = ',imdict['in_mask']+' * '+crmasks[basefits][0]
                        iraf.imcalc(imcalclist,"_tmpfile.fits","im1 * im2")

                        orig_mask = 'Orig_'+imdict['in_mask']
                        os.rename(imdict['in_mask'], orig_mask)
                        os.rename("_tmpfile.fits", imdict['in_mask'])
                        self.removeList.append(orig_mask)

                        nCRzap = crmasks[basefits][1] #WZ
                        self.logfile.write(str(nCRzap)+" cosmic ray pixels masked by pyblot in "+basefits)
                        del imcalclist, orig_mask
                        iraf.flpr('imcalc')
                        iraf.flpr()
                        self.makeEdgemask(imdict['in_mask']) #WZ
                        # self.makeIvmImage(imdict) #WZ
                    else:
                        # then we just copy crmask names to imdict
                        if imdict['in_mask']:  # this error should never happen
                            self.errorList.append((self.modName,"ERROR: bits=None, but in_mask exists!"))
                            raise Exception,"ERROR: bits=None, but in_mask exists!"
                        imdict['in_mask'] = crmasks[basefits][0]
                        self.logfile.write('Setting '+imdict['data']+' pyrizzle in_mask to '+imdict['in_mask'])
                        print 'Setting',imdict['data'],'pyrizzle in_mask to',imdict['in_mask']
                        nCRzap = crmasks[basefits][1]
                        self.logfile.write(str(nCRzap)+" cosmic ray pixels masked by pyblot in "+basefits)
                        # pdb.set_trace()
                        # self.makeIvmImage(imdict) #WZ
                del imkey,tmpstring,Xname,Xver,basefits

            parlist = PyDrOb.parlist[:]
            # print "*** Recheck parlist"
            self.parlists.append(parlist)
            self.outname=imdict['outdata']
            self.outcontext=imdict['outcontext']
            self.outweight=imdict['outweight']
            del imdict    

            if not separDriz:
                if self.noContext:
                    for junkdict in PyDrOb.parlist:
                        junkdict['outcontext'] = ''
                    del junkdict

                # then the first dict in list tells the whole data output story
                imdict = PyDrOb.parlist[0]
                self.logfile.write("All set for "+imdict['outdata']+' texptime: '+str(imdict['texptime']))

            if separDriz and (self.MatchDict[self.MatchDict.keys()[0]]['NumSci'] == 2):
                input = open('medianimages')
                L = input.readlines()
                i = 0
                for l in L:
                    L[i] = l[:-1]
                    i = i + 1
                print L
                input.close()
                RemoveList = []
                for imdict in PyDrOb.parlist:
                    print imdict['data'] 
                    if not (imdict['data'] in L):
                        RemoveList.append(imdict)
                        # pdb.set_trace()
                    for imdict in RemoveList:
                        PyDrOb.parlist.remove(imdict)
                        print "Removing... %s" % (imdict['data'])

            if separDriz:
                self.shortparlists.append(PyDrOb.parlist)

            print "Beginning pydrizzle run method . . ."
            self.logfile.write('output pscale:  '+str(round(1e8*self.asecpix)/1e8))
            print "*** Fifth step: Actual build"
            # pdb.set_trace() # _dr etc WZ
            fname=os.path.join(self.obsAlign,'alignTweak.txt')
            if (os.path.exists(fname)):
                align_name=tableio.get_str(fname,0)
                align_data=tableio.get_data(fname,(1,2,3))
                for i in range(len(PyDrOb.parlist)): # WZ Apply tweaks                  
                    for j in range(len(align_name)):
                        # rname=string.split(align_name[j],'_')[2]                      
                        if (str.find(self.parlists[0][i]['data'],align_name[j])>-1):
                            self.parlists[0][i]['xsh']=self.parlists[0][i]['xsh']+align_data[0][j]/self.parlists[0][i]['scale']
                            self.parlists[0][i]['ysh']=self.parlists[0][i]['ysh']+align_data[1][j]/self.parlists[0][i]['scale']
                            self.parlists[0][i]['rot']=self.parlists[0][i]['rot']+align_data[2][j]/self.parlists[0][i]['scale']
            # pdb.set_trace()
            PyDrOb.run(save=1, build=0)   # never 'build' ################IMPORT STEP

            if separDriz: #WZ 
                fname=os.path.join(self.obsAlign,'alignTweak.txt')
                if (not os.path.exists(fname)):
                    print "*** Sixth step: Tweak individual files"
                    parlist=[]
                    groupList=[]
                    for i in range(len(PyDrOb.parlist)): #WZ Make a list of PA
                        fname= string.split(PyDrOb.parlist[i]['outdata'],'_dz_')[1]
                        fitsfile=pyfits.open(fname)
                        v = fitsfile[0].header.get('PA_V3')
                        det = fitsfile[0].header.get('DETECTOR')
                        fitsfile.close()
                        pa = int(round(v))
                        PyDrOb.parlist[i]['pa']=pa
                        PyDrOb.parlist[i]['detector']=det
                        #pdb.set_trace()                                                                                                        
                        group=PyDrOb.parlist[i]['data'][0:6]
                        PyDrOb.parlist[i]['group']=group
                        if i==0:
                            groupList.append(group)
                        else:
                            found=0
                            for j in range(len(groupList)):
                                if (string.find(groupList[j],group)>-1):
                                # if (numpy.abs(pa-paList[j])<1): 
                                    found=1
                            if (found==0):
                                groupList.append(group)
                        parlist.append(fname)
                    for i in range(len(groupList)): # WZ Make ID by Group
                        inputString=''
                        exptime=0.  # WZ Nov 2013
                        for j in range(len(PyDrOb.parlist)):
                            if (PyDrOb.parlist[j]['group']==groupList[i]):
                                fname= PyDrOb.parlist[j]['outdata']
                                ffi      = pyfits.open(fname)
                                expt   = ffi[0].header['EXPTIME']
                                ffi.close()
                                exptime=exptime + expt
                                if len(inputString)==0:
                                    inputString=fname
                                else:
                                    inputString = inputString +','+fname
                        imcombine_input=string.split(inputString,',')
                        NumIm = len(imcombine_input)
                        # iraf.imcombine.input="@imcombine_input"
                        # iraf.imcombine.input=inputString
                        # Run imcombine
                        #iraf.imcombine.output = 'group_'+str(groupList[i])+'.fits'
                        #iraf.imcombine.sigma = ''
                        #iraf.imcombine.combine = 'median'
                        #iraf.imcombine.reject = 'minmax'
                        #iraf.imcombine.mode    = 'h'  
                        if NumIm == 1:   #1,2                                                                                         
                            nlow = 0
                            nhigh = 0
                        elif NumIm == 2:      # 2                                         
                            nlow = 0
                            nhigh = 1
                        elif NumIm == 3:  # 3                                             
                            nlow = 0
                            nhigh = 2
                        elif NumIm == 4:  # 4                                             
                            nlow = 0
                            nhigh = 3
                        elif NumIm == 5:  # 5                                             
                            nlow = 1
                            nhigh = 3
                        elif NumIm == 6:  # 6                                             
                            nlow = 1
                            nhigh = 4
                        else:
                            nlow = (NumIm+2)/4 - 1
                            nhigh = 3*NumIm/4
                        #iraf.imcombine()                                                 
                        fname='group_'+groupList[i]+'.inp'
                        inputfile = open(fname,'w')
                        print >> inputfile, NumIm, nlow, nhigh, 0
                        for j in range(len(imcombine_input)):
                            cx=string.replace(imcombine_input[j],'_dz','_cx')
                            print >> inputfile, imcombine_input[j],cx
                        inputfile.close()
                        outfile='group_'+groupList[i]+'.fits'
                        cmd="medianfilter "+fname+' '+outfile
                        # pdb.set_trace() 
                        self.logfile.write(cmd)
                        os.system(cmd)
                        fitsfile = pyfits.open(outfile, mode='update') # WZ Nov 2013
                        for ext in fitsfile:
                            ext.header.update('BUNIT', 'CPS')
                            # Update the data                                 
                            ext.data = numpy.divide(ext.data, exptime)
                        if self:
                            self.logfile.write('Converted file ' + outfile + ' from cts to cps')
                        fitsfile.close()

                        # Now make a map image                                            
                        #iraf.imcombine.reject = 'none'                                   
                        #iraf.imcombine.combine = 'sum'                                   
                        #mapString=string.replace(inputString,'_dz','_cx')                
                        #iraf.imcombine.input=mapString                                   
                        #iraf.imcombine.output = 'group_'+groupList[i]+'_map.fits'        
                        #iraf.imcombine()                                                 
                    self.alignTweak(PyDrOb.parlist,groupList)

            # delete 1st pass wght files, unless told otherwise
            if delwght and separDriz:
                for __tmpdict in PyDrOb.parlist:
                    if os.path.isfile(__tmpdict['outweight']):
                        os.remove(__tmpdict['outweight'])
                del __tmpdict

            self.wcslist.append(PyDrOb.observation.product.geometry.wcs)
            # drizzling all done.
            # Note: the WCS's of the the drizzle products change when the parlist shifts
            # change.  In particular, the CRPIX's of the images produced by drizzle 
            # with each other, and not with the WCS's that get appended to wcslist!
            # This is because drizzle sets the output image WCS, not multidrizzle.

            # Do some final fixes to drizzle products
            # if self.NumSci > 1 and float(multidrizzle.__version__[:3]) < 1.4:
            # break these conditions down....

            if not separDriz:
                # DANGER! this is a kludge to get the right count levels & ET's in
                # the simple drizzle products
                if self.NumSci > 1:
                    if units == 'counts':
                        # then drizzle doesn't get count levels right; divide by NumSci
                        _wtxt = "WARNING: dividing drizzled count levels by NumSci = "+str(self.NumSci)
                        self.logfile.write(_wtxt)
                        # self.errorList.append((self.modName,_wtxt)) #WZ
                        # divide images by NumSci
                        iraf.unlearn(iraf.imcalc)
                        # pdb.set_trace()
                        operation = "im1 / "+str(float(self.NumSci))
                        iraf.imcalc(imdict['outdata'], imdict['outdata'], operation)
                        # only the 'outdata' should be divided by NumSci -- 
                        # not 'outweight' (exposure time image) or 'outcontext'
                        del operation,_wtxt
                        iraf.flpr('imcalc')
                        iraf.flpr()

                    if float(pydriz_version[:3]) < 1.4:
                        # then drizzle gets all exposure times wrong; divide by NumSci
                        _wtxt = "WARNING: dividing drizzled exptimes by NumSci = "+str(self.NumSci)
                        self.logfile.write(_wtxt)
                        # self.errorList.append((self.modName,_wtxt)) #WZ
                        expTimeWrong = fUtil.getKeyVal(imdict['outdata'],'EXPTIME')
                        expTimeVal = [('EXPTIME',expTimeWrong/float(self.NumSci))]
                        fUtil.fixHeader(os.path.join(self.obsFits,imdict['outdata']),expTimeVal)
                        fUtil.fixHeader(os.path.join(self.obsFits,imdict['outweight']),expTimeVal)
                        if not self.noContext:
                            fUtil.fixHeader(os.path.join(self.obsFits,imdict['outcontext']),expTimeVal)
                        del expTimeWrong, expTimeVal, _wtxt

                    else:
                        # then drizzle gets ctx and wgt exposure times wrong...
                        _wtxt = "WARNING: setting context and weight exptimes to same as sci image."
                        self.logfile.write(_wtxt)
                        #self.errorList.append((self.modName,_wtxt)) #WZ
                        expTimeVal = [('EXPTIME',fUtil.getKeyVal(imdict['outdata'],'EXPTIME'))]
                        # not sci image! ctx and wgt only!
                        fUtil.fixHeader(os.path.join(self.obsFits,imdict['outweight']),expTimeVal)
                        if not self.noContext:
                            fUtil.fixHeader(os.path.join(self.obsFits,imdict['outcontext']),expTimeVal)
                        del expTimeVal, _wtxt

                # Sum subtracted sky's and write result to drizzled image header
                # if the there is a single driz im, and self.alignSky non-Null
                if self.alignSky:
                    totNcomb,totSkySub = self._sumSubSkyNcombine(PyDrOb.parlist,doSky=1)
                    fUtil.fixHeader(os.path.join(os.getcwd(),PyDrOb.parlist[0]['outdata']),\
                                    [(self.alignSky,totSkySub)])
                    self.logfile.write(self.alignSky+" = "+str(totSkySub)+" set in "+\
                                       PyDrOb.parlist[0]['outdata']+" header.")
                else:
                    totNcomb,totSkySub = self._sumSubSkyNcombine(PyDrOb.parlist,doSky=0)
                # either way, we now have totNcomb
                fUtil.fixHeader(os.path.join(os.getcwd(),PyDrOb.parlist[0]['outdata']),\
                                [("NCOMBINE",totNcomb)])
                self.logfile.write("NCOMBINE = "+str(totNcomb)+" set in "+\
                                   PyDrOb.parlist[0]['outdata']+" header.")
                del totNcomb,totSkySub
                               
            del PyDrOb

        # OK, loop over "asn" done!  All drizzling finished!
        # list the output for the mkMsg method
	#print "Please start to check the header!!__________by xingxing"
        #pdb.set_trace()
        if not separDriz:
            for list in self.parlists:
                pred_images = []
                cur_im  = ''
                next_im = ''
                sci_im = list[0]['outdata']
                if not self.noContext:
                    cxt_im = list[0]['outcontext']
                    self.contextImageList.append(cxt_im)
                wgt_im = list[0]['outweight']
                self.sciImageList.append(sci_im)
                self.weightImageList.append(wgt_im)

            for dict in list:
                next_im = dict['data'][:-7]
                if next_im != cur_im:
                    cur_im = next_im
                    pred_images.append(cur_im)
            self.outputList[sci_im] = pred_images
            if not self.noContext:
                self.outputList[cxt_im] = pred_images
            self.outputList[wgt_im] = pred_images

            # assemble all the WCS information into a list of tuples
            self.logfile.write("Making consistent WCS Info List . . .")
            #pdb.set_trace()
            wcsInfoList = self._makeWCS(template=self.sciImageList[0])
            #wcsInfoList[4][1]=self.ref.cd11 # quick fix WZ
            #wcsInfoList[5][1]=self.ref.cd12 # Keep the WCS
            #wcsInfoList[6][1]=self.ref.cd21 
            #wcsInfoList[7][1]=self.ref.cd22 
            #wcsInfoList = self._makeWCS(template=None) #WZ
            # Now we want to fix up the headers of the final drizzle products.
            # This is done by the fUtil ufunc below.
            mod_string = self.modName + ", v"+ __version__
            date = ptime()  # ptime from pUtil
            for im in self.sciImageList:
                # with all the problems with pyfits/numpy
                # we will only try to fix the header of the file
                # warn and pass if exception
                x,y = self.align.rd2xy([self.align.crval1,self.align.crval2])
                ra,dec=self.align.xy2rd([x+self.align.dx,y+self.align.dy])
                xVal =  [('CRPIX1',self.align.crpix1)]
                fUtil.fixHeader(im,xVal)
                yVal =  [('CRPIX2',self.align.crpix2)]
                fUtil.fixHeader(im,yVal)
                raVal = [('CRVAL1',self.align.crval1)]
                fUtil.fixHeader(im,raVal)
                decVal =[('CRVAL2',self.align.crval2)]
                fUtil.fixHeader(im,decVal)

                try:    fUtil.fixDrzHeader2(im,mod_string,date,addCards=None) #WZ wcsInfoList)
                except Exception, err: 
                    self.logfile.write("ERROR: Call to fixDrzHeader failed on science image, "+im+\
                                       "\n\t\t DANGER!  WCS probably WRONG!  DO NOT ARCHIVE!  DANGER!\n")
                    self.logfile.write(" . . . is this a new version of pyfits . . . ?")
                    self.errorList.append((self.modName,"Call to fixDrzHeader failed on science image, "+im))
                    self.errorList.append((self.modName,"WCS of "+im+" likely WRONG!  DO NOT ARCHIVE!"))
                    self.errorList.append((self.modName,str(err)))
                    continue
             
                
            for im in self.contextImageList:
                x,y = self.align.rd2xy([self.align.crval1,self.align.crval2])
                ra,dec=self.align.xy2rd([x+self.align.dx,y+self.align.dy])
                xVal =  [('CRPIX1',self.align.crpix1)]
                fUtil.fixHeader(im,xVal)
                yVal =  [('CRPIX2',self.align.crpix2)]
                fUtil.fixHeader(im,yVal)
                raVal = [('CRVAL1',self.align.crval1)]
                fUtil.fixHeader(im,raVal)
                decVal =[('CRVAL2',self.align.crval2)]
                fUtil.fixHeader(im,decVal)

                try: 
                    fUtil.fixDrzHeader2(im,mod_string,date,file_type="CTX",addCards=None) # wcsInfoList)
                except Exception, err: 
                    self.logfile.write("ERROR: Call to fixDrzHeader failed on context image, "+im)
                    self.logfile.write("Probabaly another flakey version of pyfits....")
                    self.errorList.append((self.modName,"Call to fixDrzHeader failed on context image, "+im))
                    self.errorList.append((self.modName,str(err)))
                    continue

            for im in self.weightImageList:
                x,y = self.align.rd2xy([self.align.crval1,self.align.crval2])
                ra,dec=self.align.xy2rd([x+self.align.dx,y+self.align.dy])
                xVal =  [('CRPIX1',self.align.crpix1)]
                fUtil.fixHeader(im,xVal)
                yVal =  [('CRPIX2',self.align.crpix2)]
                fUtil.fixHeader(im,yVal)
                raVal = [('CRVAL1',self.align.crval1)]
                fUtil.fixHeader(im,raVal)
                decVal =[('CRVAL2',self.align.crval2)]
                fUtil.fixHeader(im,decVal)
                self.fixExptime(im) #WZ    

                try:
                    fUtil.fixDrzHeader2(im,mod_string,date,file_type="WGT",addCards=None) # wcsInfoList)
                except Exception, err: 
                    self.logfile.write("ERROR: Call to fixDrzHeader failed on weight image, "+im)
                    self.logfile.write("Probabaly another flakey version of pyfits....")
                    self.errorList.append((self.modName,"Call to fixDrzHeader failed on weight image, "+im))
                    self.errorList.append((self.modName,str(err)))
                    continue

            # get the medriz_?_.fits images into the outputList. not archiveable.
            medriz_list = glob.glob("medriz_*.fits")
            if medriz_list:
                for file in medriz_list:
                    self.outputList[file] = []   # can't determine predecessors at this point but who cares?

        os.chdir(curdir)
        self.logfile.write('returning to directory '+curdir)
        return

    def fixExptime(self,image=None):
        """
        W Zheng Aug 2011.
        For UVIS and ACS images, scale exposure time in the weight image (by 0.5)
        """
        wgtfits = pyfits.open(image,'update')
        max=numpy.max(wgtfits[0].data)
        exptime = wgtfits[0].header.get('EXPTIME')
        # if (max > exptime):
        scale=exptime/max
        wgtfits[0].data =  wgtfits[0].data * scale
        wgtfits.close() # flush
        self.logfile.write('Scale the exposure by factor of '+str(scale) )
        return                                            

    def makeIvmImage(self,dic):
        """
            W Zheng Aug 2011
        """
        # pdb.set_trace()
        self.logfile.write("starting make Ivm Image: " + dic['data'])

        # reset rms image list
        #while self.rmsImageList:
        #    del self.rmsImageList[0]

        curdir = os.getcwd()
        #os.chdir(self.obsFits)
        try:
            im_crmask=dic['mask']
            im_rms=string.replace(dic['mask'],'inmask','ERR_')
            im_wgt=string.replace(dic['mask'],'inmask','wt')
            # dic['in_mask']=im_wgt #WZ still use the original definition
            maskfits = pyfits.open(im_crmask)
            flag = maskfits[0].data
            rmsfits = pyfits.open(im_rms)
            mean = numpy.mean(rmsfits[0].data)
            # numpy.max(rmsfits[0].data)
            # numpy.min(rmsfits[0].data)

            # sn0=10
            # sn = rmsfits[0].data/mean #WZ Mark bright pixels
            # idx = numpy.where(numpy.logical_and(numpy.greater(sn,sn0))
            # rmsfits[0].data[idx] = mean # numpy.abs(rmsfits[0].data[idx])

            wgtfits = pyfits.HDUList()
            wgtfits.append(pyfits.PrimaryHDU())
            wgtfits[0].header = rmsfits[0].header.copy()
            wgtfits[0].data   = maskfits[0].data / rmsfits[0].data / rmsfits[0].data

            numpy.mean(wgtfits[0].data)
            numpy.max(wgtfits[0].data)
            numpy.min(wgtfits[0].data)
            if os.path.isfile(im_wgt):
                os.remove(im_wgt)
            wgtfits.writeto(im_wgt)
            del wgtfits
        except:
            self.errorList.append((self.modName,"Cannot make a FITS object out of file "+im_wgt))
            raise Exception,"Cannot make a FITS object out of file "+im_wgt
            if len(wgtfits) > 1 or len(scifits) > 1:
                self.errorList.append((self.modName,"image file is not simple fits."))
                raise Exception,"image file is not simple fits."

            # build rms image name and open as a new file.
            self.rmsImageList.append(rmsfile)
            self.outputList[rmsfile] = [im_wgt]
            
            # make new fits obj and copy WGT/SCI hdr/data to RMS image initially
            try:
                del rmsfitsobj[0].header.ascard["EXTEND"]
            except KeyError:
                pass

            # reopen the rms image for editing.
            rmsfits = pyfits.open(rmsfile,'update')

            # ratio of default to specified output scales
            area_ratio = (self.asecpix / self.origscale)**2
            if abs(1-area_ratio) < 1e-4: area_ratio = 1
            self.logfile.write('Using area_ratio = %.6f in makeRmsImage' %(area_ratio))


            del scifits, wgtfits, im_wgt, im_sci, readVariance, totInstVar, area_ratio, sn, idx
        
            ## now fix up problem values...
            newDat = numpy.where(numpy.logical_or(numpy.greater_equal(newDat,1e38),\
                                                      numpy.less_equal(newDat,0.)),4e38,newDat) #WZ Mar 2013
            rmsfits[0].data = numpy.sqrt(newDat).astype(numpy.float32)

            # a few token updates to the header, then write it out
            rmsfits[0].header.update('FILENAME',rmsfile)
            rmsfits[0].header.update('FILETYPE','RMS')
            rmsfits.close()
            self.logfile.write('Made rms image '+rmsfile)
            del newDat, rmsfile, rmsfits
        rmsfits.close()
        maskfits.close()
        # os.chdir(curdir)
        return

    def makeEdgemask(self,im_mask): #WZ
        """
        Check the read with "width" pixels to the border. If S/N > threshold,
        make zerio in the crmask image.
            W Zheng Aug 2011
        """
        width = 15 # 20 25 10
        sn0= 7.5 # 10
        npt=0L
        # im_mask=self.parlists[0][i]['in_mask']
        # maskfits = pyfits.open(im_mask,mode='update')
        maskfits = pyfits.open(im_mask)
        im_sci=string.replace(im_mask,'inmask','SCI_')
        scifits = pyfits.open(im_sci)
        data = scifits[0].data.copy()
        naxis1 = scifits[0].header.get('NAXIS1')
        naxis2 = scifits[0].header.get('NAXIS2')
        sky = scifits[0].header.get('ALIGNSKY')
        mask = maskfits[0].data.copy()
        for j in range(0,width-1): # y
            for k in range(0,naxis1-1): # x
                if (numpy.abs(data[j,k]/sky) > sn0 and mask[j,k]==1):
                    # print j,k
                    mask[j,k]=0
                    npt = npt + 1
        #print mask[10,1000],' 10,1000'
        #print npt
        #npt=0
        for j in range(0,naxis2-1):
            for k in range(0,width-1):
                if (numpy.abs(data[j,k]/sky) > sn0 and mask[j,k]==1):
                    # print j,k
                    mask[j,k]=0
                    npt = npt + 1
        #print npt
        #print mask[1000,10],' 100,10'
        #npt=0
        for j in range(naxis2-width-1,naxis2-1):
            for k in range(0,naxis1-1):
                if (numpy.abs(data[j,k]/sky) > sn0 and mask[j,k]==1):
                    # print j,k
                    mask[j,k]=0
                    npt = npt + 1
        #print mask[2040,1000], ' 2040,1000'
        #print npt
        #npt=0
        for j in range(0,naxis2-1):
            for k in range(naxis1-width-1,naxis1-1):
                if (numpy.abs(data[j,k]/sky) > sn0 and mask[j,k]==1):
                    # print j,k
                    mask[j,k]=0
                    npt = npt + 1
        #print mask[1000,4090] ,' 1000,4090'
        #print npt
        maskfits[0].data = mask.copy()
        self.logfile.write(str(npt) + " pixels masked near the edges in image: " + im_mask)
        newfits = pyfits.HDUList()
        newfits.append(pyfits.PrimaryHDU())
        newfits[0].header = maskfits[0].header
        newfits[0].data   = mask.copy()
        # pdb.set_trace()
        scifits.close()            
        if os.path.isfile(im_mask):
            os.remove(im_mask)
        newfits.writeto(im_mask)
        # maskfits.flush()
        del npt,scifits,maskfits,newfits
        return


    def makeFlagImage(self):
        """turns the weight images produced by drizzle into flag images that
            SExtractor will use.
        """
        if not self.weightImageList:
            errtxt="No Weight Images present."
            self.errorList.append((self.modName,errtxt))
            raise Exception, errtxt
            # reset flag image list
        while self.flagImageList:
            del self.flagImageList[0]

        curdir = os.getcwd()
        os.chdir(self.obsFits)
        for im in self.weightImageList:
            try:
                wgtfits = pyfits.open(im)
            except Exception,err:
                self.errorList.append((self.modName,str(err)))
                raise Exception,err

            if len(wgtfits) > 1:
                self.errorList.append((self.modName,"image file is not simple fits."+im))
                raise Exception,"image file is not simple fits."+im

            # build flag image name
            flgfile = im.split("_drz")[0]+'_FLAG.fits'
            self.flagImageList.append(flgfile)
            self.outputList[flgfile] = [im]
            
            # create and initialize the new pyfits object
            flgfits = pyfits.HDUList()
            flgfits.append(pyfits.PrimaryHDU())
            try:
                del flgfits[0].header.ascard["EXTEND"]
            except KeyError:
                pass
            flgfits[0].header = wgtfits[0].header
            flgfits[0].data = numpy.logical_not(wgtfits[0].data).astype(numpy.int16)
            wgtfits.close()
            flgfits[0].header.update('BITPIX',16)            
            flgfits[0].header.update('FILENAME',flgfile)
            flgfits[0].header.update('FILETYPE','FLG')

            # close (write out) the flag image
            flgfits.writeto(flgfile)
            
            self.logfile.write('Made flag image '+flgfile)
            del wgtfits, flgfits

        os.chdir(curdir)
        return

    def makeRmsImage(self):
        """turns the drizzled weight images into RMS images that
            SExtractor will use.
        """
        sn0 = 20. # WZ Threshold for bright pixels
        self.logfile.write("starting makeRmsImage . . .")
        if not self.weightImageList:
            errtxt="No Weight Images present."
            self.errorList.append((self.modName,errtxt))
            raise Exception, errtxt

        # reset rms image list
        while self.rmsImageList:
            del self.rmsImageList[0]

        curdir = os.getcwd()
        os.chdir(self.obsFits)
        for im_wgt in self.weightImageList:
            #im_sci = im_wgt[:-9]+'.fits' #WZ
            #pdb.set_trace()
            im_sci = im_wgt[:-12]+'.fits' 
            if im_sci not in self.sciImageList:
                errtxt = 'makeRmsImage: '+im_sci+' not in sciImageList[]!'
                self.errorList.append((self.modName,errtxt))
                self.logfile.write(errtxt)
            try:
                wgtfits = pyfits.open(im_wgt)
                scifits = pyfits.open(im_sci)
            except:
                self.errorList.append((self.modName,"Cannot make a FITS object out of file "+im_wgt))
                raise Exception,"Cannot make a FITS object out of file "+im_wgt
            if len(wgtfits) > 1 or len(scifits) > 1:
                self.errorList.append((self.modName,"image file is not simple fits."))
                raise Exception,"image file is not simple fits."

            # build rms image name and open as a new file.
            rmsfile = im_wgt.split("_drz")[0]+'_RMS.fits'
            self.rmsImageList.append(rmsfile)
            self.outputList[rmsfile] = [im_wgt]
            
            # make new fits obj and copy WGT/SCI hdr/data to RMS image initially
            rmsfitsobj = pyfits.HDUList()
            rmsfitsobj.append(pyfits.PrimaryHDU())
            try:
                del rmsfitsobj[0].header.ascard["EXTEND"]
            except KeyError:
                pass
            rmsfitsobj[0].header = wgtfits[0].header
            rmsfitsobj[0].data   = scifits[0].data
            if os.path.isfile(rmsfile):
                os.remove(rmsfile)
            rmsfitsobj.writeto(rmsfile)
            del rmsfitsobj

            # reopen the rms image for editing.
            rmsfits = pyfits.open(rmsfile,'update')

            # ratio of default to specified output scales
            area_ratio = (self.asecpix / self.origscale)**2
            if abs(1-area_ratio) < 1e-4: area_ratio = 1
            self.logfile.write('Using area_ratio = %.6f in makeRmsImage' %(area_ratio))

            skyval  = scifits[0].header.get('ALIGNSKY')  # this rescaled below
            exptime = scifits[0].header.get('EXPTIME')
            Ncombed = scifits[0].header.get('NCOMBINE')
            if Ncombed == None:
                errtxt='Error: NCOMBINE not in '+im_sci+' header. _sumSubSkyNcombine() not run?'
                self.logfile.write(errtxt)
                self.errorList.append((self.modName,errtxt))
                raise Exception,errtxt
        
            gain,rn = pyblot._gain_rn(scifits, self.logfile, ext=0)
            # if not told to use header gain, then use 1.0 (data in electrons)
            if not self.hdrGain:
                gain = 1.54 # WZ was 1.0  # XX really?
            self.logfile.write(im_sci+":  gain,rn = "+str(gain)+","+str(rn)+\
                               "  NCOMBINE = "+str(Ncombed)+"  EXPTIME = "+str(exptime))
            if not exptime:
                raise Exception,"No EXPTIME in "+im_sci
            if (skyval == None or skyval < 1):
                warntxt = 'WARNING: found ALIGNSKY of '+str(skyval)+' in '+im_sci+\
                          ' : RMS image may be in WRONG!'
                self.logfile.write(warntxt)
                self.errorList.append((self.modName,warntxt))
                del warntxt
                if skyval == None: skyval=0

            # skyval *= area_ratio #WZ redundant
             # 1. construct variance from sky, sci[], wght[]
            # 2. clip zeros/negatives and infinities to values that will work w/sqrt(),
            #     and so that sqrt(val) = 2e30 > 1e30, => zero weight in SExtractor;
            #     have to work in Float64 to avoid Inf's.
            # 3. take sqrt() and cast as float32
            # 4. tidy header
            # 5. write it out
        
            readVariance = Ncombed*(rn/gain)*(rn/gain)
            self.logfile.write("total read variance = "+str(readVariance)+" for "+im_sci)
            if self.suppInstVar:
                # supplement factor for reference bais subtraction, etc
                # extra var per sec for dark subtraction, repaired cosmic rays, etc.
                dark = 0.00008  #WZ, UVIS dark rate
                totInstVar = readVariance + dark * exptime #WZ      
                #totInstVar =  (_rnVarSuppFac_ * readVariance) + (_exVarPerSec_ * exptime)
                self.logfile.write("adjusted instrumental variance = "+str(totInstVar)+" for "+im_sci)
            else:
                totInstVar = readVariance

            # totInstVar *= area_ratio # WZ redundant
        
            # maybe doing arithmetic in two steps will help conserve memory...
            # newDat  = ((skyval + scifits[0].data.astype(numpy.float64))/gain + totInstVar) * (exptime * area_ratio
            factor = 0.
            nullDat = factor * scifits[0].data.astype(numpy.float64) 
            mean=numpy.mean(nullDat)
            std=numpy.std(nullDat)
            min = numpy.min(nullDat)
            max = numpy.max(nullDat)
            nullDat = scifits[0].data.astype(numpy.float64)  - scifits[0].data.astype(numpy.float64) 
            mean=numpy.mean(nullDat)
            std=numpy.std(nullDat)
            min = numpy.min(nullDat)
            max = numpy.max(nullDat)
            # Set up an array that does not include source counts WZ
            # newDat  = ((skyval+ nullDat)/gain + totInstVar) * area_ratio 
            newDat  = ((skyval+ nullDat)/gain + totInstVar) * (exptime * area_ratio) # This is an array
            # newDat[] is now variance *in counts* times maximum expTime; To be divided by expTime map...
            #sn = newDat/rmsfits[0].data #WZ Mark bright pixels
            #idx = numpy.where(numpy.logical_and(numpy.greater(sn,sn0),numpy.less(wgtfits[0].data,0.)))
            #wgtfits[0].data[idx] = numpy.abs(wgtfits[0].data[idx])
            #for i in range(len(wgtfits[0].data)):
            #    for j in range(len(wgtfits[0].data[i])):
            #        if (wgtfits[0].data[i,j]< 1e-10):
            #            wgtfits[0].data[i,j]= 1e-10
            # idx = numpy.where(numpy.less(wgtfits[0].data,1e-10))
            # wgtfits[0].data[idx] = 1e-10
            # wgtfits[0].data = wgtfits[0].data + 1e-10
            newDat /= wgtfits[0].data
        
            ## now fix up problem values...
            newDat = numpy.where(numpy.logical_or(numpy.greater_equal(newDat,1e38),\
                                                      numpy.less_equal(newDat,0.)),4e38,newDat) #WZ Mar 2013
            scifits.close()
            wgtfits.close()

            rmsfits[0].data = numpy.sqrt(newDat).astype(numpy.float32)
            del scifits, wgtfits, im_wgt, im_sci, readVariance, totInstVar, area_ratio

            # a few token updates to the header, then write it out
            rmsfits[0].header.update('FILENAME',rmsfile)
            rmsfits[0].header.update('FILETYPE','RMS')
            rmsfits.close()
            self.logfile.write('Made rms image '+rmsfile)
            del newDat, rmsfile, rmsfits

        os.chdir(curdir)
        return



    def checkWeightImage(self):
        """Set negative pixels to its absolute values
            September 2011 Wei Zheng
        """

        if not self.weightImageList:
            errtxt="No Weight Images present."
            self.errorList.append((self.modName,errtxt))
            raise Exception, errtxt

        curdir = os.getcwd()
        os.chdir(self.obsFits)
        for im in self.weightImageList:
            # pdb.set_trace()
            whtfits = pyfits.open(im,'update')
            idx = numpy.where(numpy.less(whtfits[0].data,0.))
            whtfits[0].data[idx] = numpy.abs(whtfits[0].data[idx])
            whtfits.close()
        os.chdir(curdir)
        return

    def fixAstrometry(self,obs,skip):
        """Method to call the astrometer module's stuff and do the astrometric corrections.
        The astrometer constructor needs to receive the obs (i.e. DataSet) object.
        If skip=1, only one band (ref.band) will be subject to the correction #WZ
        """

        print "Now correcting astrometric zeropoint..."
        astrom=astrometer.gscMatchup(obs,skip)
        
        try:
            rval = astrom.findAstromCorrs()
        except astrometer.WebQueryError,err:
            warntxt = "Caught a WebQueryError. Astrometric matchup not successful."
            print warntxt
            self.logfile.write(warntxt)
            self.logfile.write(str(err))
            self.errorList.append((self.modName,warntxt))
            self.errorList.append((self.modName,str(err)))
            raise astrometer.WebQueryError,err
        
        if not rval:
            print "Astrometric matchup successful."
            self.logfile.write("Astrometric matchup successful.")
            self.logfile.write("Applying corrections.")
            #pdb.set_trace()
            astrom.applyCorrs()
        return
    

    def writeXml(self):
        """mark up products as xml."""

        curdir = os.getcwd()
        os.chdir(self.obsFits)

        if self.sciImageList:
            for im in self.sciImageList:
                file = xmlUtil.markupImage(im,dataset=self.obsName)
                if file not in self.outputList.keys():
                   self.outputList[file] = [im]
        if self.contextImageList:
            for im in self.contextImageList:
                file = xmlUtil.markupImage(im,dataset=self.obsName)
                if file not in self.outputList.keys():
                    self.outputList[file] = [im]
        if self.weightImageList:
            for im in self.weightImageList:
                file = xmlUtil.markupImage(im,dataset=self.obsName)
                if file not in self.outputList.keys():
                    self.outputList[file] = [im]
        if self.flagImageList:
            for im in self.flagImageList:
                file = xmlUtil.markupImage(im,dataset=self.obsName)
                if file not in self.outputList.keys():
                    self.outputList[file] = [im]
        if self.rmsImageList:
            for im in self.rmsImageList:
                file = xmlUtil.markupImage(im,dataset=self.obsName)
                if file not in self.outputList.keys():
                    self.outputList[file] = [im]
        os.chdir(curdir)
        return


    def mkMsg(self):
        """create and write module level message for this class.
        Most of this is just compiling the info. meta in a dictionary
        of lists where each list is a list of tuples describing the
        tag lines for the particular section of the message.  This tuple 
        format conforms to that used by the xmlMessage class which is
        modeled on basic python argument passing, i.e. (key,*value,**attr).
            outscale,pixfrac now converted to string - GRM following Magee
        """
        self.meta = {}
        self.meta['module']= []
        self.meta['meta']  = []
        self.meta['input'] = []
        self.meta['output']= []
        self.meta['errorlist'] = []

        self.meta['module'].append(('module','name='+self.modName,'version='+__version__,'dataset='+self.obsName))
        self.meta['module'].append(('root',self.root))
        self.meta['meta'].append(('meta',))
        self.meta['meta'].append(('configuration',))
        self.meta['meta'].append(('parameter','name=outscale',str(self.asecpix)))
        self.meta['meta'].append(('parameter','name=pixfrac',str(self.pixfrac)))
        self.meta['meta'].append(('parameter','name=kernel',self.kernel))
        self.meta['meta'].append(('depend',))
        self.meta['meta'].append(('pkg',))
        self.meta['meta'].append(('name','python'))
        self.meta['meta'].append(('version',pyversion.split()[0]))
        self.meta['meta'].append(('pkg',))
        self.meta['meta'].append(('name','xydrizzle'))
        self.meta['meta'].append(('version',pydriz_version))
        self.meta['meta'].append(('pkg',))
        self.meta['meta'].append(('name','drizzle'))
        self.meta['meta'].append(('version',drversion))
        self.meta['meta'].append(('pkg',))
        self.meta['meta'].append(('name','pyfits'))
        self.meta['meta'].append(('version',pyfits.__version__.split()[0]))
        self.meta['meta'].append(('pkg',))
        self.meta['meta'].append(('name','numpy'))
        self.meta['meta'].append(('version',numpy.__version__))
        self.meta['meta'].append(('pkg',))
        self.meta['meta'].append(('name','pyblot'))
        self.meta['meta'].append(('version',pyblot.__version__))        
        self.meta['meta'].append(('pkg',))
        self.meta['meta'].append(('name','astrometer'))
        self.meta['meta'].append(('version',astrometer.__version__))
                                 
        if self.errorList:
            self.meta['errorlist'].append(('errorlist',))
            for pkg,err in self.errorList:
                self.meta['errorlist'].append(('erroritem',err,'frompkg='+pkg))

        # input section
        self.meta['input'].append(('input',))
        for f in self.inputList:
            if string.find(f,"_asn") == -1:
                self.meta['input'].append(('file','type=image/x-fits'))
                self.meta['input'].append(('name',os.path.join("Images",f)))
            else:
                self.meta['input'].append(('file','type=image/x-fits'))
                self.meta['input'].append(('name',os.path.join("Images",f)))

        # output section
        if self.outputList:
            self.meta['output'].append(('output',))
        for f in self.outputList.keys():
            if string.find(f,".xml") == -1:
                self.meta['output'].append(('file','type=image/x-fits'))
                self.meta['output'].append(('name',os.path.join("Images",f)))
                for pred in  self.outputList[f]:
                    self.meta['output'].append(('predecessor',os.path.join("Images",pred)))
            else:
                self.meta['output'].append(('file','type=text/xml'))
                self.meta['output'].append(('name',os.path.join("Images",f)))
                for pred in  self.outputList[f]:
                    self.meta['output'].append(('predecessor',os.path.join("Images",pred)))
        

        # pass this dictionary to the class pMessage...
        msgFile = os.path.join(self.messagedir,self.modName+"_module.xml")
        mmsg = pMessage(self.meta)
        mmsg.writeMsg(msgFile)
        return


#------------------------------------ private stuff -------------------------------------#


    def _makeWCS(self, template=None):
        """ Make a list of tuples containing all the WCS and closely
        related info to be written to each output image header.
        """
        if template:  tempfits = pyfits.open(template)
        else:         tempfits = pyfits.open(self.sciImageList[0])
        temphdr = tempfits[0].header

        dec0 = temphdr['CRVAL2']
        wcstuples = [('CRPIX1', temphdr['CRPIX1']),
                     ('CRPIX2', temphdr['CRPIX2']),
                     ('CRVAL1', temphdr['CRVAL1']),
                     ('CRVAL2', temphdr['CRVAL2'])]
        tempfits.close()
        del temphdr,tempfits

        for L in self.reflogfile_add:
            self.logfile.write(L)
            
        if self.refcdmat['EXTREF']:
        	cdmat = self.refcdmat

        else:
        	PA_final = self.refotherkeys['PA_FINAL']
        	cdmat = matutil.makeCDdict(PA_final, self.asecpix)

        wcstuples.append(('CD1_1', round(1e12*cdmat['CD1_1'])/1.0e12))
        wcstuples.append(('CD1_2', round(1e12*cdmat['CD1_2'])/1.0e12))
        wcstuples.append(('CD2_1', round(1e12*cdmat['CD2_1'])/1.0e12))
        wcstuples.append(('CD2_2', round(1e12*cdmat['CD2_2'])/1.0e12))          

        for key in self.refotherkeys.keys():
            wcstuples.append((key, self.refotherkeys[key]))
    
        return wcstuples


    def _sumSubSkyNcombine(self,parlist, doSky):
        """ Sum the sky values in the input images and put the
        total in the output drizzled image header.
        """
        imSkyDict = {}
        outData = parlist[0]['outdata']        
        totNcomb = 0

        for parDict in parlist:
            if parDict['outdata'] != outData:
                raise Exception,"Images in parlist must have same outdata for _sumSubSky()."
            mefIm  = string.split(parDict['data'],'[')[0]
            NameVerString = string.split(string.split(parDict['data'],'[')[1],']')[0]
            Xname = string.split(NameVerString,',')[0]
            Xver  = int(string.split(NameVerString,',')[1])

            # get the sky and ncombine values
            ff  = pyfits.open(mefIm)
            if doSky:
                newSky = ff[Xname,Xver].header.get(self.alignSky) 
            icomb  = ff[Xname,Xver].header.get("NCOMBINE")
            if icomb == None or icomb < 1:
                self.logfile.write("WARNING: did not find NCOMBINE in "+mefIm+" ("+Xname+","+str(Xver)+").  Assuming N=1.")
                totNcomb += 1
            else:
                totNcomb += icomb
            ff.close()
            del ff,icomb

            if not doSky:
                continue
            
            if newSky == None:
                newSky = 0.0
                errtxt = "Warning: "+self.alignSky+" param not found in "+\
                         mefIm+" "+NameVerString
                self.errorList.append((self.modName,errtxt))
                self.logfile.write(errtxt)

            # append this extension sky to a list for this image
            #pdb.set_trace()
            try:
                imSkyDict[mefIm].append(newSky)
            except:
                imSkyDict[mefIm] = []
                imSkyDict[mefIm].append(newSky)
            del mefIm,NameVerString,Xname,Xver,newSky

        # end of list over input images (parDict's).
        # just normalize and return Ncombine if not doing sky
        totNcomb /= self.NumSci
        if not doSky:
            return (totNcomb,None)

        # ok, now we can just loop over mefIms, average the ext Sky's
        # for each mefIm and sum these averages for a total sky

        skyTot = 0.0
        for mefIm in imSkyDict.keys():
            imSky = 0.0
            skyList = imSkyDict[mefIm]
            for extSky in skyList:
                imSky += extSky
            imSky = imSky / len(skyList)
            skyTot = skyTot + imSky

        return (totNcomb,skyTot)


    def _rotateRect(self,_maxAng,nx,ny):
        """ rotate a rectangle and return ratio of new x,y sizes to old"""
        theta = _maxAng * math.pi/180.
        a     = ny/2.0
        b     = nx/2.0
        c     = math.sqrt(a*a + b*b)
        alpha = math.atan(b/a)

        halfwidth  = c*math.sin(alpha + theta)
        halfheight = c*math.cos(alpha - theta)

        xfac = halfwidth/b
        yfac = halfheight/a

        return (xfac,yfac)

    def printParlist(self, parlist):
        """ takes a PyDrizzle parlist and prints out the contents"""
        for i in range(len(parlist)):
            print 'Member ',i,':'
            for key in parlist[i].keys():
                print ('  %-10s :   ' % key), parlist[i][key]
            print '\n'
        return

    def help(self):
        print self.__doc__
        print self.run_drizzle.__doc__
        print self.run_all.__doc__
        return
    
    def setRef(self): #WZ, AKS
        """ Set up wcs for the output image, using the field of shift list
            new reference: center position
        saved as ref.rootname
        """

        asname=string.split(self.obsAsnDict.keys()[0],'_asn')[0]
        # pdb.set_trace()
        sfile   = os.path.join(self.obsAlign,'shifts_'+asname+'.txt')
        # sfile=self.obsAlign + '/shifts.txt' # WZ read shift par file to  find ref file     #xingxing
        fname=tableio.get_str(sfile,(0,7)) # a383_ir_F110W_asn.fits
        for i in range(len(fname[0])): # shifts_a383_ir_F110W.txt
            if (float(fname[1][i]) > 0):
                refim=fname[0][i]      # find out the reference file.   the seventh number = 1  by XX
        im  = refim + '[1]'
        ref = wcsutil.WCSObject(im)
        #ref.rootname=firstname
        ff = pyfits.open(refim)
        detector=ff[0].header.get('DETECTOR')
        # pdb.set_trace()
        if (detector == 'WFC'):
            ref.band=ff[0].header.get('FILTER1')
            if string.find(ref.band,'CLEAR') > -1:
                ref.band=ff[0].header.get('FILTER2')
        else:
            ref.band=ff[0].header.get('FILTER')
        ff.close()

        #calculate the center position
        xc = ref.naxis1 / 2.
        yc = ref.naxis2 / 2.

        #reference point at the (first) image center
        cen = [xc, yc]
        ra, dec = ref.xy2rd(cen)
        ref.crval1= ra
        ref.crval2= dec
        ref.crpix1 = xc
        ref.crpix2 = yc
        ref.rootname=refim
        # But the central position need to be updated in next step ("setRange")
        return ref

    def setgrid(self,firstname): #WZ, AKS
        """ Set up wcs for the output image, using "ref.fits" as new reference:
        center position
        saved as ref.rootname
        """
        im  = firstname + '[0]'
        #im = 'ref.fits'
        ref = wcsutil.WCSObject(im)
        #ref.rootname=firstname

        #calculate the center position
        xc = ref.naxis1 / 2.
        yc = ref.naxis2 / 2.

        #reference point at the image center
        cen = [xc, yc]
        ra, dec = ref.xy2rd(cen)
        ref.crval1= ra
        ref.crval2= dec
        ref.crpix1 = xc
        ref.crpix2 = yc
                                
        #Now the new frame with new pixel scale, PA=0. north's up
        # But the central position need to be updated in next step ("setRange")
        ref.pscale=0.065
        ref.cd11 = -ref.pscale/3600. # -1.388888E-5 # -1.80555e-5
        ref.cd12 = 0.
        ref.cd21 = 0.
        ref.cd22 = ref.pscale/3600. # 1.388888E-5 # 1.80555e-5
        ref.orient = 0.
        ref.rot=0.

        oldfits=pyfits.open(firstname)
        newfits = pyfits.HDUList()
        newfits.append(pyfits.PrimaryHDU())
        newfits[0].header.ascard = oldfits[0].header.ascard.copy()
        newfits[0].header.update('CRPIX1',ref.crpix1) 
        newfits[0].header.update('CRPIX2',ref.crpix2)
        newfits[0].header.update('CRVAL1',ref.crval1)
        newfits[0].header.update('CRVAL2',ref.crval2)
        newfits[0].data = oldfits[0].data
        #newfits.writeto("ref.fits")
        newfits.close()
        oldfits.close()
        del oldfits,newfits
        return ref

    def xyrot(self,xypos,theta=0):   # WZ, migrated from combDither.py
        """New position [Xp,Yp] after a rotation. Theta in units of deg
        """
        rad = theta * numpy.pi / 180.
        x = xypos[0]
        y = xypos[1]
        xp = x * numpy.cos(rad) - y * numpy.sin(rad)
        yp = x * numpy.sin(rad) + y * numpy.cos(rad)    
        return xp,yp


#####################################################
    
    def setRange(self): #WZ, AKS
        """ Determine the shape of output image, and update ref parameters
        """
        ra, dec = [], []
        for tab in self.obsAsnDict.keys():
            for imkey in self.MatchDict.keys():
                fitsfile = pyfits.open(imkey)
                self.wcs.naxis1 = fitsfile[1].header.get('NAXIS1')
                self.wcs.naxis2 = fitsfile[1].header.get('NAXIS2')
                self.wcs.cd11   = fitsfile[1].header.get('CD1_1')
                self.wcs.cd12   = fitsfile[1].header.get('CD1_2')
                self.wcs.cd21   = fitsfile[1].header.get('CD2_1')
                self.wcs.cd22   = fitsfile[1].header.get('CD2_2')
                self.wcs.crpix1 = fitsfile[1].header.get('CRPIX1')
                self.wcs.crpix2 = fitsfile[1].header.get('CRPIX2')
                self.wcs.crval1 = fitsfile[1].header.get('CRVAL1')
                self.wcs.crval2 = fitsfile[1].header.get('CRVAL2')
                fitsfile.close()

                xy = [0., 0.]
                r, d = self.wcs.xy2rd(xy)
                ra.append(r)
                dec.append(d)

                xy = [0., self.wcs.naxis2 + 1.]
                r, d = self.wcs.xy2rd(xy)
                ra.append(r)
                dec.append(d)

                xy = [self.wcs.naxis1 + 1., 0.]
                r, d = self.wcs.xy2rd(xy)
                ra.append(r)
                dec.append(d)

                xy = [self.wcs.naxis1 + 1., self.wcs.naxis2 + 1.]
                r, d = self.wcs.xy2rd(xy)
                ra.append(r)
                dec.append(d)        
                #print ra
                #print dec
                #print tab,imkey
                rmax, dmax = max(ra), max(dec)
                rmin, dmin = min(ra), min(dec)
                #print rmin,rmax,dmin,dmax
                #pdb.set_trace()

        #Find the maximun and minimum ra/dec values for the 
        #current set of images
        rmax, dmax = max(ra), max(dec)
        rmin, dmin = min(ra), min(dec)
        #Get image center(calculated from first image) and use the maximum
        #distance from it to the edge to calculate final image's dimensions
        xmax, ymax = self.ref.rd2xy([rmin, dmax])
        xmin, ymin = self.ref.rd2xy([rmax, dmin])
        dx = abs(xmax - xmin)
        dy = abs(ymax - ymin)
        #X = int(math.ceil(dx / 100.) * 130) #WZ
        #Y = int(math.ceil(dy / 100.) * 130)        
        X = 5000
        Y = 5000
        xc = (xmax + xmin)/2. # The new center (in old frame)
        yc = (ymax + ymin)/2.
        #pdb.set_trace()
        r, d = self.ref.xy2rd([self.ref.crpix1,self.ref.crpix2])
        r, d = self.ref.xy2rd([xc,yc])
        self.ref.crval1 = r # Update ref info
        self.ref.crval2 = d
        self.ref.crpix1 = X/2.
        self.ref.crpix2 = Y/2.
        self.ref.naxis1 = X
        self.ref.naxis2 = Y
        #dx = self.ref.crpix1 - (xmax +  xmin)/2.
        #dy = self.ref.crpix2 - (ymax +  ymin)/2.
        #xo, yo = self.ref.rd2xy([self.ref.crval1, self.ref.crval2])
        # New center, with the same RA0, Dec0
      	# dx = max(abs(xmax - xc), abs(xmin - xc)) * 2
        # dy = max(abs(ymax - yc), abs(ymin - yc)) * 2
        #r, d = self.wcs.xy2rd([(self.ref.crpix1 - xmin),(self.ref.crpix2 - ymin)])
        #r, d = self.wcs.xy2rd([(self.ref.crpix1 - xmin),(self.ref.crpix2 - ymin)])
        # Find the ceil of up to two decimal places
        #self.ref.crval1 = r
        #self.ref.crval2 = d
        ##self.ref.crpix1 = X/2. # - xmin
        ##self.ref.crpix2 = Y/2. # - ymin 
        #x0 = (xmax + xmin) / 2. # center of the output image
        #y0 = (ymax + ymin) / 2.
        #r, d = self.wcs.xy2rd([self.ref.crpix1, self.ref.crpix2])
        ## self.ref.crpix1 = xc
        ## self.ref.crpix2 = yc
        #print X,Y
        return [X, Y]

    def findMatchRef(self): #WZ

        reffile=''
        """
        for im in self.oriImageList:
            # if string.find(im,'_drz') > -1:               
            instr = drutil.getPrimaryKeyword(im+'[0]','INSTRUME')
            if instr == "ACS":
                reffile = im
        if reffile == '':
            for im in self.oriImageList:
                detector = drutil.getPrimaryKeyword(im+'[0]','DETECTOR')
                if detector == "UVIS":
                    reffile = im
        """
        refband=string.split(self.ref.rootname,'_f')[1][0:3] #WZ
        if reffile == '':
            for im in self.oriImageList:
                if string.find(im,refband) > -1:
                    reffile = im
        print "Tweak reference image: ",reffile
        #pdb.set_trace()
        return reffile

    def prepMatch(self, reffile): #WZ
        """ Calculate shift/rotation 
            new reference: center position
        """
        de=100. # distance from edge
        circle = 25. # matching range in pixel
        rmin = 1.5 # 4 FWHM minimum
        rmax = 12 # 30. # 10. # 6.
        clx = 0.2 # must be a galaxy
        fitsfile = pyfits.open(reffile)
        detector = fitsfile[0].header.get('DETECTOR')
        fitsfile.close()
        #pdb.set_trace()
        # weightfile = reffile.split(".fits")[0]+'_weight.fits'
        # cmd = 'cp '+weightfile+' temp_weight.fits'
        #sproc  = popen2.Popen3(cmd,1)
        #output = sproc.fromchild.readlines()
        #errs   = sproc.childerr.readlines()
        cmd = 'sex '+reffile+' -c '+self.pardir+'/tweak.inpar'
        sproc  = popen2.Popen3(cmd,1)
        output = sproc.fromchild.readlines()
        errs   = sproc.childerr.readlines()

        cmd = 'cp temp.cat temp0.cat'
        sproc  = popen2.Popen3(cmd,1)
        output = sproc.fromchild.readlines()
        errs   = sproc.childerr.readlines()
        x0,y0,mag0,err0,fwhm0,ra0,dec0,cls0=tableio.get_data('temp.cat',(1,2,3,4,5,6,7,8))
        xmin0=numpy.min(x0) # 2250
        xmax0=numpy.max(x0) # 4350
        ymin0=numpy.min(y0) # 1700
        ymax0=numpy.max(y0) # 4000
        # xmin0=min(x0) + de
        # xmax0=max(x0) - de
        # ymin0=min(y0) + de
        # ymax0=max(y0) - de
        tname = self.ref.tweak  # tweakfile : tweak.par under reffile/par  by XX
        # if os.path.exists(tname):
        #     twf = open(tname,'a')
        # else:
        #     twf = open(tname,'w')
        format=str('%4d %8.2f %8.2f %6.2f %6.2f %4.1f')
        # format=str('%4d %8.2f %8.2f %6.2f %6.2f %4.1f %4.0f %4.0f')
        #format2=str('%s %6.3f %5.3f %10.7f')
        #im = self.parlists[0][0]['outdata']
        im = self.outname
        fitsfile = pyfits.open(im)
        detector = fitsfile[0].header.get('DETECTOR')
        filt = fitsfile[0].header.get('FILTER')
 	fitsfile.close()

        weightfile = im.split(".fits")[0]+'_weight.fits'
        cmd = 'cp '+weightfile+' temp_weight.fits'
        sproc  = popen2.Popen3(cmd,1)
        output = sproc.fromchild.readlines()
        errs   = sproc.childerr.readlines()
        if (detector =='WFC'):
            cmd = 'sex '+im+' -c '+self.pardir+'/tweak_wfc.inpar'
        else:
            if (detector =='UVIS'):
                if string.find(filt,'275')>-1 or string.find(filt,'225')>-1:
                    cmd = 'sex '+im+' -c '+self.pardir+'/tweak_uvis_sparse.inpar'
                else:
                    cmd = 'sex '+im+' -c '+self.pardir+'/tweak_uvis.inpar'
            else:
                cmd = 'sex '+im+' -c '+self.pardir+'/tweak_ir.inpar'
        sproc  = popen2.Popen3(cmd,1)
        output = sproc.fromchild.readlines()
        errs   = sproc.childerr.readlines()

        if errs:
            print 'Program choked on '+cmd
            self.logfile.write('Program choked on '+cmd)
        else:
            n,x,y,mag,err,fwhm,ra,dec,cls=tableio.get_data('temp.cat',(0,1,2,3,4,5,6,7,8))
            flag0 = numpy.zeros(len(x0),dtype=int)
            flag = numpy.zeros(len(n),dtype=int)
            index= numpy.zeros(len(n),dtype=int)
            xmin=max(numpy.min(x),xmin0)
            xmax=min(numpy.max(x),xmax0)
            ymin=max(numpy.min(y),ymin0)
            ymax=min(numpy.max(y),ymax0)
            print 'Min/max: ', int(xmin), int(xmax), int(ymin), int(ymax)
            # pdb.set_trace()
            fname0 = im.split("_drz")[0]+'_1.cat' # This is the reference catalog  # a209_f110_1.cat
            fname1 = im.split("_drz")[0]+'_2.cat' # a209_f110_2.cat
            f=open(fname1,"w")
            k=0
            xc0 = self.ref.prod.naxis1/2.0
            yc0 = self.ref.prod.naxis2/2.0
            for i in range(len(x)):
                if (mag[i] < 90 and x[i] > xmin and x[i] < xmax and y[i] > ymin and y[i] < ymax):
                    print >> f, format % (k,x[i],y[i],mag[i],err[i],fwhm[i])
                    # print >> f, format % (k,x[i]-xc0,y[i]-yc0,mag[i],err[i],fwhm[i],x[i],y[i])
                    k=k+1
            print k, " sources used in ", im 
            f0=open(fname0,"w")
            k=0
            for j in range(len(x0)):
                if (mag0[j] < 90 and x0[j] > xmin and x0[j] < xmax and y0[j] > ymin and y0[j] < ymax):
                    print >> f0, format % (k,x0[j],y0[j],mag0[j],err0[j],fwhm0[j])
                    # print >> f0, format % (k,x0[j]-xc0,y0[j]-yc0,mag0[j],err0[j],fwhm0[j],x0[j],y0[j])
                    k=k+1
            print k, " sources used in reference " 
            f0.close()
            f.close()
        return

    def prepCheck(self, reffile): #WZ
        """ Check shift/rotation 
            new reference: center position
        """
        de=100. # distance from edge
        circle = 25. # matching range in pixel
        rmin = 3. # FWHM minimum
        rmax = 10. # 6.
        fitsfile = pyfits.open(reffile)
        detector = fitsfile[0].header.get('DETECTOR')
        fitsfile.close()
        weightfile = reffile.split(".fits")[0]+'_weight.fits'
        cmd = 'cp '+weightfile+' temp_weight.fits'
        sproc  = popen2.Popen3(cmd,1)
        output = sproc.fromchild.readlines()
        errs   = sproc.childerr.readlines()
        if (detector =='WFC'):
            cmd = 'sex '+reffile+' -c '+self.pardir+'/tweak_wfc.inpar'
        else:
            if (detector =='UVIS'):
                cmd = 'sex '+reffile+' -c '+self.pardir+'/tweak_uvis.inpar'
            else:
                cmd = 'sex '+reffile+' -c '+self.pardir+'/tweak_ir.inpar'
        sproc  = popen2.Popen3(cmd,1)
        output = sproc.fromchild.readlines()
        errs   = sproc.childerr.readlines()
        cmd = 'cp temp.cat temp1.cat'
        sproc  = popen2.Popen3(cmd,1)
        output = sproc.fromchild.readlines()
        errs   = sproc.childerr.readlines()                        
        x0,y0,mag0,err0,fwhm0,ra0,dec0,cls0=tableio.get_data('temp.cat',(1,2,3,4,5,6,7,8))
        xmin0=min(x0) + de
        xmax0=max(x0) - de
        ymin0=min(y0) + de
        ymax0=max(y0) - de
        format=str('%4d %8.2f %8.2f %6.2f %6.2f %4.1f %10.6f %10.6f %4.2f')
        #format2=str('%s %6.3f %5.3f %10.7f')
        for im in self.sciImageList:
            # print "im = ",im
            fitsfile = pyfits.open(im)
            detector = fitsfile[0].header.get('DETECTOR')
            fitsfile.close()
            weightfile = im.split(".fits")[0]+'_weight.fits'
            cmd = 'cp '+weightfile+' temp_weight.fits'
            sproc  = popen2.Popen3(cmd,1)
            output = sproc.fromchild.readlines()
            errs   = sproc.childerr.readlines()
            if (detector =='WFC'):
                cmd = 'sex '+im+' -c '+self.pardir+'/tweak_wfc.inpar'
            else:
                if (detector =='UVIS'):
                    cmd = 'sex '+im+' -c '+self.pardir+'/tweak_uvis.inpar'
                else:
                    cmd = 'sex '+im+' -c '+self.pardir+'/tweak_ir.inpar'
            sproc  = popen2.Popen3(cmd,1)
            output = sproc.fromchild.readlines()
            errs   = sproc.childerr.readlines()
            #if errs:
            #    print 'Program choked on '+cmd
            #    self.logfile.write('Program choked on '+cmd)
            n,x,y,mag,err,fwhm,ra,dec,cls=tableio.get_data('temp.cat',(0,1,2,3,4,5,6,7,8))
            #n,x,y,mag,err,fwhm = tableio.get_data('temp.cat',(0,1,2,3,4,5))
            flag0 = numpy.zeros(len(x0),dtype=int)
            flag = numpy.zeros(len(n),dtype=int)
            index= numpy.zeros(len(n),dtype=int)
            xmin=min(x) + de
            if (xmin < xmin0):
                xmin = xmin0
            xmax=max(x) - de
            if (xmax > xmax0):
                xmax = xmax0
            ymin=min(y) + de
            if (ymin < ymin0):
                ymin = ymin0
            ymax=max(y) - de
            if (ymax > ymax0):
                ymax = ymax0
            print 'Min/max: ', xmin, xmax, ymin, ymax
            fname0 = im.split("_sci")[0]+'_3.cat'
            fname1 = im.split("_sci")[0]+'_4.cat'
            f0=open(fname0,"w")
            f=open(fname1,"w")
            for i in range(len(x)):
                if (mag[i] < 90 and fwhm[i] > rmin and fwhm[i] < rmax):
                    for j in range(len(x0)):
                        if (abs(x[i]-x0[j]) < circle and abs(y[i]-y0[j]) < circle and mag0[j] < 90 and fwhm0[j] > rmin and fwhm0[j] < rmax):
                        # if (abs(x[i]-x0[j]) < circle and abs(y[i]-y0[j]) < circle and mag0[j] < 90):
                            index[i]=j
                            flag[i] = flag[i] + 1
                            flag0[j] = flag0[j] + 1
            k=0
            for i in range(len(x)):
                if (flag[i] >= 1):
                    j=index[i]
                    if (flag0[j] >= 1 and mag[i] < 90):
                        print >> f, format % (k,x[i],y[i],mag[i],err[i],fwhm[i],ra[i],dec[i],cls[i])
                        print >> f0, format % (k,x0[j],y0[j],mag0[j],err0[j],fwhm0[j],ra0[j],dec0[j],cls0[j])
                        k=k+1
            print k, " sources used in ", im
            f0.close()
            f.close()
        return

    def alignTweak(self,parList,groupList): #WZ Oct 2012                                    
        """ Check shift/rotation                                                            
            new reference: center position                                                  
        """
        de=100. # distance from edge                                                        
        circle = 25. # matching range in pixel                                              
        rmin = 3. # FWHM minimum                                                            
        rmax = 6. # 10                                                                      
        detector=parList[0]['detector']
        reffile='group_'+groupList[0]+'.fits'
        #fitsfile = pyfits.open(reffile)                                                    
        #detector = fitsfile[0].header.get('DETECTOR')                                      
        #fitsfile.close()                                                                   
        #mapfile = string.replace(reffile,'.fits','_map.fits')                              
        #cmd = 'cp '+mapfile+' temp_map.fits'                                               
        #sproc  = popen2.Popen3(cmd,1)                                                      
        #output = sproc.fromchild.readlines()                                               
        #errs   = sproc.childerr.readlines()                                                
        if (detector =='WFC'):
            cmd = 'sex '+reffile+' -c '+self.pardir+'/aligntweak_wfc.inpar'
        elif (detector =='HRC'): # WZ Apr 2013
            cmd = 'sex '+reffile+' -c '+self.pardir+'/aligntweak_hrc.inpar'
        elif (detector =='UVIS'):
            cmd = 'sex '+reffile+' -c '+self.pardir+'/aligntweak_uvis.inpar'
        else:
            cmd = 'sex '+reffile+' -c '+self.pardir+'/aligntweak_ir.inpar'
        sproc  = popen2.Popen3(cmd,1)
        output = sproc.fromchild.readlines()
        errs   = sproc.childerr.readlines()
        cmd = 'cp temp.cat temp0.cat'
        sproc  = popen2.Popen3(cmd,1)
        output = sproc.fromchild.readlines()
        errs   = sproc.childerr.readlines()
        format=str('%4d %8.2f %8.2f %6.2f %6.2f %4.1f %10.6f %10.6f %4.2f')
        #format2=str('%s %6.3f %5.3f %10.7f')                                           print 'Tweak groups. This might take a bit' # WZ Apr 2013            
        imList=[]
        for i in range(len(parList)):
            fname= string.split(parList[i]['outdata'],'_dz_')[1]
            imList.append(fname)
        for i in range(len(groupList)):
            imfile = 'group_'+groupList[i]+'.fits'
            ## WZ: 2013
            if (imfile == reffile):
                tx=0.
                ty=0.
                tr=0.
            else:
                #mapfile = string.replace(imfile,'.fits','_map.fits')                           
                #cmd = 'cp '+mapfile+' temp_map.fits'                                           
                #sproc  = popen2.Popen3(cmd,1)                                                  
                #output = sproc.fromchild.readlines()                                           
                #errs   = sproc.childerr.readlines()                     
                if (detector =='WFC'):
                    cmd = 'sex '+imfile+' -c '+self.pardir+'/aligntweak_wfc.inpar'
                elif (detector =='HRC'): # WZ Apr 2013
                    cmd = 'sex '+imfile+' -c '+self.pardir+'/aligntweak_hrc.inpar'
                elif (detector =='UVIS'):
                    cmd = 'sex '+imfile+' -c '+self.pardir+'/aligntweak_uvis.inpar'
                else:
                    cmd = 'sex '+imfile+' -c '+self.pardir+'/aligntweak_ir.inpar'
                sproc  = popen2.Popen3(cmd,1)
                output = sproc.fromchild.readlines()
                errs   = sproc.childerr.readlines()
                #if errs:                                                                       
                #    print 'Program choked on '+cmd                                             
                #    self.logfile.write('Program choked on '+cmd)                               
                n,x,y,mag,err,fwhm,ra,dec,cls=tableio.get_data('temp.cat',(0,1,2,3,4,5,6,7,8))
                #n,x,y,mag,err,fwhm = tableio.get_data('temp.cat',(0,1,2,3,4,5))                
                cmd = 'tweak_fine temp0.cat temp.cat'
                print 'Tweak_fine  ',imfile,' ',reffile # WZ April 2013
                # pdb.set_trace()                                                               
                sproc  = popen2.Popen3(cmd,1)
                matchout = str(sproc.fromchild.readlines())
                errs   = sproc.childerr.readlines()
                if errs:
                    self.logfile.write('Program choked on '+cmd)
                    match = -1
                else:
                    tx=float(matchout[2:].split()[0]) # Tweak value                             
                    ty=float(matchout[2:].split()[1])
                    tr=float(matchout[2:].split()[2]) # /numpy.pi*180. #WZ rad > degree         
            for j in range(len(parList)):
                if (groupList[i]==parList[j]['group']):
                    parList[j]['tx']=tx
                    parList[j]['ty']=ty
                    parList[j]['tr']=tr
        fname=os.path.join(self.obsAlign,'alignTweak.txt')
        fw = open(fname,'w')
        fmt=str('%s\t%6.2f %6.2f %10.7f %3d %s')
        for j in range(len(parList)):
            if (string.find(parList[j]['outdata'],'SCI_1')>-1):
                rname=string.split(parList[j]['outdata'],'_')[2]
                print >> fw, fmt % (rname,parList[j]['tx'],parList[j]['ty'],\
                    parList[j]['tr'],parList[j]['pa'],parList[j]['group'])
        fw.close()
        return


    def run_Match(self,reffile,tweak): #WZ
        """ Check shift/rotation 
            Done in two steps:
            (1) "match": If first round, this is the coarse result
            (2) "tweak": for later rounds, fine tune the values 
        """
        circle=25 # 15
        match= 0
        clx = 0.2
        im = self.outname
        f1=str('%s %6.2f %6.2f %10.7f')
        f2=str('%s %6.2f %6.2f %10.7f %s %s')
        # pdb.set_trace()
        xc0 = self.ref.prod.naxis1/2.0
        yc0 = self.ref.prod.naxis2/2.0
        target=string.split(self.outname,'_')[0]
        tmp = string.replace(self.root,'datasets','ingest')
        tmp1=string.split(tmp,target)[0]
        rname = os.path.join(tmp1,target,'par','ref.cat')
        #rname = string.replace(self.root,'datasets','ingest')+'/ref.cat'
        if (os.path.isfile(rname) and self.detector=='UVIS'):
            #fname0 = rname
            #WZ
            #xc0 = self.ref.prod.naxis1/2.0
            #yc0 = self.ref.prod.naxis2/2.0
            # rdata=tableio.get_str(rname,(0,1,2,3)) 
            rdata=tableio.get_data(rname,(1,2,3)) #WZ Dec 2012 
            xx=rdata[0] # -xc0
            yy=rdata[1] # -yc0
            mm=rdata[2]
            fname0 = im.split("_drz")[0]+'_1.cat'
            f=open(fname0,"w")
            format=str('%4d %8.2f %8.2f %6.2f')
            for i in range(len(rdata)):
                print >> f, format % (i,xx[i],yy[i],mm[i])
            self.logfile.write('Use '+rname)
            del rdata,rname,xx,yy,mm
        else:
            fname0 = im.split("_drz")[0]+'_1.cat' # This is the reference catalog
        fname1 = im.split("_drz")[0]+'_2.cat'
        inpname=string.split(self.outname,'_drz')[0]
        if os.path.exists('tweak.log'):
            firsttry=0
            x,y,r=tableio.get_data('tweak.log',(1,2,3))             
            j=len(x)
            x0=x[j-1]
            y0=y[j-1]
            r0=r[j-1]
            del x,y,r
            # Occasionally, residual offsets are large and beyond the tweak 
            # range. Therefore an added run of match for additional offsets:
            #cmd = 'match '+fname1+' 1 2 3 '+fname0+' 1 2 3 scale=1.0 matchrad='+str(circle)
            # pdb.set_trace()
            #sproc  = popen2.Popen3(cmd,1)
            #matchout = str(sproc.fromchild.readlines())
            #errs   = sproc.childerr.readlines()
            #if errs:
            #    self.logfile.write('Initial guess not found. Check images and make a guess')
            #    exit()
            #else:
            #    print matchout
            #    xoff=float(string.split((matchout[2:].split()[1]),'=')[1])
            #    yoff=float(string.split((matchout[2:].split()[4]),'=')[1])
            #    roff=float(string.split((matchout[2:].split()[3]),'=')[1])
            # x, y, r are the additional offsets
            # cmd ='tweak '+fname1+' '+fname0+' '+str(xoff)+' '+str(yoff)+' '+str(roff)
            if (tweak>1 and self.detector=='UVIS'):
                cmd = 'tweak ' + fname1 +' '+ fname0 + ' 5'
            else:
                cmd = 'tweak ' + fname1 +' '+ fname0  # + ' 0 0 ' + str(xc0) + ' ' + str(yc0)  # Assuming small shifts
            # pdb.set_trace()
            sproc  = popen2.Popen3(cmd,1)
            matchout = str(sproc.fromchild.readlines())
            # pdb.set_trace()
            errs   = sproc.childerr.readlines()
            if errs:
                self.logfile.write('Program choked on '+cmd)
                match = -1
            else:
                # dx=float(string.split((matchout[2:].split()[1]),'=')[1])
                # dy=float(string.split((matchout[2:].split()[4]),'=')[1])
                # dr=float(string.split((matchout[2:].split()[3]),'=')[1])
                chi2=float(string.split(string.split(matchout,'=')[1])[0])
                if (chi2>8): # No good. Try a maller tweak box
                    matchout2=matchout
                    chi2a=chi2
                    cmd = 'tweak ' + fname1 +' '+ fname0 + ' 10' 
                    sproc  = popen2.Popen3(cmd,1)
                    matchout = str(sproc.fromchild.readlines())
                    errs   = sproc.childerr.readlines()
                    chi2=float(string.split(string.split(matchout,'=')[1])[0])
                    if chi2>chi2a: # No help. Use the old one
                        matchout=matchout2
                        chi2=chi2a
                if (os.path.isfile(rname) and self.detector=='UVIS'):  # Try simple tweak
                    matchout2=matchout
                    chi2a=chi2
                    cmd = 'tweak_easy ' + fname1 +' '+ fname0 + ' 10' 
                    sproc  = popen2.Popen3(cmd,1)
                    matchout = str(sproc.fromchild.readlines())
                    errs   = sproc.childerr.readlines()
                    chi2=float(string.split(string.split(matchout,'=')[1])[0])
                    if chi2>chi2a: # No help. Use the old one
                        matchout=matchout2
                dx=float(matchout[2:].split()[0])
                dy=float(matchout[2:].split()[1])
                dr=float(matchout[2:].split()[2])
                x = x0 - dx
                y = y0 - dy
                r = r0 + dr
                substr=string.split(matchout,'Chi^2=')[1]
                substr2=string.split(substr,')')[0]
                chi2=float(string.split(substr2)[0])
                date = ptime()  # ptime from pUtil
                print "Tweak.log updated"
                self.logfile.write('Tweak: dX= '+str(x)+' dY= '+str(y)+' dR= '+str(r))
                twf = open('tweak.log','a')
                print >> twf, f2 % (inpname, x,y,r,matchout[2:-4],date)
                # del x,y,r,dx,dy,dr,xoff,yoff,roff,date
                del x,y,r,dx,dy,dr,date
        else: # This is the first try
            firsttry=1
            chi2=999.
            if os.path.exists(self.ref.tweak):
                x,y,r=tableio.get_data(self.ref.tweak,(1,2,3))             
                j=len(x)
                x0=x[j-1]
                y0=y[j-1]
                r0=r[j-1]
            else:
                cmd = 'match '+fname1+' 1 2 3 '+fname0+' 1 2 3 scale=1.0 matchrad='+str(circle)
                sproc  = popen2.Popen3(cmd,1)
                matchout = str(sproc.fromchild.readlines())
                errs   = sproc.childerr.readlines()
                if errs:
                    self.logfile.write('Initial guess not found. Check images and make a guess')
                    exit()
                else:
                    print matchout
                    x0=float(string.split((matchout[2:].split()[1]),'=')[1])
                    y0=float(string.split((matchout[2:].split()[4]),'=')[1])
                    r0=float(string.split((matchout[2:].split()[3]),'=')[1])
            # direct output
            self.logfile.write('First shifts: dX= '+str(x0)+' dY= '+str(y0)+' dR= '+str(r0))
            twf = open('tweak.log','w')
            print >> twf, f1 % (inpname,x0,y0,r0) # dX, dY, rot
            print "Tweak.log generated"
        twf.close()
        return chi2

    def checkFinalMatch(self,reffile,tweak=0): #WZ
        """ Calculate shift/rotation 
            new reference: center position
        """
        circle=25 # 15
        match=0
        clx = 0.2
        #format=str('%4d %8.2f %8.2f %6.3f %6.3f')
        tname = self.root + '/tweak.log'
        if os.path.exists(tname):
            twf = open(tname,'a')
        else:
            twf = open(tname,'w')
        ipf = open("inp.lis","w")
        format=str('%s %6.2f %5.2f %10.7f')
        # for item in self.parlists: #WZ
        im = self.parlist[0]['outdata']
        rootname=string.split(im,'_sci')[0]
        #pdb.set_trace()
        print >> ipf, rootname
        if string.find(reffile,im) > -1: 
            dx = dy = dr = dX = dY = rot = 0.0 # This is the file self
        else:
            dr=dx=dy=dX=dY=rot=0.                
            fname0 = im.split("_sci")[0]+'_3.cat'
            fname1 = im.split("_sci")[0]+'_4.cat'
            #pdb.set_trace()
            cmd = 'match '+fname1+' 1 2 3 '+fname0+' 1 2 3 scale=1.0 matchrad='+str(circle)
            sproc  = popen2.Popen3(cmd,1)
            matchout = str(sproc.fromchild.readlines())
            errs   = sproc.childerr.readlines()
            if errs:
                print 'Program choked on '+cmd
                self.logfile.write('Program choked on '+cmd)
                match=-1
            else:
                if (tweak):
                    cmd = 'tweak '+' '+fname1+' '+fname0
                    print "Calculating tweaking parameters for "+im+"  ..."
                    # pdb.set_trace()
                    sproc  = popen2.Popen3(cmd,1)
                    matchout = str(sproc.fromchild.readlines())
                    errs   = sproc.childerr.readlines()
                    if errs:
                        print 'Program choked on '+cmd
                        self.logfile.write('Program choked on '+cmd)
                        match=-1
                    else:
                        print matchout
                        rot=float(matchout[2:].split()[0])
                        dX=float(matchout[2:].split()[1])
                        dY=float(matchout[2:].split()[2])
                        rms=matchout.split('rms:')[1]
                        print 'rms in alignment: ',rms[:-5],im
                        self.logfile.write('rms in alignment: '+' '+rms[:-5]+' '+im)
                else:
                    print matchout
                    dX=float(string.split((matchout[2:].split()[1]),'=')[1])
                    dY=float(string.split((matchout[2:].split()[4]),'=')[1])
                    rot=float(string.split((matchout[2:].split()[3]),'=')[1])
                wcs = self.ref.copy()
                if (string.find(item[0]['data'],'[') > -1):
                    inpname=string.split(item[0]['data'],'[')[0]
                else:
                    inpname=item[0]['data']
                fitsfile = pyfits.open(inpname)
                wcs.naxis1 = fitsfile[1].header.get('NAXIS1')
                wcs.naxis2 = fitsfile[1].header.get('NAXIS2')
                wcs.cd11 = fitsfile[1].header.get('CD1_1')
                wcs.cd12 = fitsfile[1].header.get('CD1_2')
                wcs.cd21 = fitsfile[1].header.get('CD2_1')
                wcs.cd22 = fitsfile[1].header.get('CD2_2')
                wcs.crpix1 = fitsfile[1].header.get('CRPIX1')
                wcs.crpix2 = fitsfile[1].header.get('CRPIX2')
                wcs.crval1 = fitsfile[1].header.get('CRVAL1')
                wcs.crval2 = fitsfile[1].header.get('CRVAL2')
                fitsfile.close()
                #xc = wcs.naxis1 / 2 + 1  # old frame
                #yc = wcs.naxis2 / 2 + 1
                #print "Center: ",xc,yc
                #rac, decc = wcs.xy2rd([xc,yc])
                #x,y = self.ref.rd2xy([rac,decc]) # on the output grid
                ## print x,y," on the common grid"  

            self.logfile.write('DQ: '+im+' dX= '+str(dX)+' dY= '+str(dY)+' dR= '+str(rot))
            print >> twf, format % (string.split(im,'.fits')[0],dX,dY,rot) # image_name, dX, dY, rot
        twf.close()
        ipf.close()
        return match

    def updateTweak(self): #WZ                                                          
        """                                                                             
           Find the best fit from 'tweak.log'                                           
           copy the best 'bck' file into 'drz'                                          
           Update the final result into a file in general output directory              
        """
        #pdb.set_trace()
        j=string.find(self.root,self.obsName)
        tmp=string.split(self.root,'/')
        tmp1='run_'+tmp[-2]+'.log'
        pname=os.path.join(self.root[0:j-1],tmp[-2],tmp1)
        if os.path.exists(pname):
            twf = open(pname,'a')
        else:
            twf = open(pname,'w')
            txt='     Name  \t    dX       dY    Rot (rad) ResX ResY Npt   Date'
            print >> twf, txt
            del txt
        f = open('tweak.log','r')
        i=0
        min=99999.
        for line in f.readlines():
            if (string.find(line,'Chi^2')>-1):
                tmp=string.split(line,'Chi^2=')[1]
                chi2=float(string.split(tmp)[0])
                if (chi2<min):
                    min=chi2
                    id=i
                    best=line[:-1]
            i=i+1
        f.close()
        tmp=string.split(best)
        name=tmp[0]
        dx=float(tmp[1])
        dy=float(tmp[2])
        dr=float(tmp[3])
        chi2=float(string.split(string.split(best,'=')[1],' ')[1])
        rex=float(string.split(string.split(best,'rms:')[1])[0])
        rey=float(string.split(string.split(best,'dy')[1])[0])
        npt=int(string.split(string.split(best,'Use')[1])[0])
        date=(string.split(best,'pairs')[1])[:-1]
        format=str('%s\t%8.2f %8.2f %10.7f %4.2f %4.2f %3d %s')
        print >> twf, format % (name,dx,dy,dr,rex,rey,npt,date)
        twf.close()
        # pdb.set_trace()                                                               
        if (id < i-1):
            bck='bck'+str(id)
            back=string.replace(self.outname,'drz',bck)
            os.rename(back,self.outname)
            back=string.replace(self.outcontext,'drz',bck)
            os.rename(back,string.replace(back,bck,'drz'))
            back=string.replace(self.outweight,'drz',bck)
            os.rename(back,string.replace(back,bck,'drz'))
            del best,line,chi2,dx,dy,rex,rey,npt,date,tmp,bck,back
        return


def read_sextractor(filter = 'f555w'):
	'''
	CREATED: 08/11/11 (Viana)
	UPDATED: 08/15/11 (Viana)
	'''
	
	# Define the output path.
	output_path = os.environ['DATASETS']
	output_path += '/' + filter + '/'
	
	# Gather the files.
	query = output_path + '*.sex'
	catalog_list_1 = glob.glob(query)

	# Zero the counters.	


#############WZZZZZZZZZZ
