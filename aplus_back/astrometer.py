#!/usr/bin/env python

# $Id: astrometer.py,v 1.8 2006/05/22 17:38:53 meurer Exp $
# ---------------------------------------------------------------------

__version__      = '$Revision: 1.8 $ '[11:-3]
__version_date__ = '$Date: 2006/05/22 17:38:53 $ '[7:-3]
__author__       = 'J. Blakeslee, R. White, K. Anderson'

import sys,os,string,popen2,cStringIO,math,types,urllib
import wcsclass,pyfits
from   msg   import pMessage
import xmlUtil
import time
import subprocess
pyversion = sys.version
# Host and URL information for search will need to be changed when
# new catalog version is released.

_HOST_ = "gsss.stsci.edu"
_URL_ = "/webservices/vo/ConeSearch.aspx"
_CATALOG_ = "GSC23"

_minSig_ = 0.2        # don't believe rms less than this (arcsec)

# Sets maximum number of queries to try before bailing out.
MAXWQtry=10

class WebQueryError(IOError):
    """Provides an exception class for web query errors and timeouts."""
    

class gscMatchup:
    """ This class provides various methods for searching and matching
    against extended version of GSC2 catalog, determining the pointing
    error, and correcting the image header CRVAL1,CRVAL2 keywords.
    The methods have mostly been cobbled together from Rick's scripts and
    wfp align module functions.

    Pipeline should run it like this:

      mm=astrometer.gscMatchup(obs)
      if (mm.findAstromCorrs() != 0):
          # then it didn't work!
          # note this in long and move on to next module
      else:
          mm.applyCorrs()

    The mkMsg method is remains callable, though deprecated for pipeline
    purposes as this module and class is used by the combDither module rather
    than called directly by the pipeline.
        
    """
    def __init__(self,obs):
        """ constructor gets Observation object (obs) from pipeline. """
        self.modName        = string.split(string.split(str(self))[0],'.')[0][1:]
        self.root           = obs.newobspath      # root path of the observation dir
        self.obsName        = obs.newobs
        self.pardir         = obs.newpardir
        self.Imagedir       = obs.fitsdir         # path to fits Images/ dir
        self.Imagedir       = obs.fitsdir         # path to fits Images/ dir
        self.messagedir     = obs.messagedir      # where the module message will go.
        self.configdir      = obs.configdir
        self.logfile        = obs.logfile
        self.paramFile      = os.path.join(self.pardir,self.modName + '.param')
        self.inparFile      = os.path.join(self.pardir,self.modName + '.inpar')
        self.sciImlist      = obs.sciImageList
        self.ctxImlist      = obs.contextImageList
        self.wgtImlist      = obs.weightImageList
        self.rmsImlist      = obs.rmsImageList
        # make astrom working dir
        self.astromdir      = obs.astromdir = os.path.join(self.root,'astrom')
        if not os.path.isdir(self.astromdir):
            os.mkdir(self.astromdir)

        self.errorList      = []
        self.inputList      = []
        self.outputList     = {}
        self.logfile.write('****************************************')
        self.logfile.write(self.modName+": Instantiated astrometer object.")

        # basic input check
        self.Nim=Nim = len(self.sciImlist)
        if Nim < 1 or Nim!=len(self.rmsImlist):
            errtxt="Image Lists of zero or differing lenths."
            self.errorList.append((self.modName,errtxt))
            self.logfile.write('Error: '+errtxt)
            raise IOError, errtxt
        # init the input list
        for i in self.sciImlist:
            self.inputList.append(i)


    def applyCorrs(self):
        """ Apply the astrometric corrections:
        update sci,context,weight,RMS image header CRVAL1,CRVAL2 values.
        """
        curdir = os.getcwd()
        os.chdir(self.Imagedir)
        
        for imlist in [self.sciImlist, self.ctxImlist, self.wgtImlist, self.rmsImlist]:
            for im in imlist:
                self.logfile.write("Applying header updates: "+im)
                ff = pyfits.open(im,'update')
                ff[0].header.update('CRVAL1',self.CRVAL1new)
                ff[0].header.update('CRVAL2',self.CRVAL2new)
                ff[0].header.update('AMDRA',self.bestrans['dx'])
                ff[0].header.update('AMDDEC',self.bestrans['dy'])
                ff[0].header.update('AMNTCH',self.bestrans['Nmed'])
                ff[0].header.update('AMSGRA',self.bestrans['sigx'])
                ff[0].header.update('AMSGDEC',self.bestrans['sigy'])
                
                # document with comments
                ff[0].header.add_comment("astrom: CRVAL1,CRVAL2 corrected by wfp astrometer module.")
                ff[0].header.add_comment("astrom: Old CRVAL1,CRVAL2: %.8f %.8f" \
                                         %(self.crval1in,self.crval2in))
                ff[0].header.add_comment("astrom: Applied corrections (arcsec) dRA,dDec:  %.4f %.4f"\
                                         %(self.bestrans['dx'],self.bestrans['dy']))
                ff[0].header.add_comment("astrom: Correction errors (arcsec) e_dRA,e_dDec: %.4f %.4f"\
                                         %(self.bestrans['dx_err'],self.bestrans['dy_err']))
                ff[0].header.add_comment("astrom: Nmatch,sigx,sigy (arcsec):  %d %.4f %.4f"\
                                         %(self.bestrans['Nmed'],self.bestrans['sigx'],self.bestrans['sigy']))
                ff[0].header.add_comment("APLUS: Properietary data by CLASH pipeline")
                ff.close()
                self.logfile.write("Image: " + im)
                self.logfile.write("********************")
                self.logfile.write("Old CRVAL1: " + str(self.crval1in) + " Old CRVAL2: " + str(self.crval2in))
                self.logfile.write("New CRVAL1: " + str(self.CRVAL1new) + " New CRVAL2: " + str(self.CRVAL2new))
                self.logfile.write("RA Shift : " + str(self.bestrans['dx']) + " Dec Shift : " + str(self.bestrans['dy']))
                self.logfile.write("Sigma RA  : " + str(self.bestrans['sigx']) + " Sigma Dec : " + str(self.bestrans['sigy']))
                self.logfile.write("No. match : " + str(self.bestrans['Nmed']))

        os.chdir(curdir)
        self.logfile.write("Updated image headers with new CRVAL's.")
        self.logfile.write("****************************************")

        return
        

    def findAstromCorrs(self):
        """ This 'meta-method'  runs the methods necessary for determining
        astrometric zeropoint corrections.
        """
        self.logfile.write("Entered findAstromCorrs - will run: "+\
                           "makeGSCcat(), makeObjcats(), doMatching().")

        if self.makeGSCcat() != 0:
            return -1
        if self.makeObjcats()!= 0:
            return -1
        if self.doMatching() != 0:
            # here we want to remake the GSCcat using a "chopSpur" flag,
            # if the cat has a goodly number of objects
            if(self.Ncull >= 10):
                print "Retrying matchup with only GSC objs detected in 2 bands..."
                self.logfile.write("Retrying matchup with only GSC objs detected in 2 bands...")
                if self.makeGSCcat(chopSpur=1) != 0:
                    return -1
                if self.makeObjcats()!= 0:
                    return -1
                if self.doMatching() != 0:
                    return -1

        return 0
            

    def doMatching(self):
        """ Match object .proj catalogs against GSC proj catalog. """

        self.alltrans = []
        curdir = os.getcwd()
        os.chdir(self.astromdir)

        for cat in self.projcats:
            self.logfile.write("Beginning match attempt for %s" %cat)
            nobj = max(80,self.Nrefobjs)
            base = 'match %s 1 2 3 %s 1 2 3 identity recalc nobj=%d medtf medsigclip=2.7 '\
                   %(cat,self.GSCmatchin,nobj)
            converged=0
            mrad=4.5
            xsh=0.
            ysh=0.
            Retry=0
            while not converged and mrad<10 and Retry<11:
                ferr,mdict=self._runmatch(base,mrad,xsh,ysh)
                if ferr:
                    # match not found, make radius bigger and retry
                    mrad  += 0.5
                    Retry += 1
                    self.logfile.write("Retry: %d ... set mrad: %s" %(Retry,mrad))
                    continue

                dx   = mdict['dx']
                dy   = mdict['dy']
                sigx = mdict['sigx']
                sigy = mdict['sigy']

                # if the extra shift is too big compared to the
                # match radius, update input shift and retry
                if max(abs(dx-xsh),abs(dy-ysh)) > min(0.4, mrad/10.):
                    xsh=dx
                    ysh=dy
                    Retry += 1
                    continue

                # if sigma seems too big, shrink mrad
                if (max(sigx,sigy) > 0.7 and mrad>2.9):
                    mrad -= 0.5
                    Retry += 1
                    continue

                # otherwise, it looks like we have a good fit
                # but we want to tune it up with smaller matchrad
                # Don't increment retry in this final case.
                if mrad > 3.9 or (mrad>2.9 and max(sigx,sigy)>0.4):
                    mrad = max(mrad-1,2.5)
                    continue

                self.alltrans.append(mdict)
                converged=1
                self.logfile.write("""%s trans:
                mdx,mdy   =  %.4f %.4f
                edx,edy   =  %.4f %.4f
                sigx,sigy =  %.4f %.4f
                """%(cat, mdict['dx'],mdict['dy'],mdict['dx_err'],mdict['dy_err'],\
                     mdict['sigx'],mdict['sigy']))
        if len(self.alltrans)<1:
            errtxt = "WARNING: No Image successfully matched against extended GSC!"
            self.logfile.write(errtxt)
            self.errorList.append((self.modName,errtxt))
            return -1

        os.chdir(curdir)
        # finally, choose best match transform and calculate
        # the new CRVAL's
        self._pickBestrans()

        return 0


    def _pickBestrans(self):
        """ Choose the best from among available offsets. """
        ###
        ### CHECK OVER THIS PLEASE!!!!
        ###
        min_uncty = 999.
        for trans in self.alltrans:
            # maxerr = max(trans['dx_err'],trans['dy_err'])
            rerr = math.sqrt(trans['dx_err']**2 + trans['dy_err']**2)
            if rerr < min_uncty:
                min_uncty = rerr
                self.bestrans = trans
            del rerr

        self.logfile.write("Pick  dx,dy:  %.4f %.4f   edx,edy:  %.4f %.4f   sigx,sigy:  %.4f %.4f  rerr: %.4f"\
                           %(self.bestrans['dx'],self.bestrans['dy'],self.bestrans['dx_err'],\
                             self.bestrans['dy_err'],self.bestrans['sigx'],self.bestrans['sigy'],
                             min_uncty))

        crval1_corr = (self.bestrans['dx']/3600.) / math.cos(self.crval2in * math.pi/180.)
        crval2_corr = (self.bestrans['dy']/3600.)

        self.CRVAL1new = self.crval1in + crval1_corr
        self.CRVAL2new = self.crval2in + crval2_corr

        self.logfile.write('CRVAL1new,CRVAL2new: '+str((self.CRVAL1new,self.CRVAL2new)))

        return


    def _runmatch(self,base,mrad,xsh,ysh):
        """Run match and parse result. """
        cmd = base+'matchrad='+str(mrad)+' xsh=%.2f ysh=%.2f' %(xsh,ysh)
        sproc  = popen2.Popen3(cmd,1)
        output = sproc.fromchild.readlines()
        errs   = sproc.childerr.readlines()
        sproc.fromchild.close()
        sproc.childerr.close()
     
        if errs and not (len(errs) == 1 and errs[0][0:27] == 'atFindMedtf: RMS <= 0.0, so'):
            self.logfile.write('match produced the following message(s): ')
            for f in errs:
                self.logfile.write(string.strip(f))
                self.errorList.append(('match',string.strip(f)))
        matchline = ''
        medianline = ''
        for line in output:
            if line[0:6] == 'TRANS:':
                matchline = line
            elif line[0:6] == 'MEDTF:':
                medianline = line
            else:
                self.logfile.write("ERROR:. Unexpected output from match program:\n\t\t\t"+line)
                return (-1,None)
        if (not matchline) or (not medianline):
            self.logfile.write("Match program failed.")
            return (-1,None)
        
        mdict = self._process_median(medianline,mrad,xsh,ysh)
        if not mdict:
            return (-1,None)
        return (0,mdict)


    def _process_median(self, medianline,mrad,xsh,ysh):
        """Parse, process and write out match MEDTF output line."""
        fields = string.split(medianline)
        mdict={}
        if len(fields) != 8 or fields[0] != 'MEDTF:':
            errtxt = "ERROR:  parsing match TRANS output: \n"+str(fields)
            self.logfile.write(errtxt)
            return None

        # read the fit parameters into variables
        mdx   = float(fields[1].split('=')[1])
        mdy   = float(fields[2].split('=')[1])
        adx   = float(fields[3].split('=')[1])
        ady   = float(fields[4].split('=')[1])
        sigx  = float(fields[5].split('=')[1])
        sigy  = float(fields[6].split('=')[1])
        Nmed  =   int(fields[7].split('=')[1])

        # check to see if it's a reasonable transformation
        if max(sigx,sigy) > 1.2:
            errtxt = 'Poor median trans: x,y rms = '+str(sigx)+' '+str(sigy)
            self.logfile.write(errtxt)
            self.errorList.append((self.modName,errtxt))
            return None

        # finally, populate MatchDict median keywords:
        mdict['Nmed']  = Nmed
        mdict['sigx']  = sigx
        mdict['sigy']  = sigy
        mdict['dx']    = mdx
        mdict['dy']    = mdy
        mdict['dx_err'] = max(sigx,_minSig_)/math.sqrt(Nmed-1.0)
        mdict['dy_err'] = max(sigy,_minSig_)/math.sqrt(Nmed-1.0)
        self.logfile.write(("rad,xsh,ysh: (%.1f,%.2f,%.2f)"%(mrad,xsh,ysh))+\
                           "  med x,y shifts for "+str(Nmed)+" objects; "+\
                           str(mdx)+" "+str(mdy)+"  rms: "+str(sigx)+" "+str(sigy))
        return mdict

    def makeObjcats(self):
        """ Make object catalogs for all input science images. """
        curdir = os.getcwd()
        os.chdir(self.astromdir)
        # make the sxtr paramfiles and run with them
        self._setupSxtrFiles()

        self.objcats = []
        self.projcats = []
        for imsci in self.sciImlist:
            rmsfile = imsci.split("_sci")[0]+'_RMS.fits'
            if rmsfile not in self.rmsImlist:
                errtxt = "Expected file %s not in rmsImlist!" %rmsfile
                self.logfile.write(errtxt)
                self.errorList.append((self.modName,errtxt))
                return -1
            tmpcat = imsci.split("_sci")[0]+'.tmp'
            cmd    = 'sex %s -c %s -CATALOG_NAME %s -WEIGHT_IMAGE %s' \
                     %(os.path.join(self.Imagedir,imsci), self.inparFile, \
                       tmpcat, os.path.join(self.Imagedir,rmsfile))
            self.logfile.write(' '+cmd)
            sproc   = popen2.Popen4(cmd, 1)
            errs    = sproc.fromchild.readlines()
            sproc.fromchild.close()

            if errs:
                self.logfile.write("Sxtr seems to have choked a bit:")
                for line in errs:
                    self.logfile.write(line)
                    self.errorList.append((self.modName,line))

            try:
                sxtlines = open(tmpcat).readlines()
            except:
                sxtlines = []
            del cmd,sproc,rmsfile,errs,tmpcat
            # refine that catalog, please
            objcat = imsci.split("_sci")[0]+'.obj'
            projcat = imsci.split("_sci")[0]+'.proj'
            fp = open(objcat,'w')
            ngood=0

            # the following seems very contrived, but helps in some cases
            if len(sxtlines)<250: bigMagerr = 0.4
            else: bigMagerr = 0.35
            self.logfile.write('cat: %s  Magerr_lim: %.2f' %(objcat,bigMagerr))
            for line in sxtlines:
                flds = line.split()
                if len(flds)<9 or flds[0][0] == '#':
                    continue
                mag    = float(flds[3])
                magerr = float(flds[4])
                ellip  = float(flds[5])
                fwhm   = float(flds[6])
                a_im   = float(flds[7])
                b_im   = float(flds[8])
                if(magerr>bigMagerr or ellip>0.7 or b_im < 0.6 or max(a_im,fwhm/2)>400):
                    continue
                fp.write(line)
                ngood += 1
            fp.close()
            if ngood>5:
                self.objcats.append(objcat)
            else:
                continue

            # finally, write the file with projected coord offsets
            cmd = 'project_coords %s 1 2 %.6f %.6f outfile=%s asec' \
                  %(objcat,self.RActr,self.Decctr,projcat)
            self.logfile.write("**>  "+cmd)
            projcoords = popen2.Popen4(cmd)
            _outlines = projcoords.fromchild.readlines()
            projcoords.fromchild.close()
            if len(_outlines) > 0:
                # I've never seen this happen, but better check...
                errtxt = "ERROR: project_coords mysteriously failed!"
                self.logfile.write(errtxt)
                self.errorList.append((self.modName,errtxt))
                for _line in _outlines:
                    self.logfile.write(_line)
                    print _line
                return -1
            self.projcats.append(projcat)
            self.logfile.write("Object astrometric catalog %s constructed."%projcat)
            del projcoords,cmd,objcat,projcat
            # ok, all done making catalogs for this image
            
	os.chdir(curdir)
        if len(self.objcats) < 1:
            self.logfile.write("makeObjCats: No catalogs made for matching?")
            return -1
        self.logfile.write("Made %d object catalogs for matching."%len(self.objcats))
        return 0

    def _setupSxtrFiles(self):
        "Make basic input files for running sextractor."
        # first the output parameter list file
        fp = open(self.paramFile,'w')
        fp.write("NUMBER \nALPHA_J2000 \nDELTA_J2000 \nMAG_ISO \nMAGERR_ISO\n"+\
                 "ELLIPTICITY \nFWHM_IMAGE \nA_IMAGE \nB_IMAGE \n"+\
                 "MAG_AUTO \nMAGERR_AUTO \nMAG_APER(1) \nMAGERR_APER")
        fp.close()
        # now for the input configuration file
        fp = open(self.inparFile,'w')
        fp.write("CATALOG_NAME    ztemp.cat\n")
        fp.write("CATALOG_TYPE    ASCII_HEAD\n")
        fp.write("PARAMETERS_NAME %s\n" %self.paramFile)
        fp.write("DETECT_TYPE     CCD\n")
        fp.write("DETECT_MINAREA  9\n")
        fp.write("DETECT_THRESH   6.5\n")
        fp.write("ANALYSIS_THRESH 6.5\n")
        fp.write("FILTER          N\n")
        fp.write("FILTER_NAME     $PIPELINE/configs/gauss_1.5_3x3.conv \n")
        fp.write("DEBLEND_NTHRESH 12\n")
        fp.write("DEBLEND_MINCONT 0.05\n")
        fp.write("CLEAN           Y\n")
        fp.write("CLEAN_PARAM     1.2\n")
        fp.write("MASK_TYPE       CORRECT\n")
        fp.write("PHOT_APERTURES  10\n")
        fp.write("PHOT_AUTOPARAMS 2.5,3.5\n")
        fp.write("SATUR_LEVEL     128000.0\n")
        fp.write("MAG_ZEROPOINT   0.0\n")
        fp.write("MAG_GAMMA       4.0\n")
        fp.write("GAIN            1.0\n")
        fp.write("PIXEL_SCALE     1.0\n")
        fp.write("SEEING_FWHM     2.0\n")
        fp.write("STARNNW_NAME    $PIPELINE/configs/default.nnw\n")
        fp.write("BACK_SIZE       128\n")
        fp.write("BACK_FILTERSIZE 3\n")
        fp.write("BACKPHOTO_TYPE  LOCAL\n")
        fp.write("BACKPHOTO_THICK 26\n")
        fp.write("CHECKIMAGE_TYPE NONE    \n")
        fp.write("CHECKIMAGE_NAME NONE\n")
        fp.write("MEMORY_OBJSTACK 10000\n")
        fp.write("MEMORY_PIXSTACK 2600000         \n")
        fp.write("MEMORY_BUFSIZE  2048\n")
        fp.write("VERBOSE_TYPE    QUIET\n")
        fp.write("WEIGHT_TYPE     MAP_RMS\n")
        fp.write("WEIGHT_IMAGE    rms.fits\n")
        fp.write("WEIGHT_THRESH   0,9.9e29\n")
        fp.write("INTERP_TYPE     NONE\n")


    def makeGSCcat(self,chopSpur=0):
        """ Constructs the GSC catalog file to be used in astrometric matching
        for this field.  Name of catalog stored as parameter 'self.GSCmatchin'.
        Returns 0 if successful, -1 if not.
        """
        curdir = os.getcwd()
        os.chdir(self.astromdir)
        self.logfile.write("Entered 'makeGSCcat' to make astrometric ref catalog...")
        # get input wcs information from first image in sci list
        ff = pyfits.open(os.path.join(self.Imagedir,self.sciImlist[0]))
        inwcs = wcsclass.BasicWCS(ff[0].header)
        ff.close()
        del ff
        
        self.crpix1,self.crpix2  = (inwcs.wcs['CRPIX1'],inwcs.wcs['CRPIX2'])
        self.crval1in = inwcs.wcs['CRVAL1']
        self.crval2in = inwcs.wcs['CRVAL2']
        NX=self.NX = inwcs.wcs['NAXIS1']
        NY=self.NY = inwcs.wcs['NAXIS2']

        # get RA,Dec of central pixel and the search radius
        Xctr = (self.NX+1)/2
        Yctr = (self.NY+1)/2
        self.RActr,self.Decctr = inwcs.xy2rd((Xctr,Yctr))     # both in deg
        
        ((rahh,ramm,rass), (decdd,decmm,decss)) = self._deg2hms(self.RActr,self.Decctr)
        rad_amin = 0.5*math.sqrt(NX*NX + NY*NY) * inwcs.wcs['PSCALE']/60.
        rad_amin = round(1000*(0.1+rad_amin))/1000.

        self.logfile.write('Input WCS: CRVAL1,2: '+str((self.crval1in,self.crval2in))+\
                           ' crpix1,2: '+str((self.crpix1,self.crpix2)))
        self.logfile.write('Xctr,Yctr: '+str((Xctr,Yctr))+'  RA_ctr,Dec_ctr: '+str((self.RActr,self.Decctr)))
        self.logfile.write('Making query to '+_HOST_+'   RA,Dec: '+\
                           str(((rahh,ramm,rass),(decdd,decmm,decss)))+' rad_amin = '+str(rad_amin))
        
        fh = cStringIO.StringIO()
        # the STScI query fails a lot, randomly. This loop actually helps sometimes.
        for iwqtry in range(MAXWQtry):
            try:
                self.webquery(host=_HOST_, url=_URL_, method="GET", file=fh,
                    RA = "%s" % self.RActr,
                    DEC = "%s" % self.Decctr,
                    SR=rad_amin/60.0,
                    CAT=_CATALOG_,
                    FORMAT="CSV")
            except IOError,err:
                errtxt = str(err)
                self.logfile.write("WARNING: webquery "+str(iwqtry+1)+" failed...\n ")
                self.logfile.write(errtxt)
                self.errorList.append((self.modName,errtxt))
                time.sleep(2)
                if (iwqtry<MAXWQtry-1):
                    sys.stderr.write(" trying again...\n")
                else:
                    sys.stderr.write(" sorry, no luck.\n")
                    print "Web qeury on ",_HOST_,"failed."
                    print errtxt
                    raise WebQueryError,err
                continue
            break
        # read and format the output, first line is a header so we will ignore it
        output = string.split(fh.getvalue(),'\n')[2:]
        fh.close()
        gsclines = [i.replace(',', ' ') for i in output if i != '']
        self.Nrefobjs = len(gsclines)
        #
        # AT THIS POINT WE NEED TO BAIL OUT IF THERE AREN'T ANY GSC OBJS FOUND!
        #
        if(self.Nrefobjs<6):
            errtxt = "WARNING: Too few (%d) GSC objects: no astrometric recalibration possible."%(self.Nrefobjs)
            self.logfile.write("WARNING: NOT ENUFF GSC OBJS TO CONTINUE!")
            self.logfile.write(errtxt)
            self.errorList.append((self.modName,errtxt))
            return -1
        
        # write stuff to data files
        fdump = open(os.path.split(_URL_)[1].split('.')[0]+'.scat', 'w')
        culledfile = os.path.split(_URL_)[1].split('.')[0]+'.cull'
        fcull = open(culledfile, 'w')
        self.GSCmatchin = os.path.split(_URL_)[1].split('.')[0]+'.proj'
        
        self.Ncull=0
        for line in gsclines:
            fdump.write(line+'\n')
            flds = line.split()
            # check Ra
            if flds[1][0] not in "1234567890.":
                continue
            # FpgMag, JpgMag, Vmag
            # FpgMagErr, JpgMagErr, VmagErr
            _mlist = [float(flds[7]),float(flds[8]),float(flds[12])]
            _elist = [float(flds[23]),float(flds[24]),float(flds[28])]
            _mlist.sort()
            _elist.sort()
            if _mlist[0] > 88.:  continue
            if chopSpur:
                if _mlist[1] > 88.:  continue
                
            mag = _mlist[0]
            magerr = _elist[0]
            del _mlist,_elist
            
            self.Ncull += 1
            # hstID, ra, dec, mag, raEpsilon, decEpsilon, magerr, epoch
            oline = '%-14s  %s  %s %6.2f   %s %s %5.2f %s\n' \
                    %(flds[0],flds[1],flds[2],mag,flds[4],flds[5],magerr,flds[6])
            fcull.write(oline)
        fdump.close()
        fcull.close()
        del fdump,fcull,line,oline,flds,mag,magerr

        self.logfile.write("Culling: kept %d GSC objects."%self.Ncull)
        print "Culling: kept %d GSC objects."%self.Ncull
        
        # finally, write the file with projected coord offsets
        cmd = 'project_coords %s 1 2 %.6f %.6f outfile=%s asec' \
              %(culledfile,self.RActr,self.Decctr,self.GSCmatchin)
        self.logfile.write("**>  "+cmd)

        projcoords = popen2.Popen4(cmd)
        _outlines = projcoords.fromchild.readlines()
        projcoords.fromchild.close()
        if len(_outlines) > 0:
            # I've never seen this happen, but better check...
            errtxt = "Error: project_coords mysteriously failed!"
            self.logfile.write(errtxt)
            self.errorList.append((self.modName,errtxt))
            for _line in _outlines:
                self.logfile.write(_line)
                print _line
            return -1
        del projcoords

        self.logfile.write("Astrometric reference catalog %s constructed."%self.GSCmatchin)
	os.chdir(curdir)
        return 0


    def _gscReform(self,lines):
        """Reformat GSC2 output"""
        newlines =[]
        Nobj=0
        if lines[0][:4] == "HTTP":
            # remove HTTP header
            for i in xrange(len(lines)):
                line = string.rstrip(lines[i])
                if line == '': break
            else:
                # no null line to end HTTP header -- just print it
                i = -1
            lines = lines[i+1:]
        hdr = []
        for line in lines[1:]:
            line = line.rstrip()
            if line and (line[:5] != "[EOD]"):
                fields = line.split(',')
                if len(fields) > 10:
                    newlines.append(line)
                    Nobj+=1
        if Nobj == 0:
            self.logfile.write("WARNING: No GSC objects returned!")
        self.logfile.write("GSC query returned %d objects."%Nobj)
        print "GSC query returned %d objects."%Nobj
        self.Nrefobjs = Nobj
        return newlines


    def _deg2hms(self, ra, dec):
        try:
            ra = float(ra)
            dec = float(dec)
        except ValueError:
            raise ValueError("deg2hrms: BAD input: %s %s" % (ra,dec))
        # convert to hms, dms format
        if dec > 90 or dec < -90:
            raise ValueError("deg2hrms: BAD input ra,dec: %s %s" % (ra,dec))
        ra = ra/15
        while ra>24:
            ra = ra-24
        while ra<0:
            ra = ra+24
        h = int(ra)
        m = int(60*(ra-h))
        s = 60*(60*(ra-h)-m)
        if s+0.0005 > 60:
            s = 0.0
            m = m+1
            if m==60:
                m = 0
                h = h+1
                if h==24:
                    h = 0
        h = "%02d" % h
        m = "%02d" % m
        s = "%06.3f" % s

        if dec<0:
            dsign = "-"
            dec = -dec
        else:
            dsign = "+"
        dd = int(dec)
        dm = int(60*(dec-dd))
        ds = 60*(60*(dec-dd)-dm)
        if ds+0.0005 > 60:
            ds = 0.0
            dm = dm+1
            if dm==60:
                dm = 0
                dd = dd+1
        dd = "%s%02d" % (dsign,dd)
        dm = "%02d" % dm
        ds = "%06.3f" % ds
        return ((h,m,s),(dd,dm,ds))



    def webquery(self, args=(), **kw):
        """Write output of a specified URL to stdout or file.
        Keywords for query may be specified as a sequence of pairs in
        args or as keywords.  Special keywords that define the URL include:
        host (default 'localhost'), url (default null), method (default 'POST'),
        and port (default 80).  The file keyword specifies an output filename
        or file handle (default sys.stdout).  Additional keywords are passed
        as parameters to the query.
            This method swiped in toto from Rick.
        """
    
        args = list(args)
        for key, value in kw.items():
        	args.append((key,value))
        port = 80
        method = "POST"
        url = ""
        host = urllib.localhost()
        outfile = sys.stdout
        query = []
        for key, value in args:
            if key == "port":
                port = int(value)
            elif key == "method":
                method = value.upper()
            elif key == "url":
                url = value
            elif key == "host":
                host = value
            elif key == "file":
                outfile = value
            elif value is None:
                query.append(urllib.quote(key))
            else:
                query.append('%s=%s' % (urllib.quote(key),urllib.quote_plus(str(value))))
        query = '&'.join(query)
        if isinstance(outfile, types.StringType):
            outfile = open(outfile,"w")
        if url[:1] == "/":
            # don't add an extra slash (purely for aesthetic purposes)
            url = "http://%s:%d%s" % (host,port,url)
        else:
            url = "http://%s:%d/%s" % (host,port,url)
        if not query:
            query = None
        elif method == "GET":
            url = "%s?%s" % (url,query)
            query = None
        inurl = urllib.urlopen(url,query)
        print url,query
        s = inurl.read(102400)
        while s:
            outfile.write(s)
            s = inurl.read(102400)


    def writeXml(self):
        """This module must remark the image headers with WFP xml.
        """
        curdir = os.getcwd()
        os.chdir(self.Imagedir)
        allImageLists = [self.sciImlist, self.ctxImlist, self.wgtImlist, self.rmsImlist]
        
        for imlist in allImageLists:
            for im in imlist:
                file = xmlUtil.markupImage(im,dataset=self.obsName)
                
                # Don't write these images as output of this module, which
                # really doesn't have any.
                
                #if file not in self.outputList.keys():
                #   self.outputList[file] = [im]
                
        os.chdir(curdir)
        return

    def mkMsg(self):
        """Create and write module level message for this class.
        Most of this is just compiling the info. meta is a dictionary
        of lists where each list is a list of tuples describing the
        tag lines for the particular section of the message.  This tuple 
        format conforms to that used by the xmlMessage class which is
        modeled on basic python argument passing, i.e. (key,*value,**attr).
        """
        # getting the version of project_coords
        project_coords_cmd = 'project_coords --version'
        outp = popen2.Popen4(project_coords_cmd)
        outpline = outp.fromchild.readlines()
        pcoorVer = outpline[0].split()[-1]
        
        self.meta = {}
        self.meta['module']= []
        self.meta['meta']  = []
        self.meta['input'] = []
        self.meta['output']= []
        self.meta['errorlist'] = []
    
        self.meta['module'].append(('module','name='+self.modName,'version='+__version__,'dataset='+self.obsName))
        self.meta['module'].append(('root',self.root))
        self.meta['meta'].append(('meta',))
        self.meta['meta'].append(('depend',))
        self.meta['meta'].append(('pkg',))
        self.meta['meta'].append(('name','python'))
        self.meta['meta'].append(('version',pyversion.split()[0]))
        self.meta['meta'].append(('pkg',))
        self.meta['meta'].append(('name','pyfits'))
        self.meta['meta'].append(('version',pyfits.__version__.split()[0]))
        self.meta['meta'].append(('pkg',))
        self.meta['meta'].append(('name','project_coords'))
        self.meta['meta'].append(('version',pcoorVer))
        self.meta['meta'].append(('pkg',))
        self.meta['meta'].append(('name','Guide Star Catalog'))
        self.meta['meta'].append(('version',_URL_.split("/")[-1].split("q")[0]))
    
        # SExtractor info
        sub  = subprocess.Popen(['sex', '--version'], stdin=subprocess.PIPE, stdout=subprocess.PIPE, close_fds=True)
        outp = sub.stdout.readlines()
        name = outp[0].split()[0]
        ver  = outp[0].split()[2]
        self.meta['meta'].append(('pkg',))
        self.meta['meta'].append(('name',name))
        self.meta['meta'].append(('version',ver))
        cmdline1 = 'sex fitsfile -c self.InParFileName'
        self.meta['meta'].append(('commandline',cmdline1))
        del outp,sub,name,ver
    
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

