#!/usr/bin/env python
# $Id: combFilter.py,v 2.20 2011/09/12 20:15:45 zheng Exp $
# ---------------------------------------------------------------------

__version__      = '$Revision: 2.20 $ '[11:-3]
__version_date__ = '$Date: 2011/09/12 $ '[7:-3]
__author__       = 'Wei Zheng'

import os, string
import pdb
import popen2
import pUtil,xmlUtil,fUtil
import pyfits
import numpy
import subprocess
from   pyraf import iraf
from   msg   import pMessage
from   sys   import version

pyversion = version
_delFiltsDef = ['G800L']

class detectionImageError(Exception):
    """Exception raised in making detection image."""
    def __init__(self,value):
        self.value = value
    def __repr__(self):
        return repr(self.value)


class detectionImage:
    """ Class provides a detectionImage object with methods to create an
    optimally combined (chi_square?) image from an observation's fits files.
    The detection image is written to the observations Images directory.
    """

    def __init__(self,obs,Band=None,excludeFilters=None,noContext=None):
        """ constructor gets an Observation object (obs) from the pipeline.  """
        self.modName          = string.split(string.split(str(self))[0],'.')[0][1:]
        self.root             = obs.newobspath                      # root path of the observation dir
        self.obsName          = obs.newobs
        self.obsPars          = obs.newpardir
        self.obsFits          = obs.fitsdir                         # path to fits files
        self.messagedir       = obs.messagedir                      # where the module message will go.
        self.configdir        = obs.configdir
        self.logfile          = obs.logfile       
        self.defaultInParFile = os.path.join(obs.pardir,self.modName + '_'+ Band + '.inpar') #WZ
        self.defaultParamFile = os.path.join(obs.pardir,self.modName + '.param')
        self.inParFileName    = os.path.join(self.obsPars, 'noisy.inpar')
        self.inpars           = pUtil.readInParFile(self.defaultInParFile)      # readInParFile returns a dict  
        # obs.dethresh          =self.inpar['DETECT_THRESH']
        # obs.anathresh         =self.inpar['ANALYSIS_THRESH']
        # self.detImName        = 'detectionImage.fits'
        # self.detWgtName       = 'detectionWeight.fits'
        self.edgeMaskName     = obs.edgeMaskName = 'Edgemask.fits' #XX
        self.errorList        = []
        self.inputList        = []
        self.outputList       = {}
        if Band=='nir': #WZ
            self.detImName        = 'detection_nir.fits'
            self.detWgtName       = 'detection_nir_wht.fits'
            Detector='IR'
        elif Band=='red':   # only choose the last infared filter F160W F140W  by XX
            self.detImName        = 'detection_red.fits'
            self.detWgtName       = 'detection_red_wht.fits'
            Detector='IR'
        elif Band=='opt':
            self.detImName        = 'detection_opt.fits'
            self.detWgtName       = 'detection_opt_wht.fits'
            Detector='WFC'
        elif Band=='uv':
            self.detImName        = 'detection_uv.fits'
            self.detWgtName       = 'detection_uv_wht.fits'
            Detector='UVIS'
        else:
            self.logfile.write("Incorrect band name: "+Band)
            exit()
        self.logfile.write(self.modName+": Instantiated detectionImage object.")
        # excludeFilters must be a list of filter names
        # i.e. --excludefilt=['F850LP','F225W']
        # in a pipeline context, this form will be passed to the constructor.
    
        self.excludeFilters = excludeFilters
        self.det = Detector
        if noContext:
            checkLists = [obs.sciImageList,obs.weightImageList,obs.flagImageList,obs.rmsImageList]
        else:
            checkLists = [obs.sciImageList,obs.contextImageList,obs.weightImageList,obs.flagImageList,obs.rmsImageList]
         
        # remove drizzled grism images and other oddities from various lists and further processing
        self.logfile.write("Checking for input images not with detector: "+str(Detector))
        for _obsList in checkLists:
            for ii in range(len(_obsList)-1,-1,-1):
                im = _obsList[ii]
                imfilt = fUtil.getFilter(self.obsFits+'/'+im)
                if imfilt == _delFiltsDef:
                    self.logfile.write(im+" ("+imfilt+") removed from further obs processing.")
                    del _obsList[ii]   # now it belongs to the ages
        
	#pdb.set_trace() #XX
        for _obsList in checkLists:
            if len(_obsList) != len(obs.sciImageList):
                errtxt = "Error: inconsistency in Image Lists after removing unwanted detector images."
                self.logfile.write(errtxt)
                raise Exception,errtxt

        self.sciImageList = obs.sciImageList
        for i in self.sciImageList:
            self.inputList.append(i)

        if len(self.sciImageList) == 0:
            raise detectionImageError("Image list has no images to make detection image.")
        
    def getinpar(self, param):
        """Returns value of input parameter given a parameter name"""
        return self.inpars[param]

    def setinpar(self, param, value):
        """Sets value of input parameter given a parameter name and value"""
        self.inpars[param] = value
        return
    
    def setpars(self):
        """ This method in this class generates a phoney (unimportant) .inpar file
        to run sextractor and get out the background and noise stats. """

        self.logfile.write('Generating parameter sets for noise and background measurements...')

        # flag       = self.getinpar('FLAG_IMAGE')
        filter     = self.getinpar('FILTER_NAME')
        starnnw    = self.getinpar('STARNNW_NAME')

        self.setinpar('PARAMETERS_NAME', self.defaultParamFile)
        self.setinpar('FILTER_NAME'    , os.path.join(self.configdir, filter ))     
        # self.setinpar('FLAG_IMAGE'     , os.path.join(self.configdir, flag   ))
        self.setinpar('STARNNW_NAME'   , os.path.join(self.configdir, starnnw))

        # Following parameters are for the phoney .inpar file. 

        self.setinpar('CATALOG_NAME'   , '/dev/null')
        self.setinpar('CATALOG_TYPE'   , 'NONE')
        self.setinpar('DETECT_THRESH'  , '100000')
        self.setinpar('CHECKIMAGE_TYPE', 'NONE')
        self.setinpar('CHECKIMAGE_NAME', 'NONE')
        self.setinpar('VERBOSE_TYPE'   , 'NORMAL')
            # setinpar gets us the noise and background numbers

        inparfile  = open(self.inParFileName, 'w')
        
        for param, value in self.inpars.items():
            inparfile.write("%-18s   %-20s\n" % (param, value))            
        inparfile.close()

        self.logfile.write('phoney parameter file '+os.path.basename(self.inParFileName)+' written.')
        
        return

    def getStats(self,Band):
        """Runs SExtractor with phoney pars to get noise and background stats. 
       Now implement this using the Popen4 class, which sends both stdout
       and stderr to the fromchild file object.  We suppress stdout from 
       SExtractor by shunting stuff to /dev/null via the par file, so we
       will only get stderr anyway, at least, in theory.  The Popen3 class
       or rather, the factory function is uses, popen3, has proved
       to be buggy at times, so this is an attempt to move away from it.
        """ 

        self.logfile.write('Generating background and noise statistics for fits images')
        self.statsList = []

        # pdb.set_trace()
        for fitsfile in self.sciImageList:
            # Use user-defined detector for the detection image
            det = fUtil.getDetector(os.path.join(self.obsFits,fitsfile))
            if det == self.det:
                if (Band =='red'):     
                    if string.find(fitsfile,'160') > -1 or string.find(fitsfile,'140') > -1:
                        yes = 1
                    else:
                        self.logfile.write('Not including '+fitsfile+' in detection image.')
                        yes = 0
                else:
                    yes = 1
            else:
                self.logfile.write('Not including '+fitsfile+' in detection image.')
                yes = 0

            if self.excludeFilters:
                # pass on user excluded filters for the detection image
                filt = fUtil.getFilter2(os.path.join(self.obsFits,fitsfile))
                if filt in self.excludeFilters:
                    self.logfile.write('User requested exclude filter found in data.')
                    self.logfile.write('Not including'+filt+' in detection image.')
                    yes = 0

            if (yes == 1):
                self.logfile.write('running SExtractor on ' + fitsfile)
                cmd     = 'sex ' +  os.path.join(self.obsFits, fitsfile) + ' -c ' + self.inParFileName
                subproc = popen2.Popen4(cmd)
                lines   = subproc.fromchild.readlines()
                self.logfile.write('logging background and noise stats.')
        
                for line in lines:
                    if string.find(line,'(M+D) Background:') != -1:
                        # pdb.set_trace()
                        linef = string.split(line)
                        self.logfile.write(fitsfile+': '+linef[1]+linef[2]+'  '+linef[3]+linef[4])
                        self.statsList.append((fitsfile,linef[2],linef[4]))          # linef[2] is background, linef[4] is noise.    ############ attention here by XX
                    elif string.find(line,'*ERROR*') != -1:
                        self.errorList.append((self.modName,line))
                        self.logfile.write('Error in running SExtractor: '+ line)
        return 


    def mkIrafIm(self):
        """ make the detection image.  This method creates the detection 
        Image using the CHI^2 algorithm.
        """ 

        self.logfile.write('Generating detection Image...')

        cwd = os.getcwd()
        os.chdir(self.obsFits)                  # cd to the observations Images dir
        self.logfile.write('cd to observation Images dir ' + self.obsFits)

        for i in range(len(self.statsList)):
            image,background,noise = self.statsList[i]
            expr="((a - %s)/%s)**2" % (background,noise)
            self.logfile.write('Generating background subtracted, squared image: '+'subtracted_'+image)
            iraf.imexpr(expr,'subtracted_'+image,a=image)

        inputString = string.join(map(lambda x: 'subtracted_'+x,self.sciImageList),',')

        iraf.unlearn(iraf.imsum)
        iraf.imsum(input=inputString, output='dummy.fits', option='sum')
        iraf.unlearn(iraf.imexpr)

        self.logfile.write('Generating detectionImage.fits for observation...')
        iraf.imexpr('sqrt(a)',self.detImName,a="dummy.fits")

        self.logfile.write('Removing dummy image...')
        iraf.imdelete('dummy.fits',verify='no')
        self.logfile.write('Removing subtraction images...')
        iraf.imdelete(inputString ,verify='no')

        # the header of the image produced by iraf is junk. fix it up
    
        self._fixIrafHeader()
        os.chdir(cwd)                           # return to orig dir.
        return


    def mkChi2Im(self):
        """ 
        Make the detection image using Numeric package, not IRAF.
        This method creates the detection Image using the CHI^2 algorithm.
        """ 
        cwd = os.getcwd()
        os.chdir(self.obsFits)            # cd to the observations Images dir
        self.logfile.write('cd to '+self.obsFits+' to generate detection image...')

        # some initializations for the big loop
        firstpass = 1
        exptime = 0
        pred_images = []
        for im in self.statsList:
            imName   = im[0]
            sky      = float(im[1])
            noise    = float(im[2])
            ffi      = pyfits.open(im[0])
            bitpix   = ffi[0].header['BITPIX']
            naxis1   = ffi[0].header['NAXIS1']
            naxis2   = ffi[0].header['NAXIS2']
            ######## do a few tests for sanity
            if noise < 1e-30:
                self.errorList.append((self.modName,im+': Zero noise!'))
                raise ValueError, 'Zero noise!'
            print imName,'sky,noise =',sky,noise
            #### warn if it's not a 32-bit float
            if bitpix != -32:
                warnTXT = 'WARNING: BITPIX of image '+str(imName)+' not -32.'
                self.errorList.append((self.modName,warnTXT))
                self.logfile.write(warnTXT)
            #### increment exposure times
            if not ffi[0].header.has_key('EXPTIME'):
                self.errorList.append((self.modName,im+': Image header lacks EXPTIME key!'))
                raise KeyError,'Image header lacks EXPTIME key!'
            exptime += ffi[0].header['EXPTIME']
            #########################################
            pred_images.append(imName)  # This is for the predecessor images list 
            if firstpass:
                firstpass=0
                ### open det image and give it a header from the first image
                self.detFF = pyfits.HDUList()
                self.detFF.append(pyfits.PrimaryHDU())
                try:
                    del self.detFF[0].header.ascard["EXTEND"]
                except KeyError:
                    pass
                #self.detFF = fits.FITS(self.detImName,'w')
                self.detFF[0].header = ffi[0].header
                NX = naxis1
                NY = naxis2

                self.logfile.write('Generating chi^2 detectionImage.fits for observation...')
                #### initialize the data array of detection image
                if ffi[0].header['BITPIX'] != -32:
                    ffi[0].header['BITPIX'] = -32
                xmask = numpy.not_equal(ffi[0].data,0).astype(numpy.float32)
                x = (ffi[0].data  - sky)/noise * xmask #old WZ
                x = x*x
                x = ffi[0].data
                #self.detFF[0].data = x.astype(numpy.float32)
                #Sum  = x - x
                Norm = Sum
                print "Norm: ", Norm, imName
                err=ffe[0].data                
                Var = err * err
                #var = err * err
                Wht = 1./Var
                Sum = ffi[0].data * Wht
                #for i in range(ffi[0].shape[0]):
                #    for j in range(ffi[0].shape[1]):
                #        w = 1./err[i,j]/err[i,j]
                #        wht[i,j] += w
                #        sum[i,j] += (ffi[0].data[i,j] * w)
            else:
                # check sizes of subsequent images
                if(naxis1 != NX or naxis2 != NY):
                    warnTXT = 'WARNING: Image '+'imName'+ \
                              ' is a different size than the detection Image.'
                    self.errorList.append((self.modName,warnTXT))
                    self.logfile.write(warnTXT)
                    raise ValueError,warnTXT

                # NOTE! do not use the "**" operator in the next step, or the
                # astype gets changed from float32 to Float64 

                # jpb 11/Jan -- just a temporary fix, until we do it right
                xmask = numpy.not_equal(ffi[0].data,0).astype(numpy.float32)
                x = ffi[0].data 
                x = (ffi[0].data - sky)/noise * xmask
                self.detFF[0].data += x*x
            ffi.close()
            del ffi
        # end of loop, detection image made -- take sqr root
        # update exposure time -- Drizzle doesn't update TEXPTIME, but we can
        #self.detFF[0].data = numpy.sqrt(self.detFF[0].data)
        self.detFF[0].header.update('EXPTIME', exptime)
        self.detFF[0].header.update('TEXPTIME',exptime)
        self.detFF[0].header.update('ORIGIN','Properietary data by CLASH pipeline:'+ostring) #WZ
        # sew it up
        self.detFF.writeto(self.detImName)
        del self.detFF

        self.logfile.write('detectionImage.fits created; fixing header...')
        
        # the header of the image needs to be fixed up
        self._fixHeader()
        self.logfile.write('done. Returning to original directory.')
        self.outputList[self.detImName] = pred_images
        os.chdir(cwd)      # ride off into sunset
        
        return


    def mkInvarIm(self):
        """ 
        Make the detection image using Numeric package, not IRAF.
        This method creates a detection image using the inverse variance
        algorithm.
        """ 
        cwd = os.getcwd()
        os.chdir(self.obsFits)            # cd to the observations Images dir
        self.logfile.write('cd to '+self.obsFits+' to generate detection image...')

        # some initializations for the big loop
        firstpass = 1
        exptime = 0
        pred_images = []
        print "Checking detection Image"
        # pdb.set_trace()
        for im in self.statsList:
            imName   = im[0]
            sky      = float(im[1])
            noise    = float(im[2])
            ffi      = pyfits.open(im[0])
            bitpix   = ffi[0].header['BITPIX']
            naxis1   = ffi[0].header['NAXIS1']
            naxis2   = ffi[0].header['NAXIS2']
            cenwave  = float(ffi[0].header['PHOTPLAM'])/1e4  #WZ
            #pdb.set_trace() #XX
 	    print "Make detection image using the files like :  _drz.fits  _RMS.fits  _drz_context.fits by xingxing"
            #ename = string.split(imName,'_sci')[0]+'_RMS.fits'    #by xingxing
	    ename = string.split(imName,'_drz')[0]+'_RMS.fits'
            ffe   = pyfits.open(ename)
            #cname = string.split(imName,'_sci')[0]+'_sci_context.fits'      #by xingxing
	    cname = string.split(imName,'_drz')[0]+'_drz_context.fits'
            ffc   = pyfits.open(cname)
            ######## do a few tests for sanity
            if noise < 1e-30:
                self.errorList.append((self.modName,im+': Zero noise!'))
                raise ValueError, 'Zero noise!'
            print imName,'sky,noise =',sky,noise
            #### warn if it's not a 32-bit float
            if bitpix != -32:
                warnTXT = 'WARNING: BITPIX of image '+str(imName)+' not -32.'
                self.errorList.append((self.modName,warnTXT))
                self.logfile.write(warnTXT)
            #### increment exposure times
            if not ffi[0].header.has_key('EXPTIME'):
                self.errorList.append((self.modName,im+': Image header lacks EXPTIME key!'))
                raise KeyError,'Image header lacks EXPTIME key!'
            exptime += ffi[0].header['EXPTIME']
            #########################################
            pred_images.append(imName)  # This is for the predecessor images list
            detector=ffi[0].header['DETECTOR']
            cd11=ffi[0].header['CD1_1']
            cd12=ffi[0].header['CD1_2']
            pscale= numpy.sqrt(cd11*cd11+cd12*cd12)*3600.
            if firstpass:
                firstpass=0
                ### open det image and give it a header from the first image
                self.detFF = pyfits.HDUList()
                self.detFF.append(pyfits.PrimaryHDU())
                self.wgtFF = pyfits.HDUList()
                self.wgtFF.append(pyfits.PrimaryHDU())
                try:
                    del self.detFF[0].header.ascard["EXTEND"]
                except KeyError:
                    pass
                #self.detFF = fits.FITS(self.detImName,'w')
                self.detFF[0].header = ffi[0].header
                NX = naxis1
                NY = naxis2

                self.logfile.write('Generating an inverse variance detectionImage.fits for observation...')
                #### initialize the data array of detection image
                if ffi[0].header['BITPIX'] != -32:
                    ffi[0].header['BITPIX'] = -32
                sens0=ffi[0].header['PHOTFLAM']
                # pdb.set_trace()
                if (detector == 'WFC'):
                    wtf= 0. # WZ even weighting
                    ratio = numpy.square(pscale/0.05) # 1.69
                else:
                    if detector == 'IR':
                        wtf= 0. # WZ even weighting
                        # wtf=2.0 # square-wv weighting
                        ratio = numpy.square(pscale/0.128) # 0.26 
                    else:
                        if detector == 'UVIS':
                            wtf= 0. # WZ even weighting            #added by xingxing
                            ratio = numpy.square(pscale/0.0396) # 2.69
		# pdb.set_trace()##
                scale= ratio * numpy.power(cenwave,wtf)
                print im[0],' with scale of ',scale
                xmask = numpy.not_equal(ffi[0].data,0).astype(numpy.float32)
                #x     = ((ffi[0].data - sky)/noise/noise) * xmask
                #self.detFF[0].data = x.astype(numpy.float32)
                #weight = 1.0/noise/noise
                #del x
                #x = (ffi[0].data  - sky)/noise * xmask #old WZ
                #x = x*x
                x = ffi[0].data   #WZ
                self.detFF[0].data = (x-x).astype(numpy.float32)
                Sum  = x - x
                Wht = x - x
                # wht=ffc[0].data
                err=ffe[0].data
                Var = err * err
                #var = err * err
                # Wht = wht * scale
                Wht = scale / Var #WZ
                Sum = ffi[0].data * Wht
                print "",'Count: ',ffi[0].data[2921,2834]," error: ",err[2921,2834],' Weight: ', Wht[2921,2834],' Sum: ', Sum[2921,2834]
                #for i in range(ffi[0].shape[0]):
                #    for j in range(ffi[0].shape[1]):
                #        w = 1./err[i,j]/err[i,j]
                #        wht[i,j] += w
                #        sum[i,j] += (ffi[0].data[i,j] * w)
            else:
                # check sizes of subsequent images
                if(naxis1 != NX or naxis2 != NY):
                    warnTXT = 'WARNING: Image '+'imName'+ \
                              ' is a different size than the detection Image.'
                    self.errorList.append((self.modName,warnTXT))
                    self.logfile.write(warnTXT)
                    raise ValueError,warnTXT

                # NOTE! do not use the "**" operator in the next step, or the
                # astype gets changed from float32 to Float64 

                sens=ffi[0].header['PHOTFLAM']
                detector=ffi[0].header['DETECTOR']
                # pdb.set_trace()
                if (detector == 'WFC'):
                        wtf=0.0 # even weighting
                        ratio = numpy.square(pscale/0.05) # 1.69
                else:
                    if detector == 'IR':
                        wtf=0 # square-wv weighting
                        ratio = numpy.square(pscale/0.128) # 0.26 
                    else:
                        if detector == 'UVIS':
                            wtf= 0
                            ratio = numpy.square(pscale/0.0396) # 2.69
                scale= sens0 / sens * ratio * numpy.power(cenwave,wtf)
                print im[0],' with scale of ',scale
                xmask = numpy.not_equal(ffi[0].data,0).astype(numpy.float32)
                #x = ((ffi[0].data - sky)/noise/noise) * xmask
                #self.detFF[0].data += x.astype(numpy.float32)
                #weight += 1.0/noise/noise
                #del x
                #Sum  = x - x
                #wht=ffc[0].data
                err=ffe[0].data                
                var = err * err
                wht = scale / var
                # wht = wht * scale
                Sum += ffi[0].data * wht
                Wht += wht
                print "",'Count: ',ffi[0].data[2921,2834]," error: ",err[2921,2834],' Weight: ', wht[2921,2834],' Sum: ', Sum[2921,2834],' Weight: ',Wht[2921,2834]
                #for i in range(ffi[0].shape[0]):
                #    for j in range(ffi[0].shape[1]):
                #        w = 1./err[i,j]/err[i,j]
                #        wht[i,j] += w
                #        sum[i,j] += (ffi[0].data[i,j] * w)
        Norm = numpy.sum(Wht)
        print "Norm: ", Norm
        mx=numpy.max(Wht)
        mn=numpy.min(Wht)
        av=numpy.mean(Wht)
        sd=numpy.std(Wht)
        Sum /= Wht
        Var = 1./Wht
        #pdb.set_trace()
        #Err = numpy.sqrt(Var)
        # detWt  = Wht / Norm
        self.detFF[0].data = Sum
        self.wgtFF[0].data = Wht / Norm
        #detFF = sum                
        ffi.close()
        del ffi
        ffe.close()
        del ffe
        # end of loop, detection image made, divide by weight.
        #self.detFF[0].data = self.detFF[0].data/weight

        #The detection image created by using inverse_variance should have 
        #zero background and a gaussian noise with rms=1. 
        #We apply to it algorithms as the FDR 
        
        # update exposure time -- Drizzle doesn't update TEXPTIME, but we can
        self.detFF[0].header.update('EXPTIME', exptime)
        self.detFF[0].header.update('TEXPTIME',exptime)
        # sew it up
        #pdb.set_trace()
        self.detFF.writeto(self.detImName)
        del self.detFF
        self.wgtFF.writeto(self.detWgtName) #WZ
        del self.wgtFF
        self.logfile.write('detectionImage.fits created; fixing header...')
        
        # the header of the image needs to be fixed up
        #pdb.set_trace()
        # self._fixHeader()
        self.logfile.write('done. Returning to original directory.')
        self.outputList[self.detImName] = pred_images
        os.chdir(cwd)      # ride off into sunset
        
        return


    ############## function to make weight image for detection image ###############
    def mkWghtIm(self, doMedFilt=0, medfiltsize=5):
        """ Makes the weight image to be used by Sextractor when finding
        sources in the (inverse variance-weighted) detection image.
        The method for constructing the weight image follows analogously
        to that used in making the (inverse-variance) detection image.
        """ 
        cwd = os.getcwd()
        os.chdir(self.obsFits)            # cd to the observations Images dir
        self.logfile.write('cd to '+self.obsFits+' to generate detection image...')
        # some initializations for the big loop
        firstpass = 1
        exptime = 0
        pred_images = []
        #pdb.set_trace()
        for im in self.statsList:
            wgtName  = im[0][:-5]+'_weight.fits'
            # we don't use 'sky' . . . safest not to define
            noise    = float(im[2])
            
            # now open the drizzled weight image for this filter
            ffi      = pyfits.open(wgtName)
            bitpix   = ffi[0].header['BITPIX']
            naxis1   = ffi[0].header['NAXIS1']
            naxis2   = ffi[0].header['NAXIS2']
            ######## do a few tests for sanity
            if noise < 1e-30:
                self.errorList.append((self.modName,im+': Zero noise!'))
                raise ValueError, 'Zero noise!'
            print wgtName,'noise =',noise
            #### warn if it's not a 32-bit float
            if bitpix != -32:
                warnTXT = 'WARNING: BITPIX of image '+wgtName+' not -32.'
                self.errorList.append((self.modName,warnTXT))
                self.logfile.write(warnTXT)
            #### increment exposure times
            if not ffi[0].header.has_key('EXPTIME'):
                self.errorList.append((self.modName,wgtName+': Image header lacks EXPTIME key!'))
                raise KeyError, wgtName+': Image header lacks EXPTIME key!'
            exptime += ffi[0].header['EXPTIME']
            #########################################
            pred_images.append(wgtName)  # This is for the predecessor images list 
            if firstpass:
                firstpass=0
                NX = naxis1
                NY = naxis2
                ### open det image and give it a header from the first image
                wgtFF = pyfits.HDUList()
                wgtFF.append(pyfits.PrimaryHDU())
                try:
                    del wgtFF[0].header.ascard["EXTEND"]
                except KeyError:
                    pass
                #wgtFF = fits.FITS(self.detWgtName,'w')
                wgtFF[0].header = ffi[0].header
                if wgtFF[0].header['BITPIX'] != -32:
                    wgtFF[0].header['BITPIX'] = -32
                self.logfile.write('Generating '+self.detWgtName+' for observation...')

                # initialize the data array of detection image
                # weight the weights according to inverse var...
                x = (ffi[0].data)/noise/noise
                wgtFF[0].data = x.astype(numpy.float32)
                weightNorm = 1.0/noise/noise
                del x
            else:
                # check sizes of subsequent images
                if(naxis1 != NX or naxis2 != NY):
                    warnTXT = 'WARNING: Image '+wgtName+ \
                              ' is a different size than the detection Image.'
                    self.errorList.append((self.modName,warnTXT))
                    self.logfile.write(warnTXT)
                    raise ValueError,warnTXT

                x = (ffi[0].data)/noise/noise
                wgtFF[0].data += x.astype(numpy.float32)
                weightNorm += 1.0/noise/noise
                del x
            ffi.close()
            del ffi,noise,wgtName,bitpix,naxis1,naxis2
        # end of loop, detection Weight image made, but divide by weights.
        wgtFF[0].data = wgtFF[0].data/weightNorm
        
        # and now multiply by edge mask, if available
        if os.path.isfile(self.edgeMaskName):
            self.logfile.write('Multiplying detection wgt image by edgemask.')
            _edgeFF = pyfits.open(self.edgeMaskName)
            if _edgeFF[0].header['NAXIS1'] != NX or _edgeFF[0].header['NAXIS2'] != NY:
                self.logfile.write('ERROR: '+self.edgeMaskName+' size different from detection image?!')
                self.logfile.write(' skipping edge mask multiplication!')
            else:
                wgtFF[0].data = wgtFF[0].data * _edgeFF[0].data
            _edgeFF.close()
            del _edgeFF
        else:
            self.logfile.write('No edgemask found for detection wgt image.')
        
        # update exposure time -- Drizzle doesn't update TEXPTIME, but we can
        wgtFF[0].header.update('EXPTIME', exptime)
        wgtFF[0].header.update('TEXPTIME',exptime)
        # sew it up
        #wgtFF.writeto(self.detWgtName) #WZ
        self.logfile.write('detectionWeight.fits created; fixing header...')

        # the header of the image needs to be fixed up
        self._fixHeader(fitsfile=self.detWgtName)
        self.logfile.write('done. Returning to original directory.')
        self.outputList[self.detWgtName] = pred_images

        if doMedFilt:
            self.logfile.write('median filtering detection weight image with filtersize = '+str(medfiltsize))
            self._medFilterMaskedWgtIm(self.detWgtName, medfiltsize=medfiltsize)
        os.chdir(cwd)      # ride off into sunset
        del wgtFF
        return


    def _medFilterMaskedWgtIm(self, imname, medfiltsize=5):
        """ Method to median filter an image, intended to be the
        detection weight image.  Default filter size is 5 pix.
        It first makes a mask from the input image so that no input
        zero pixels (which are specially flagged as having no weight)
        become nonzero after filtering.
        """
        tmpMed  = '_tmpMED.fits'
        tmpMask = '_tmpMask.fits'
        iraf.flpr()
        iraf.flpr(iraf.median)
        iraf.flpr(iraf.imcalc)

        iraf.unlearn(iraf.median)
        iraf.median.input   = imname
        iraf.median.output  = tmpMed
        iraf.median.xwindow = medfiltsize
        iraf.median.ywindow = medfiltsize
        iraf.median.mode = 'h'
        iraf.median()

        iraf.unlearn(iraf.imcalc)
        iraf.imcalc(imname, tmpMask, "if im1 .eq. 0 then 0 else 1")
        _instring = tmpMed+','+tmpMask
        iraf.imcalc(_instring, tmpMed, "im1 * im2")

        os.rename(tmpMed,imname)
        os.remove(tmpMask)
        iraf.flpr(iraf.median)
        iraf.flpr(iraf.imcalc)
        iraf.flpr()

        return


    def writeXml(self):
        """Mark up the detection Image produced by this module."""
        cwd = os.getcwd()
        os.chdir(self.obsFits)
        for im in [self.detImName,self.detWgtName]:
            ofile = xmlUtil.markupImage(im,dataset=self.obsName)
            self.outputList[ofile] = [im]
        os.chdir(cwd)                                           # return to orig dir.
        return

    def mkMsg(self):
        """create and write module level message for this class.
        Most of this is just compiling the info. meta is a dictionary
        of lists where each list is a list of tuples describing the
        tag lines for the particular section of the message.  This tuple 
        format conforms to that used by the xmlMessage class which is
        modeled on basic python argument passing, i.e. (key,*value,**attr).
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
        self.meta['meta'].append(('depend',))
        self.meta['meta'].append(('pkg',))
        self.meta['meta'].append(('name','python'))
        self.meta['meta'].append(('version',pyversion.split()[0]))
        self.meta['meta'].append(('pkg',))
        self.meta['meta'].append(('name','pyfits'))
        self.meta['meta'].append(('version',pyfits.__version__.split()[0]))
        self.meta['meta'].append(('pkg',))
        self.meta['meta'].append(('name','numpy'))
        self.meta['meta'].append(('version',numpy.__version__))

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

    def _fixHeader(self, fitsfile=None):
        """  Fix up the header of the detection Image (default).
            optional parameter 'fitsfile' can be used to specify some other
            image (notably, the detectionWeight image) for header fix.
            """
        # the list of keywords that will be copied out of what the drizzle
        # task produced.  ASCards for these keywords copied in this order.
        # NOTE: I'm throwing EXPTIME and TEXPTIME into this lot (jpb, 4/Oct/O1).
        if fitsfile == None:
            fitsfile = self.detImName

        keylist = [ 'SIMPLE', 
                'BITPIX', 
                'NAXIS' , 
                'NAXIS1',  
                'NAXIS2', 
                'TELESCOP', 
                'INSTRUME', 
                'DETECTOR',
                'EXPTIME',
                'TEXPTIME',
                'EQUINOX', 
                'CRPIX1', 
                'CRPIX2', 
                'CRVAL1', 
                'CRVAL2', 
                'CTYPE1', 
                'CTYPE2', 
                'CD1_1', 
                'CD1_2', 
                'CD2_1', 
                'CD2_2', 
                'LTV1', 
                'LTV2', 
                'LTM1_1', 
                'LTM2_2',
                #'PA_V3',
                #'PA_FINAL',
                'AMDRA',
                'AMDDEC',
                'AMNTCH',
                'AMSGRA',
                'AMSGDEC'
                ]

        oldfits = pyfits.open(fitsfile,"update")
        oldfits_headerKeys = oldfits[0].header.ascard.keys()

        for key in oldfits_headerKeys:
                if key not in keylist:
                    del oldfits[0].header.ascard[key]

        oldfits[0].header.update('FILENAME',fitsfile)
        oldfits[0].header.update('DATASET',self.obsName)
        ostring = self.modName + ' Ver. '+__version__
        oldfits[0].header.update('ORIGIN','Properietary data by CLASH pipeline:'+ostring) #WZ
        oldfits[0].header.update('DATE',pUtil.ptime())
        try:
            oldfits[0].header.update('OBJECT',oldfits[0].header['OBJECT'])
        except:
            pass
        oldfits.close()
        del oldfits
        return


    def _fixIrafHeader(self, fitsfile=None):
        """fix up the header of the detection Image."""

        # This is the list of keywords that will be copied out of what
        # the iraf drizzle task produced. Ascards for these keywords
        # will be laid down in this order.
        # New keywords added for astrometrically corrected CRVALs.
        # See Bugzilla bug #
        if fitsfile == None:
            fitsfile = self.detImName

        keylist = [ 'SIMPLE', 
                'BITPIX', 
                'NAXIS' , 
                'NAXIS1',  
                'NAXIS2', 
                'TELESCOP', 
                'INSTRUME', 
                'DETECTOR', 
                'EQUINOX', 
                'CRPIX1', 
                'CRPIX2', 
                'CRVAL1', 
                'CRVAL2', 
                'CTYPE1', 
                'CTYPE2', 
                'CD1_1', 
                'CD1_2', 
                'CD2_1', 
                'CD2_2', 
                'LTV1', 
                'LTV2', 
                'LTM1_1', 
                'LTM2_2',
                'PA_V3',
                'PA_FINAL'
                'AMDRA'
                'AMDDEC'
                'AMNTCH'
                'AMSGRA'
                'AMSGDEC'
              ]

        oldfits = pyfits.open(fitsfile,"update")
        oldfits_headerKeys = oldfits[0].header.ascard.keys()

        for key in oldfits_headerKeys:
                if key not in keylist:
                    del oldfits[0].header.ascard[key]

        oldfits[0].header.update('FILENAME',fitsfile)
        oldfits[0].header.update('DATASET',self.obsName)
        ostring = self.modName + ' Ver. '+__version__
        oldfits[0].header.update('ORIGIN','Properietary data by CLASH pipeline:'+ostring) #WZ
        #oldfits[0].header.update('ORIGIN','ACS Science Data Pipeline:'+ostring)
        oldfits[0].header.update('DATE',pUtil.ptime())
        try:
            oldfits[0].header.update('OBJECT',oldfits[0].header['OBJECT'])
        except:
            pass

        # Now figure out the total exposure time (sum) of the images in the detection image.

        exptime = 0
        for im in self.sciImageList:
            fo = pyfits.open(im)
            exptime += fo[0].header['EXPTIME']
            fo.close()
            del fo
        oldfits[0].header.update('EXPTIME',exptime)
        oldfits[0].header.update('TEXPTIME',exptime)
        oldfits.close()
        del oldfits
    #   os.remove(self.detImName)
    #   os.rename('temp.fits',self.detImName)
        return

