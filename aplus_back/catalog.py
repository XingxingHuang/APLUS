#!/usr/bin/env python

# $Id: catalog.py,v 1.14 2006/01/30 21:22:36 anderson Exp $
# ---------------------------------------------------------------------

__version__      = '$Revision: 1.14 $ '[11:-3]
__version_date__ = '$Date: 2006/01/30 21:22:36 $ '[7:-3]
__author__       = "Ken Anderson, anderson@pha.jhu.edu"


import os
import string
import popen2
import pdb
import pUtil,xmlUtil
import extinction
import fUtil
import subprocess
from   msg   import pMessage
from   sys   import version
pyversion = version 

class Catalog: 
    """ 
    Base class for the filterCatalog and detectionCatalog classes. 
    (See those modules for more detailed descriptions of how to use them).
    Provides methods for setting and writing the required parameter
    files used by SExtractor and for executing SExtractor within
    the pipeline environment. Those classes override these methods as needed.
    The richField argument determines which parameter set will be use for UVIS
    datasets.  If true, the richfield parameter set is specified.

    NOTES:
    - use of the run() method for the filterCatalog instance is deprecated.
      in this version.  All filter catalogs should be produced by dual image
      execution, i.e. the run2() method should be called.
      Only when one filter's data is available should run() be called.
    - input parameters are those in the .inpar file
    - default input parameters are read from  the SExtractor.inpar file
      in the pipeline pars directory ($PIPELINE/pars).
    - setinpar method is available.  The pipeline will call this
      method to adjust output parameters for a particular observation's
      data and environment.
    - output parameters are read from the SExtractor.param file
      in the pipeline pars directory ($PIPELINE/pars).
    - the setoutpar method is NOT available. output parameters should be
      adjusted in the default .param this module reads.
    - cross-correlation using ASSOC* input parameters is not supported yet
    - writeXmlCatalog method now write a .xml file from the .cat file
    """ 
    def __init__(self, obs, richField=None, excludeList=None): 
        self.modName       = string.split(string.split(str(self))[0],'.')[0][1:]    # module name
        self.root          = obs.newobspath         # root path of the observation dir
        self.configdir     = obs.configdir
        self.obsName       = obs.newobs
        self.obsCats       = obs.catdir
        self.obsPars       = obs.newpardir
        self.obsFits       = obs.fitsdir            # path to fits files
        self.messagedir    = obs.messagedir         # where the module message will go.
        self.obsFitsList   = obs.fitslist           # list of fits files in newfits
        self.sciImageList  = obs.sciImageList       # list of dither-combined filter images
        self.flagImageList = obs.flagImageList
        self.rmsImageList  = obs.rmsImageList
        self.logfile       = obs.logfile            # append log entries to this file
        self.errorList     = []
        self.inputList     = []
        self.outputList    = {}
        self.checkimagesList = []                   # used for the markup
        self.inParFileList   = []                   # used for the markup
        self.excludeList     = excludeList          # used for the markup
        
        self.filterCorrections = {}    # dictionary of image names/filter extinction corrections

        self.logfile.write(self.modName+" object instantiated.")

        # full path to detection Image
        self.detectionImage = os.path.join(self.obsFits,'detection_'+obs.band+'.fits')
        # if obs.detector=='WFC':
        #    self.detectionImage = os.path.join(self.obsFits,'detection_opt.fits')
        # elif obs.detector=='IR':
        #    self.detectionImage = os.path.join(self.obsFits,'detection_nir.fits')

        # Now extract the RA and Dec positions (CRVAL1 and CRVAL2 from the header) from the
        # detection Image and determine E(B-V) at that position.
        
        ra  = float(fUtil.getKeyVal(self.detectionImage,"CRVAL1"))
        dec = float(fUtil.getKeyVal(self.detectionImage,"CRVAL2"))
        
        # getExtinction function returns a tuple of strings
        # (extinction_correction, data_quality_flags)
        
        self.logfile.write("Determining E(B-V) for detection Image coords: "+str(ra)+", "+str(dec))
        #pdb.set_trace()
        extCorr,dqFlags = extinction.getEBV(ra,dec)
        
        # convert the extinction correction, extCorr to a float for use in arithmetic
        
        self.logfile.write("E(B-V) determined: "+extCorr)
        self.eBV = float(extCorr)
        self.eBV_error = self.eBV * 0.16   # This is the expected uncertainty in E(B-V): 16%
        
        # make a list of the data quality flags
        
        self.dqList = dqFlags.strip().split()[2:]

        # Define different input parameters for different detectors
        # This requirement is described in Buzilla bug #1188
        # WZ: This has been changed, as the selection is independent
        # of detector type...

        # self.defaultInParFile = os.path.join(obs.pardir,self.modName+'.inpar')
        if (obs.detector == 'WFC'):
            self.defaultInParFile = os.path.join(obs.pardir,self.modName+'_opt.inpar')
        elif (obs.detector == 'IR'):
                self.defaultInParFile = os.path.join(obs.pardir,self.modName+'_nir.inpar')
        elif (obs.detector == 'UVIS'):
            self.defaultInParFile = os.path.join(obs.pardir,self.modName+'_uv.inpar')
        # self.defaultInParFile = os.path.join(obs.pardir,self.modName+'.inpar')
        self.logfile.write("Selected input parameter set."+self.defaultInParFile)

        #output parameters are the same for any detector.
        self.defaultParamFile = os.path.join(obs.pardir,self.modName+'.param') # default .param file
        self.inpars  = pUtil.readInParFile(self.defaultInParFile)              # readInParFile returns a dict
        self.outpars = pUtil.readParamFile(self.defaultParamFile)    # readParamFile returns a list of tuples

        # Some checking of images before proceeding.

        if not self.sciImageList:
            errtxt = "No Science Images listed in attribute sciImageList."
            self.errorList.append((self.modName,errtxt))
            raise AttributeError, errtxt
        if not self.flagImageList:
            errtxt = "No Flag Images listed in attribute flagImageList."
            self.errorList.append((self.modName,))
            raise AttributeError, errtxt
        if len(self.sciImageList) != len(self.flagImageList):
            errtxt = "Number of Flag Images does not match the number of Science Images."
            self.errorList.append((self.modName,errtxt))
            raise AttributeError, errtxt

        for im in self.sciImageList:
            self.inputList.append(os.path.basename(im)) 
        for im in self.flagImageList:
            self.inputList.append(os.path.basename(im))

    def setpars(self,fitsfile,band):
        """ generate the .inpar and .param files for a passed fitsfile."""

        self.logfile.write('Generating parameter sets for observations...')

        # Here we have to deal with the fact that the CHECKIMAGE_TYPE value could be a comma
        # seperated string of several actual image types, eg., SEGMENTATION,OBJECTS,APERTURES, etc.
        # We have to split up this field and then create CHECKIMAGE_NAMES for each one.
    
        checkImageTypes = string.split(self.getinpar('CHECKIMAGE_TYPE'),",") # a list of image types  
        # substring = string.split(self.getinpar('CHECKIMAGE_NAME'),"BACK") #WZ 
        # checkImageName = substring[0] + '_' + band + 'BACK' + substring[1]
        # Set up file specific parameters for writing the par file.

        fitsname        = string.split(fitsfile,'.')[0]
        bandname        = string.replace(fitsname,'sci',band)
        flag            = string.split(fitsfile,'_sci')[0]+'_FLAG.fits' #WZ
        #flag            = string.split(fitsfile,'_drz')[0]+'_FLAG.fits'
        catalogFileName = os.path.join(self.obsCats, bandname + '.cat')
        inParFileName   = os.path.join(self.obsPars, bandname + '.inpar')
        paramFileName   = os.path.join(self.obsPars, bandname + '.param')
        filter          = self.getinpar('FILTER_NAME')
        starnnw         = self.getinpar('STARNNW_NAME')

        #pdb.set_trace()
        if not os.path.isfile(os.path.join(self.obsFits,flag)):   
            raise RuntimeError,"Flag image missing for "+fitsfile  

        self.setinpar('CATALOG_NAME'   , catalogFileName)  
        self.setinpar('PARAMETERS_NAME', paramFileName)
        # self.setinpar('CHECKIMAGE_NAME', checkImageName)
        self.setinpar('FILTER_NAME'    , os.path.join(self.configdir, filter ))
        # self.setinpar('FLAG_IMAGE'     , os.path.join(self.obsFits,   flag   )) #WZ
        self.setinpar('STARNNW_NAME'   , os.path.join(self.configdir, starnnw))

        

        #---------------------- Begin MAG_ZEROPOINT ----------------------#
        
        # Determine the zeropoint of the fitsfile and set the MAG_ZEROPOINT parameter in the inpar file
        # We are now passing the zeropoint of the fits images to SExtractor so that is reports
        # magnitudes as actual magnitudes rather than instrumental one.
        # Bugzilla bug # 2729 now requires that the MAG_ZEROPOINT for each filter also have the filter
        # extinction correction subtracted as well.
       
        self.logfile.write("Calculating zero point for "+fitsfile)

        #pdb.set_trace()
        fitsImage = os.path.join(self.obsFits,fitsfile)
        try:
            imfilter= fUtil.filterResponseName(fitsImage)
        except fUtil.filterError:
            imfilter = fUtil.twoFilterName(fitsImage)
        try:
            filterXCorr       = extinction.filterFactor(imfilter) * self.eBV
            filterXCorrError  = filterXCorr*0.16
        except fUtil.filterError,err:
            self.errorList.append(("extintion.py","filterError: "+str(err)))
            self.logfile.write("filterError: "+str(err))
            filterXCorr       = 0.0
            filterXCorrError  = 0.0
        
        self.filterCorrections[fitsfile]=(filterXCorr,filterXCorrError)

        zpoint = fUtil.zeroPoint(fitsImage)
        self.logfile.write(fitsfile+" Header Zero point: "+str(zpoint))
        self.logfile.write("correcting zeropoint by the filter extinction factor,"+str(filterXCorr))
        # Now the zpoint is going to be corrected by the extinction factor, filterXCorr
        zpoint = zpoint - filterXCorr
        self.logfile.write(fitsfile+" Zero point is now: "+str(zpoint))
        self.setinpar('MAG_ZEROPOINT' , str(zpoint))

        #---------------------- End MAG_ZEROPOINT ----------------------#

    
        # Here we create a string of new output checkimage file names.
        # Add each one of these new name to the images list and the output list.
        checkimage_names = ''
        for im in checkImageTypes:
            imFileName = bandname+'_'+im[:4]+'.fits' 
            if os.path.join(self.obsFits,imFileName) not in self.checkimagesList: 
                checkimage_names += imFileName+','
                self.checkimagesList.append(os.path.join(self.obsFits,imFileName))
                self.outputList[os.path.basename(imFileName)] = [fitsfile]

        self.setinpar('CHECKIMAGE_NAME', checkimage_names[:-1])               # don't write a trailing comma
        inparfile  = open(inParFileName, 'w')
        outparfile = open(paramFileName, 'w')

        for param, value in self.inpars.items():
            inparfile.write("%-18s   %-20s\n" % (param, value))            
        for f in range(len(self.outpars)):
            param,value = self.outpars[f]
            if value:
                outparfile.write(param + '\n')                
        inparfile.close()
        outparfile.close()
        if inParFileName not in self.inParFileList:
            self.inParFileList.append(inParFileName)
        self.logfile.write('parameter files for ' + fitsfile + ' written.')
        return    
        
    def getinpar(self, param):
        """Returns value of input parameter given a parameter name"""
        return self.inpars[param]
    
    def getoutpar(self, param):
        """Returns value of input parameter given a parameter name"""
        return self.outpars[param]
    
    def setinpar(self, param, value):
        """Sets value of input parameter given a parameter name and value"""
        self.inpars[param] = value
        return

    def run(self):
        """Runs SExtractor on the list of fits files in sciImageList in single image
        mode.
        """
        curdir = os.getcwd()
        os.chdir(self.obsFits)
        self.catalogList = []
        for fitsfile in self.sciImageList:
            self.logfile.write('running SExtractor on ' + fitsfile)
            fitsname        = string.split(fitsfile,'.')[0]
            inParFileName   = os.path.join(self.obsPars, fitsname + '.inpar')
            catalogFileName = os.path.join(self.obsCats, fitsname + '.cat')     
            cmd = 'sex ' + os.path.join(self.obsFits, fitsfile) + ' -c ' + inParFileName
            try:
                subproc = popen2.Popen4(cmd)
                stderr_lines = subproc.fromchild.readlines()
                # call the errorSearch ufunc and trim out the junk.
                # errorSearch returns a dictionary of unique error strings and the
                # number of occurrences found in the passed list (errs).
                foundErrs = pUtil.errorSearch(stderr_lines)
                if foundErrs:
                    self.logfile.write('SExtractor produced the following message(s) on ' + fitsfile)
                    for f in foundErrs.keys():
                        self.logfile.write(f+": occurred "+str(foundErrs[f])+" times.")
                        self.errorList.append((self.modName,f+": occurred "+str(foundErrs[f])+" times."))
            except Exception, err:
                self.logfile.write('SExtractor failed:')
                self.errorList.append((self.modName,str(err)))
                raise(Exception, err)
            if os.path.isfile(catalogFileName):
                self.catalogList.append(catalogFileName)
                self.outputList[os.path.basename(catalogFileName)] = [fitsfile]
        os.chdir(curdir)
        return

    def writeXml(self):

        """
        writeXml calls the xmlUtil catalog functions on a list of source catalogs
        produced by SExtractor of the observation object and also calls the
        markupImage function (xmlUtil) on the segmentation images produced
        by SExtractor. New requirements for the markup now mean the pipeline
        has to include the filter name in the catalog tag of the xml file.
        filterName func returns the _useful_ filter name, i.e. the one that is
        not clear.  This is now part of the argument list that goes to the 
        xmlStartCat function. A new requirement has been placed upon
        this writeXml method to allow certain fields of a catalog file
        to be excluded from the xml markup.  This is done to avoid redundant
        fields for the same objects being sent to the SDA.  See Bugzilla
        bug #1436.

            """
        self.logfile.write('Generating the XML catalog(s) from SExtractor output.')

        lcount = 0
        for catfile in self.catalogList:
            base     = string.split(catfile,'.')[0]
            # the following is done to just get the filter name out of the 
            # corresponding image.
            dir,file = os.path.split(catfile)
            image    = os.path.splitext(file)[0]+'.fits'
            fullfile = os.path.join(self.obsFits,image)
            if self.excludeList:
                self.logfile.write("N.B. The following fields are being excluded from the xml output of "+file+":")
                for field in self.excludeList:
                    self.logfile.write(field)
            xmlCatalogName  =  base +'.xml'
            pardict = pUtil.readInParFile(self.inParFileList[lcount])
            xmlUtil.xmlStartCat(xmlCatalogName,self.obsName,imgfile=fullfile)
            try:
                xmlUtil.xmlStartConfig(xmlCatalogName)
            except Exception,err:
                self.errorList.append((self.modName,str(err)))
                self.logfile.write('function xmlStartConfig encountered a problem.')
                raise Exception(err)
            try:
                xmlUtil.xmlPars(pardict,xmlCatalogName)
                self.logfile.write('Input Paramaters prepended to xml catalog file.')
            except Exception,err:
                self.errorList.append((self.modName,str(err)))
                self.logfile.write('function xmlPars encountered a problem.')
                raise Exception(err)
            
            # Now we make another pardict for the extinction stuff for each of the
            # catalogs/images for the filter catalogs.  Don't do this for the detection
            # catalog, as the magnitudes are bogus anyway.
            if self.modName != "detectionCatalog":
                extParDict = {}
                extParDict["EBV"]                     = self.eBV
                extParDict["EBV_Error"]               = self.eBV_error
                extParDict["Extinction_Hcons_Flag"]   = self.dqList[0]
                extParDict["Extinction_Asteroid_Flag"]= self.dqList[1]
                extParDict["Extinction_Glitch_Flag"]  = self.dqList[2]
                extParDict["Extinction_Source_Flag"]  = self.dqList[3]
                extParDict["EXtinction_NoList_Flag"]  = self.dqList[4]
                extParDict["Extinction_BigObj_Flag"]  = self.dqList[5]
                extParDict["Extinction_NoIRAS_Flag"]  = self.dqList[6]
                extParDict["Extinction"]              = self.filterCorrections[image][0]
                extParDict["Extinction_Error"]        = self.filterCorrections[image][1]
                # and send this to the xmlPars function
                try:
                    xmlUtil.xmlPars(extParDict,xmlCatalogName)
                    self.logfile.write('Extinction Correction Paramaters prepended to xml catalog file.')
                except Exception,err:
                    self.errorList.append((self.modName,str(err)))
                    self.logfile.write('function xmlPars encountered a problem with extinction parameters')
                    raise Exception(err)
            
            try:
                xmlUtil.xmlEndConfig(xmlCatalogName)
            except Exception,err:
                self.errorList.append((self.modName,str(err)))
                self.logfile.write('function xmlEndConfig encountered a problem')
                raise Exception(err)
            
            # mark up the data
            try:
                xmlUtil.xmlData(catfile,xmlCatalogName,excludeList=self.excludeList)
                self.logfile.write('XML catalog written for ' + os.path.basename(catfile))
            except Exception,err:
                self.errorList.append((self.modName, str(err)))     
                self.logfile.write('function xmlify encountered a problem.')
                raise Exception(err)
            xmlUtil.xmlEndCat(xmlCatalogName)
            self.outputList[os.path.basename(xmlCatalogName)] = [os.path.basename(catfile)]
            lcount = lcount +1
        # Now markup the checkimages generated by SExtractor.
        for im in self.checkimagesList:
            try:
                opath = xmlUtil.markupImage(im,dataset=self.obsName)
                self.outputList[os.path.basename(opath)] = [os.path.basename(im)]
            except Exception,err:
                print "xmlUtil.markupImage threw an exception: "+str(err)
                print "New pyfits?"
                print "No markup on checkimages done."
                self.logfile.write("Warning: xmlUtil.markupImage threw an exception.")
                self.logfile.write(str(err))
                self.logfile.write("Suspected problem in pyfits.")
                self.errorList.append((self.modName,str(err)))
        return


    def mkMsg(self):
        """create and write module level message for this class.
        Most of this is just compiling the info. meta in a dictionary
        of lists where each list is a list of tuples describing the
        tag lines for the particular section of the message.  This tuple 
        format conforms to that used by the xmlMessage class which is
        modeled on basic python argument passing, i.e. (key,*value,**attr).
        """
        self.meta = {}
        self.meta['module']    = []
        self.meta['meta']      = []
        self.meta['input']     = []
        self.meta['output']    = []
        self.meta['errorlist'] = []

        self.meta['module'].append(('module','name='+self.modName,'version='+__version__,'dataset='+self.obsName))
        #instname = string.split(string.split(str(self))[0],'.')[1]
        self.meta['module'].append(('root',self.root))
        #self.meta['module'].append(('instance', instname))
        self.meta['meta'].append(('meta',))
        self.meta['meta'].append(('depend',))
        self.meta['meta'].append(('pkg',))
        self.meta['meta'].append(('name','python'))
        self.meta['meta'].append(('version',pyversion.split()[0]))

        # SExtractor info
        sub  = subprocess.Popen(['sex', '--version'], stdin=subprocess.PIPE, stdout=subprocess.PIPE, close_fds=True)
        outp = sub.stdout.readlines()
        name = outp[0].split()[0]
        ver  = outp[0].split()[2]
        self.meta['meta'].append(('pkg',))
        self.meta['meta'].append(('name',name))
        self.meta['meta'].append(('version',ver))
        cmdline1 = 'sex fitsfile -c inParFileName '
        self.meta['meta'].append(('commandline',cmdline1))
        del outp,sub,name,ver


        if self.errorList:
            self.meta['errorlist'].append(('errorlist',))
            for pkg,err in self.errorList:
                self.meta['errorlist'].append(('erroritem',err,'frompkg='+pkg))

        # input section
        self.meta['input'].append(('input',))
        for f in self.inputList:
            self.meta['input'].append(('file','type=image/x-fits'))
            self.meta['input'].append(('name',os.path.join("Images",f)))

        # output section
        if self.outputList:
            self.meta['output'].append(('output',))
        for f in self.outputList.keys():
            if string.find(f,"fits.xml") != -1:
                self.meta['output'].append(('file','type=text/xml'))
                self.meta['output'].append(('name',os.path.join("Images",f)))
                for pred in  self.outputList[f]:
                    self.meta['output'].append(('predecessor',os.path.join("Images",pred)))
            elif string.find(f,".fits") != -1:
                self.meta['output'].append(('file','type=image/x-fits'))
                self.meta['output'].append(('name',os.path.join("Images",f)))
                for pred in  self.outputList[f]:
                    self.meta['output'].append(('predecessor',os.path.join("Images",pred)))
            elif string.find(f,".xml") != -1:
                self.meta['output'].append(('file','type=text/xml'))
                self.meta['output'].append(('name',os.path.join("Catalogs",f)))
                for pred in  self.outputList[f]:
                    self.meta['output'].append(('predecessor',os.path.join("Catalogs",pred)))
            else:
                self.meta['output'].append(('file','type=text/ascii'))
                self.meta['output'].append(('name',os.path.join("Catalogs",f)))
                for pred in  self.outputList[f]:
                    self.meta['output'].append(('predecessor',os.path.join("Images",pred)))
        

        # pass this dictionary to the class pMessage...
        msgFile = os.path.join(self.messagedir,self.modName+"_module.xml")
        mmsg = pMessage(self.meta)
        mmsg.writeMsg(msgFile)
        return 
