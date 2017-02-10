#!/usr/bin/env python

# $Id: colorCatalog.py,v 2.00 2012/06/18 23:14:13 Zheng Exp $
# ---------------------------------------------------------------------

__version__      = '$Revision: 2.00 $ '[11:-3]
__version_date__ = '$Date: 2012/06/18 23:14:13 $ '[7:-3]
__author__       = "Wei Zheng, zheng@pha.jhu.edu"

0
## This is the systematic error in the zero point calibration
## From Benitez email, Fri, 24 May 2002 19:04:55
## Update when the zeropoints are improved
## Add band option. 
zp_error=0.03 # was 0.01 WZ

import os,string,glob
import pdb
import pUtil, xmlUtil, fUtil, extinction
import numpy.oldnumeric as Numeric
import numpy.oldnumeric.mlab as MLab
import math
import pyfits
numver = Numeric.__version__
import tableio
from   msg import pMessage
from   sys import version
pyversion = version

class colorCatalog:
    """
    This module provides the ColorCatalog class which has only
    the run method to create a multicolor catalog.  Constructor is passed
    an WFC3 pipeline DataSet object and looks for a detectionImage catalog and 
    throws an exception if this catalog does not exist.  This detectionImage
    catalog must be created by a call to the detectionCatalog class of the
    catalogs module. See that module for those details.

    usage:

    mcat = colorCatalog.colorCatalog(obs)
    mcat.run()

    Look for the mulitcolor.cat catalog and .columns file in the particular 
    DataSet's catalog directory.

    """

    def __init__(self,obs,catalogList,band,excludeList=["X_IMAGE", "Y_IMAGE", "FWHM_IMAGE",'CLASS_STAR','ALPHA_J2000','DELTA_J2000']): # WZ
        self.modName    = string.split(string.split(str(self))[0],'.')[0][1:] # module name
        self.aperData   = obs.base.pipeline+"/aperData"                       # aperture data files
        self.root       = obs.newobspath                                      # root path observation dir
        self.obsName    = obs.newobs
        self.obsCats    = obs.catdir
        self.obsPars    = obs.newpardir
        self.obsFits    = obs.fitsdir            # path to fits files
        self.messagedir = obs.messagedir         # where the module message will go.
        self.logfile    = obs.logfile            # append log entries to this file
        self.errorList  = []
        self.inputList  = []
        self.catalogList= []
        self.outputList = {}
        self.filterPars = {}                     # this is a dictionary of filter parameter dictionaries
                                                 # like {"filter1":{"par1":val1,par2:val2,...},"filter2":...}
        self.excludeList= excludeList
        self.filterNames    = {}                 # dict of { filter1:fitsimage1, filter2:fitsimage2,...}
        self.columnsPopList = []
        dname='detection_'+band+'.fits' #WZ
        self.detectionImage = os.path.join(obs.fitsdir,dname)
        dname='detection_'+band+'.cat'

        for cat in catalogList:                 # explicit list of filter catalogs passed
            self.catalogList.append(os.path.basename(cat))

        self.logfile.write('Instantiating colorCatalog object...')

        if not os.path.exists(self.detectionImage):
            errtxt = "Instance of colorCatalog cannot find a detection image in this dataset."
            self.errorList.append((self.modName,errtxt))
            self.logfile.write(errtxt)
            raise IOError,errtxt

        # Now extract the RA and Dec positions (CRVAL1 and CRVAL2 from the header) from the
        # detection Image for use later with the extinction correction.
        ra  = float(fUtil.getKeyVal(self.detectionImage,"CRVAL1"))
        dec = float(fUtil.getKeyVal(self.detectionImage,"CRVAL2"))

        # getExtinction function returns a tuple of strings (extinction_correction, data_quality_flags)
        self.logfile.write("Determining E(B-V) for detection Image coords: "+str(ra)+", "+str(dec))
        extCorr,self.DQflags = extinction.getEBV(ra,dec)
        # convert the extinction correction, extCorr to a float for use in arithmetic
        self.logfile.write("E(B-V) determined: "+extCorr)
        self.eBV = float(extCorr)
        self.eBV_error = self.eBV * 0.16   # This is the expected uncertainty in E(B-V): 16%

        if not os.path.exists(os.path.join(self.obsCats,dname)):
            errtxt = "Instance of colorCatalog cannot find a detection image catalog in this dataset."
            self.errorList.append((self.modName,errtxt))
            self.logfile.write(errtxt)
            raise IOError,errtxt

        for f in self.catalogList:
            self.inputList.append(f)             # record all catalogs for this observation as input.

        # This is a list of initial output parameters for the final multicolor catalog
        # final magnitudes are added to the catalog after these appear.

        # 0self.outColumns = ['NUMBER','X_IMAGE','Y_IMAGE','FWHM_IMAGE']
        self.outColumns = ['NUMBER','X_IMAGE','Y_IMAGE','FWHM_IMAGE','CLASS_STAR','ALPHA_J2000','DELTA_J2000'] # WZ

        # Now read the detection catalog and create the dictionary of headers.
        # selectSet is a dictionary where keys are the column names and values are
        # the column numbers, starting at 1.  We get the columns above from the
        # detection catalog.  All other columns will come from the individual
        # filter catalogs.

        self.detCatalog = os.path.join(self.obsCats,dname)
        catalog         = open(self.detCatalog).readlines()
        selectSet       = pUtil.makeHeaderDict(catalog)
        detectionList   = []
        for key in self.outColumns:
            detectionList.append(selectSet[key] - 1)
        self.detectionColumns = tuple(detectionList)     # the get_data function requires a tuple

        # call the private _apcorrSetup to generate the object radii list.
        # This will be used to then create the aperture corrections for the objects
        # in the various filters of the dataset.
        
        self.iradList = self._apcorrSetup(catalog,selectSet)

    def run(self,band): 
        """Make the multicolor catalog."""
        # define the zero point error
        zpoint_err = 0.01

        self.logfile.write('Processing catalogs...')

        detection_variables = tableio.get_data(self.detCatalog,self.detectionColumns)
        nsources = len(detection_variables[0])
        # WZ Should be a dimension check here with the filter catalogs before crashing
        ncats    = len(self.catalogList)

        #Create a series of arrays to keep the relevant data
        flux    = Numeric.zeros((ncats,nsources))*1.0          # ISO flux ... for now.
        fluxerr = Numeric.zeros((ncats,nsources))*1.0          # ISO flux error... for now.
        m       = Numeric.zeros((ncats,nsources))*1.0          # MAG_ISO
        em      = Numeric.zeros((ncats,nsources))*1.0          # MAGERR_ISO
        m_corr  = Numeric.zeros((ncats,nsources))*1.0          # extinction corrected magnitude
        em_corr = Numeric.zeros((ncats,nsources))*1.0          # extinction corrected magnitude error
        ap_corr = Numeric.zeros((ncats,nsources))*1.0          # aperture   correction
        m_bpz   = Numeric.zeros((ncats,nsources))*1.0          # MAG_BPZ    extinction/aperture corrected
        em_bpz  = Numeric.zeros((ncats,nsources))*1.0          # MAGERR_BPZ extinction/aperture corrected


        preds = [os.path.basename(self.detCatalog)]

        # this is a list of those indices of the catalogList which raised the filterError exception
        # i think i need to know this, but i'm not sure why yet....
        self.badcats = []
        for i in range(ncats):
            selectSet = {}
            pardict = {}
            imname  =  os.path.splitext(self.catalogList[i])[0]
            j=string.rfind(imname,'_')
            basefits= imname[0:j] + '_sci.fits'    
            fitsfile= os.path.join(self.obsFits,basefits)
            detector= fUtil.getDetector(fitsfile)
            try:
                imfilter= fUtil.filterResponseName(fitsfile)

            except fUtil.filterError:
                self.badcats.append(i)
                imfilter = fUtil.twoFilterName(fitsfile)

            try:
                filterXCorr       = extinction.filterFactor(imfilter) * self.eBV
                filterXCorrError  = filterXCorr*0.16
            except fUtil.filterError,err:
                self.errorList.append(("extintion.py","filterError: "+str(err)))
                self.logfile.write("filterError: "+str(err))
                filterXCorr       = 0.0
                filterXCorrError  = 0.0

            try:
                f = pyfits.open(fitsfile) # WZ 
                instr = f[0].header.get('INSTRUME')
                f.close()
                del f
                if instr == 'ACS':
                    #filter = 'F814W'   #WZ
                    filter = fUtil.filterNameACS(fitsfile)
                elif instr == 'WFC3':
                    #filter = fUtil.filterName(fitsfile)
                    filter = fUtil.filterNameWFC3(fitsfile)
                if filter:
                    print "Find filter ", filter
                else:
                    raise filterError,"No filter found for fits image: "+fitsfile
                #filter = fUtil.filterName(fitsfile)
                ap_corr[i,:] = self._apcorr(detector,filter)
            except KeyError,err:
                print str(err)
                self.logfile.write(str(err))
                # if we can't do it, well, we can't do it.
                pass
            except fUtil.filterError,err:
                print str(err)
                self.logfile.write(str(err))
                # No filter aperture data, no aperture correction.  move on.
                pass
            
            
            dqlist = self.DQflags.strip().split()[2:]

            pardict["Extinction_Hcons_Flag"]   = dqlist[0]
            pardict["Extinction_Asteroid_Flag"]= dqlist[1]
            pardict["Extinction_Glitch_Flag"]  = dqlist[2]
            pardict["Extinction_Source_Flag"]  = dqlist[3]
            pardict["EXtinction_NoList_Flag"]  = dqlist[4]
            pardict["Extinction_BigObj_Flag"]  = dqlist[5]
            pardict["Extinction_NoIRAS_Flag"]  = dqlist[6]
            pardict["Extinction"]              = filterXCorr
            pardict["Extinction_Error"]        = filterXCorrError

            pardict["EBV"]        = self.eBV
            pardict["EBV_Error"]  = self.eBV_error

            self.filterNames[imfilter]= basefits
            catalog = ''
            catalog = os.path.join(self.obsCats,self.catalogList[i])
            preds.append(self.catalogList[i])                # predecessors list

            # Now we search through each filter's catalog for the keywords below
            # and determine which columns those values are in.
            # Subtract 1 from the column numbers to use as array indices
            # i.e the catalog column numbers start at 1 not zero.

            catalogObj = open(catalog).readlines()
            selectSet  = pUtil.makeHeaderDict(catalogObj)

            #Info for flux columns
            fluxList = []
            fluxList.append(selectSet['FLUX_ISO']   - 1)
            fluxList.append(selectSet['FLUXERR_ISO']- 1)
            self.fluxColumns = tuple(fluxList)    # the get_data function interface requires a tuple

            # Build the various columns arrays with the get_data function.
            # We read raw fluxes and errors into the flux,fluxerr arrays.
            # They are afterwards transformed to magnitudes
            t1,t2 = tableio.get_data(catalog,self.fluxColumns)
            if (len(t1) != nsources):
                self.logfile.write("Catalog dimension mismatch: ", str(len(t1)),' ', nsources)
                self.logfile.write("Check patrameters in detectionCatalog.inpar and filterCatalog.inpar")
            #flux[i,:],fluxerr[i,:] = tableio.get_data(catalog,self.fluxColumns)
            flux[i,:],fluxerr[i,:] = t1,t2
            flux[i,:] = pUtil.deNAN(flux[i,:])

            # Those objects with flux equal or less than 0 are assigned a magnitude of 99
            # and a limiting magnitude equal to their SExtractor photometric error. This
            # is interpreted by BPZ as a nondetection with zero flux and 1-sigma error
            # equal to the limiting magnitude

            nondetected = Numeric.less_equal(flux[i,:],0.0)*Numeric.greater(fluxerr[i,:],0.0)

            # Those objects with error flux and flux equal to 0 are assigned a magnitude of -99
            # and a flux of 0, which is interpreted by SExtractor as a non-observed object

            nonobserved = Numeric.less_equal(fluxerr[i,:],0.0)

            # When flux error > 100*(flux), mark as nonobserved (Benitez, 24-Oct-03).
            
            nonobserved = Numeric.where(fluxerr[i,:] > 100*(abs(flux[i,:])),1.0,nonobserved[:])
            
            detected    = Numeric.logical_not(nonobserved+nondetected)

            # Get the zero point for the final magnitudes

            zpoint  = fUtil.zeroPoint(fitsfile)       # pass the fits file to zeroPoint func
            pardict["Zeropoint"] = str(zpoint)
            pardict["Zeropoint_Error"] = zpoint_err
            self.logfile.write("Photometric ZeroPoint of "+
                                   os.path.basename(fitsfile)+": "+str(zpoint))

            #---------------- Temporary zero point correction --------------------#
            ## Commented 24-Sep-2002 as per Bugzilla bug #1800
            ##
            ## zpointCor=fUtil.zeroPointCorrection(imfilter)
            ## zpoint+= zpointCor
            ## self.logfile.write("ZeroPoint %s corrected by %.2f mag" % (imfilter,zpointCor))
            ## print "ZeroPoint %s corrected by %.2f mag" % (imfilter,zpointCor)
            ##
            #---------------- End temporary zero point correction ----------------#                   
            
            flux[i,:]=Numeric.clip(flux[i,:],1e-100,1e100)
            m[i,:]  = Numeric.where(detected,-2.5*Numeric.log10(abs(flux[i,:]))+zpoint,m[i,:])
            m[i,:]  = Numeric.where(nondetected,99.0,m[i,:])
            m[i,:]  = Numeric.where(nonobserved,-99.0,m[i,:])
            # the filter specific extinction correction is applied.
            m_corr[i,:] = Numeric.where(nondetected,99.0,m[i,:]-filterXCorr)
            m_corr[i,:] = Numeric.where(nonobserved,-99,m_corr[i,:])

            m_bpz[i,:] = Numeric.where(nondetected,99.0,m_corr[i,:] + ap_corr[i,:])
            m_bpz[i,:] = Numeric.where(nonobserved,-99,m_bpz[i,:])

            # clip values from being too small or large, i.e. 0 or inf.
            fluxerr[i,:]=Numeric.clip(fluxerr[i,:],1e-100,1e100)
            em[i,:]  = Numeric.where(detected,2.5*Numeric.log10(1.0+abs(fluxerr[i,:]/flux[i,:])),em[i,:])
            em[i,:]  = Numeric.where(nondetected,2.5*Numeric.log10(abs(fluxerr[i,:]))-zpoint,em[i,:])
            em[i,:]  = Numeric.where(nonobserved,0.0,em[i,:])
            # will need to determine how the mag error changes with the app of the extinction
            # correction. For now, the error in the corrected mags is the same.
            em_corr[i,:] = Numeric.where(nondetected,em[i,:]+filterXCorr,em[i,:])

            # no spec to handle errors for the bpz magnitudes.
            # make them the same as the corrected errors.
            em_bpz[i,:]  = em_corr[i,:]
            
            #-------------------- Limiting Magnitude Section --------------------#
            # N.B. This section merely determines the limiting magnitude, if it can.
            # The limMag value is currently not used.  Hence the try/except clause.
            #To define the typical 1-sigma limiting magnitude (which now changes from object to object),
            #we take the median of the 5% faintest 1-sigma fluxerror_auto fluxes. This is obviously 
            #biased, but shows which is the depth in the deepest part of the images with the current 
            #SExtractor parameters, etc.

            # one_sigma_mags = Numeric.compress(Numeric.less_equal(abs(em[i,:]-0.7526),0.02),m[i,:])
            # n_one  = len(one_sigma_mags)
            # test snippet:

            dm = 0.02
            n_one = 0
            while n_one < 21:
                one_sigma_mags = Numeric.compress(Numeric.less_equal(abs(em[i,:]-0.7526),dm),m[i,:])
                n_one  = len(one_sigma_mags)
                dm += 0.01
                if dm > 0.03 and dm < 1.1:
                    message="Warning: not enough 1-sigma objects in the catalog. Using dm=+-%.2f" % dm
                    self.logfile.write(message)
                if dm >= 1.1:
                    message="Warning: Stopped searching for 1-sigma objects at dm=+-%.2f" % dm
                    self.logfile.write(message)
                    break

            try:
                limMag = MLab.median(Numeric.sort(one_sigma_mags)[-20:])
                self.logfile.write("Limiting Mag  "+ basefits + ":"+str(limMag))
                print "Limiting Mag  "+ basefits + ":"+str(limMag)     
            except Exception,err:
                print "Exception raised: could not determine limiting magnitude."
                print err

            #---------------------- End  Limiting Magnitude Section ----------------------#
            
            #self.outColumns.append(imfilter +'_MAG_ISO')
            #self.outColumns.append(imfilter +'_MAGERR_ISO')
            self.outColumns.append(imfilter +'_MAGCORR_ISO')
            self.outColumns.append(imfilter +'_MAGERRCORR_ISO')
            self.outColumns.append(imfilter +'_APER_CORR')
            self.outColumns.append(imfilter +'_MAG_BPZ')
            self.outColumns.append(imfilter +'_MAGERR_BPZ')
        
            # save the filter specific pardict into the filterPars dictionary
            self.filterPars[imfilter] = pardict

        #Save the results to file
        #self.colorcat = os.path.join(self.obsCats,'multicolor_',band,'.cat')
        self.colorcat = self.obsCats+'/multicolor_'+band+'.cat' #WZ
        # Build the header as a list of tuples then write the header into the file.
        self.header = []
        for key in self.outColumns:
            self.header.append((key,len(self.header)+1))
            
        self._writeHeader(self.colorcat,self.header)
        self._writeColumns(self.header)     # write the .columns file for BPZ
        
        # OK. header is now written
        vars = list(detection_variables)
        for i in range(ncats): 
            #vars.append((m[i,:]))
            #vars.append((em[i,:]))
            vars.append((m_corr[i,:]))
            vars.append((em_corr[i,:]))
            vars.append((ap_corr[i,:]))
            vars.append((m_bpz[i,:]))
            vars.append((em_bpz[i,:]))
        

        variables = tuple(vars)
        format = '%i\t'+'%4f  '*(len(variables)-1)
        self.logfile.write('Writing data to multicolor catalog...')
        tableio.put_data(self.colorcat,variables,format=format,append='yes')
        self.outputList[os.path.basename(self.colorcat)] = preds
        self.logfile.write('Multicolor catalog complete.')
        return

    def writeXml(self):

        """
        writeXml marks the multicolor catalog with the pipeline protocol markup. 
        A new requirement has been placed upon this writeXml method to allow certain 
        fields of the catalog file to be excluded from the xml markup.  This is done 
        to avoid redundant fields for the same objects being sent to the SDA. 
        See Bugzilla bug #1436.
    
        New functionality in this method now has each filter producing it's own photometry
        catalog.  This functionality is implement via the excludeList which allows the
        method to select different fields to be excluded from the catalog markup.  This method
        now uses the columnPopList to remove a filter from the list of filters, i.e. the APER_CORR,
            MAG_BPZ and MAGERR_BPZ of that filter will not be excluded.  The filter is then re-inserted
        into the list and iterates to the next filter.
        See Bugzilla bug #1599.

        The aperture correction specified in bug # 2708, is now part of both the multicolor
        and xml catalogs.  The field is included in the markup as the new final magnitudes being
        sent to bpz,
        _APER_CORR
        _MAG_BPZ
        _MAGERR_BPZ

        """
        self.logfile.write("Generating the XML catalogs from the multicolor catalog.")

        self.logfile.write("N.B. The following fields are being excluded from the xml output: ")
        if self.excludeList:
            for field in self.excludeList:
                self.logfile.write(field)

        self.columnsPopList = self.filterNames.keys()
        for current_filter in self.filterNames.keys():
            i = self.columnsPopList.index(current_filter)
            self.columnsPopList.remove(current_filter)
            xmlCatalogName =  os.path.join(self.obsCats,current_filter+'_photometry.xml')
            imgfile = os.path.join(self.obsFits,self.filterNames[current_filter])
            xmlUtil.xmlStartCat(xmlCatalogName,self.obsName,imgfile=imgfile)
            xmlUtil.xmlStartConfig(xmlCatalogName)
            xmlUtil.xmlPars(self.filterPars[current_filter],xmlCatalogName)
            xmlUtil.xmlEndConfig(xmlCatalogName)
            self.logfile.write('Input Paramaters prepended to xml photometry catalog file.')

            newXList = []
            newXList.extend(self.excludeList)
            # Append the current filter fields below to the exclude list for the markup.
            # Only the aperture correction and the bpz mag and err go into the xml catalog.
            # See Bugzilla bug # 2708.
            newXList.append(current_filter+"_MAG_ISO")
            newXList.append(current_filter+"_MAGERR_ISO")
            newXList.append(current_filter+"_MAGCORR_ISO")
            newXList.append(current_filter+"_MAGERRCORR_ISO")
            for item in self.columnsPopList:
                newXList.append(item+"_MAG_ISO")
                newXList.append(item+"_MAGERR_ISO")
                newXList.append(item+"_MAGCORR_ISO")
                newXList.append(item+"_MAGERRCORR_ISO")
                newXList.append(item+"_APER_CORR")
                newXList.append(item+"_MAG_BPZ")
                newXList.append(item+"_MAGERR_BPZ")

            self.columnsPopList.insert(i,current_filter)
            xmlUtil.xmlColorData(self.colorcat,xmlCatalogName,newXList)
            self.logfile.write("XML catalog written: "+current_filter+'_photometry.xml')

            xmlUtil.xmlEndCat(xmlCatalogName)
            self.outputList[os.path.basename(xmlCatalogName)] = [os.path.basename(self.colorcat)]
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
        self.meta['meta'].append(('name','numpy.oldnumeric'))
        self.meta['meta'].append(('version',numver))
        self.meta['meta'].append(('pkg',))
        self.meta['meta'].append(('name','tableio.py'))
        self.meta['meta'].append(('version',tableio.__version__))


        if self.errorList:
            self.meta['errorlist'].append(('errorlist',))
            for pkg,err in self.errorList:
                self.meta['errorlist'].append(('erroritem',err,'frompkg='+pkg))

        # input section
        self.meta['input'].append(('input',))
        for f in self.inputList:
            self.meta['input'].append(('file','type=text/ascii'))
            self.meta['input'].append(('name',os.path.join("Catalogs",f)))

        # output section
        if self.outputList:
            self.meta['output'].append(('output',))
        for f in self.outputList.keys():
            if string.find(f,".xml") != -1:
                self.meta['output'].append(('file','type=text/xml'))
                self.meta['output'].append(('name',os.path.join("Catalogs",f)))
                for pred in  self.outputList[f]:
                    self.meta['output'].append(('predecessor',os.path.join("Catalogs",pred)))
            else:
                self.meta['output'].append(('file','type=text/ascii'))
                self.meta['output'].append(('name',os.path.join("Catalogs",f)))
                for pred in  self.outputList[f]:
                    self.meta['output'].append(('predecessor',os.path.join("Catalogs",pred)))

        # pass this dictionary to the class pMessage...
        msgFile = os.path.join(self.messagedir,self.modName+"_module.xml")
        mmsg = pMessage(self.meta)
        mmsg.writeMsg(msgFile)
        return 


#----------------------------------------------------------------------------------------------------#
#---------------------------------- "private" helper methods ----------------------------------------#
#----------------------------------------------------------------------------------------------------#

    def _apcorrSetup(self,dcat,detectionSelectSet):
        """private method called by the constructor to help initialise aperture correction stuff.
        This method reads the detection catalog and greps the ISOAREA_IMAGE data out of that catalog.
        This method then converts the areas for the objects in the catalog into radii and returns
        a Numeric array to the constructor.  These ISO areas are the same for all filters so this
        only has to happen once, hence the constructor calls this.  The catalog arg  passed is the
        detection catalog which is just a readlines() list, the second arg is the selectSet for the
        detection catalog header.
        """
        iradList = []
        isoAreaIndex = detectionSelectSet["ISOAREA_IMAGE"] - 1 # index is one less than the col no.
        isoAreaData = tableio.get_data(self.detCatalog,isoAreaIndex)
        for area in isoAreaData:
            iradList.append(max(2.0 ,(math.sqrt(area/math.pi))))
        return iradList
    
    def _apcorr(self,detector,filter):
        """Calculate the aperture corrections for this dataset for the given detector and filter.
        This is a factor to be added to the magnitudes calculated above array [m].
        This "apcor" is to be in magnitudes and apcor must be <= 0.  This new magnitude or
        array of magnitudes is then written into the multicolor.cat file and presented to
        BPZ as the magnitude to use, as indicated in the multicolor.columns file.
        Function returns a numeric array of aperture corrections.  The fiducial radius
        specified herein is specified in Bugzilla bug #2708, as are the other specifications
        for this method.
        """
        apcorr = []
        apcorrArray = None
        irad_fiducial = 14
        aperFiles = {
            "UVIS" : os.path.join(self.aperData,"newUVIS_2011.ee_new_csky.dat"),
            "IR"   : os.path.join(self.aperData,"newIR_1210.ee_new_csky.dat"  ),
            "WFC" : os.path.join(self.aperData,"newWFC_0803.ee_new_csky.dat"),
            #{}"HRC" : os.path.join(self.aperData,"newHRC_0803.ee_new_csky.dat")
            }
            
        try:
            aperDataFile = aperFiles[detector]
        except KeyError:
            raise KeyError,"No known aperture data for detector, "+detector

        # now get the filter aperture data for the indicated detector.
        # raise a filterError if the lookup fails.
        
        aperCatalog = open(aperDataFile).readlines()
        selectSet   = pUtil.makeHeaderDict(aperCatalog)
        colName     = "EE_"+filter
        try:
            aperIndex   = selectSet[colName]-1
        except KeyError:
            raise fUtil.filterError,"No aperture data for filter "+filter
        
        aperData    = tableio.get_data(aperDataFile,aperIndex)

        # the aperData contains the encircled energy as a function of radius in
        # pixels.  Conveniently, the pixel radius = index + 1 of the value we
        # want from the array for that radius.
        
        eeFiducial  = aperData[irad_fiducial -1]

        for irad in self.iradList:
            if irad > irad_fiducial:
                apcorr.append(0)
                continue
            
            # irad1, irad2 are the adjacent integers of the radius of the object.
            # used to do a linear interpolation on the encircled energy curve.

            irad1    = int(irad)
            irad2    = irad1 + 1
            iradFrac = irad - irad1
            ee_irad1 = aperData[irad1-1]
            ee_irad2 = aperData[irad2-1]
            k = (ee_irad2 - ee_irad1)
            ee_irad = k*iradFrac + ee_irad1
            apcorr.append(min(0,+2.5*math.log10(ee_irad/eeFiducial)))
        apcorrArray = Numeric.array(apcorr)
        return apcorrArray

        
    def _writeHeader(self,cat,header):
        """Write the header to the color catalog file.
        """
        self.logfile.write('Writing multicolor catalog header...')
        hd = open(cat,'w')
        hd.write('## ' + pUtil.ptime() + '\n')
        hd.write('## ' + self.modName + ' Catalog file for Observation: ' + self.obsName + '\n')
        hd.write('## This proprietary file was by the CLASH Pipeline.\n##\n')
        for pair in header:
            name,col = pair
            hd.write('# '+ str(col) + '\t' + name + '\n')
        hd.close()
        return

    def _writeColumns(self,header):
        """Write the .columns file for BPZ. will match the data in the header of the
        multicolor catalog.
        """
        xFilts = []
        badFilts = []
        ignoreFilters = ['FR388N','FR423N','FR462N' ,'FR505N','FR551N','FR601N' ,
                          'FR656N','FR716N','FR782N' ,'FR853N','FR931N', 'FR1016N',
                          'FR459M','FR647' ,'FR914M']
        
        # columnsFile = os.path.join(self.obsCats,'multicolor.columns')
        columnsFile = string.replace(self.colorcat,'.cat','.columns') #WZ
        cfile = open(columnsFile,'w')
        cfile.write('## ' + pUtil.ptime() + '\n')
        cfile.write('## ' + self.modName + ' .columns file for Observation: ' + self.obsName + '\n')
        cfile.write('## (This file was generated automatically by the APLUS Pipeline.  Do not edit.)\n##\n')
        cfile.write('## N.B. The columns listed are the extinction corrected/aperture corrected\n')
        cfile.write('## magnitude columns which appear as MAG_BPZ and MAGERR_BPZ of the \n')
        cfile.write('## multicolor catalog file. These columns are now sent to BPZ.\n')
        cfile.write('## (See Bugzilla bug #2708.)\n##\n')
        i=0
        for i in range(len(header)):
            name,col=header[i]
            line=name.split('_')
            if len(line) < 3:
                continue
            elif line[-1] == "CORR":
                continue
            elif len(line) > 6:
                if (line[3]+"_"+line[4]) not in xFilts:
                    xFilts.append(line[3]+"_"+line[4])
                continue
            elif line[3] in ignoreFilters:
                if line[3] not in badFilts:
                    badFilts.append(line[3])
                continue
            elif line[-2] != "MAG":
                continue
            elif line[-1] != "BPZ":
                continue
            else:
                # colname=line[3]+'_'+line[2] # WZ  Only runs if F814W or F850Lp images are present
                colname=line[0]+'_'+line[1]+'_'+line[2]+'_'+line[3]
                if colname=='HST_ACS_WFC_F814W' or colname=='HST_ACS_WFC_F850LP':
                    n_mo=str(col)

                cfile.write(colname+'\t'+str(col)+','+str(col+1)+'\tAB\t%.2f\t0.0\n' % zp_error)

        for filt in xFilts:
            self.logfile.write("Excluding complex filter, "+filt+" from bpz processing")
            self.errorList.append((self.modName,"Complex Filter "+filt+" excluded from bpz processing"))
        for badfilt in badFilts:
            self.logfile.write("Excluding catalog from bpz interface for filter, "+badfilt)
            self.errorList.append((self.modName,"Filter "+badfilt+" excluded from bpz processing"))

        ## From Benitez email, Fri, 24 May 2002 19:04:55
        ## Hack to add the F814W column
        try: 
            cfile.write('M_0\t%s\nID\t1\n'%n_mo) #WZ
        except: pass
        cfile.close()
        self.outputList[os.path.basename(columnsFile)] = [os.path.basename(self.colorcat)]
        return

##################################################################
#WZ: planted from utils/fUtil.py

def temfilterNameACS(fitsfile):
    """returns a the filter name of the real (i.e. not clear) filter used
    in the image.  Function determines this like the above filterResponseName func
    but only returns the filter name.
    """
    f     = pyfits.open(fitsfile)
    filt1 = f[0].header.get('FILTER1')
    filt2 = f[0].header.get('FILTER2')
    f.close()
    del f

    if not filt1 and not filt2:
        return None

    acs = filters.ACSFilterInfo()
    if filt1 in acs.clear and filt2 in acs.clear:
        raise filterError,"WARNING: Both filter keywords in image header are clear."

    if filt1 not in acs.clear and filt2 not in acs.clear:
        raise filterError,"WARNING: Neither filter keyword in image header is clear."

    if filt1 not in acs.clear:
        filter = filt1
    else:
        filter = filt2
    del acs
    return filter

def filterNameWFC3(fitsfile):
    """returns a the filter name of the real (i.e. not clear) filter used
    in the image.  Function determines this like the above filterResponseName func
    but only returns the filter name.
    """
    f     = pyfits.open(fitsfile)
    filt = f[0].header.get('FILTER')
    f.close()
    del f

    if not filt: # and not filt2:
        return None

    wfc3 = filters.WFC3FilterInfo()
    #pdb.set_trace()
    #if filt in wfc3.clear: # and filt2 in wfc3.clear: # WZ
    #    raise filterError,"WARNING: Both filter keywords in image header are clear."

    #if filt not in wfc3.clear: # and filt2 not in wfc3.clear: # WZ
    #    raise filterError,"WARNING: Neither filter keyword in image header is clear."

    #if filt not in wfc3.clear:
    if filt not in wfc3.filters:
        filter = filt
    filter = filt # WZ
    del wfc3
    return filter
