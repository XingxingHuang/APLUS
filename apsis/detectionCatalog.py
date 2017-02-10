#!/usr/bin/env python

# $Id: detectionCatalog.py,v 1.5 2003/10/02 20:43:40 anderson Exp $
# ---------------------------------------------------------------------
# this module defines the detectionCatalog class which subclassed from
# the Catalog class (see the catalog module in this package).
# Overrides the inherited setpars() and run() methods.
# Function main() can be called directly with the caller
# passing a DataSet object and a weight image file.
# In the pipeline context, the weight image will usually be the
# detectionWeight image which the combFilter module has produced.
# K Anderson
# 9-11-01

__version__      = '$Revision: 1.5 $ '[11:-3]
__version_date__ = '$Date: 2003/10/02 20:43:40 $ '[7:-3]
__author__       = "K Anderson <anderson@pha.jhu.edu"


import os
import string
import popen2
import pUtil,xmlUtil
from   catalog import Catalog
import pdb

class detectionCatalog(Catalog):
    """
    Inherits from the base class Catalog.  
    Overrides the setpars and run methods

    Example usage:

    Detection Catalog creation:
    ===========================
    import detectionCatalog
    dcat = detectionCatalog.detectionCatalog(obs)     # Create an instance
                              # obs is the observation object (instance of ObsDir)
    dcat.setpars()
    dcat.run()                                        # make the detection catalog

    """
    def setpars(self):
        """override the setpars method and set the par files for the detection Image.
        """

        self.fitsname = string.split(os.path.basename(self.detectionImage),'.')[0]
        #self.inpars   = pUtil.readInParFile(self.defaultInParFile)             
        self.logfile.write('Generating parameter sets for detection Image...')
        self.inParFileList = []

        # This section is removing Flag Image parameters from the output parameter list (outpars) and the 
        # FLAG_IMAGE file name parameter in the input paramaters list (inpars).

        try:
            del self.inpars['FLAG_IMAGE']
        except:
            pass

        # Here we have to deal with the fact that the CHECKIMAGE_TYPE value could be a comma
        # seperated string of several actual image types, eg., SEGMENTATION,OBJECTS,APERTURES, etc.
        # We have to split up this field and then create CHECKIMAGE_NAMES for each one.

        checkImageTypes      = string.split(self.getinpar('CHECKIMAGE_TYPE'),",")   # a list of image types
        self.checkimagesList = []

        filter             = self.getinpar('FILTER_NAME')
        starnnw            = self.getinpar('STARNNW_NAME')
        catalogFileName    = os.path.join(self.obsCats, self.fitsname + '.cat')
        self.inParFileName = os.path.join(self.obsPars, self.fitsname + '.inpar')
        paramFileName      = os.path.join(self.obsPars, self.fitsname + '.param')

        self.setinpar('CATALOG_NAME'   , catalogFileName)
        self.setinpar('PARAMETERS_NAME', paramFileName)
        self.setinpar('FILTER_NAME'    , os.path.join(self.configdir, filter ))
        self.setinpar('STARNNW_NAME'   , os.path.join(self.configdir, starnnw))

        # Here we create a string of new output checkimage file names.
        # Add each one of these new name to the images list and the output list.
        checkimage_names = ''
        for im in checkImageTypes:
            imFileName        = self.fitsname + '_' + im[:4] + '.fits'
            checkimage_names += imFileName+','
            self.checkimagesList.append(os.path.join(self.obsFits,imFileName))
            self.outputList[os.path.basename(imFileName)] = [os.path.basename(self.detectionImage)]
        self.setinpar('CHECKIMAGE_NAME', checkimage_names[:-1])                  # don't write the trailing comma

        inparfile  = open(self.inParFileName, 'w')
        outparfile = open(paramFileName, 'w')
        
        for param, value in self.inpars.items():
            inparfile.write("%-18s   %-20s\n" % (param, value))            
        for f in range(len(self.outpars)):
            param,value = self.outpars[f]
            if value:
                outparfile.write(param + '\n')
                
        inparfile.close()
        outparfile.close()
        self.inParFileList.append(self.inParFileName)
        self.logfile.write('parameter files for ' + self.fitsname + ' written.')
        return

    def run(self):
        """override the run method of the Catalog class to run the detection image only."""
        #pdb.set_trace()
        curdir = os.getcwd()
        os.chdir(self.obsFits)
        self.catalogList = []
        # reset the input list which gets set in the superclass constructor
        # to the the list of science and flag images.
        self.inputList = []
        catalogFileName  = os.path.join(self.obsCats, self.fitsname + '.cat')
        self.logfile.write('running SExtractor on ' + self.fitsname)
        cmd = 'sex ' + self.detectionImage +' -c ' + self.inParFileName
        self.inputList.append(os.path.basename(self.detectionImage))
        try:
            subproc = popen2.Popen4(cmd)
            stderr_lines = subproc.fromchild.readlines() 
            # call the errorSearch ufunc and trim out the junk.
            # errorSearch returns a dictionary of unique error strings and the
            # number of occurrences found in the passed list (errs).
            foundErrs = pUtil.errorSearch(stderr_lines)
            if foundErrs:
                self.logfile.write('SExtractor produced the following message(s) on ' + self.fitsname)
                for f in foundErrs.keys():
                    self.logfile.write(f+": occurred "+str(foundErrs[f])+" times")
                    self.errorList.append((self.modName,f+": occurred "+str(foundErrs[f])+" times"))
        except Exception, err:
            self.logfile.write('SExtractor failed:')
            self.errorList.append((self.modName,str(err)))
            raise Exception, err
        if os.path.isfile(catalogFileName):
            self.catalogList.append(catalogFileName)
            self.outputList[os.path.basename(catalogFileName)] = [os.path.basename(self.detectionImage)]
        self._hackit(catalogFileName)
        os.chdir(curdir)
        return

    def _hackit(self,cat,keep_apertures=[1,2,3]):
        """hack the detectionCatalog.cat file to take out a bunch of the aperture data. 
        Default is only to keep the first three apertures in the final catalog but caller
        can change this by passing a keep_apertures list (this is aperture number and *not*
        the radius).  This will hack the columns indicated by 

        MAG_APER
        MAGERR_APER
        FLUX_APER
        FLUXERR_APER
        """
        dir,old_file = os.path.split(cat)
        headerList = []
        headerList = pUtil.fillHeader(cat)       # this returns a list of the catalog header.

        # go through the header and find the columns to keep.  We are looking for
        #'FLUX_APER', 1)
        #'FLUX_APER', 2)
        #'FLUX_APER', 3)
        #'FLUXERR_APER', 1)
        #'FLUXERR_APER', 2)
        #'FLUXERR_APER', 3)
        #'MAG_APER', 1)
        #'MAG_APER', 2)
        #'MAG_APER', 3)
        #'MAGERR_APER', 1)
        #'MAGERR_APER', 2)
        #'MAGERR_APER', 3)

        newheader = []
        for i in headerList:
            if len(i) == 2:
                newheader.append(i)
            else:
                if i[2] not in keep_apertures:
                    continue
                else:
                    newheader.append(i)

        #return newheader

        cols = []
        for i in newheader:
            cols.append(i[0]-1)

        new_rows = []
        for row in open(cat).readlines():
            if '#' in row: continue
            fields = row.split()
            arow = ''
            for column in cols:
                arow += '  '+fields[column]
            new_rows.append(arow)

        # OK, we have the newheader and the new data . 
        # We need to renumber the columns in the header.  the newheader
        # list has the old catalog's column identifiers and that needs to get fixed.

        new_newheader = []
        for i in range(len(newheader)):
            if len(newheader[i]) == 2:
                new_newheader.append((i+1,newheader[i][1]))
            else:
                new_newheader.append((i+1,newheader[i][1],newheader[i][2]))

        # Now we are just going to overwrite the original detectionImage.cat file
        # (well, at least whatever was passed to this func, anyway)

        file = open(cat,'w')
        self.logfile.write("private method _hackit, trimming aperture parameters")
        self.logfile.write("_hackit overwriting detectionImage.cat file.")
        file.write("## Date: "+pUtil.ptime()+"\n")
        file.write("## This file has been modified from its original form by the WFP Pipeline.\n")
        file.write("## Some aperture fields have been removed.\n")
        file.write("## This file written by the WFP Pipeline.  Do not edit.\n")
        for item in new_newheader:
            file.write('# '+str(item[0])+'\t'+str(item[1])+'\n')
        for row in new_rows:
            file.write(row+'\n')
        file.close()
        return

def main(obs,wtfile):
    """Just run it all and set the appropriate parameters for weighting image source detection."""

    dcat = detectionCatalog(obs)
    pdb.set_trace()
    ### The responsibility for running SExtractor using detection weight images
    ### falls here...The pipeline could just call this main function rather
    ### doing what it does now, which is just what is below.
    ### i.e. the same effect can be produced by the pipeline by this
    ### import detectionCatalog
    ### detectionCatalog.main(obs,some_detectionweight_image)
    if not os.path.isfile(wtfile):
        wgtfile = os.path.join(dcat.obsFits,wtfile)
    else:
        wgtfile = wtfile
    dcat.setinpar('WEIGHT_TYPE','MAP_WEIGHT')
    dcat.setinpar('WEIGHT_IMAGE',wgtfile)
    dcat.setinpar('INTERP_TYPE','NONE')
    dcat.setpars()
    dcat.run()
    dcat.writeXml()
    dcat.mkMsg()           # write the module message 
    return
