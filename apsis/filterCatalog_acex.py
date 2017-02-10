#!/usr/bin/env python

# $Id: filterCatalog.py,v 1.11 2005/02/01 21:36:56 anderson Exp $
# ---------------------------------------------------------------------
# subclass of catalog.Catalog class.  Provides the run2() method to 
# allow the caller to run SExtractor in dual image mode.  Also
# provides a main function to run the whole thing and set up the
# parameters which the ACS group has identified as the nominal
# input parameter set by using provided weight images and RMS maps.
# K Anderson 

__version__      = '$Revision: 1.11 $ '[11:-3]
__version_date__ = '$Date: 2005/02/01 21:36:56 $ '[7:-3]
__author__       = "Ken Anderson, anderson@pha.jhu.edu"


import os
import string
import popen2
import numpy.oldnumeric as Numeric
import tableio
from   catalog import Catalog
from   pUtil   import errorSearch,fillHeader,ptime
from   fUtil   import zeroPoint

class filterCatalog(Catalog):
    """
    Example usage:

    Filter Catalog creation:
    ========================

    se = filterCatalog.filterCatalog(obs)             # Create an instance
                              # obs is the DataSet object (instance of DataSet)
    se.setinpar('FILTER_NAME', 'gauss_2.0_3x3.conv')  # change a parameter from .inpar
    se.setpars()                                      # generates the .inpar and .param files in obsdir/par
    se.run2()                                         # Run source extractor--dual image mode. 
    se.writeXmlCatalog()                              # writes .xml files from .cat files in obsdir/cat
    se.logfile.write(string)                          # write an entry to the logfile


    NOTES:
    - the _magFix private function is in this module and called within the run
      method.  It is only available for the filterCatalogs which is why it is here.
    - use of the inherited run() method for the filterCatalog instance is deprecated.
      All filter catalogs should be produced by dual image
      execution. Only when one filter's data is available should run() be called.
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
    - writeXml method now writes a .xml file from the .cat file

    """

    def __init__(self,obs,sparseField=None,excludeList=["X_IMAGE", "Y_IMAGE"]):
        """adjusted constructor to allow a passed excludeList for the 
        exclusion of a passed list of fields from the xml markup.
        """

        Catalog.__init__(self,obs,sparseField,excludeList)
    

    def run2(self):
        """Run SExtractor on the sciImageList in dual image mode (i.e. detection image vs filter image)"""
        if os.path.isfile(self.detectionImage):
            pass
        else:
            errtxt="This method is not appropriate. No Detection Image found for this observation."
            self.errorList.append((self.modName,errtxt))
            raise IOError,errtxt
        curdir = os.getcwd()
        os.chdir(self.obsFits)
        self.catalogList = []   
        for fitsfile in self.sciImageList:
            print "Generating source catalog for",fitsfile
            self.logfile.write('running SExtractor on ' + fitsfile)
            fitsname        = string.split(fitsfile,'.')[0]
            inParFileName   = os.path.join(self.obsPars, fitsname + '.inpar')
            catalogFileName = os.path.join(self.obsCats, fitsname + '.cat')
            cmd = 'sex ' + self.detectionImage + ' ' + os.path.join(self.obsFits, fitsfile) \
                      + ' -c ' + inParFileName
            try:
                subproc = popen2.Popen4(cmd)
                stderr_lines = subproc.fromchild.readlines() 
                # call the errorSearch ufunc and trim out the junk.
                # errorSearch returns a dictionary of unique error strings and the
                # number of occurrences found in the passed list (errs).

                foundErrs = errorSearch(stderr_lines)  # ufunc errorSearch from pUtil
                if foundErrs:
                    self.logfile.write('SExtractor produced the following message(s) on ' + fitsfile)
                    for f in foundErrs.keys():
                        self.logfile.write(f+": occurred "+str(foundErrs[f])+" times")
                        self.errorList.append((self.modName,f+": occurred "+str(foundErrs[f])+" times"))
            except Exception, err:
                self.logfile.write('run2 method encountered a problem: ')
                self.logfile.write(str(err))
                self.errorList.append((self.modName,str(err)))
                raise Exception, err
            if not os.path.isfile(catalogFileName):
                    raise RuntimeError,catalogFileName+" was not written!"
            self.logfile.write("Calling _magFix method")
            self.logfile.write("Fixing bad mags in "+os.path.split(catalogFileName)[1])
            self._magFix(catalogFileName)
            self.catalogList.append(catalogFileName)
            self.outputList[os.path.basename(catalogFileName)] = [fitsfile]
        os.chdir(curdir)
        return

    def _magFix(self,catalogFile):
        """This private method receives a path to a catalog file and sifts through the
        MAGERR field looking for values > 10.  It sets the corresponding MAG field = -99 and
        sets that object's MAGERR field to 0.0.  catalogFile is a path not a file object."""

        # fillHeader will return a list of tuples where which looks like
        # 
        # [(1, 'NUMBER'),
        # (2, 'X_IMAGE'),
        # (3, 'Y_IMAGE'),
        # ...
        # (12, 'MAG_ISOCOR'),
        # (13, 'MAGERR_ISOCOR'),
        # (14, 'FLUX_APER', 1)
        # (15, 'FLUX_APER', 2),
        # (16, 'FLUX_APER', 3),
        # ...
        # ]
        #
        # The tuples are either of length 2 or 3.  If len is 3, the 3rd item of the
        # tuple is the nth occurance of that column identifier.  This occurs on those
        # columns of MAGs and MAGERRs for a series of increasingly larger apertures. 

        # newFieldList will be a list of Numeric arrays containing the columns of the catalogs.
        # This list will contain fields which have not been altered, i.e. all fields other than
        # MAG_* and MAGERR_*, and the new MAG and MAGERR fields which have been corrected.
        # Once the list is complete, it is tuple-ized and send to the tableio pu_data function.

        newFieldList     = []
        newMagsList      = []
        newMagErrsList   = []
        newMagHeaders    = []
        newMagErrHeaders = []
        newHeaders       = []
        magCols          = []
        magErrCols       = []
        selectSet        = fillHeader(catalogFile)

        print "Searching catalog for required columns, MAG, MAGERR"
        for i in range(len(selectSet)):
            if len(selectSet[i]) == 2:
                column,name = selectSet[i]
                paramNames = name.split("_")
                if "MAG" in paramNames:
                    magCols.append((column,name))
                elif "MAGERR" in paramNames:
                    magErrCols.append((column,name))
                else:
                    oldField = tableio.get_data(catalogFile,(column-1))
                    newFieldList.append(oldField)
                    newHeaders.append(name)
                    continue
            else:
                column,name,id = selectSet[i]
                paramNames = name.split("_")
                if "MAG" in paramNames:
                    magCols.append((column,name,id))
                elif  "MAGERR" in paramNames:
                    magErrCols.append((column,name,id))
                else:
                    oldField = tableio.get_data(catalogFile,(column-1))
                    newFieldList.append(oldField)
                    newHeaders.append(name)
                    continue

        # We now have
        #  catalog field  --> list
        # --------------------------------
        #        MAG_*    --> magCols
        #     MAGERR_*    --> magErrCols
        #
        # The algorithm will be to step through the magErrCols columns, extracting those fields
        # via get_data and getting Numeric arrays.  The matching mag columns are slurped as well.
        # We search the magErrCols arrays looking for >= 10 values and then marking the those mags
        # as -99.0 and the matching magerrs as 0.0
        # See Bugzilla bug #2700

        for item in magErrCols:
            magErrAperId  = None
            # item may be of len 2 or 3
            if len(item) == 2:
                magErrColId,magErrColName = item
            else:
                magErrColId,magErrColName,magErrAperId = item

            magErrKind = magErrColName.split("_")[1]           # ISO, ISOCORR, etc.

            print "\n\nMAG type:",magErrKind
            if magErrAperId: print magErrColName,"Aper id is",magErrAperId
            print "Getting\t",magErrColName,"\tfield",magErrColId

            # MAGERR array:
            magErrs = tableio.get_data(catalogFile,magErrColId -1)

            matchingMagColName    = None
            matchingMagColId      = None

            #----------------------- Search for matching MAG_* field -----------------------#

            for magitems in magCols:                

                # We know that the magErrColName is MAGERR and if magErrNameId is true then
                # the tuple is of len 3, i.e. a MAGERR_APER field.  We look for the matching
                # MAG_APER field id, 1, 2, 3... etc.

                if len(magitems) == 3:
                    magColId,magColName,magAperId = magitems
                    if magColName == "MAG_"+magErrKind:
                        matchingMagColName = magColName
                        #print "Found matching field type:",magColName,"in field",magColId
                        if magAperId == magErrAperId:
                            print "Found matching aperture id."
                            print "MAG_APER id: ",magAperId,"MAGERR_APER id: ",magErrAperId
                            matchingMagColId = magColId
                            matchingMags = tableio.get_data(catalogFile,magColId-1)
                            break
                    else:
                        continue
                else:
                    magColId,magColName = magitems
                    if magColName == "MAG_"+magErrKind:
                        print "Found matching field type:",magColName,"in field",magColId
                        matchingMagColName = magColName
                        matchingMagColId   = magColId
                        matchingMags = tableio.get_data(catalogFile,magColId-1)
                        break
                    else:
                        continue

            #--------------------------------------------------------------------------------#


            print " MAG err field:",magErrColName,magErrColId
            print "     Mag field:",matchingMagColName,matchingMagColId

            # Now the grunt work on the arrays,
            # magErrs, matchingMags
            #
            # update: flagging all MAGs as -99 when the corresponding MAGERR > 10
            # introduced a bug which unintentionally reset the magnitudes
            # SExtractor had flagged with a MAG = 99.0 and a MAGERR = 99.0
            # This now checks for a MAGERR of 99 and does not reset the MAG value
            # if MAGERR = 99.0 but does for all other MAGERRS > 10.0

            badMagErrs1 = Numeric.where(magErrs >= 10,1,0)
            badMagErrs2 = Numeric.where(magErrs != 99.0,1,0)
            badMagErrs  = badMagErrs1 * badMagErrs2
            del badMagErrs1,badMagErrs2
            newMags     = Numeric.where(badMagErrs,-99.0,matchingMags)
            newMagErrs  = Numeric.where(badMagErrs,0.0,magErrs)

            newMagsList.append(newMags)
            newMagHeaders.append(matchingMagColName)
            newMagErrsList.append(newMagErrs)
            newMagErrHeaders.append(magErrColName)

        # concatenate the lists.  This is done to preserve the MAG_APER and MAGERR_APER
        # grouping of the original SExtractor catalog.

        newFieldList = newFieldList + newMagsList
        newFieldList = newFieldList + newMagErrsList
        newHeaders   = newHeaders   + newMagHeaders
        newHeaders   = newHeaders   + newMagErrHeaders

        newVariables = tuple(newFieldList)

        # rename the old catalog file as catalogFile.old
        os.rename(catalogFile,catalogFile+".old")
        self.outputList[os.path.basename(catalogFile)+".old"] = [os.path.basename(catalogFile)]
        fob = open(catalogFile,'w')
        fob.write("## " + ptime() + "\n")
        fob.write("## " + self.modName + " catalog regenerated by _magFix method.\n")
        fob.write('## (This file was generated automatically by the ACS Pipeline.)\n##\n')
        fob.write("## This catalog has been photometrically corrected to remove\n")
        fob.write("## 'bad' magnitude values.\n")
        fob.write("##\n")
        for i in range(len(newHeaders)):
            fob.write("# "+str(i+1)+"\t"+newHeaders[i]+"\n")
        fob.close()
        tableio.put_data(catalogFile,newVariables,append="yes")

        return


def main(obs,wtfile):
    """Just cobbles all the method calls and runs it all with one call.
    """

    fcat = filterCatalog(obs)

    ### The responsibility for running SExtractor using the rms and weight images
    ### falls here...
    ### Here we set the parameters to use the weight map of the detection image
    ### detIm.detWgtName and the RMS map of the filters which came out of 
    ### combDither (see obs.rmsImageList).  Each filter has an associated
    ### rms map.  We have to set each filter's parameter file individually.
    ### This is using the reconfigured setpars method which now takes a fitsfile
    ### from the caller and writes the appropriate par file.

    fcat.setinpar('WEIGHT_TYPE','MAP_WEIGHT,MAP_RMS')
    fcat.setinpar('INTERP_TYPE','NONE')
    wgtfile = wtfile

    for i in range(len(fcat.sciImageList)):
        im      = fcat.sciImageList[i]
        rms_im  = fcat.rmsImageList[i]
        rmsfile = rms_im
        fcat.setinpar('WEIGHT_IMAGE',wgtfile+','+rmsfile)
        fcat.setpars(im)
    fcat.run2()
    fcat.writeXml()
    fcat.mkMsg()           # write the module message 
    return
