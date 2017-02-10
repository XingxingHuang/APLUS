#!/usr/bin/env python

# $Id: fUtil.py,v 2.10 2012/07/01 21:30:04 zheng Exp $
# ---------------------------------------------------------------------
# ACS pipeline generic interface to some kind of image package which
# is pyfits right now.
#
__version__      = '$Revision: 2.10 $ '[11:-3]
__version_date__ = '$Date: 2012/07/01 21:30:04 $ '[7:-3]
__author__       = "Wei Zheng, zheng@pha.jhu.edu"

import os,re,types
import glob
import string
import math
import pyfits,numpy      # This is pyfits in pipeline/lib/python
import filters
import pdb

class filterError(Exception):
    """raise a filter error exception."""

def getNXNY(fitsfile):
    """Return the (NAXIS1,NAXIS2) keyword values as a tuple from FITS header.
       Uses ext=0 for simple fits, ext=1 for mef"""

    f  = pyfits.open(fitsfile)
    if len(f) == 1:
        ext=0
    else:
        ext=1
    NX = f[ext].header.get('NAXIS1')
    NY = f[ext].header.get('NAXIS2')
    f.close()
    del f
    return (NX,NY)


def getNX(fitsfile,ext=0):
    """Return the NAXIS1 keyword value from FITS header.
       Defaults to first extension."""

    f  = pyfits.open(fitsfile)
    NX = f[ext].header.get('NAXIS1')
    f.close()
    del f
    return NX


def getNY(fitsfile,ext=0):
    """return the NAXIS2 keyword value from FITS header.
       Defaults to first extension."""

    f  = pyfits.open(fitsfile)
    NY = f[ext].header.get('NAXIS2')
    f.close()
    del f
    return NY


def getTel(fitsfile,ext=0):
    """returns the value of TELESCOP keyword in a fits header."""
    f   = pyfits.open(fitsfile)
    tel = f[ext].header.get('TELESCOP')
    f.close()
    del f
    return tel
    

def getInstr(fitsfile,ext=0):
    """returns the value of INSTRUME keyword in a fits header."""
    f   = pyfits.open(fitsfile)
    instr = f[ext].header.get('INSTRUME')
    f.close()
    del f
    return instr
    

def getDetector(fitsfile,ext=0):
    """returns the value of DETECTOR keyword in a fits header."""
    f   = pyfits.open(fitsfile)
    det = f[ext].header.get('DETECTOR')
    f.close()
    del f
    return det

def getTarget(fitsfile,ext=0):
    """returns the value of TARGNAME keyword in a fits header."""
    f   = pyfits.open(fitsfile)
    det = f[ext].header.get('TARGNAME')
    f.close()
    del f
    return det
    
def getFilter(fitsfile,ext=0):
    """returns the value of FILTER keyword in a fits header."""
    f   = pyfits.open(fitsfile)
    filt= f[ext].header.get('FILTER')
    f.close()
    del f
    return filt

def getFilter1(fitsfile,ext=0):
    """returns the value of FILTER1 keyword in a fits header."""
    f   = pyfits.open(fitsfile)
    filt1= f[ext].header.get('FILTER1')
    f.close()
    del f
    return filt1


def getFilter2(fitsfile,ext=0):
    """returns the value of FILTER2 keyword in a fits header."""
    f   = pyfits.open(fitsfile)
    filt2= f[ext].header.get('FILTER2')
    f.close()
    del f
    return filt2
    

def getIdcTab(fitsfile,ext=0):
    """return a full path to the IDCTAB file specified in the header of the
    passed fitsfile.  Will interpret jref$ iraf env var if it is present in 
    the value of the IDCTAB keyword.  This keyword must be in the primary header
    or else raises a FITS_SevereError exception.
    """
    f   = pyfits.open(fitsfile)
    tab = f[ext].header.get('IDCTAB')
    if tab[0:5] == 'jref$':
        jref = os.environ['jref']
        file = tab[5:]
        idcfile = jref + file
    elif tab[0:5] == 'iref$':
        iref = os.environ['iref']
        file = tab[5:]
        idcfile = iref + file
    elif os.path.isabs(tab):
        idcfile = tab
    else:
        raise pyfits.FITS_SevereError,"Cannot locate IDCTAB"
    if not os.path.isfile(idcfile):
        raise IOError,"Specified IDCTAB file, "+idcfile+" does not exist."
    else:
        return idcfile

def getChipId(fitsfile,ext,chipname='CCDCHIP'):
    """return the chip id found for a keyword which identifies the
    chip id in an extension header. For ACS this is CCDCHIP but for WFPC2 
    it is DETCHIP.  This function requires the extension number.
    returns None if not found.
    """
    f = pyfits.open(fitsfile)
    val = f[ext].header.get(chipname)
    f.close()
    del f
    return val

def getKeyVal(fitsfile,key,ext=0):
    """get any keyword value from a specified (default = 0) extension of
    a fits file."""
    f   = pyfits.open(fitsfile)
    val = f[ext].header.get(key)
    f.close()
    del f
    return val
    

def pixelScale(fitsfile,ext=0):
    "Return the pixel scale from the WCS info of the passed fitsfile."
    f = pyfits.open(fitsfile)
    pixscale = None
    # make sure we understand the units first
    ctype1 = f[ext].header.get('CTYPE1')
    ctype2 = f[ext].header.get('CTYPE2')
    if(ctype1 == None or ctype2 == None):
        raise KeyError,'Could not find CTYPE Keys in '+fitsfile
    ctype1.strip()
    ctype2.strip()
    if(ctype1 == 'RA---TAN' and ctype2 == 'DEC--TAN'):
        cd_11 = float(f[ext].header.get('CD1_1'))
        cd_12 = float(f[ext].header.get('CD1_2'))
        cd_21 = float(f[ext].header.get('CD2_1'))
        cd_22 = float(f[ext].header.get('CD2_2'))
        pixscale = math.sqrt((cd_11**2 + cd_12**2 + cd_21**2 + cd_22**2)/2.0)
        pixscale *= 3600.
    del f
    return pixscale

def filterResponseName(fitsfile,ext=0):
    """returns a generalised filter name of the form 'telescope_instrument_detector_filter'
    Filter will generally be one of FILTER1 or FILTER2 when one of the filter keywords is 
    populated with a 'CLEARxxx' value. This employs the ACSFilterInfo class from the filters
    module. This does not employ other fUtil functions to avoid opening and closing fits objects
    for the same fitsfile multiple times.
    """
    f     = pyfits.open(fitsfile)
    tel   = f[ext].header.get('TELESCOP')
    instr = f[ext].header.get('INSTRUME')
    det   = f[ext].header.get('DETECTOR')
    filter = filterName(fitsfile)
    #print "Please check the filter for the following commands!  I think it will be OK! ___by xingxing"
    #pdb.set_trace()
    acs = filters.ACSFilterInfo()
    if filter:
        if filter in acs.ramp:
            lrfwave = int(f[ext].header.get('LRFWAVE'))
            fname = tel+'_'+instr+'_'+det+'_'+filter+'_'+str(lrfwave)
        else:
            fname = tel+'_'+instr+'_'+det+'_'+filter
    else:
        raise filterError,"No filter found for fits image: "+os.path.basename(fitsfile)
    return fname

def filterName(fitsfile,ext=0):             #XX   change this module to be able to read both ACS and WFC3 filters.
    """returns a the filter name of the real (i.e. not clear) filter used
    in the image.  Function determines this like the above filterResponseName func
    but only returns the filter name.
    """
    f     = pyfits.open(fitsfile)
    filt1 = f[ext].header.get('FILTER1')
    filt2 = f[ext].header.get('FILTER2')

    if not filt1 and not filt2:
        # print "You are runing with the WFC3 filter! (filterName in fUtil.py) xingxing"#by xingxing
        #pdb.set_trace()
	filter12 = f[ext].header.get('FILTER')
        f.close()
	del f
        return filter12                #by xingxing

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

def twoFilterName(fitsfile,ext=0):
    """returns a full filter name of the real filter used in the image.
    This function should be used in conjunction with the filterName fucntion
    when it determines that both filters are populated by non-clear filter
    values.  This will then return a filter name like
    HST_ACS_HRC_POL60UV_F435W
    and it is up to the user to use this as they need. Note, there is no checking
    of the filter values here.  If one filter is clear, the caller will get it
    back, like
        HST_ACS_HRC_CLEAR1S_F250W
    """
    f     = pyfits.open(fitsfile)
    #print "Please check the FILTER:_________by xingxing"
    #pdb.set_trace()  # by XX   different for ACS and WFC3
    if f[ext].header.get('INSTRUME')=="ACS":
        filt1 = f[ext].header.get('FILTER1')
        filt2 = f[ext].header.get('FILTER2')
        tel   = f[ext].header.get('TELESCOP')
        instr = f[ext].header.get('INSTRUME')
        det   = f[ext].header.get('DETECTOR')
        fname = tel+'_'+instr+'_'+det+'_'+filt1+'_'+filt2
    elif  f[ext].header.get('INSTRUME') == "WFC3":
	print "Using the WFC3 (twoFilterName in fUtil.py)! by xingxing" 
        filt = f[ext].header.get('FILTER')
        tel   = f[ext].header.get('TELESCOP')
        instr = f[ext].header.get('INSTRUME')
        det   = f[ext].header.get('DETECTOR')
        fname = tel+'_'+instr+'_'+det+'_'+filt
    else :
        print "ERROR: Please check the instrument you are using!!!(in fUtil.py)   ==ERROR==by xingxing"
        pdb.set_trace()  # by XX
    f.close()
    del f
    return fname

def zeroPoint(fitsfile):
    """extract zero point information from header of a fitsfile. return the zero point.  
    Function will determine whether the image is a COUNT image or a COUNTRATE image from 
    the BUNIT header keyword.
    """
    f = pyfits.open(fitsfile)
    flt = filterName(fitsfile)
    
    # the PHOTZPT, PHTOPLAM, and PHOTFLAM keywords
    # This will only be done for simple fits images right now.

    if len(f) > 1:
        print "WARNING: This is not a simple fits file.\n Defaulting to 1st extension of "+fitsfile
        try:    pzpt  = f[1].header.get('PHOTZPT')
        except: raise pyfits.FITS_SevereError,"ERROR, can't get PHOTZPT keyword value from "+fitsfile

        try:    pflam = f[1].header.get('PHOTFLAM')
        except: raise pyfits.FITS_SevereError, "ERROR, can't get PHOTFLAM keyword value from "+fitsfile

        try:    pplam = f[1].header.get('PHOTPLAM')
        except: raise pyfits.FITS_SevereError,"ERROR, can't get PHOTPLAM keyword value from "+fitsfile

        try:    bunit = string.strip(f[1].header.get('BUNIT'))
        except: raise pyfits.FITS_SevereError,"ERROR, can't get BUNIT keyword value from "+fitsfile

        try:    exptime = f[1].header.get('EXPTIME')
        except: raise pyfits.FITS_SevereError,"ERROR, can't get EXPTIME keyword value from "+fitsfile

        
    else:
        try:    pzpt  = f[0].header.get('PHOTZPT')
        except: raise pyfits.FITS_SevereError,"ERROR, can't get PHOTZPT keyword value from "+fitsfile

        try:    pflam = f[0].header.get('PHOTFLAM')
        except: raise pyfits.FITS_SevereError,"ERROR, can't get PHOTFLAM keyword value from "+fitsfile

        try:    pplam = f[0].header.get('PHOTPLAM')
        except: raise pyfits.FITS_SevereError,"ERROR, can't get PHOTPLAM keyword value from "+fitsfile

        try:    bunit = string.strip(f[0].header.get('BUNIT'))
        except: raise pyfits.FITS_SevereError,"ERROR, can't get BUNIT keyword value from "+fitsfile

        try:    exptime = f[0].header.get('EXPTIME')
        except: raise pyfits.FITS_SevereError,"ERROR, can't get EXPTIME keyword value from "+fitsfile


    # This calculation is defined by Meurer email of 3-Oct-2001
    # at http://acs.pha.jhu.edu/general/team-only/email/acs-pipeline-dev/0075.shtml)
    # calculation revised according to new spec in Bugzilla bug #921.
    #
    # specification of count units and count rate units which can be used for
    # this zero point calculation comes from Bugzilla bug #1085.

	#by __________-xingxing  Add some elements!
    count_units     = ['DN','COUNTS','COUNT','ELECTRON','ELECTRONS','PHOTON','PHOTONS','electron','electrons']
    countrate_units = ['DN/S','COUNTS/S','COUNT/S','ELECTRON/S','ELECTRONS/S','PHOTON/S','PHOTONS/S','electron/s','electrons/s']

#     if bunit in count_units:
#   try:
#       zpt = pzpt - 2.5*math.log10(pflam/exptime) - 5.0*math.log10(pplam/5475.4)
#   except:
#       print "WARNING: Could not determine zeropoint. Returning default value."
#       zpt = 30
#     elif bunit in countrate_units:
#   try:
#       zpt = pzpt - 2.5*math.log10(pflam) - 5.0*math.log10(pplam/5475.4)
#   except:
#       print "WARNING: Could not determine zeropoint. Returning default value."
#       zpt = 30
#     else:
#   raise pyfits.FITS_SevereError,"ERROR: unrecognized BUNIT keyword value in  "+fitsfile+": "+bunit

#     return zpt


# Below is a modified bunit condidtional.
# old version is commented above
# -  f. menanteau

    if bunit in count_units:
        try:
            zpt = pzpt - 2.5*math.log10(pflam/exptime) - 5.0*math.log10(pplam/5475.4)
        except:
            print "WARNING: Could not determine zeropoint. Trying the table."
            zpt = getTableZeroPoint(fitsfile) + 2.5*math.log10(exptime)
            print "WARNING: Using zeropoint from table for %s, set to: %f" % (flt,zpt)
        #else:
     #   print "WARNING: Could not determine zeropoint. Returning default value."
     #   zpt = 30
    elif bunit in countrate_units:
        try:
            zpt = pzpt - 2.5*math.log10(pflam) - 5.0*math.log10(pplam/5475.4)
        except:
            print "WARNING: Could not determine zeropoint. Trying the table."
            zpt = getTableZeroPoint(fitsfile)
            print "WARNING: Using zeropoint from table for %s, set to: %f" % (flt,zpt)
        #else:
        #   print "WARNING: Could not determine zeropoint. Returning default value."
        #   zpt = 30
    else:
        raise pyfits.FITS_SevereError,"ERROR: unrecognized BUNIT keyword value in  "+fitsfile+": "+bunit

    return zpt

def zeroPointCorrection(filter):
    """Temporary corrections to the zeropoints until new ones are provided by STScI"""
    corr={
    'F220W':       1.255,   # this nothing but a rough extrapolation
    'F250W':       1.250,   # this nothing but a rough extrapolation
    'F330W':       1.245,   # this nothing but a rough extrapolation
    'F435W':       1.230,
    'F475W':       1.226,
        'F502N':       1.195,   # this nothing but a rough interpolation
    'F555W':       1.155,
    'F606W':       1.132,
    'F625W':       1.132,
        'F658N':       1.122,   # this nothing but a rough interpolation
    'F775W':       1.113,
        'F800L':       1.100,   # this nothing but a rough interpolation
    'F814W':       1.089,
    'F850LP':      1.016
    }
    zptc=corr[filter[12:]]
    return 2.5*math.log10(zptc)

def limitingMag(fitsfile):
    """determine the limiting magnitude of the passed fitsfile.
    """
    ### ??? don't know how to do this yet. ###
    ### limit = ???
    limit = 25.0
    return limit


def makeAsnDict(dir):
    """Returns a dictionary of asn table names with an embedded list of associated files.
       i.e. {'jsimh400t_asn.fits' : [ 'jsimh430t_crj.fits','jsimh420t_crj.fits','jsimh410t_crj.fits','jsimh440t_crj.fits'],
             'jsimh500t_asn.fits' : [ 'jsimh530t_crj.fits','jsimh520t_crj.fits','jsimh510t_crj.fits','jsimh540t_crj.fits'],
         'jsimh600t_asn.fits' : [ 'jsimh630t_crj.fits','jsimh620t_crj.fits','jsimh610t_crj.fits','jsimh640t_crj.fits'] }

    Determined from reading the asn files for an observation. arg is the directory where
    all files, asn and data files are located.
    """
    curdir = os.getcwd()
    os.chdir(dir)
    asnlist = glob.glob('*_asn.fits')
    asnlist.sort() #WZ
    obsdict = {}
    for table in asnlist:
        imlist = []
        tab = pyfits.open(table)
        for i in range(len(tab[1].data)):
            mtype = tab[1].data[i].field('MEMTYPE')
            if string.find(mtype,'PROD-CR') == 0:
                suffix = '_crj.fits'
                prefix = tab[1].data[i].field('MEMNAME')
                imlist.append(prefix+suffix)
            elif string.find(mtype,'EXP-DTH') == 0:              # These matching string is unknown at the moment
                suffix = '_flt.fits'
                prefix = tab[1].data[i].field('MEMNAME')
                imlist.append(prefix+suffix)
        obsdict[table] = imlist
    os.chdir(curdir)

    return obsdict

def makeAsnDict2(dir,filters):
    """Returns a dictionary of asn table names with an embedded list of associated files.
       i.e. {'jsimh400t_asn.fits' : [ 'jsimh430t_crj.fits','jsimh420t_crj.fits','jsimh410t_crj.fits','jsimh440t_crj.fits'],
             'jsimh500t_asn.fits' : [ 'jsimh530t_crj.fits','jsimh520t_crj.fits','jsimh510t_crj.fits','jsimh540t_crj.fits'],
         'jsimh600t_asn.fits' : [ 'jsimh630t_crj.fits','jsimh620t_crj.fits','jsimh610t_crj.fits','jsimh640t_crj.fits'] }

    Determined from reading the asn files for an observation. arg is the directory where
    all files, asn and data files are located.
    """
    curdir = os.getcwd()
    os.chdir(dir)
    asnlist = glob.glob('*'+filters+'_asn.fits')
    obsdict = {}
    for table in asnlist:
        imlist = []
        tab = pyfits.open(table)
        for i in range(len(tab[1].data)):
            mtype = tab[1].data[i].field('MEMTYPE')
            if string.find(mtype,'PROD-CR') == 0:
                suffix = '_crj.fits'
                prefix = tab[1].data[i].field('MEMNAME')
                imlist.append(prefix+suffix)
            elif string.find(mtype,'EXP-DTH') == 0:              # These matching string is unknown at the moment
                suffix = '_flt.fits'
                prefix = tab[1].data[i].field('MEMNAME')
                imlist.append(prefix+suffix)
        obsdict[table] = imlist
    os.chdir(curdir)

    return obsdict


def convertExtension(fitsfile,extension=1):
    """Takes a full path to a fits file, breaks out the passed extension into simple fits
    file and returns the full path name of this file to the caller.  This file will be written
    in the same dir as the passed fitsfile.  Returns None if the passed fitsfile does not have
    multiple extensions.
    """
    # check that the file name delivered is a full path

    if not os.path.isabs(fitsfile):
        raise IOError, "Path to image file is not a full path."

    warningsList = []

    oldfits = pyfits.open(fitsfile)
    if len(oldfits) == 1:
        warningsList.append("WARNING: File "+fitsfile+" is not multi-extension. No extension to convert.")
        print "File "+fname+" is not multi-extension. No extension to convert."
        return None
    dir,fname = os.path.split(fitsfile)
    base,ext  = os.path.splitext(fname)
    xname = oldfits[extension].header['EXTNAME'].upper()
    try:
        xver = str(oldfits[extension].header['EXTVER'])
    except:
        xver = ''
    if xver:
        newfile = base+"_"+xname+"_"+xver + ext
    else:
        newfile = base+"_"+xname + ext

    tfile   = os.path.join(dir,newfile)
    newfits = pyfits.HDUList()
    newfits.append(pyfits.PrimaryHDU())
    newfits[0].header.ascard = oldfits[0].header.ascard.copy()
    try:
        del newfits[0].header.ascard['BSCALE']
        warningsList.append("WARNING: BSCALE keyword found in primary header: "+fname)
        print "WARNING: BSCALE keyword found in primary header: "+fname
    except: pass
    try:
        del newfits[0].header.ascard['BZERO']
        warningsList.append("WARNING: BZERO keyword found in primary header: "+fname)
        print "WARNING: BZERO keyword found in primary header: "+fname
    except: pass

    # These are primary header keywords that need to disappear.

    try:    del newfits[0].header.ascard['EXTEND']
    except: pass
    try:    del newfits[0].header.ascard['NEXTEND']
    except: pass

    for card in newfits[0].header.ascardlist():
        if card.key == "HISTORY":
            newfits[0].header.ascard.remove(card)
            continue
    # delete some of the extension-only keywords

    try:    del oldfits[extension].header.ascard['XTENSION']
    except: pass
    try:    del oldfits[extension].header.ascard['INHERIT']
    except: pass
    try:    del oldfits[extension].header.ascard['EXTNAME']
    except: pass
    try:    del oldfits[extension].header.ascard['EXTVER']
    except: pass
    try:    del oldfits[extension].header.ascard['PCOUNT']
    except: pass
    try:    del oldfits[extension].header.ascard['GCOUNT']
    except: pass

    print "breaking out extension ",extension," of file ",fname

    for card in oldfits[extension].header.ascard:
        key = card.key
        val = card.value
        if not key:
            continue
        if key == "HISTORY" or key == "COMMENT":
            continue

        # override a primary keyword value with the extension value

        if newfits[0].header.has_key(key):
            newfits[0].header.update(key,oldfits[extension].header[key])
            continue
        if key == 'NAXIS1':
            try:   inaxis = newfits[0].header.ascard.index_of('NAXIS')
            except KeyError,err:
                raise pyfits.FITS_FatalError,"Cannot find NAXIS keyword!"
            newfits[0].header.ascard.insert(inaxis+1,card)
            continue
        elif key == 'NAXIS2':
            try:       inaxis2 = newfits[0].header.ascard.index_of('NAXIS1') + 1
            except:
                try:   inaxis2 = newfits[0].header.ascard.index_of('NAXIS') + 1
                except KeyError,err:
                    raise pyfits.FITS_FatalError,"Cannot find NAXIS keyword!"
            newfits[0].header.ascard.insert(inaxis2,card)
            continue
        else:
            newfits[0].header.ascard.append(card)

    newfits[0].data = oldfits[extension].data
    newfits[0].header.update('FILENAME',newfile)
    newfits.writeto(tfile)
    oldfits.close()
    del oldfits,newfits
    return tfile,warningsList


def convertFile(fitsfile):
    """Takes a full path to a fits file and converts it to simple fits format if it can.
    Names the converted extensions to be like
       
    j1570d10t_crj.fits ===> NEXTEND number of files named like 

          j1570d10t_crj_EXTNAME_EXTVER.fits

    """
    warningsList = []

    if not os.path.isabs(fitsfile):
        raise IOError, "Path to image file is not a full path."
    dir,fname = os.path.split(fitsfile)
    base,ext  = os.path.splitext(fname)
    
    fo = pyfits.open(fitsfile)

    if len(fo) == 1:
        warningsList.append("WARNING: File "+fname+" is not multi-extension. Not converted.")
        print "File "+fname+" is not multi-extension. Not converted."
        return None

    # first check that BSCALE and BZERO are *not* in the primary header.
    # remove them if they are and record a WARNING

    try:
        del fo[0].header.ascard['BSCALE']
        warningsList.append("WARNING: BSCALE keyword found in primary header: "+fname)
        print "WARNING: BSCALE keyword found in primary header: "+fname
    except:
        pass
    try:
        del fo[0].header.ascard['BZERO']
        warningsList.append("WARNING: BZERO keyword found in primary header: "+fname)
        print "WARNING: BZERO keyword found in primary header: "+fname
    except:
        pass    

    convertList = []
    for i in range(1,len(fo)):
        print "breaking out extension ",i," of file ",fname
        xname = fo[i].header['EXTNAME'].upper()
        try:
            xver = str(fo[i].header['EXTVER'])
        except:
            xver = ''
        if xver:
            newfile = base+"_"+xname+"_"+xver + ext
        else:
            newfile = base+"_"+xname + ext

        nfits = os.path.join(dir,newfile)
        convertList.append(nfits)

        #f1 = pyfits.open(nfits,"w") # old fits syntax

        f1=pyfits.HDUList()
        f1.append(pyfits.PrimaryHDU())

        # Adjust header keywords for simple fits format.
        # Delete some keywords which must go away in simple fits format
        # First, get a copy of the Primary header into the new file

        f1[0].header.ascard = fo[0].header.ascard.copy()

        # These are primary header keywords that need to disappear.

        try:    del f1[0].header.ascard['EXTEND']
        except: pass
        try:    del f1[0].header.ascard['NEXTEND']
        except: pass

        # delete some of the extension-only keywords

        try:    del fo[i].header.ascard['XTENSION']
        except: pass
        try:    del fo[i].header.ascard['INHERIT']
        except: pass
        try:    del fo[i].header.ascard['EXTNAME']
        except: pass
        try:    del fo[i].header.ascard['EXTVER']
        except: pass
        try:    del fo[i].header.ascard['PCOUNT']
        except: pass
        try:    del fo[i].header.ascard['GCOUNT']
        except: pass


        # get the index of first HISTORY comment in the header we're updating
        # remember to increment this in the loop as we insert keys
        try:   ihist = f1[0].header.ascard.index_of('HISTORY')
        except KeyError:
            print 'Primary header has no HISTORY keyword.'
            print 'Appending extension keywords...'
            ihist = len(f1[0].header.ascard)

        # Now go through all keywords of extension and pass them
        # into the new simple fits header.
        # This is how I'd like to do it, but it doesn't work this way
        # for key,item in fo[i].header.ascard.keys(),fo[i].header.ascard:
        # keylist = fo[i].header.ascard.keys()

        # Below, an ascard object has attributes key,value,comment which allows 
        # direct access to these values.

    for item in fo[i].header.ascard:
        key = item.key
        val = item.value
        if not key and not val:
            continue
        elif not key and val:
            f1[0].header.ascard.insert(ihist,item)  
            ihist += 1
            continue

        #print "setting keyword,",key

        # if the 'key' is already in the header because it was copied
        # from the primary, we just want to update it with the value
        # of the keyword in the extension; no need to update ihist

        if (f1[0].header.has_key(key) and key != 'COMMENT' and key != 'HISTORY'):
            try:
                f1[0].header.update(key,fo[i].header[key])
            except pyfits.FITS_SevereError,err:
                warningsList.append("WARNING: FITS Error encountered in header,"+str(err))
                print "FITS Error encountered in header:",err
            continue

        # if it's NAXIS1, insert after NAXIS

        if key == 'NAXIS1':
            try:
                inaxis = f1[0].header.ascard.index_of('NAXIS')
            except KeyError,err:
                raise pyfits.FITS_FatalError,"Cannot find NAXIS keyword!"

            # NOTE: dereferencing 'key' to get item at this point 
            # was no good.  For repeated 'keys', as for example with
            # 'COMMENT' lines, this just gave back the first one.
            # item = fo[i].header.ascard[key]

            f1[0].header.ascard.insert(inaxis+1,item)
            ihist += 1

        # if it's NAXIS2, insert after NAXIS1 if it's already there,
        # otherwise, insert after NAXIS

        elif key == 'NAXIS2':
            try:
                inaxis2 = f1[0].header.ascard.index_of('NAXIS1') + 1
            except:
                try:
                    inaxis2 = f1[0].header.ascard.index_of('NAXIS') + 1
                except KeyError,err:
                    raise pyfits.FITS_FatalError,"Cannot find NAXIS keyword!"

                f1[0].header.ascard.insert(inaxis2,item)
                ihist += 1

        # if key is not NAXIS1 or NAXIS2, just insert before HISTORY
        else:
            f1[0].header.ascard.insert(ihist,item)
            ihist += 1

    # Fix the filename keyword
    f1[0].header.update('FILENAME',newfile)
    # And finally, get the data.
    f1[0].data = fo[i].data
    f1.writeto(nfits)
    #f1.close() # not needed anymore
    del f1

    fo.close()
    return convertList,warningsList

def fixHeader(fpath,list):
    """Adjust the keyword value pairs of a header of an image.
    This function does this 'in place', by just overwriting the
    old file with the new one. list is a list of tuples which are
    (keyword,value) pairs.
    """
    curdir = os.getcwd()
    path,file = os.path.split(fpath)
    if path:
        os.chdir(path)

    oldfits = pyfits.open(file,"update")
    if len(oldfits) > 1:
        raise IOError,"file is not a simple fits file."

    for keyword,value in list:      
        oldfits[0].header.update(keyword,value)

    oldfits[0].header.update('COMMENT','fUtil function call to adjust keyword values.')
    index = oldfits[0].header.ascard.index_of('COMMENT')

    for keyword,value in list:
        index += 1
        card = pyfits.Card('COMMENT',value=keyword+' '+str(value))
        oldfits[0].header.ascard.insert(index,card)
    oldfits.close()
    os.chdir(curdir)
    return

def fixDrzHeader(fpath,module,date,file_type="SCI",addCards=None):
    """Adjust the header of a drizzled image produced by pydrizzle/drizzle.  Mostly,
    this means removing a lot of keywords which have no meaning in the context of
    an image which is created from a number of different images from potentially any
    number of sources.  This list of keywords can be found in Bugzilla bug #1333.
    The addCards option allows the caller to add new ascards to the header.
    addCards is a list of tuples like
    [("newkeyword1","value1"),
     ("newkeyword2","value2"),
     ]
    """
    keepem = ["SIMPLE",
          "BITPIX",
          "NAXIS",
          "NAXIS1",
          "NAXIS2",
          #"NAXIS3",
          "TELESCOP",
          "INSTRUME",
          "DETECTOR",
          "FILTER", 
          "FILTER1",
          "FILTER2",
          "LRFWAVE",
          "BUNIT",
          "CRPIX1",
          "CRPIX2",
          "CRVAL1",
          "CRVAL2",
          "CTYPE1",
          "CTYPE2",
          "CD1_1",
          "CD1_2",
          "CD2_1",
          "CD2_2",
          "LTV1",
          "LTV2",
          "LTM1_1",
          "LTM2_2",
          "IDCTAB",
          "PHOTMODE",
          "PHOTFLAM",
          "PHOTZPT",
          "PHOTPLAM",
          "PHOTBW",
          "ORIGIN",
          "DATE",
          "FILENAME",
          "EXPTIME",
          "ALIGNSKY",
          "NCOMBINE",
          "NDRIZIM",
          "PA_V3",
          ]
    curdir = os.getcwd()
    path,file = os.path.split(fpath)
    if path:
        os.chdir(path)

    oldfits = pyfits.open(file)
    if len(oldfits) > 1:
        raise IOError,"file is not a simple fits file."

    newfits = pyfits.HDUList()
    newfits.append(pyfits.PrimaryHDU())

    # Tried using the pyfits "update" mode on the open call but this is
    # logically problematic when you start deleting ascards which do not
    # have a keyword (blank cards, section header cards, etc). And we noticed
    # that pyfits winds up doing exactly what we do below when the size of the
    # header changes dramatically, i.e. writes a temp file and does a copy.
    # Pyfits automatically attaches an EXTEND ascard onto the header of the
    # appended Primary HDU with the value T.  This seems like a bug since pyfits
    # has no idea whether the file to be made is an mef file or not.
    # this card is deleted and not required for a simple fits file.

    try:
        del newfits[0].header.ascard["EXTEND"]
    except:
        pass

    # the first regex pattern will match all drizzle keywords. according to the spec in 
    # Bugzilla bug # 1333, the drizzle specific keywords should only be kept in the SCI -type data 
    # (i.e. the _drz.fits image).  The second pattern will match all context specific keywords 
    # (CONnnnn) and should only be kept in the context image (i.e. _drz_context.fits) we need 
    # to keep the two patterns seperate and use passed file_type to figure out which to keep for the 
    # passed image file.
    # set the matching pattern based on file_type.

    if file_type == "SCI":
        eng_pars = ["CCDAMP",
                "CCDGAIN",
                "CCDOFSTA",
                "CCDOFSTB",
                "CCDOFSTC",
                "CCDOFSTD",
                "ATODGNA",
                "ATODGNB",
                "ATODGNC",
                "ATODGND",
                "READNSEA",
                "READNSEB",
                "READNSEC",
                "READNSED"
                ]
        for i in eng_pars:
            keepem.append(i)
        pattern = '^D[D0-9][0-9]{2}[a-zA-Z0-9]{2}.?.?'
    elif file_type == "CTX":
        pattern = '^CON[0-9][0-9][0-9][0-9]'
    else:
        pattern = None

    for ascard in oldfits[0].header.ascard:
        if pattern:
            result = re.match(pattern,ascard.key)
        else: 
            result = None
        if not ascard.key:
            continue
        elif ascard.key in keepem:
            newfits[0].header.update(ascard.key,ascard.value,ascard.comment)
            continue
        elif result:
            newfits[0].header.update(ascard.key,ascard.value,ascard.comment)
            continue

    newfits[0].header.update("FILENAME",file)
    newfits[0].header.update("FILETYPE",file_type,after="FILENAME")
    newfits[0].header.update("ORIGIN",module)
    newfits[0].header.update("DATE",date)
    
    if addCards:
        for card in addCards:
            if newfits[0].header.has_key(card[0]):
                keyind  = newfits[0].header.ascard.index_of(card[0])
                # prevkey = newfits[0].header.ascard[keyind-1].key
                try:
                    comment = newfits[0].header.ascard[keyind].comment
                except:
                    comment = ''
                del newfits[0].header[card[0]]

                # The format, gentlemen, please!
                try:
                    format = card[2]
                except:
                    if isinstance(card[1],types.FloatType):
                        format='%20.12G / %s'
                    else:
                        format=None

                #newfits[0].header.update(card[0],card[1],after=prevkey)
                newCard = pyfits.Card(card[0],card[1],comment,format=format)
                
                newCard = pyfits.Card(card[0],card[1],comment)
                newfits[0].header.ascard.insert(keyind,newCard)
                del keyind,comment,newCard
            else:
                newfits[0].header.update(card[0],card[1])
            
    newfits[0].data = oldfits[0].data
    newfits.writeto("temp.fits")
    newfits.close()
    oldfits.close()
    del newfits,oldfits
    os.remove(file)
    os.rename('temp.fits', file)
    os.chdir(curdir)
    return

def fixDrzHeader2(fpath,module,date,file_type="SCI",addCards=None):
    """Adjust the header of a drizzled image produced by pydrizzle/drizzle.  Mostly,
    this means removing a lot of keywords which have no meaning in the context of
    an image which is created from a number of different images from potentially any
    number of sources.  This list of keywords can be found in Bugzilla bug #1333.
    The addCards option allows the caller to add new ascards to the header.
    addCards is a list of tuples like
    [("newkeyword1","value1"),
     ("newkeyword2","value2"),
     ]
    """
    keepem = ["SIMPLE",
          "BITPIX",
          "NAXIS",
          "NAXIS1",
          "NAXIS2",
          #"NAXIS3",
          "TELESCOP",
          "INSTRUME",
          "DETECTOR",
          "FILTER",      #XX by ____xingxing
          "FILTER1",
          "FILTER2",
          "LRFWAVE",
          "BUNIT",
          "CRPIX1",
          "CRPIX2",
          "CRVAL1",
          "CRVAL2",
          "CTYPE1",
          "CTYPE2",
          "CD1_1",
          "CD1_2",
          "CD2_1",
          "CD2_2",
          "LTV1",
          "LTV2",
          "LTM1_1",
          "LTM2_2",
          "IDCTAB",
          "PHOTMODE",
          "PHOTFLAM",
          "PHOTZPT",
          "PHOTPLAM",
          "PHOTBW",
          "ORIGIN",
          "DATE",
          "FILENAME",
          "EXPTIME",
          "ALIGNSKY",
          "NCOMBINE",
          "NDRIZIM",
          "PA_V3",
          ]
    curdir = os.getcwd()
    path,file = os.path.split(fpath)
    if path:
        os.chdir(path)

    oldfits = pyfits.open(file)
    if len(oldfits) > 1:
        raise IOError,"file is not a simple fits file."

    newfits = pyfits.HDUList()
    newfits.append(pyfits.PrimaryHDU())

    # Tried using the pyfits "update" mode on the open call but this is
    # logically problematic when you start deleting ascards which do not
    # have a keyword (blank cards, section header cards, etc). And we noticed
    # that pyfits winds up doing exactly what we do below when the size of the
    # header changes dramatically, i.e. writes a temp file and does a copy.
    # Pyfits automatically attaches an EXTEND ascard onto the header of the
    # appended Primary HDU with the value T.  This seems like a bug since pyfits
    # has no idea whether the file to be made is an mef file or not.
    # this card is deleted and not required for a simple fits file.

    try:
        del newfits[0].header.ascard["EXTEND"]
    except:
        pass

    # the first regex pattern will match all drizzle keywords. according to the spec in 
    # Bugzilla bug # 1333, the drizzle specific keywords should only be kept in the SCI -type data 
    # (i.e. the _drz.fits image).  The second pattern will match all context specific keywords 
    # (CONnnnn) and should only be kept in the context image (i.e. _drz_context.fits) we need 
    # to keep the two patterns seperate and use passed file_type to figure out which to keep for the 
    # passed image file.
    # set the matching pattern based on file_type.

    if file_type == "SCI":
        eng_pars = ["CCDAMP",
                "CCDGAIN",
                "CCDOFSTA",
                "CCDOFSTB",
                "CCDOFSTC",
                "CCDOFSTD",
                "ATODGNA",
                "ATODGNB",
                "ATODGNC",
                "ATODGND",
                "READNSEA",
                "READNSEB",
                "READNSEC",
                "READNSED"
                ]
        for i in eng_pars:
            keepem.append(i)
        pattern = '^D[D0-9][0-9]{2}[a-zA-Z0-9]{2}.?.?'
    elif file_type == "CTX":
        pattern = '^CON[0-9][0-9][0-9][0-9]'
    else:
        pattern = None

    for ascard in oldfits[0].header.ascard:
        if pattern:
            result = re.match(pattern,ascard.key)
        else: 
            result = None
        if not ascard.key:
            continue
        elif ascard.key in keepem:
            newfits[0].header.update(ascard.key,ascard.value,ascard.comment)
            continue
        elif result:
            newfits[0].header.update(ascard.key,ascard.value,ascard.comment)
            continue

    newfits[0].header.update("FILENAME",file)
    newfits[0].header.update("FILETYPE",file_type,after="FILENAME")
    newfits[0].header.update("ORIGIN",module)
    newfits[0].header.update("DATE",date)
    
    if addCards:
        for card in addCards:
            if newfits[0].header.has_key(card[0]):
                keyind  = newfits[0].header.ascard.index_of(card[0])
                try:
                    comment = newfits[0].header.ascard[keyind].comment
                except:
                    comment = ''
                del newfits[0].header[card[0]]
                newCard = pyfits.Card(card[0],card[1],comment)
                newfits[0].header.ascard.insert(keyind,newCard)
                del keyind,comment,newCard
            else:
                newfits[0].header.update(card[0],card[1])
    newfits[0].data = oldfits[0].data
    newfits.writeto("temp.fits")
    newfits.close()
    oldfits.close()
    del newfits,oldfits
    os.remove(file)
    os.rename('temp.fits', file)
    os.chdir(curdir)
    return

def printTab(fitsfile):
    """print out the data table of an association table (or any fits really).
    where the association table data is in extension 1.
    """
    fo = pyfits.open(fitsfile)
    if fo[0].header['FILETYPE'] != 'ASN_TABLE':
        raise TypeError,"fits file is not an association table."

    print "\nAsn table data for ",fitsfile,":"
    for i in fo[1].data:
        print i
    return


def createfits(newfile,NX,NY,bitpix):
    """
    Create a new simple FITS disk file of specified size and type;
    return the file ojbect in update mode.  The header is minimal,
    and the data array is initialized to all zeros.
    """
    if bitpix == 8:
        type = numpy.uint8
    elif bitpix == 16:
        type = numpy.int16
    elif bitpix == -16:
        type = numpy.uint16
    elif bitpix == 32:
        type = numpy.int32
    elif bitpix == -32:
        type = numpy.float32
    elif bitpix == -64:
        type = numpy.float64
    else:
        raise pyfits.FITS_SevereError, "Invalid BITPIX specification."

    # create object template in memory
    fobj = pyfits.HDUList()
    fobj.append(pyfits.PrimaryHDU())

    # make data array -- NOTE! the numpy/pyfits indexing is backwards!!
    fobj[0].data = numpy.zeros((NY,NX), dtype=type)

    # fix up the header -- 
    # actually, this no longer seems necessary with pyfits, but best safe
    fobj[0].header['BITPIX'] = bitpix
    fobj[0].header['NAXIS'] = 2  
    fobj[0].header.update('NAXIS1',NX,after="NAXIS")
    fobj[0].header.update('NAXIS2',NY,after="NAXIS1")
    if fobj[0].header.has_key('EXTEND'):
       del fobj[0].header.ascard['EXTEND']

    # write the file to disk
    if os.path.isfile(newfile):
        os.remove(newfile)
    fobj.writeto(newfile)
    del fobj

    # reopen for updating, and return
    ff = pyfits.open(newfile,'update')
    if ff[0].header['NAXIS1'] != NX or ff[0].header['NAXIS2'] != NY:
        raise pyfits.FITS_SevereError, "blech! nasty naxis bug!"
    return ff

def getTableZeroPoint(fitsfile):
    """The ACS/WFC Zeropoints table
    (http://www.stsci.edu/hst/acs/analysis/zeropoints/#tablestart)
    in dictionary form.  We use this to determine the filter dependant
    zeropoint for normalised data.

    Updated to be consistent with the PHOTLAM
    keys in the GOODS catalogs. F. Menanteau Mar 2005"""

    ABMagZeroPoints = { 
        "F225W":   24.06, #UVIS
        "F275W":   24.14, #UVIS
        "F350LP":  26.95, #UVIS
        "F435W":   25.6528602704, # was 25.673,
        "F475W":   26.068,
        "F502N":   22.271,
        "F550M":   24.892,
        "F555W":   25.718,
        "F606W":   26.4933852945, # was 26.486,
        "F625W":   25.898,
        "F658N":   22.747,
        "F660N":   21.676,
        "F775W":   25.6405052886, # was 25.654,
        "F814W":   25.937,
        "F850LP":  24.8431369524, # was 24.862,
        "F892N":   22.371,
        "F105W":   26.27, #IR
        "F125W":   26.24, #IR
        "F140W":   26.46, #IR
        "F160W":   25.95  #IR

        }

    filter = filterName(fitsfile)
    zpt = ABMagZeroPoints[filter]
    return zpt
