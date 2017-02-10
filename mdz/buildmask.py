"""
Functions to build mask files for PyDrizzle.
    - buildMaskImage(rootname,bitvalue,extname='DQ',extver=None):
        This function takes the DQ array(or array from extension given)
        and builds a bit mask as an input weighting mask for 'drizzle'.
        
    - buildShadowMaskImage(rootname,detnum,replace=no):
        This function builds a weighting image for use with drizzle from 
        the WFPC2 shadow mask functions derived from 'wmosaic'.
"""

#
# Revision History:
#
#   14-Nov-2001 - Revised to produce proper simple FITS image starting
#                   with 'SIMPLE=T' keyword instead of 'XTENSION'. WJH
#
#   25-Jan-2002 - Converted to work with PyFITS V0.6.1 and numpy. WJH
#
#   12-June-2002 - Added check to remove old versions of or existing mask 
#                  files before building a new one. 
#
import string,os,types
from pyraf import iraf 

import fileutil
import drutil

import pyfits
from numpy import *

yes = iraf.yes
no = iraf.no

def buildMask(dqarr,bitvalue):
    """ Builds a bit-mask from an input DQ array and a bitvalue flag"""
    _maskarr = bitwise_or(dqarr,array([bitvalue]))
    return choose(greater(_maskarr,bitvalue),(1,0)).astype(int16)

def buildMaskImage(rootname,bitvalue,extname='DQ',extver=1):
    """ Builds mask image from rootname's DQ array 
        If there is no valid 'DQ' array in image, then return
        an empty string.
    """
    # If no bitvalue is set, assume no mask is desired
    if bitvalue == None:        
        return None

    # build output name
    _indx = string.find(rootname,'.fits')
    if _indx > 0:
        maskname = rootname[:_indx]+'_inmask'+repr(extver)+'.fits'
        whtname = rootname[:_indx]+'_wtmask'+repr(extver)+'.fits' #WZ
    else:
        maskname = rootname+'_inmask'+repr(extver)+'.fits'
        whtname = rootname+'_wtmask'+repr(extver)+'.fits'    #WZ
    # If an old version of the maskfile was present, remove it and rebuild it.
    if fileutil.findFile(maskname) == yes:
        os.remove(maskname)
    if fileutil.findFile(whtname) == yes: #WZ
        os.remove(whtname)
        
    # Open input file with DQ array
    fdq = pyfits.open(rootname)
    try:
        # Read in DQ array
        dqarr = fdq[extname,extver].data

        # Build mask array from DQ array
        maskarr = buildMask(dqarr,bitvalue)

        #Write out the mask file as simple FITS file
        fmask = pyfits.open(maskname,'append')
        maskhdu = pyfits.PrimaryHDU(data=maskarr)
        fmask.append(maskhdu)

        #Close files
        fmask.close()
        del fmask
    except:
        return None

    fdq.close()
    del fdq
            
    # Return the name of the mask image written out
    return maskname

# Utility functions for generating X and Y arrays for
# creating the WFPC2 shadow masks.
# Coefficients are from WMOSAIC function 'get_edge2.x'
"""
_coeffs = {'1':[[52.20921,0.009072887,0.009072887],[42.62779,0.009122855,-1.106709E-5]],
          '2':[[21.77283,0.01842164,-1.398300E-5],[47.68184,0.00265608,-1.468158E-5]],
          '3':[[44.18944,0.0138938,-1.412296E-5],[30.23626,0.008156041,-1.491324E-5]],
          '4':[[44.56632,0.003509023,-1.278723E-5],[40.91462,0.01273679,-1.063462E-5]]}
"""
def _func_Shadow_WF1y(x,y):
    return y + 0.5 - (42.62779 + 0.009122855*y - 1.106709E-5*y**2)
def _func_Shadow_WF1x(x,y):
    return x + 0.5 - (52.20921 + 0.009072887*x - 9.941337e-6*x**2)

def _func_Shadow_WF2y(x,y):
    return y + 0.5 - (47.68184 + 0.00265608*y - 1.468158E-5*y**2)
def _func_Shadow_WF2x(x,y):
    return x + 0.5 - (21.77283 + 0.01842164*x - 1.398300E-5*x**2)

def _func_Shadow_WF3y(x,y):
    return y + 0.5 - (30.23626 + 0.008156041*y - 1.491324E-5*y**2 )
def _func_Shadow_WF3x(x,y):
    return x + 0.5 - (44.18944 + 0.0138938*x - 1.412296E-5*x**2)

def _func_Shadow_WF4y(x,y):
    return y + 0.5 - (40.91462 + 0.01273679*y - 1.063462E-5*y**2)
def _func_Shadow_WF4x(x,y):
    return x + 0.5 - (44.56632 + 0.003509023*x - 1.278723E-5*x**2)

# Function for creating the weighting image for WFPC2 drizzling
# from the 'wmosaic' shadow mask polynomials.     
def buildShadowMaskImage(rootname,detnum,replace=no):
    """ Builds mask image from WFPC2 shadow calibrations.
      detnum - string value for 'DETECTOR' detector  
    """  
    # insure detnum is a string
    if not isinstance(detnum,types.StringType):
        detnum = repr(detnum)
    
    _funcroot = '_func_Shadow_WF'  
    # build template shadow mask's filename
    _mask = 'wfpc2_inmask'+detnum+'.fits'
    
    if rootname != None:
        # build output name by stripping off old extension from rootname...
        _name = drutil.buildNewRootname(rootname)
        #print 'found rootname of ',_name
        # ... now add our own syntax and extension to the output name.
        _indx = string.find(_name,'.')
        if _indx > 0:
            maskname = _name[:_indx]+'_inmask'+detnum+'.fits'
        else:
            maskname = _name+'_inmask'+detnum+'.fits'
    else:
        maskname = None
    # Check to see if file exists...
    if fileutil.findFile(_mask) == no or replace:  
    # If not, create the file.
    # This takes a long time to run, so it should be done
    # only when absolutely necessary...
        try:        
            _funcx = _funcroot+detnum+'x'
            _funcy = _funcroot+detnum+'y'
            
            _xarr = clip(fromfunction(eval(_funcx),(800,800)),0.0,1.0).astype(int16)
            _yarr = clip(fromfunction(eval(_funcy),(800,800)),0.0,1.0).astype(int16)
            maskarr = _xarr * _yarr
            
            #Write out the mask file as simple FITS file
            fmask = pyfits.open(_mask,'append')
            maskhdu = pyfits.PrimaryHDU(data=maskarr)
            fmask.append(maskhdu)
            
            #Close files
            fmask.close()
            del fmask
        except:
            return None

    # Now, copy template mask file to output file, if necessary
    if maskname != None:
        fileutil.copyFile(_mask,maskname,replace=no)
    # Return the name of the mask image written out
    return maskname
    
