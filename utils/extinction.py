#!/usr/bin/env python

# $Id: extinction.py,v 1.2 2003/02/27 23:59:20 anderson Exp $ 
# ---------------------------------------------------------------------
# This module provides functionality to calculate the extinction correction
# for a given position on the sky.  This is done using the 
# Schlegel, Finkbeiner & Davis (SFD) dust maps (ApJ, 1998, 500, 525)
# and the interface to these maps provided by SFD.
#
# RA & DEC coordinates are passed by the caller and converted herein
# to the required galactic coords, l, b. via the iraf astutil task
# galactic.  The dust_getval code is then run twice, once to get the
# extinction correction and another to get the data quality flags for
# the sky position.  See the README.C in SFD for a description of the
# mask values returned.


__version__      = '$Revision: 1.2 $ '[11:-3]
__version_date__ = '$Date: 2003/02/27 23:59:20 $ '[7:-3]
__author__       = "Ken Anderson, anderson@pha.jhu.edu"

import popen2,fUtil
from pyraf import iraf
from iraf  import astutil
import pdb

def getEBV(ra,dec):
    """function recieves a coordinate pair, RA and DEC from the caller, 
    converts to galactic coords and runs the dust_getval code, installed under the pipeline 
    environment.  Returns an extinction correction in magnitudes and an error object 
    (a list of strings) of possible reported problems with the region of the sky.
    """
    # convert ra and dec to l,b using astutil.  ra should be in hours (freaking iraf).
    # we emulate pipes to Stdin and Stdout so we don't write no stinking files
    raanddec = [str(ra/15) +" "+str(dec)+" 2000"]
    conversion = astutil.galactic(Stdin=raanddec,print_c="no",Stdout=1)
    # conversion is a list of strings, this has only one element:
    # eg. ['     227.5430   46.1912']
    # which is l and b
    gall = conversion[0].split()[0]
    galb = conversion[0].split()[1]
    
    #pdb.set_trace()
    # ok, we have l and b. now onto the extinction stuff. build the dust_val command line
    cmd    = "dust_getval "+ gall+ " "+ galb+" interp=y verbose=n"
    output = _calc_ebv(cmd)
    # output is a list of strings, only one element in this case. looks like
    # [' 227.543  46.191      0.03452\n']
    # dust_getval returns the original coords and the extinction correction in mags
    # which is the last piece of that string
    eBV = output[0].split()[2]
    # next run dust_getval with the mask option to look for anomolies in the maps
    cmd  = "dust_getval "+ gall+ " "+ galb+" map=mask verbose=n"
    mask = _calc_ebv(cmd)
    # quality is a string of data quality flags returned by dust_getval when map=mask.
    # looks like
    # ' 227.543  46.191  3hcons OK      OK      OK      OK      OK      OK     \n'
    quality = mask[1]
    # return this with the extinction.
    return eBV,quality


def filterFactor(filter):
    """caller passes an ACS filter of the form "DET_FILTER" and function returns the
    extinction correction factor, a float, for that filter.  Now this function defines a dictionary
    of extinction correction factors directly, but this also exists as a file in 
    $PIPELINE/maps.  It is anticipated that these values will not change often, if at all,
    hence, the dictionary is defined here rather than created on the fly from the file, but
    that could be changed if it is anticipated that these numbers might change a lot.

    eg.
    >>>filterFactor("HRC_F555W")
    '3.24695724147'
    """
    ffactors = {
     # Updated, WZ Oct 2010
    "HST_WFC3_UVIS_F225W": 7.474, # 8.8,
    "HST_WFC3_UVIS_F275W": 6.14,
    "HST_WFC3_UVIS_F336W": 5.09,
    "HST_WFC3_UVIS_F390W": 4.514,
    "HST_WFC3_UVIS_F350LP": 2.335,
    "HST_WFC3_UVIS_F606W": 2.928500386,
    "HST_WFC3_UVIS_F814W": 1.84740717341,
    "HST_WFC3_IR_F098M": 1.3, #WZ estimate 1.015,
    "HST_WFC3_IR_F105W": 1.015,
    "HST_WFC3_IR_F110W": 0.876,
    "HST_WFC3_IR_F125W": 0.757,
    "HST_WFC3_IR_F140W": 0.609,
    "HST_WFC3_IR_F160W": 0.470,
    "HST_ACS_WFC_F435W": 4.11697363315,
    "HST_ACS_WFC_F555W": 3.24230342654,
    "HST_ACS_WFC_F606W": 2.928500386,
    "HST_ACS_WFC_F814W": 1.84740717341,
    "HST_ACS_WFC_F475W": 3.74714182372,
    "HST_ACS_WFC_F625W": 2.67121669327,
    "HST_ACS_WFC_F775W": 2.01774028108,
    "HST_ACS_WFC_F850LP": 1.47335876958,
     # things above are added by xingxing from aplus  
    "HST_ACS_WFC_F502N": 3.52366637215,
    "HST_ACS_WFC_F892N": 1.51713294198,
    "HST_ACS_WFC_F658N": 2.52193964747,
    "HST_ACS_WFC_F550M": 3.05672088958,
    "HST_ACS_HRC_F435W": 4.11227228078,
    "HST_ACS_HRC_F555W": 3.24695724147,
    "HST_ACS_HRC_F606W": 2.94741773243,
    "HST_ACS_HRC_F814W": 1.82249557542,
    "HST_ACS_HRC_F475W": 3.724544959,
    "HST_ACS_HRC_F625W": 2.67859346295,
    "HST_ACS_HRC_F775W": 2.02818977096,
    "HST_ACS_HRC_F850LP": 1.44407689634,
    "HST_ACS_HRC_F344N":  5.10086305785,
    "HST_ACS_HRC_F502N":  3.52410034736,
    "HST_ACS_HRC_F892N":  1.5170227107,
    "HST_ACS_HRC_F658N":  2.52245685895,
    "HST_ACS_HRC_F220W":  8.81083859281,
    "HST_ACS_HRC_F250W":  6.52297815722,
    "HST_ACS_HRC_F330W":  5.17376227866,
    "HST_ACS_HRC_F550M":  3.05823161415,
    }
    try:
        factor = ffactors[filter]
    except KeyError,err:
        raise fUtil.filterError,"No correction factor for filter "+filter

    return factor

#############################################################################################

def _calc_ebv(cmd):
    sproc  = popen2.Popen3(cmd,1)
    output = sproc.fromchild.readlines()
    errs   = sproc.childerr.readlines()
    return output
