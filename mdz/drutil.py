"""
Utility functions for PyDrizzle that rely on PyRAF's interface to IRAF 
tasks.  
"""
#
# Revision History:
#   Nov 2001: Added function (getLTVOffsets) for reading subarray offsets 
#               from extension header's LTV keywords.  WJH
#   Mar 2002: Restructured to only contain IRAF based functions.
#  

import string,os, pdb
from pyraf import iraf

import numpy as N

import fileutil

# Convenience definitions
DEGTORAD = fileutil.DEGTORAD

no = iraf.no
yes = iraf.yes

# Constants
IDCTAB	= 1
DRIZZLE = 2
TRAUGER = 3


# List of supported default file types
# It will look for these file types by default
# when trying to recognize input rootnames.
EXTLIST =  ['_crj.fits','_flt.fits','_sfl.fits','_raw.fits','.c0h','.hhh','.fits']

#################
#
#
#		Utility Functions based on IRAF
#
#
#################		
	
def findNumExt(filename):
	# Open the file given by 'rootname' and return the 
	# number of extensions written out for the observation.
	# It will be up to the Exposure classes to determine how
	# many IMSETS they have to work with.
	#
	# Only look in the primary extension or first group.
	iraf.keypar(filename,'NEXTEND',silent='yes')
	
	# This may need to be changed to support simple image without
	# extensions (such as simple FITS images).
	if iraf.keypar.value == 0 or iraf.keypar.value == '':
		raise ValueError,"There are NO extensions to be read in this image!"
	
	return iraf.keypar.value


# General method for returning what type of detector was
# used for the observation.  This would rely on the class
# variable 'PyDrizzle.instrument' to act as the catalog of
# supported instruments. 
def getPrimaryKeyword(filename,keyword):

    iraf.keypar(filename,keyword)
    value =  iraf.keypar.value
    if value == '':
        value = None
    
    return value

def updateKeyword(filename,keyword,value):
    """ Silently update/add keyword to header with given value.
        Used for adding DEXPTIME to extension headers.
    """
    # Remember how hedit was set to start with...
    _add = iraf.hedit.add
    _update = iraf.hedit.update
    _verify = iraf.hedit.verify
    _show = iraf.hedit.show
    _delete = iraf.hedit.delete
    
    iraf.hedit.add=yes
    iraf.hedit.update=yes
    iraf.hedit.verify=no
    iraf.hedit.show=no
    iraf.hedit.delete=no
    iraf.hedit(filename,keyword,value)
    # Reset optional parameters to previous setting
    iraf.hedit.add=_add
    iraf.hedit.update=_update
    iraf.hedit.verify=_verify
    iraf.hedit.show=_show
    iraf.hedit.delete=_delete   

# Utility function to extract subarray offsets
# from LTV keywords.  Could be expanded to use
# different keywords for some cases.
# Added: 7-Nov-2001 WJH
def getLTVOffsets(rootname):

    _ltv1 = getPrimaryKeyword(rootname,'LTV1')
    if _ltv1 == None: 
        _ltv1 = 0.
        print 'No LTV keywords found... Defaulting to 0. for offsets.'
        
    _ltv2 = getPrimaryKeyword(rootname,'LTV2')
    if _ltv2 == None: _ltv2 = 0.

    return _ltv1,_ltv2
       
def getChipId(rootname):
	
	chip = 1
	_found = yes

	try:
		iraf.keypar(rootname,"CCDCHIP",silent='yes')
		chip = int(iraf.keypar.value)	
	except ValueError:
		chip = 1
		_found = no
	
	if _found == no:
		try:
			iraf.keypar(rootname,"DETECTOR",silent='yes')
			chip = int(iraf.keypar.value)
		except ValueError:
			chip = 1

	return chip

 
def getIDCFile(filename,keyword=None,dir=None):
    # Open the primary header of the file and read the name of 
    # the IDCTAB.
    # Parameters:
    #   filename -  filename of primary extension to read IDCTAB info
    #
    #   keyword(optional) - header keyword with name of IDCTAB 
    #            --OR-- 'HEADER' if need to build IDCTAB name from scratch
    #           (default value: 'IDCTAB')
    #   dir(optional) - directory with default drizzle coeffs tables
    #                   (default value: 'drizzle$coeffs'
    #
    # This function needs to be generalized to support
    # keyword='HEADER', as would be the case for WFPC2 data.
    #

    if string.lower(keyword) == 'idctab':
        # keyword specifies header keyword with IDCTAB name
        try:
	        iraf.keypar(filename,keyword,silent='yes')
	        idcfile = iraf.keypar.value        			
        except:
	        print 'Warning: No IDCTAB specified in header!'
	        idcfile = None       
    else:
        # Need to build IDCTAB filename from scratch
        if dir == None:
            default_dir = 'drizzle$coeffs/'
        else:
            default_dir = dir

        instrument = getPrimaryKeyword(filename,'INSTRUME')
        detector = getPrimaryKeyword(filename,'DETECTOR')
        
        if instrument == 'WFPC2':
            if detector == '1':
                detname = 'pc'
            else:
                detname = 'wf'
            idcfile = default_dir+detname+detector+'-'+string.lower(keyword)

        elif instrument == 'STIS':
            idcfile = default_dir+'stis-'+string.lower(detector)

        elif instrument == 'NICMOS':
            camera = getPrimaryKeyword(filename,'CAMERA')
            if camera != None:
                idcfile = default_dir+'nic-'+repr(camera)
            else:
                idcfile = None
        else:
            idcfile = None           
    
    # Account for possible absence of IDCTAB name in header
    if idcfile == 'N/A': 
        idcfile = None

    # Now we need to recursively expand any IRAF symbols to full paths...
    indx = 1
    while indx > 0:
        indx = string.find(idcfile,'$')
        if indx > 0:
            fullpath = iraf.envget(idcfile[:indx])
            idcfile = fullpath+idcfile[indx+1:]
    
    return idcfile

#
# Function to read in coefficients from ASCII Cubic Drizzle
# Coefficients table.
# 
    
def readCubicTable(idcfile):
    # Assumption: this will only be used for cubic file...
    order = 3
    # Also, this function does NOT perform any scaling on
    # the coefficients, it simply passes along what is found
    # in the file as is...

    print "Read x0, y0"
    # pdb.set_trace()
    ifile = open(idcfile,'r')
    # Search for the first line of the coefficients
    _line = fileutil.rAsciiLine(ifile)
    while string.lower(string.strip(_line)) != 'cubic':
        _line = fileutil.rAsciiLine(ifile)        
    # Read in each row of coefficients, without line breaks or newlines
    # split them into their values, and create a list for A coefficients
    # and another list for the B coefficients
    _line = fileutil.rAsciiLine(ifile)
    a_coeffs = string.split(_line)
    x0 = N.float(a_coeffs[0]) 
    _line = fileutil.rAsciiLine(ifile)
    a_coeffs[len(a_coeffs):] = string.split(_line)
    # Scale coefficients for use within PyDrizzle
    for i in range(len(a_coeffs)):
        a_coeffs[i] = N.float(a_coeffs[i])
    
    _line = fileutil.rAsciiLine(ifile)
    b_coeffs = string.split(_line)
    y0 = N.float(b_coeffs[0])
    _line = fileutil.rAsciiLine(ifile)
    b_coeffs[len(b_coeffs):] = string.split(_line)
    # Scale coefficients for use within PyDrizzle
    for i in range(len(b_coeffs)):
       b_coeffs[i] = N.float(b_coeffs[i]) 

    ifile.close()
    del ifile
    # Now, convert the coefficients into a Numeric array
    # with the right coefficients in the right place.
    # Populate output values now...
    fx = N.zeros(shape=(order+1,order+1),dtype=N.float32)
    fy = N.zeros(shape=(order+1,order+1),dtype=N.float32)
    # Assign the coefficients to their array positions
    fx[0,0] = 0.
    fx[1] = N.array([a_coeffs[2],a_coeffs[1],0.,0.],dtype=N.float32)    
    fx[2] = N.array([a_coeffs[5],a_coeffs[4],a_coeffs[3],0.],dtype=N.float32)
    fx[3] = N.array([a_coeffs[9],a_coeffs[8],a_coeffs[7],a_coeffs[6]],dtype=N.float32)
    fy[0,0] = 0.
    fy[1] = N.array([b_coeffs[2],b_coeffs[1],0.,0.],dtype=N.float32)    
    fy[2] = N.array([b_coeffs[5],b_coeffs[4],b_coeffs[3],0.],dtype=N.float32)
    fy[3] = N.array([b_coeffs[9],b_coeffs[8],b_coeffs[7],b_coeffs[6]],dtype=N.float32)
    
    # Used in Pattern.computeOffsets()
    refpix = {}
    refpix['XREF'] = None 
    refpix['YREF'] = None
    refpix['V2REF'] = x0 
    refpix['V3REF'] = y0 
    refpix['XDELTA'] = 0.
    refpix['YDELTA'] = 0.
    refpix['PSCALE'] = None
    refpix['centered'] = yes
    
    return fx,fy,refpix,order

# Function to compute the index of refraction for MgF2 at
# the specified wavelength for use with Trauger coefficients
def _MgF2(lam):
    _sig = pow((1.0e7/lam),2)
    return N.sqrt(1.0 + 2.590355e10/(5.312993e10-_sig) +
        4.4543708e9/(11.17083e9-_sig) + 4.0838897e5/(1.766361e5-_sig))
        
# Function to read Trauger ASCII file and return cubic coefficients    
def readTraugerTable(idcfile,wavelength):
    # Trauger coefficients only result in a cubic file...
    order = 3
    numco = 10
    a_coeffs = [0] * numco
    b_coeffs = [0] * numco
    indx = _MgF2(wavelength)
    
    ifile = open(idcfile,'r')
    # Search for the first line of the coefficients
    _line = fileutil.rAsciiLine(ifile)
    while string.lower(string.strip(_line)) != 'trauger':
        _line = fileutil.rAsciiLine(ifile)        
    # Read in each row of coefficients,split them into their values, 
    # and convert them into cubic coefficients based on 
    # index of refraction value for the given wavelength
    # Build X coefficients from first 10 rows of Trauger coefficients
    j = 0
    while j < 20:
        _line = fileutil.rAsciiLine(ifile)
        _lc = string.split(_line)
        if j < 10:
            a_coeffs[j] = float(_lc[0])+float(_lc[1])*(indx-1.5)+float(_lc[2])*(indx-1.5)**2
        else:
            b_coeffs[j-10] = float(_lc[0])+float(_lc[1])*(indx-1.5)+float(_lc[2])*(indx-1.5)**2
        j = j + 1
                
    ifile.close()
    del ifile
    
    # Now, convert the coefficients into a Numeric array
    # with the right coefficients in the right place.
    # Populate output values now...
    fx = N.zeros(shape=(order+1,order+1),dtype=N.float32)
    fy = N.zeros(shape=(order+1,order+1),dtype=N.float32)
    # Assign the coefficients to their array positions
    fx[0,0] = 0.
    fx[1] = N.array([a_coeffs[2],a_coeffs[1],0.,0.])    
    fx[2] = N.array([a_coeffs[5],a_coeffs[4],a_coeffs[3],0.])
    fx[3] = N.array([a_coeffs[9],a_coeffs[8],a_coeffs[7],a_coeffs[6]])
    fy[0,0] = 0.
    fy[1] = N.array([b_coeffs[2],b_coeffs[1],0.,0.])    
    fy[2] = N.array([b_coeffs[5],b_coeffs[4],b_coeffs[3],0.])
    fy[3] = N.array([b_coeffs[9],b_coeffs[8],b_coeffs[7],b_coeffs[6]])
    
    # Used in Pattern.computeOffsets()
    refpix = {}
    refpix['XREF'] = None 
    refpix['YREF'] = None
    refpix['V2REF'] = None 
    refpix['V3REF'] = None
    refpix['XDELTA'] = 0.
    refpix['YDELTA'] = 0.
    refpix['PSCALE'] = None
    refpix['centered'] = yes
    
    return fx,fy,refpix,order

def rotateCubic(fxy,theta):
    # This function transforms cubic coefficients so that
    # they calculate pixel positions oriented by theta (the same 
    # orientation as the PC).
    # Parameters: fxy - cubic-coefficients Numeric array
    #               theta - angle to rotate coefficients
    # Returns new array with same order as 'fxy'
    #
    # Set up some simplifications    
    newf = fxy * 0.
    cost = N.cos(DEGTORAD(theta))
    sint = N.sin(DEGTORAD(theta))
    cos2t = pow(cost,2)
    sin2t = pow(sint,2)
    cos3t = pow(cost,3)
    sin3t = pow(sint,3)
    
    # Now compute the new coefficients
    newf[1][1] = fxy[1][1] * cost - fxy[1][0] * sint
    newf[1][0] = fxy[1][1] * sint + fxy[1][0] * cost
    
    newf[2][2] = fxy[2][2] * cos2t - fxy[2][1] * cost * sint + fxy[2][0] * sin2t
    newf[2][1] = fxy[2][2] * 2 * cost * sint + fxy[2][1] * (cos2t - sin2t) + fxy[2][0] * 2 * sint * cost
    newf[2][0] = fxy[2][2] * sin2t + fxy[2][1] * cost * sint + fxy[2][0] * cos2t

    newf[3][3] = fxy[3][3] * cos3t - fxy[3][2] * sint * cos2t + fxy[3][1] * sin2t * cost - fxy[3][0] * sin3t
    newf[3][2] = fxy[3][3] * 3. * cos2t * sint + fxy[3][2]* (cos3t - 2. * sin2t * cost) +fxy[3][1] * (sin3t + 2 * sint * cos2t) - fxy[3][0] * sin2t * cost
    newf[3][1] = fxy[3][3] * 3. * cost * sin2t + fxy[3][2] *(2.*cos2t*sint - sin3t) + fxy[3][1] * (2 * sin2t * cost + cos3t) + fxy[3][0] * sint * cos2t
    newf[3][0] = fxy[3][3] * sin3t + fxy[3][2] * sin2t * cost + fxy[3][1] * sint * cos2t + fxy[3][0] * cos3t
    
    return newf

		
def buildRotMatrix(theta):
	_theta = DEGTORAD(theta)
	_mrot = N.zeros(shape=(2,2),dtype=N.float32)
	_mrot[0] = (N.cos(_theta),N.sin(_theta))
	_mrot[1] = (-N.sin(_theta),N.cos(_theta))

	return _mrot

def rotatePos(pos, theta,offset=None,scale=None):

    if scale == None:
        scale = 1.
        
    if offset == None:
        offset = N.array([0.,0.])
    mrot = buildRotMatrix(theta)
    xr = ((pos[0] * mrot[0][0]) + (pos[1]*mrot[0][1]) )/ scale + offset[0]
    yr = ((pos[0] * mrot[1][0]) + (pos[1]*mrot[1][1]) )/ scale + offset[1]

    return xr,yr        
        
# Function for determining the positions of the image
# corners after applying the geometry model.
# Returns a dictionary with the limits and size of the 
# image.
def getRange(members,ref_wcs,verbose=None):
    xma,yma = [],[]
    xmi,ymi = [],[]
    #print "WZ: running multidrizzle and get range"
    #pdb.set_trace()
    # Compute corrected positions of each chip's common point	
    crpix = (ref_wcs.crpix1,ref_wcs.crpix2)
    ref_rot = ref_wcs.orient
    _rot = ref_wcs.orient - members[0].geometry.wcslin.orient

    for member in members:	
        # Need to make sure this is populated with proper defaults
        # for ALL exposures!
        _model = member.geometry.model
        _wcs = member.geometry.wcs
        _wcslin = member.geometry.wcslin
        _theta = _wcslin.orient - ref_rot
        _ratio = ref_wcs.pscale / _wcslin.pscale
        _mod_ratio = _model.pscale / _wcslin.pscale
        _pscale = ref_wcs.pscale * _mod_ratio
        
        if _model.refpix != None:
            xdelta = _model.refpix['XDELTA'] / _ratio
            ydelta = _model.refpix['YDELTA'] / _ratio
        else:
            xdelta = 0.0
            ydelta = 0.0
        
        # Compute the corrected corner positions for each chip
        xypos = member.calcNewEdges(delta=(xdelta,ydelta),pscale=_pscale)
        
        if _theta != 0.0:
            #rotate coordinates to match output orientation 
            # Now, rotate new coords
            _mrot = buildRotMatrix(_theta)
            xypos = N.dot(xypos,_mrot)

        _oxmax = N.maximum.reduce(xypos[:,0])
        _oymax = N.maximum.reduce(xypos[:,1])
        _oxmin = N.minimum.reduce(xypos[:,0])
        _oymin = N.minimum.reduce(xypos[:,1])  

        xma.append(_oxmax)
        yma.append(_oymax)
        xmi.append(_oxmin)
        ymi.append(_oymin)

    # Determine the full size of the metachip
    xmax = N.maximum.reduce(xma)
    ymax = N.maximum.reduce(yma)
    ymin = N.minimum.reduce(ymi)
    xmin = N.minimum.reduce(xmi)
    
    # Compute offset from center that distortion correction shifts the image.
    # This accounts for the fact that the output is no longer symmetric around
    # the reference position...
    # We are adding one here to account for integerization of values.
    nref = ( abs(xmin + xmax + 0.5), abs(ymin + ymax + 0.5) )
    
    # Now, compute overall size based on range of pixels and offset from center.
    xsize = int(xmax - xmin + nref[0]+ 0.5) 
    ysize = int(ymax - ymin + nref[1]+ 0.5) 
        
    meta_range = {}
    meta_range = {'xmin':xmin,'xmax':xmax,'ymin':ymin,'ymax':ymax,'nref':nref}	
    meta_range['xsize'] = xsize
    meta_range['ysize'] = ysize 
    
    if verbose == yes:
        print 'Meta_WCS:' 
        print '    NREF         :',nref
        print '    X range      :',xmin,xmax
        print '    Y range      :',ymin,ymax
        print '    Computed Size: ',xsize,ysize

    return meta_range

def buildNewRootname(filename,extn=None):
    """ Build rootname for a new file.
        Use 'extn' for new filename if given, does NOT
        append a suffix/extension at all.
        
        Does NOT check to see if it exists already.  
        Will ALWAYS return a new filename.
    """
    # Search known suffixes to replace ('_crj.fits',...)
    extlist = EXTLIST
    # Also, add a default where '_dth.fits' replaces
    # whatever extension was there ('.fits','.c1h',...)
    extlist.append('.')
    for suffix in extlist:
        _indx = string.find(filename,suffix)
        if _indx > 0: break

    if _indx < 0:
         # default to entire rootname
        _indx = len(filename)

    if extn == None: extn = ''
    
    return filename[:_indx]+extn
    
def buildRootname(filename,ext=None):
    """
    Built a rootname for an existing file and given extension.
    Any user supplied extensions to use for searching for file
    need to be provided as a list of extensions.  
    
    Usage:
        rootname = buildRootname(filename,ext=['_dth.fits'])
        
    """
    # Get complete list of filenames from current directory
    flist = os.listdir(os.curdir)
    #First, assume given filename is complete and verify
    # it exists...
    #pdb.set_trace()
    rootname = None
    extlist = []
    for name in flist:
        if name == filename:
            rootname = filename
    
    # If we have an incomplete filename, try building a default
    # name and seeing if it exists...
    if rootname == None:
        # Set up default list of extensions to try...
        if ext == None:
            extlist = EXTLIST
        else:
            for i in ext:
                extlist.append(i)
        
        for extn in extlist:
            rname = string.lower(filename) + extn
            for name in flist:
                if rname == name:
                    rootname = name
                    break
            if rootname != None:
                break

    # If we still haven't found the file, see if we have the
    # info to build one...
    if rootname == None and ext != None:
        # Check to see if we have a full filename to start with...
        _indx = string.find(filename,'.')
        if _indx > 0: 
            rootname = filename[:_indx]+ext[0]
        else:
            rootname = filename + ext[0]

    # It will be up to the calling routine to verify 
    # that a valid rootname, rather than 'None', was returned.
    return rootname



