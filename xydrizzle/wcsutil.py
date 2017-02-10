import string, copy, os

import pyfits
import numpy as N

import xydrizzle, drutil, fileutil
import pdb

from pyraf import iraf

# Convenience definitions...
yes = iraf.yes
no = iraf.no
DEGTORAD = fileutil.DEGTORAD
RADTODEG = fileutil.RADTODEG
DIVMOD = fileutil.DIVMOD

#
# History
#
# 30-Mar-2002 WJH: Added separate XYtoSky interface function.
# 19-Apr-2002 WJH: Corrected 'ddtohms' for error in converting neg. dec.
#

#################
#
#
#		Coordinate Transformation Functions
#
#
#################		

def XYtoSky(input, pos, idckey='IDCTAB', linear=yes, verbose=no):
    """ Convert input pixel position(s) into RA/Dec position(s).
        Output will be either an (ra,dec) pair or a 'list' of (ra,dec) 
        pairs, not a numpy, to correspond with the input position(s).
        
        Parameter:
            input - Filename with extension specification of image
            pos   - Either a single [x,y] pair or a list of [x,y] pairs.
            idckey - Keyword which points to the IDC table to be used.
            linear - If no, apply distortion correction for image.
             
    """
    
    # Start by making sure we have a valid extension specification.
    _insplit = string.split(input,'[')
    if len(_insplit) == 1: 
        raise IOError, 'No extension specified for input image!'
    
    # Now we need to insure that the input is an array:
    if not isinstance(pos,N.ndarray):
        if N.array(pos).getrank() > 1:
            pos = N.array(pos)
       
    # Set up Exposure object
    _exposure = xydrizzle.Exposure(input,idckey=idckey)

    ra,dec = _exposure.geometry.XYtoSky(pos,linear=linear,verbose=verbose)
    
    if not isinstance(ra,N.ndarray):
        # We are working with a single input, return single values
        return ra,dec
    else:
        # We are working with arrays, so we need to convert them 
        # from 2 arrays with RA in one and Dec in the other to 1 array
        # with pairs of RA/Dec values.
        _radec = N.zeros(shape=(len(ra),2),dtype=ra.dtype)
        _radec[:,0] = ra
        _radec[:,1] = dec
        
        return _radec.tolist()


def ddtohms(xsky,ysky,verbose=no):

    """ Convert sky position(s) from decimal degrees to HMS format."""

    xskyh = xsky /15.
    xskym = (xskyh - N.floor(xskyh)) * 60.
    xskys = (xskym - N.floor(xskym)) * 60.

    yskym = (N.abs(ysky) - N.floor(N.abs(ysky))) * 60.
    yskys = (yskym - N.floor(yskym)) * 60.

    if isinstance(xskyh,N.ndarray):
        rah,dech = [],[]
        for i in xrange(len(xskyh)):
            rastr = repr(int(xskyh[i]))+':'+repr(int(xskym[i]))+':'+repr(xskys[i])
            decstr = repr(int(ysky[i]))+':'+repr(int(yskym[i]))+':'+repr(yskys[i])
            rah.append(rastr)
            dech.append(decstr)
            if verbose:
                print 'RA = ',rastr,', Dec = ',decstr
    else:
        rastr = repr(int(xskyh))+':'+repr(int(xskym))+':'+repr(xskys)
        decstr = repr(int(ysky))+':'+repr(int(yskym))+':'+repr(yskys)
        rah = rastr
        dech = decstr
        if verbose:
            print 'RA = ',rastr,', Dec = ',decstr
            
    return rah,dech

#################
#
#
#		Coordinate System Class
#
#
#################		

class WCSObject:
    """ This class should contain the WCS information from the 
	    input exposure's header and provide conversion functionality 
	    from pixels to RA/Dec and back. 
        
        It knows how to update the CD matrix for the new solution.

        rootname: this needs to be in a format supported by IRAF 
        and 'keypar', specifically: 
                filename.hhh[group] or filename.fits[ext]
    """
    def __init__(self, rootname,shape=None):
        # Initialize wcs dictionaries: 
        #   wcsdef - default values for new images
        #   wcstrans - translation table from header keyword to attribute 
        self.wcsdef = {'crpix1':0.0,'crpix2':0.0,'crval1':0.0,'crval2':0.0,'cd11':1.0,
	        'cd12':1.0,'cd21':1.0,'cd22':1.0,'orient':1.0,'naxis1':0,'naxis2':0,'pscale':1.0}
        self.wcstrans = {'CRPIX1':'crpix1','CRPIX2':'crpix2','CRVAL1':'crval1','CRVAL2':'crval2',
            'CD1_1':'cd11','CD1_2':'cd12','CD2_1':'cd21','CD2_2':'cd22',
            'ORIENTAT':'orient', 'NAXIS1':'naxis1','NAXIS2':'naxis2','pixel scale':'pscale'}
        # Now, read in the CRPIX1/2, CRVAL1/2, CD1/2_1/2 keywords.
        # Simplistic, but easy to understand what you are asking for. 

        _exists = yes
        if rootname != None:
            self.rootname = rootname
        else:
            self.rootname = 'New'
        # Look for extension specification in rootname
        _indx = string.find(self.rootname,'[')
        # If none are found, use entire rootname
        if _indx < 0:
            _indx = len(self.rootname)
        
        #pdb.set_trace()
        # Determine whether we are working with a new image or not.
        _filename = self.rootname[:_indx]
        _dirindx = string.rfind(_filename,os.sep)
        _logindx = string.rfind(_filename,'$')
        # Check to see if there is a path pre-pended to filename
        if _dirindx < 0 and _logindx < 0: 
            # No path...
            _dir = None
            _rootname = _filename

        else:
            # Path specified is pulled out, rootname found seperately.
            if _dirindx < 0: _dirindx = len(_filename)
            _dir = _filename[:_dirindx]
            _logindx = string.rfind(_dir,'$')
            if _logindx > 0:
                # We have an IRAF logical given in path that needs expanding
                _dir = iraf.osfn(_dir)

            _rootname = _filename[_dirindx+1:]
            _filename = _dir+os.sep+_rootname
            
        _exists = fileutil.checkFileExists(_rootname,directory=_dir)

        if string.find(_filename,'.fits') > 0:
            _fits = yes
        else:
            _fits = no

        if _exists == yes and _fits == yes:
            _fimg = pyfits.open(_filename)
            # Extract extnum
            if _indx > 0:
                _indx2 = string.find(self.rootname[_indx+1:],']')
                if self.rootname[_indx+1:_indx+_indx2+1] in list(string.digits):
                    # We were given an extension number directly
                    extnum = int(self.rootname[_indx+1:_indx+_indx2+1])
                else:
                    # Look for extension with given name
                    _indxc = string.find(self.rootname[_indx+1:],',')
                    if _indxc < 0:
                        # No EXTVER number given, so use first occurance of 
                        # extension with specified EXTNAME
                        extnum = self.rootname[_indx+1:_indx+_indx2+1]
                    else:
                        # extract EXTNAME,EXTVER number and return extension number
                        _extname = self.rootname[_indx+1:_indx+_indxc+1]
                        _extn = int(self.rootname[_indx+_indxc+2:_indx+_indx2+1])
                        
                        extnum = fileutil.findExtname(_fimg,_extname,extver=_extn)
            
            # Initialize WCS object with keyword values...  
            try:
                for key in self.wcstrans.keys():
                    _dkey = self.wcstrans[key]
                    if _dkey != 'pscale':
                        self.__dict__[_dkey] = _fimg[extnum].header[key]

                self.new = no

            except:
                raise IOError,'Image %s does not contain valid WCS keywords!' % self.rootname
            _fimg.close()
            del _fimg
        elif _exists == yes and _fits == no:
                        # Initialize WCS object with keyword values...
            try:
                iraf.keypar(self.rootname,'i_naxis1',silent='yes')
                self.naxis1 = int(iraf.keypar.value)
                iraf.keypar(self.rootname,'i_naxis2',silent='yes')
                self.naxis2 = int(iraf.keypar.value)
                
                for key in self.wcstrans.keys():
                    _dkey = self.wcstrans[key]
                    if _dkey != 'pscale' and _dkey != 'naxis1' and _dkey != 'naxis2':
                        iraf.keypar(self.rootname,key,silent='yes')
                        self.__dict__[_dkey] = float(iraf.keypar.value) 

                self.new = no

            except:
                raise IOError,'Image %s does not contain valid WCS keywords!' % self.rootname

        else:
            # or set default values...
            self.new = yes
            for key in self.wcsdef.keys():
                self.__dict__[key] = self.wcsdef[key]
            if shape != None:
                # ... and update with user values.
                self.naxis1 = int(shape[0])
                self.naxis2 = int(shape[1])
                self.pscale = float(shape[2])

        # Set up a true CD array which can be used numerically
        #self.cd = array([[self.cd11,self.cd12],[self.cd21,self.cd22]])
        
        if shape == None:
            self.pscale = N.sqrt(N.power(self.cd11,2)+N.power(self.cd21,2)) * 3600.
        # Establish an attribute for the linearized orient
        # defined as the orientation of the CD after applying the default
        # distortion correction.
        self._orient_lin = 0.

    # You never know when you want to print out the WCS keywords...
    def __str__(self):
        block = 'WCS Keywords for ' + self.rootname + ': \n'
        for key in self.wcstrans.keys():
            _dkey = self.wcstrans[key]
            strn = string.upper(key) + " = " + repr(self.__dict__[_dkey]) + '\n'
            block = block + strn

        return block

    def updateWCS(self, pixel_scale=None, orient=None,refpos=None,refval=None,size=None):
        """
        Create a new CD Matrix from the absolute pixel scale 
        and reference image orientation.
        """
        # Set up parameters necessary for updating WCS
        # Check to see if new value is provided,
        # If not, fall back on old value as the default

        _updateCD = no
        if orient != None:
            pa = DEGTORAD(orient)
            self.orient = orient
            self._orient_lin = orient
            _updateCD = yes
        else:
            # In case only pixel_scale was specified
            pa = DEGTORAD(self.orient)
        
        if pixel_scale != None:
            self.pscale = pixel_scale
            _updateCD = yes
        else:
            # In case, only orient was specified
            pixel_scale = self.pscale
            
        if size != None:
	        self.naxis1 = size[0]
	        self.naxis2 = size[1]

        if refpos != None:
	        self.crpix1 = refpos[0]
	        self.crpix2 = refpos[1]
        if self.crpix1 == None:
	        self.crpix1 = self.naxis1/2.
	        self.crpix2 = self.naxis2/2.

        if refval != None:
	        self.crval1 = refval[0]
	        self.crval2 = refval[1]

        # Reset WCS info now...
        if _updateCD:
            # Only update this should the pscale or orientation change...
            pscale = pixel_scale / 3600.

            self.cd11 = -pscale * N.cos(pa)
            self.cd12 = pscale * N.sin(pa)
            self.cd21 = pscale * N.sin(pa)
            self.cd22 = pscale * N.cos(pa)

		
    def xy2rd(self,pos):
        """
        This method would apply the WCS keywords to a position to
        generate a new sky position. 

        The algorithm comes directly from 'imgtools.xy2rd'

        translate (x,y) to (ra, dec)
        """
        if isinstance(pos,N.ndarray):
            # If we are working with an array of positions, 
            # point to just X and Y values
            posx = pos[:,0]
            posy = pos[:,1]
        else:
            # Otherwise, we are working with a single X,Y tuple
            posx = pos[0]
            posy = pos[1]
            
        xi = self.cd11 * (posx - self.crpix1) + self.cd12 * (posy - self.crpix2)
        eta = self.cd21 * (posx - self.crpix1) + self.cd22 * (posy - self.crpix2)

        xi = DEGTORAD(xi)
        eta = DEGTORAD(eta)
        ra0 = DEGTORAD(self.crval1)
        dec0 = DEGTORAD(self.crval2)

        ra = N.arctan((xi / (N.cos(dec0)-eta*N.sin(dec0)))) + ra0
        dec = N.arctan( ((eta*N.cos(dec0)+N.sin(dec0)) / 
                (N.sqrt((N.cos(dec0)-eta*N.sin(dec0))**2 + xi**2))) )

        ra = RADTODEG(ra)
        dec = RADTODEG(dec)
        ra = DIVMOD(ra, 360.)
        
        # Otherwise, just return the RA,Dec tuple.
        return ra,dec


    def rd2xy(self,skypos,hour=no):
        """
        This method would use the WCS keywords to compute the XY position
        from a given RA/Dec tuple (in deg).
        """
        det = self.cd11*self.cd22 - self.cd12*self.cd21

        if det == 0.0: 
	        raise ArithmeticError,"singular CD matrix!"

        cdinv11 = self.cd22 / det
        cdinv12 = -self.cd12 / det
        cdinv21 = -self.cd21 / det
        cdinv22 = self.cd11 / det

        # translate (ra, dec) to (x, y)

        ra0 = DEGTORAD(self.crval1)
        dec0 = DEGTORAD(self.crval2)
        if hour == yes:
            skypos[0] = skypos[0] * 15.
        ra = DEGTORAD(skypos[0])
        dec = DEGTORAD(skypos[1])

        bottom = N.sin(dec)*N.sin(dec0) + N.cos(dec)*N.cos(dec0)*N.cos(ra-ra0)
        if bottom == 0.0: 
            raise ArithmeticError,"Unreasonable RA/Dec range!"

        xi = N.cos(dec) * N.sin(ra-ra0) / bottom
        eta = (N.sin(dec)*N.cos(dec0) - N.cos(dec)*N.sin(dec0)*N.cos(ra-ra0)) / bottom
        xi = RADTODEG(xi)
        eta = RADTODEG(eta)

        x = cdinv11 * xi + cdinv12 * eta + self.crpix1
        y = cdinv21 * xi + cdinv22 * eta + self.crpix2

        return x,y

    def rotateCD(self,orient):
        """ Rotates WCS CD matrix to new orientation given by 'orient'
        """
        # Determine where member CRVAL position falls in ref frame            
        # Find out whether this needs to be rotated to align with
        # reference frame.

        _delta = self.orient - orient
        if _delta == 0.: 
            return

        # Start by building the rotation matrix...
        _rot = drutil.buildRotMatrix(_delta)
        # ...then, rotate the CD matrix and update the values...
        _cd = N.array([[self.cd11,self.cd12],[self.cd21,self.cd22]])
        _cdrot = N.dot(_cd,_rot)
        self.cd11 = _cdrot[0][0]
        self.cd12 = _cdrot[0][1]
        self.cd21 = _cdrot[1][0]
        self.cd22 = _cdrot[1][1]
        self.orient = orient


    def copy(self,deep=yes):
        """ Makes a (deep)copy of this object for use by other objects.
        """
        if deep:
            return copy.deepcopy(self)   
        else:
            return copy.copy(self)
