# Wei Zheng February 2011
# This version is more advanced as it contain
import types,string,os,copy

from pyraf import iraf
from iraf import stsdas, dither
import pdb
# Import PyDrizzle utility modules
import buildmask, fileutil, wcsutil, drutil, matutil, pyfits #WZ
import numpy as N

#Add buildasn to namespace for use by other programs
import buildasn

yes = iraf.yes
no = iraf.no

# List of supported instruments/detectors	
INSTRUMENT = ["ACS","WFPC2","STIS","WFC3"]
DEXPTIME = 'EXPTIME'
GLOBALfilts = ['NULL','NULL']

# Version 
version = "2.6.2 (14-Jul-2006)"

# For History of changes and updates, see 'History'

def _ptime():
    import time

    # Format time values for keywords IRAF-TLM, and DATE 
    _ltime = time.localtime(time.time())
    tlm_str = time.strftime('%H:%M:%S (%d/%m/%Y)',_ltime)
    #date_str = time.strftime('%Y-%m-%dT%H:%M:%S',_ltime)
    return tlm_str



#################
#
#
#		Geometry/Distortion Classes
#
#
#################		

class GeometryModel:
    """
    Base class for Distortion model. 
    There will be a separate class for each type of 
    model/filetype used with drizzle, i.e., IDCModel and 
    DrizzleModel. 
    
    Each class will know how to apply the distortion to a
    single point and how to convert coefficients to an input table
    suitable for the drizzle task.    

    Coefficients will be stored in CX,CY arrays.
    """
    # 	
    # 	
    # 	
    # 	 
    # 	
    # 	
    # 	
    NORDER = 3

    def __init__(self):
        "  This will open the given file and determine its type and norder."

        # 	Method to read in coefficients from given table and 
        # 	populate the n arrays 'cx' and 'cy'. 
        # 	This will be different for each type of input file,
        # 	IDCTAB vs. drizzle table.		

        # Set these up here for all sub-classes to use...
        # But, calculate norder and cx,cy arrays in detector specific classes.
        self.cx = None
        self.cy = None
        self.refpix = None
        self.norder = self.NORDER

        # default values for these attributes
        self.direction = 'forward'

        self.pscale = 1.0

    def convert(self, tmpname):
        """
         Open up an ASCII file, output coefficients in drizzle
          format after converting them as necessary.
        First, normalize these coefficients to what drizzle expects
        Normalize the coefficients by the MODEL/output plate scale.
        
        16-May-2002:
        Revised to work with higher order polynomials by John Blakeslee.
        """
        cx = self.cx / self.pscale
        cy = self.cy / self.pscale
        x0 = self.refpix['XDELTA']
        y0 = self.refpix['YDELTA']

        # Now, write out the coefficients into an ASCII 
        # file in 'drizzle' format.
        lines = []

        lines.append('# ACS polynomial distortion coefficients\n')
        lines.append('# Extracted from "%s" \n'%self.name)
        if self.norder==3:
            lines.append('cubic\n')
        elif self.norder==4:
            lines.append('quartic\n')
        elif self.norder==5:
            lines.append('quintic\n')
        else:
            raise ValueError, "Drizzle cannot handle poly distortions of order %d"%self.norder
    
        str = '%16.8f %16.8g %16.8g %16.8g %16.8g \n'% (x0,cx[1,1],cx[1,0],cx[2,2],cx[2,1])
        lines.append(str)
        str = '%16.8g %16.8g %16.8g %16.8g %16.8g \n'% (cx[2,0],cx[3,3],cx[3,2],cx[3,1],cx[3,0])
        lines.append(str)
        if self.norder>3:
            str = '%16.8g %16.8g %16.8g %16.8g %16.8g \n'% (cx[4,4],cx[4,3],cx[4,2],cx[4,1],cx[4,0])
            lines.append(str)
        if self.norder>4:
            str = '%16.8g %16.8g %16.8g %16.8g %16.8g %16.8g \n'% (cx[5,5],cx[5,4],cx[5,3],cx[5,2],cx[5,1],cx[5,0])
            lines.append(str)
        lines.append("\n")
        
        str = '%16.8f %16.8g %16.8g %16.8g %16.8g \n'% (y0,cy[1,1],cy[1,0],cy[2,2],cy[2,1])
        lines.append(str)
        str = '%16.8g %16.8g %16.8g %16.8g %16.8g \n'% (cy[2,0],cy[3,3],cy[3,2],cy[3,1],cy[3,0])
        lines.append(str)
        if self.norder>3:
            str = '%16.8g %16.8g %16.8g %16.8g %16.8g \n'% (cy[4,4],cy[4,3],cy[4,2],cy[4,1],cy[4,0])
            lines.append(str)
        if self.norder>4:
            str = '%16.8g %16.8g %16.8g %16.8g %16.8g %16.8g \n'% (cy[5,5],cy[5,4],cy[5,3],cy[5,2],cy[5,1],cy[5,0])
            lines.append(str)

        output = open(tmpname,'w')
        output.writelines(lines)
        output.close()

    def apply(self, pixpos,scale=1.0):
        """
         Apply coefficients to a pixel position or a list of positions. 
          This should be the same for all coefficients tables.
        Return the geometrically-adjusted position 
        in arcseconds from the reference position as a tuple (x,y).
              
        Compute delta from reference position
        """     
        if self.cx == None:
            return pixpos
        
        _cx = self.cx / scale
        _cy = self.cy / scale
        _convert = no 
        _p = pixpos
        
        if not isinstance(_p, N.ndarray):
            _p = N.array(_p)
            _convert = yes
      
        dxy = _p - (self.refpix['XREF'],self.refpix['YREF'])
                
        # Apply coefficients from distortion model here...
        c = _p * 0.
        for i in range(self.norder+1):
            for j in range(i+1):
                c[:,0] = c[:,0] + _cx[i][j] * pow(dxy[:,0],j) * pow(dxy[:,1],(i-j))
                c[:,1] = c[:,1] + _cy[i][j] * pow(dxy[:,0],j) * pow(dxy[:,1],(i-j))
        xc = c[:,0]
        yc = c[:,1]
        
        # Convert results back to same form as original input 
        if _convert:
            xc = xc.tolist()
            yc = yc.tolist()        
            # If a single tuple was input, return just a single tuple
            if len(xc) == 1:
                xc = xc[0]
                yc = yc[0]
                
        return xc,yc

class IDCModel(GeometryModel):
    """  
    This class will open the IDCTAB, select proper row based on 
    chip/direction and populate cx,cy arrays.
    We also need to read in SCALE, XCOM,YCOM, XREF,YREF as well.
    """
    def __init__(self, idcfile, chip=1, direction='forward'):
        GeometryModel.__init__(self)
        #
        # Norder must be derived from the coeffs file itself, 
        # then the arrays can be setup. Thus, it needs to be 
        # done in the sub-class, not in the base class.
        # Read in table. 
        # Populate cx,cy,scale, and other variables here.
        #
        self.name = idcfile

        #print "Check IDC table" # WZ
        f=pyfits.open(idcfile)
        det=f[0].header.get('DETECTOR')
        f.close()
        if (det == 'WFC'):
            # pdb.set_trace()
            self.cx,self.cy,self.refpix,self.norder = matutil._readACSIDCtab(idcfile,
		        chip=chip, direction=direction, verbose=0,
                        filter1=GLOBALfilts[0], filter2=GLOBALfilts[1])
        else:
        # self.cx,self.cy,self.refpix,self.norder = fileutil.readIDCtab(idcfile,
            self.cx,self.cy,self.refpix,self.norder = matutil._readIDCtab(idcfile,
		        chip=chip, direction=direction, verbose=0,
                        filter=GLOBALfilts[0])
        self.pscale = self.refpix['PSCALE']
		
class DrizzleModel(GeometryModel):
    """  
    This class will read in an ASCII Cubic
    drizzle coeffs file and populate the cx,cy arrays. 
    """

    def __init__(self, idcfile):	
        GeometryModel.__init__(self)
        #
        # We now need to read in the file, populate cx,cy, and
        # other variables as necessary.	
        #
        self.name = idcfile
        self.cx,self.cy,self.refpix,self.norder = drutil.readCubicTable(idcfile)			
        self.pscale = self.refpix['PSCALE']        		

           
class TraugerModel(GeometryModel):
    """
    This class will read in the ASCII Trauger coeffs 
    file, convert them to SIAF coefficients, then populate
    the cx,cy arrays. 
    """
    NORDER = 3

    def __init__(self, idcfile,lam):	
        GeometryModel.__init__(self)		
        self.name = idcfile
        self.cx,self.cy,self.refpix,self.norder = drutil.readTraugerTable(idcfile,lam)
        self.pscale = self.refpix['PSCALE']
        #
        # Read in file here.
        # Populate cx,cy, and other variables. 
        #

	
class ObsGeometry:
    """
	Base class for Observation's geometric information.

	This class must know how to recognize the different
	types of distortion coefficients tables and instantiate
	the correct class for it. 
    """

    def __init__(self, rootname, idcfile, idckey=None, chip=1, direction="forward",new=None):
        """
         We need to setup the proper object for the GeometryModel
          based on the format of the provided idctab. 
        
         We need to trap attempts to address values of 'chip' that
         are not part of the observation; such as in sub-arrays.
        """
        self.idcfile = idcfile
        self.direction = direction	
        self.ikey = None
        
        if new == no:            
            self.wcs = wcsutil.WCSObject(rootname)
            self.wcslin = wcsutil.WCSObject(rootname)
            
            # Based on the filetype, open the correct geometry model
            ikey = string.lower(idckey)
            self.ikey = ikey
                   
            if ikey == 'idctab':
                # pdb.set_trace()
                self.model =  IDCModel(self.idcfile, 
                    chip=chip, direction=self.direction)
            elif ikey == 'cubic':
                self.model = DrizzleModel(self.idcfile)
            elif ikey == 'trauger':
                _lam = drutil.getPrimaryKeyword(rootname,'PHOTPLAM')
                self.model = TraugerModel(self.idcfile,float(_lam)/10.)
            else:
                raise ValueError, "Unknown type of coefficients table %s"%idcfile

            # Insure that a default model pscale has been set
            if self.model.refpix['PSCALE'] == None:
                self.model.pscale = 1.0
                self.model.refpix['PSCALE']  = self.model.pscale
            if self.model.refpix['XREF'] == None:
                self.model.refpix['XREF'] = self.wcs.crpix1
                self.model.refpix['YREF'] = self.wcs.crpix2
            # Generate linear WCS to linear CD matrix
            self.undistortWCS()
        else:
            # For new images with no distortion, CD matrix is sufficient.
            self.wcs = wcsutil.WCSObject(None)
            self.wcslin = wcsutil.WCSObject(None)
            self.model = GeometryModel()
        
    def apply(self, pixpos,delta=None,pscale=None,verbose=no):
        """
         This method applies the model to a pixel position 
          to calculate the new position. 
          Depending on the value of direction, this could mean
          going from distorted/raw to a corrected positions or
          the other way around.
          
          If a specific pixel scale is provided, this will be 
          used to determine the final output position.
        
        """
        if delta == None:
            # Use default from table.
            deltax = 0.
            deltay = 0.
        else:
            deltax = delta[0]
            deltay = delta[1]

        
        if pscale == None:
            pscale = self.model.pscale

        v2,v3 = self.model.apply(pixpos)

        # If there was no distortion applied to
        # the pixel position, simply shift by new
        # reference point.
        
        if self.model.cx == None:
            if self.model.refpix != None:
                if self.model.refpix['XREF'] == None:
                    refpos = (self.wcs.crpix1,self.wcs.crpix2)	
                else:
                    refpos = (self.model.refpix['XREF'],self.model.refpix['YREF'])
            else:
                refpos = (self.wcs.crpix1,self.wcs.crpix2)

            xpos = v2 - refpos[0] + deltax
            ypos = v3 - refpos[1] + deltay
        else:
            # For sub-arrays, we need to account for the offset
            # between the CRPIX of the sub-array and the model reference
            # position in the full frame coordinate system.
            # For full images, this correction will be ZERO.
            # This offset, though, has to be scaled by the relative plate-scales.
            _ratio = pscale / self.model.pscale
            _delta_refx = (self.wcs.crpix1 - self.model.refpix['XREF'])/_ratio
            _delta_refy = (self.wcs.crpix2 - self.model.refpix['YREF'])/_ratio
            
            xpos = ( v2 / pscale ) + deltax + _delta_refx
            ypos = ( v3 / pscale ) + deltay + _delta_refy
            
        # Return the geometrically-adjusted position as a tuple (x,y)
        return xpos,ypos	

    def undistortWCS(self):
        """
        This method applies the distortion to a 1x1 box at the reference
        position and updates the WCS based on the results.  
        
        This method is based directly on the 'drizzle' subroutine 'UPWCS'
        written by R. Hook.
        """
        # Check to see if we have a valid model to apply to the WCS
        if self.model.cx == None:
            # We do not have a model to apply, so simply copy the 
            # original WCS
            self.wcslin = self.wcs.copy()
            return
            
        # define the reference points
        _cpix1 = self.wcs.crpix1
        _cpix2 = self.wcs.crpix2
        _xy = N.array([(_cpix1,_cpix2),(_cpix1+1.,_cpix2),(_cpix1,_cpix2+1.)])

        # apply the distortion to them
        #_xc,_yc = self.model.apply(_xy,scale=self.model.pscale)
        _xc,_yc = self.apply(_xy)
        
        # Now, work out the effective CD matrix of the transformation
        _am = _xc[1] - _xc[0]
        _bm = _xc[2] - _xc[0]
        _cm = _yc[1] - _yc[0]
        _dm = _yc[2] - _yc[0]
        
        # Check the determinant for singularity
        _det = (_am * _dm) - (_bm * _cm)        
        if ( _det == 0.0):
            print 'Matrix is singular! Can NOT update WCS.'
            return 
        
        _a = _dm / _det
        _b = -_bm / _det
        _c = -_cm / _det
        _d = _am / _det
                
        self.wcslin.cd11 = _a * self.wcs.cd11 + _c * self.wcs.cd12
        self.wcslin.cd21 = _a * self.wcs.cd21 + _c * self.wcs.cd22
        self.wcslin.cd12 = _b * self.wcs.cd11 + _d * self.wcs.cd12
        self.wcslin.cd22 = _b * self.wcs.cd21 + _d * self.wcs.cd22
        self.wcslin.orient = N.arctan2(self.wcslin.cd12,self.wcslin.cd22) * 180./N.pi
        self.wcslin.pscale = N.sqrt(N.power(self.wcslin.cd11,2)+N.power(self.wcslin.cd21,2))*3600.
            
    def XYtoSky(self, pos,verbose=no,linear=no):
        """
        This method applies the distortion model to a pixel position
         and calculates the sky position in RA/Dec.
        """
        if linear == no and self.model.refpix != None:
            # Perform full solution including distortion
            dcx,dcy = self.model.apply(pos)
            
            if dcx != None:
                _crpix = (self.model.refpix['XREF'],self.model.refpix['YREF'])
                _pscale = self.model.pscale
                _cpos = (dcx/_pscale + _crpix[0], dcy/_pscale + _crpix[1])
                
                # Now apply linear CD matrix appropriate to model
                xsky,ysky = self.wcslin.xy2rd(_cpos)
            else:
                xsky,ysky = self.wcs.xy2rd(pos)

        else:
            xsky,ysky = self.wcs.xy2rd(pos)
            
        # Format the results for easy understanding, if desired...
        if verbose == yes:
            rastr,decstr = wcsutil.ddtohms(xsky,ysky,verbose=verbose)
            
        # Return the skypos as a tuple (x,y)
        return xsky,ysky

    def SkytoXY(self, skypos, verbose=no, hour=no):
        """
         This method applies the model to an RA/Dec 
          and calculates the pixel position.
          RA and Dec need to be in decimal form!
        
        This needs to be expanded to include full distortion model
        as well, i.e. inverse distortion solution.
        """
        x,y = self.wcs.rd2xy(skypos,hour=hour)

        if verbose == yes:
            print 'X = ',x,' Y = ',y

        # Return the pixel position as a tuple (x,y)
        # return pos
        return x,y

#################
#
#
#		Exposure Classes 
#
#
#################		
class Exposure:
    """  
    This class will provide the basic functionality for keeping 
    track of an exposure's parameters, including distortion model,
    WCS information, and metachip shape. 
    """

    def __init__(self,expname, idckey=None,new=no,wcs=None,mask=None,wht=None,nchips=None): #WZ

        # This name should be formatted for use in image I/O
        self.name = expname

        # Keep track of any associated mask file created for
        # this exposure from its DQ file, or other mask file.
        #pdb.set_trace()
        self.maskname = mask
        self.whtname = wht #WZ

        if new == no:		
            _chip = drutil.getChipId(expname)
            _extver = int(drutil.getPrimaryKeyword(expname,'EXTVER'))
            if _extver == None: _extver = 1
            # Determine whether there is any offset for the image
            # as in the case of subarrays (based on LTV keywords)
            _ltv1,_ltv2 = drutil.getLTVOffsets(expname)
        else:
            _chip = 1
            _ltv1,_ltv2 = 0.,0.
            _extver = 1

        self.chip = repr(_chip)		
        # We need the pixel value of the starting pixel
        # NOT the pixel value of the full image 0,0 pixel.
        self.offset = (-float(_ltv1),-float(_ltv2))
                
        # Read the name of idcfile from image header if not explicitly
        # provided by user.
        if idckey != None:
            idcfile = drutil.getIDCFile(expname,keyword=idckey)
        else:
            idcfile = None

        # Get distortion model and WCS info.
        self.geometry = ObsGeometry(expname, idcfile, idckey=idckey, 
                chip=_chip, new=new)
        
        # Remember the name of the IDC file used...
        self.idcfile = idcfile

        # Define shape here...
        # nx,ny,pixel scale
        #
        if wcs != None:
            # We have been passed a WCS to use
            self.geometry.wcs = wcs
            self.geometry.model.pscale = wcs.pscale
            if expname != None:
                self.geometry.wcs.rootname = expname

        self.naxis1 = self.geometry.wcs.naxis1
        self.naxis2 = self.geometry.wcs.naxis2
        self.pscale = self.geometry.wcs.pscale
        self.shape = (self.naxis1,self.naxis2,self.pscale)
               
        # Generate BLOT output name specific to this Exposure
        _blot_extn = '_blot_sci'+repr(_extver)+'.fits'
        self.outblot = drutil.buildNewRootname(self.name,extn=_blot_extn)

    def getShape(self):
        """
        This method gets the shape after opening and reading 
        the input image. This version will be the default
        way of doing this, but each instrument may override
        this method with one specific to their data format.
        """
        return self.shape

    def setShape(self,size,pscale):
        """
        This method will set the shape for a new file if
        there is no image header information. 
        
        Size will be defined as (nx,ny) and pixel size
        """
        self.shape = (size[0],size[1],pscale)

    def calcNewCorners(self,delta=None,pscale=None):
        """
        This method will compute a new shape based on the positions of 
        the corners AFTER applying the geometry model.
        
        Parameter: delta - offset from nominal crpix center position
        
        These new position for each corner should be calculated by calling
        self.geometry.apply() on each corner position.
        This should also take into account the output scale as well.
        """
        corners = N.zeros(shape=(4,2),dtype=N.float32)
        xin = [0] * 4
        yin = [0] * 4

        xin[0]=self.offset[0]
        xin[1]=self.offset[0]
        xin[2]=self.shape[0] + self.offset[0]
        xin[3]=self.shape[0] + self.offset[0]
        yin[0]=self.offset[1]
        yin[1]=self.shape[1] + self.offset[1]
        yin[2]=self.shape[1] + self.offset[1]
        yin[3]=self.offset[1]
        
        _scale = self.geometry.model.pscale
        _ratio = pscale / _scale
        
        for i in range(len(xin)):
            xout,yout = self.geometry.model.apply((xin[i],yin[i]),scale=_scale)
            corners[i] = ( (xout*_ratio) + delta[0] - xin[0], (yout*_ratio) + delta[1] - yin[0])
        
        return corners
    def calcNewEdges(self,delta=None,pscale=None):
        """
        This method will compute arrays for all the pixels around
        the edge of an image AFTER applying the geometry model.
        
        Parameter: delta - offset from nominal crpix center position
        
        Output:   xsides - array which contains the new positions for
                          all pixels along the LEFT and RIGHT edges
                  ysides - array which contains the new positions for
                          all pixels along the TOP and BOTTOM edges
        The new position for each pixel is calculated by calling
        self.geometry.apply() on each position.
        """
        #print 'Start calcNewEdges at ',_ptime()
        # build up output arrays

        # build up arrays for pixel positions for the edges
        # These arrays need to be: array([(x,y),(x1,y1),...])
        numpix = self.naxis1*2 + self.naxis2 * 2
        border = N.zeros(shape=(numpix,2),dtype=N.float32)

        # Now determine the appropriate values for this array
        # We also need to account for any subarray offsets
        xmin = self.offset[0]
        xmax = self.naxis1 + self.offset[0]
        ymin = self.offset[1]
        ymax = self.naxis2 + self.offset[1]
        
        # Build range of pixel values for each side
        # Add 1 to make them consistent with pixel numbering in IRAF
        # Also include the LTV offsets to represent position in full chip
        #   since the model works relative to full chip positions.
        xside = N.arange(self.naxis1) + 1 + xmin
        yside = N.arange(self.naxis2) + 1 + ymin
        
        #Now apply them to the array to generate the appropriate tuples
        #bottom
        _range0 = 0
        _range1 = self.naxis1
        border[_range0:_range1,0] = xside
        border[_range0:_range1,1] = ymin
        #top
        _range0 = _range1
        _range1 = _range0 + self.naxis1
        border[_range0:_range1,0] = xside
        border[_range0:_range1,1] = ymax
        #left
        _range0 = _range1
        _range1 = _range0 + self.naxis2
        border[_range0:_range1,0] = xmin
        border[_range0:_range1,1] = yside
        #right
        _range0 = _range1
        _range1 = _range0 + self.naxis2
        border[_range0:_range1,0] = xmax
        border[_range0:_range1,1] = yside

        
        # calculate new edge positions
        border[:,0],border[:,1] = self.geometry.apply(border,pscale=pscale)        

        border = border + delta
        
        #print 'Calculated corrected border positions at ',_ptime()

        return border

    def getWCS(self):
        return self.geometry.wcs
    def showWCS(self):
        print self.geometry.wcs

            
class Pattern:
    """ 
     Set default values for these to be overridden by 
     instrument specific class variables as necessary.
    """
    IDCKEY = 'IDCTAB'
    PARITY = {'detector':[[1.0,0.0],[0.0,1.0]]}
    REFDATA = {'detector':[[1.,1.],[1.,1.]]}
    SEPERATOR = ','             # Used for parsing EXTN syntax
    
    SCI_EXTN = '[sci,1]'        # default specification for SCI extension/group
    NUM_IMSET = 3               # Number of extensions in an IMSET
    EXT0 = '[0]'
    
    def __init__(self, pars=None):
        # Set these up for use...
        self.members = []
        self.pars = pars

        # Set IDCKEY, if specified by user...
        if self.pars['idckey'] == None:
            self.idckey = self.IDCKEY
        else:
            self.idckey = self.pars['idckey']
            
    def setFileFormat(self,filename):
        """ This method examines the image to determine whether
            it is a simple FITS file without extensions.  If so, the
            it resets the SCI_EXTN,EXT0,and SEPERATOR attributes to
            work correctly with the file format. 
        """
        
    def addMembers(self,filename):
        """ Build rootname for each SCI extension, and 
            create the mask image from the DQ extension. 
            It would then append a new Exposure object to 'members' 
            list for each extension. 
        """
        # Build rootname here for each SCI extension...
        #print "sci input images"
        #pdb.set_trace()
        for i in range(self.nmembers):
            extname = filename+self._setSciExtn(extn=i+1)
            maskname = buildmask.buildMaskImage(filename,self.pars['bitvalue'],extver=i+1)
            #self.members.append(Exposure(extname, idckey=self.idckey, mask=maskname, nchips=self.nmembers)) #WZ
            if (maskname != None):
                #pdb.set_trace()
                whtname = string.replace(maskname,'_inmask','_wt')
            else:
                whtname=None
            # self.members.append(Exposure(extname, idckey=self.idckey, mask=maskname, wht=whtname, nchips=self.nmembers)) #WZ
            self.members.append(Exposure(extname, idckey=self.idckey, mask=maskname, wht=maskname, nchips=self.nmembers)) #WZ
                    
    def _setSciExtn(self,extn=None):
        """ Builds extension specification for accessing SCI extension/group.
        """
        if extn != None:
            _lensep = len(self.SEPERATOR)
            _indx = string.find(self.SCI_EXTN,self.SEPERATOR) + _lensep
            return self.SCI_EXTN[:_indx]+repr(extn)+self.SCI_EXTN[_indx+1:]
        else:
            return self.SCI_EXTN        
    
    def buildProduct(self,filename,output,ref): #WZ
        """
        Create Exposure object for meta-chip product after applying 
        distortion model to input members.
        """
        # Now, set up output product based on info from input files
        _psize = self.pars['psize']        
        _rot = self.pars['rot']

        #pdb.set_trace()
        output_wcs = self.buildMetachip(ref,psize=_psize,orient=_rot) #WZ
        
        # Update the final output size with the meta-chip size
        self.size = (output_wcs.naxis1,output_wcs.naxis2)
        
        # Need to compute shape for dither product here based on output_shapes
        # and reference pixels for each observation.	
        print "build product for ", filename
        #pdb.set_trace()
        self.product = Exposure(self.outdata,wcs=output_wcs, new=yes)
        self.product.exptime = self.exptime
        self.ref=ref #WZ
        
    def buildMetachip(self,ref,psize=None,orient=None): #WZ
        """ Build up the new metashape based on the 
         corrected size and position for each Exposure.
         (Pattern method)
        """    
        #pdb.set_trace()
        # Build an initial value for size of meta-chip
        _geo_ref = self.members[0].geometry
        self.ref = ref.copy() #WZ
        #wcs_ref = self.ref.copy() #WZ temporarily commented out, Feb 18 
        _wcs_ref = _geo_ref.wcslin
        _model_ref = _geo_ref.model

        # Get pscale from model
        if psize == None:
            pscale = _model_ref.pscale
        else: 
            # User specified a different pixel scale to use...
            pscale = psize

        print "buildMetachip"
        #pdb.set_trace() #WZ
        _shape = (_wcs_ref.naxis1,_wcs_ref.naxis2,pscale)	
        # Build new WCS for output metachip here       
        # It will be based on the undistorted version of members[0] WCS
        meta_wcs = _geo_ref.wcslin.copy() # WZ test Mar 2011: Ste back
        #meta_wcs = self.ref.copy() #WZ Here is the difference prod?

        # Use the undistorted WCS to update the linearized orientation
        # for all the input chips...
        for member in self.members:
            member.geometry.wcs._orient_lin = meta_wcs.orient
                    
        # Now compute the CRPIX for the new output frame
        # for the CRVAL from the reference chip
        _crval = (_wcs_ref.crval1,_wcs_ref.crval2)
        _cpix = (_wcs_ref.crpix1, _wcs_ref.crpix2)
    
        # Update meta_wcs with initial guess for CRPIX and CRVAL
        meta_wcs.updateWCS(refval=_crval,refpos=_cpix,pixel_scale=pscale)
        if orient != None:
            meta_wcs.orient = orient
        #Determine range of pixel values for corrected image
        # Using verbose=yes will return additional info on range calculation
        meta_range = drutil.getRange(self.members,meta_wcs,verbose=no)
        if orient != None:
            meta_wcs.orient = _geo_ref.wcslin.orient # set back WZ
            #meta_wcs.orient = self.prod.orient #WZ
        
        # Update WCS based on new size
        xsize = int(meta_range['xsize'] + 1.0)
        ysize = int(meta_range['ysize'] + 1.0)
        meta_wcs.naxis1 = xsize
        meta_wcs.naxis2 = ysize
        cen = ((xsize/2.),(ysize/2.))
        if orient != None:
            # Do final rotation of WCS, if user specified one
            meta_wcs.rotateCD(orient)
                
        # Shifts position of CRPIX to reflect new size
        # Instead of being centered on (0.,0.) like the original guess.
        # CRVAL for this position remains the same, though, as it is still
        # the same point in the sky/image.
        #        
        _ratio = _model_ref.refpix['PSCALE'] / pscale
        _rpix = ( _model_ref.refpix['XDELTA'] * _ratio, _model_ref.refpix['YDELTA'] * _ratio)
        if orient != None:
            # rotate XDELTA/YDELTA for computation of rotated crpix position
            _rotm = drutil.buildRotMatrix(_geo_ref.wcslin.orient - orient )
            _rpix = N.dot(N.array(_rpix),_rotm)
            
        meta_wcs.crpix1 = cen[0] + _rpix[0]
        meta_wcs.crpix2 = cen[1] + _rpix[1]
        
        # reset the CRPIX and CRVAL to the image center now.
        meta_wcs.crval1,meta_wcs.crval2 = meta_wcs.xy2rd(cen)
        meta_wcs.crpix1 = cen[0]
        meta_wcs.crpix2 = cen[1]
                
        return meta_wcs

    # This method would use information from the product class and exposure class
    # to build the complete parameter list for running the 'drizzle' task
    def buildPars(self,ref=None,def_wcs=None,tweak=0):
        """This method would build a list of parameters to run 'drizzle'
         one a single input image.
            The reference image info will be passed as a SkyField object.
            The default output reference frame will be passed as 'def_wcs'
            for comparison to the user's selected object 'ref'.
        """
        # pars contains the drizzle parameters for each input(in order):
        #  data,outdata,outnx,outny,scale,xsh,ysh,rot
        parlist = []

        # Extract offsets from parameters dictionary
        # There will only be values when multiple observations are combined.
        print "single Buildpars"
        #pdb.set_trace() #WZ
        offsets = {}
        if self.pars.has_key('xoff') == 1:
            offsets['xsh'] = self.pars['xoff']
        else:
            offsets['xsh'] = 0.0

        if self.pars.has_key('yoff') == 1:
            offsets['ysh'] = self.pars['yoff']
        else:
            offsets['ysh'] = 0.0

        # Setup the parameters for the output field
        # User-specified parameters will come from SkyField object
        # Default values will be set up as a SkyField object
        # 
        # _delta contains the offset from the default product's
        # CRVAL position.
        if ref != None:
            # Extract the total exptime for this output object
            if ref.exptime == None:
                _texptime = self.exptime
            else:
                _texptime = ref.exptime
                
            # Now define the default WCS for this product
            if def_wcs == None:
                #pdb.set_trace()
                def_wcs = sef.ref # self.product.geometry.wcs # this is a key: # self.ref
            #
            # Initialize the SkyField Object with the 
            # user settings.
            _field = ref

            # Check to make sure we have a complete WCS 
            # if not, fill in using the product's default WCS
            _field.mergeWCS(def_wcs)
            _field.wcs.rootname=def_wcs.rootname
            #pdb.set_trace() #Z
            ref_wcs = _field.wcs
            
            # Account for differences in input size relative to 
            # output user-specified reference frame.
            _dcen = ((def_wcs.naxis1 - ref_wcs.naxis1)/2.,(def_wcs.naxis2 - ref_wcs.naxis2)/2.)       

            # Compute delta crvals now...
            # We are only interested in any deltas induced by
            # the user in the input SkyField object, so we have to compare
            # the user's WCS (from ref) to the default situation.
            ref_crv = (ref_wcs.crval1,ref_wcs.crval2)
            _wcs_x,_wcs_y = def_wcs.rd2xy(ref_crv)
            _dx,_dy = (ref_wcs.crpix1 - _wcs_x),(ref_wcs.crpix2 - _wcs_y)
            
            _delta_x = _dx + _dcen[0]
            _delta_y = _dy + _dcen[1]            
            print "DELTA: ", _delta_x, _delta_y
            
            # Update output size based on this final shift
            #   without overriding any user-specified shape.
            # When we update the output WCS, we also need to update the
            #   product WCS to reflect the new values as well...
            #_prod_wcs = self.product.geometry.wcs #WZ
            #print "Test point"
            #pdb.set_trace()
            _prod_wcs = self.ref.model.copy() #WZ
            if (ref.psize == None):
                _prod_wcs.pscale = ref.wcs.pscale #WZ for the first loop
                #_prod_wcs.cd11 = -ref.wcs.pscale/3600.
                #_prod_wcs.cd12 = 0.
                #_prod_wcs.cd21 = 0.
                #_prod_wcs.cd22 = ref.wcs.pscale/3600.
            """ WZ commented out
            if _field.shape == None: 
                ref_wcs.naxis1 = ref_wcs.naxis1 + int(abs(_dx * 2)+0.5)
                ref_wcs.naxis2 = ref_wcs.naxis2 + int(abs(_dy * 2)+0.5)            
                _prod_wcs.naxis1 = _prod_wcs.naxis1 + int(abs(_dx * 2)+0.5)
                _prod_wcs.naxis2 = _prod_wcs.naxis2 + int(abs(_dy * 2)+0.5)
                _prod_wcs.crpix1 = _prod_wcs.crpix1 + _dx
                _prod_wcs.crpix2 = _prod_wcs.crpix2 + _dy
            else:
                _prod_wcs.naxis1 = ref_wcs.naxis1
                _prod_wcs.naxis2 = ref_wcs.naxis2
                _prod_wcs.crpix1 = ref_wcs.crpix1
                _prod_wcs.crpix2 = ref_wcs.crpix2
            """                
        else:
            _delta_x = _delta_y = 0.
            # Define a default WCS
            ref_wcs = self.product.geometry.wcs
            # Extract the total exposure time
            _texptime = self.product.exptime           

                    
        # Build up name for output weighting and context images
        #_outdata = self.product.name
        #indx = string.find(_outdata,'.')
        #if indx > 0:
        #    outweight = _outdata[:indx]+'_weight'+_outdata[indx:]
        #    outcontext = _outdata[:indx]+'_context'+_outdata[indx:]
        #else:
        #    outweight = _outdata+'_weight.fits'
        #    outcontext = _outdata+'_context.fits'

        i=0
        # pdb.set_trace()
        if (tweak['val']>0):
            tweak_data=tableio.get_str(tweak['name'],(0,1,2,3))            
        for member in self.members:
            in_wcs = member.geometry.wcslin
            in_refpix = member.geometry.model.refpix
            #x0 = self.ref.prod.crpix1 # supposed to be the image center
            #y0 = self.ref.prod.crpix2 # Better use _prod_wcs
            x0 = _prod_wcs.crpix1 # supposed to be the image center
            y0 = _prod_wcs.crpix2 
            parameters = {}
            parameters['data'] = member.name
            wcs = ref_wcs.copy()
            ## pdb.set_trace()
            fitsfile = pyfits.open(string.split(member.name,'[')[0],"update") #WZ
            wcs.naxis1 = fitsfile[1].header.get('NAXIS1')
            wcs.naxis2 = fitsfile[1].header.get('NAXIS2')
            wcs.cd11 = fitsfile[1].header.get('CD1_1')
            wcs.cd12 = fitsfile[1].header.get('CD1_2')
            wcs.cd21 = fitsfile[1].header.get('CD2_1')
            wcs.cd22 = fitsfile[1].header.get('CD2_2')
            wcs.crpix1 = fitsfile[1].header.get('CRPIX1')
            wcs.crpix2 = fitsfile[1].header.get('CRPIX2')
            wcs.crval1 = fitsfile[1].header.get('CRVAL1')
            wcs.crval2 = fitsfile[1].header.get('CRVAL2')
            wcs.filter = fitsfile[0].header.get('FILTER')
            """
            wcs.filter1 = fitsfile[0].header.get('FILTER1')
            wcs.filter2 = fitsfile[0].header.get('FILTER2')
            if (wcs.filter1 != 'None'):
                member.filter = string.lower(wcs.filter1)
            else:
                member.filter = string.lower(wcs.filter2)
            """
            xc = wcs.naxis1 / 2 # + 1  + in_refpix['XDELTA']# old frame
            yc = wcs.naxis2 / 2 # + 1  + in_refpix['YDELTA']
            fitsfile.close()
            #print "Old center: ",xc,yc
            pscale =  fitsfile[1].header.get('IDCSCALE') #WZ
            #pscale = N.sqrt(N.power(wcs.cd11,2)+N.power(wcs.cd21,2))*3600.
            rac, decc = wcs.xy2rd([xc,yc]) #WZ
            #print rac, decc # This is the WCS for the center of a new image
            #pdb.set_trace()
            #x,y = self.ref.prod.rd2xy([rac,decc]) # on the output grid
            x,y = _prod_wcs.rd2xy([rac,decc]) # on the output grid
            # ratio = _prod_wcs.pscale / self.ref.pscale # WZ out/inp
            # ratio = ref.psize / self.ref.pscale # WZ out/inp
            ratio = ref_wcs.pscale / member.geometry.model.pscale
            #x = x * ratio
            #y = y * ratio # input                        
            x = x + in_refpix['XDELTA'] # WZ old frame
            # actually it lands with an offset because of distortion
            y = y + in_refpix['YDELTA'] 
            x = x - x0 # with respect to the new center
            y = y - y0
            x1,y1 =self.xyrot([x,y],in_wcs.orient) # after rotation
            #print x1,y1," on the common grid"
            # offsets are important as they reflect matching results
            parameters['xsh']=  offsets['xsh'] + _delta_x # x1 should not apply twice
            parameters['ysh']=  offsets['ysh'] + _delta_y # y1  
            parameters['rot'] = in_wcs.orient # - ref_wcs.orient #WZ
            print "Rot: ",parameters['rot']," dX: ",parameters['xsh']," dY: ",parameters['ysh']
            #pdb.set_trace()

            # Check to see if a mask was created at all...
            # if not, set it to ''
            #pdb.set_trace()
            if member.maskname == None:
                _maskname = ''
            else:
                _maskname = member.maskname
                if (os.path.exists(member.whtname)):
                    #member.whtname = string.replace(_maskname,'inmask','wht'); #WZ
                    parameters['mask'] = member.maskname
                    parameters['in_mask'] = member.whtname
                #if member.whtname ==None:
                else:
                    parameters['in_mask'] = _maskname
                    parameters['mask'] = _maskname

            # Setup parameters for special cases here...
            parameters['outsingle'] = self.outsingle
            parameters['outsweight'] = self.outsweight
            parameters['outscontext'] = self.outscontext
            parameters['outblot'] = member.outblot
            parameters['blotnx'] = member.naxis1
            parameters['blotny'] = member.naxis2
            
            #Setup parameters for normal operations
            parameters['outdata'] = self.outdata
            parameters['outweight'] = self.outweight
            parameters['outcontext'] = self.outcontext
            parameters['outnx'] = ref_wcs.naxis1
            parameters['outny'] = ref_wcs.naxis2

            dx=dy=dr=0.
            if (tweak['val']>0): #WZ
                for i in range(len(tweak_data[0])):
                    rootname = string.split(tweak_data[0][i],'_ini')[0]
                    if (string.find(member.name,rootname) > -1):
                        dx = float(tweak_data[1][i]) # in output grid
                        dy = float(tweak_data[2][i])
                        dr = float(tweak_data[3][i])
                        print "Round ",tweak['val'],": Shifts ",dx,dy,"; Rotation ", \
                              dr,"; ",member.name 
            # Calculate any rotation relative to the orientation
            # AFTER applying ONLY the distortion coefficients without
            # applying any additional rotation 
            # to make PA=0.0 ... #WZ
            parameters['rot']=in_wcs.orient + dr*180./N.pi #WZ
            print "Rot: ",parameters['rot']," dX: ",parameters['xsh']," dY: ", \
              parameters['ysh']
            #pdb.set_trace()                                                   
            # parameters['rot'] = in_wcs.orient # - ref_wcs.orient #WZ
            if i == 0:
                ref.wcs.orient = in_wcs.orient #WZ
            parameters['texptime'] = _texptime[0]
            parameters['expstart'] = _texptime[1]
            parameters['expend'] = _texptime[2]
            # Need to revise how this gets computed...
            # The pixel scale of the product corresponds to the 
            # desired output pixel scale, and the model pscale for
            # the member represents the un-scaled pixel size for the input.
            # parameters['xsh'] = offsets['xsh'] + _delta_x + dx
            # parameters['ysh'] = offsets['ysh'] + _delta_y + dy
            parameters['scale'] = ref_wcs.pscale / member.geometry.model.pscale
            # pdb.set_trace() #WZ
            
            # Parameters only used by 'wdrizzle'
            parameters['geomode'] = 'user'

            # Set up the idcfile for use by 'drizzle'
            indx = string.find(member.name,'.')
            coeffs = member.name[:indx]+'_coeffs'+member.chip+'.dat'
            member.geometry.model.convert(coeffs)
            parameters['coeffs'] = coeffs

            # Now pass along the remainder of the user specified parameters
            if self.pars['units'] != None:
                parameters['units'] = self.pars['units']
            else:
                parameters['units'] = 'cps'

            if self.pars['pixfrac'] != None:	
                parameters['pixfrac'] = self.pars['pixfrac']
            else:
                parameters['pixfrac'] = 1.0

            if self.pars['kernel'] != None:
                parameters['kernel'] = self.pars['kernel']
            else:
                parameters['kernel'] = 'square'
            i = i + 1
            parlist.append(parameters)

        # Now, combine them for a complete set of pars for all exposures
        # in this pattern/observation.
        #
        #print "Finish single builPars"
        #pdb.set_trace()
        return parlist

    def computeCubicCoeffs(self):
        """
        Method for converting cubic and Trauger coefficients tables
        into a usable form.  It also replaces 'computeOffsets' for
        those tables as well.
        """
        # For each chip in the observation...
        _pscale1 = None
        for img in self.members:
            _chip = img.chip
            # scale all chips to first chip plate scale...
            if _pscale1 == None or img.chip == '1':
                _pscale1 = self.REFDATA[_chip]['psize']
                _reftheta = self.REFDATA[_chip]['theta']
        
        for img in self.members:
            # ... get the model and type of coefficients table used...
            _model = img.geometry.model
            _ikey = img.geometry.ikey
            _chip = img.chip
            _refdata = self.REFDATA[_chip]

            # ... determine the plate scale and scaling factor between chips...
            if img.chip == '1':
                _pscale = _refdata['psize']
                _ratio = 1.0
            else:
                _pscale = _refdata['psize']
                _ratio = _refdata['psize'] / _pscale1
            # Record the plate scale for each chip's model that was used
            # to compute the coefficients, not the plate scale from the
            # image header.
            _model.pscale = _pscale
            _model.refpix['PSCALE'] = _pscale

            # Correct the coefficients for the differences in plate scales
            _model.cx = _model.cx * N.array([_model.pscale/_ratio])
            _model.cy = _model.cy * N.array([_model.pscale/_ratio])
            _model.refpix['XREF'] = self.REFPIX['x']
            _model.refpix['YREF'] = self.REFPIX['y']

            # Correct the offsets for the plate scales as well...
            _model.refpix['XDELTA'] = _model.refpix['V2REF'] / _ratio
            _model.refpix['YDELTA'] = _model.refpix['V3REF'] / _ratio
            _model.refpix['centered'] = yes
            

    def computeOffsets(self,parity=None):
        """
        This version of 'computeOffsets' calculates the zero-point
        shifts to be included in the distortion coefficients table
        used by 'drizzle'. 
        It REQUIRES a parity matrix to convert from 
        V2/V3 coordinates into detector image X/Y coordinates. This
        matrix will be specific to each detector.
        """
        vref = []

        # print "Check parity"
        #pdb.set_trace() #WZ
        # Set up the parity matrix here for a SINGLE chip
        if parity == None:
            # Use class defined dictionary as default
            parity = self.PARITY
        
        nmembers = len(self.members)

        v2sum,v3sum = 0.,0.
        for member in self.members:
            in_model = member.geometry.model
            pscale = in_model.pscale
            
            detector = drutil.getPrimaryKeyword(member.name,'DETECTOR')
            pmat = parity[detector]
            
            v2 = in_model.refpix['V2REF']
            v3 = in_model.refpix['V3REF']

            if (detector == 'UVIS'):
                deg = -135 # UVIS specific
                angle = (deg + 90.0) * 3.141519 / 180.0
                (v2, v3) = (v2 * N.cos(angle) + v3 * N.sin(angle),
                            -v2 * N.sin(angle) + v3 * N.cos(angle))
                ## if member.chip == '2':  #AKS - Default:1
                                                
            #vpos = N.matrixmultiply((v2,v3),pmat) / pscale
            vpos = N.dot((v2,v3),pmat) / pscale
            vref.append(vpos)
            v2sum = v2sum + vpos[0]
            v3sum = v3sum + vpos[1]
            
        v2com = v2sum / nmembers
        v3com = v3sum / nmembers
               
        # Compute position of each chip's common point relative
        # to the output chip's reference position.
        i = 0
        #pdb.set_trace()
        for chip in self.members:
            # Update member's geometry model with computed
            # reference position...
            _refp = chip.geometry.model.refpix
            _refp['XDELTA'] = vref[i][0] - v2com 
            _refp['YDELTA'] = vref[i][1] - v3com 
            _refp['centered'] = no
            i = i + 1

    def setNames(self,filename,output):
        """ 
        Define standard name attibutes:
                ext0name    - Primary extension of input file
                outname     - Default final output name
                outdata     - Name for drizzle science output 
                outsingle   - Name for output for single image
        """
        self.ext0name = filename + self.EXT0
        self.outname = output
        
        # Define output file names for final case
        self.outdata = drutil.buildNewRootname(output,extn='.fits') #WZ
        self.outweight = drutil.buildNewRootname(output,extn='_weight.fits')
        self.outcontext = drutil.buildNewRootname(output,extn='_context.fits')
        #self.outdata = drutil.buildNewRootname(output,extn='_sci.fits')
        #self.outweight = drutil.buildNewRootname(output,extn='_sci_weight.fits')
        #self.outcontext = drutil.buildNewRootname(output,extn='_sci_context.fits')
        
        # Define output file names for separate output for each input
        self.outsingle = drutil.buildNewRootname(filename,extn='_single.fits') #WZ
        self.outsweight = drutil.buildNewRootname(filename,extn='_single_weight.fits')
        self.outscontext = drutil.buildNewRootname(filename,extn='_single_context.fits')
        #self.outsingle = drutil.buildNewRootname(filename,extn='_single_sci.fits')
        #self.outsweight = drutil.buildNewRootname(filename,extn='_single_sci_weight.fits')
        #self.outscontext = drutil.buildNewRootname(filename,extn='_single_sci_context.fits')
        
    ########
    #
    # User interface methods
    #
    ########
    def getWCS(self):
        return self.members[0].getWCS()
    
    def getExptime(self,filename):
        _extname = filename+self.EXT0
        _exptime = float(drutil.getPrimaryKeyword(_extname,'EXPTIME'))
        
        _expstart = float(drutil.getPrimaryKeyword(_extname,'EXPSTART'))
        _expend = float(drutil.getPrimaryKeyword(_extname,'EXPEND'))
         
        return (_exptime,_expstart,_expend)
        
    def DeltaXYtoOffset(self,delta):
        """
        Converts provided delta(x,y) pixel offset into a 
        delta(RA,Dec) offset in arcseconds. 
        """
        _wcs = self.product.getWCS()
        _geom = self.product.geometry
        
        new_rd = _geom.XYtoSky((_wcs.crpix1 - delta[0],_wcs.crpix2 - delta[1]))
        delta_ra = (_wcs.crval1 - new_rd[0]) * 3600.
        delta_dec = (_wcs.crval2 - new_rd[1]) * 3600.
        
        return (delta_ra,delta_dec)
        
    def xyrot(self,xypos,theta=0):   # WZ, migrated from combDither.py
        """New position [Xp,Yp] after a rotation. Theta in units of deg
        """
        rad = theta * N.pi / 180.
        x = xypos[0]
        y = xypos[1]
        xp = x * N.cos(rad) - y * N.sin(rad)
        yp = x * N.sin(rad) + y * N.cos(rad)    
        return xp,yp
        

class WFC3Observation(Pattern):
    """This class defines an observation with information specific
       to WFC3 exposures, including knowledge of how to mosaic both
       chips."""

    # Define a class variable for the gap between the chips
    PARITY = {'UVIS':[[-1.0,0.0],[0.0,1.0]],'IR':[[-1.0,0.0],[0.0,1.0]]}

    def __init__(self, filename, output, ref, pars=None):      #WZ

        # Now, initialize Pattern with all member Exposures...
        Pattern.__init__(self, pars=pars)
        
        # build required name attributes:
        #   ext0name, outname, output, outsingle
        self.setNames(filename,output)

        # Determine how many 'chips' make up the observation	
        self.nmembers = int(drutil.findNumExt(self.ext0name)) / self.NUM_IMSET
        
        # Set EXPTIME for exposure
        self.exptime = self.getExptime(filename)
        
        # Build up list of chips in observation
        self.addMembers(filename)
                
        # Only need to worry about IDC tables for WFC3...
        self.computeOffsets()

        # Set up the input members and create the product meta-chip
        self.buildProduct(filename, output, ref)

class ACSObservation(Pattern):
    """This class defines an observation with information specific
       to ACS WFC exposures, including knowledge of how to mosaic both
       chips."""

    # Define a class variable for the gap between the chips
    PARITY = {'WFC':[[1.0,0.0],[0.0,-1.0]],'HRC':[[-1.0,0.0],[0.0,1.0]],'SBC':[[-1.0,0.0],[0.0,1.0]]}

    def __init__(self, filename, output, ref, pars=None): #WZ

        # Now, initialize Pattern with all member Exposures...
        Pattern.__init__(self, pars=pars)
        
        # build required name attributes:
        #   ext0name, outname, output, outsingle
        self.setNames(filename,output)

        # Determine how many 'chips' make up the observation	
        self.nmembers = int(drutil.findNumExt(self.ext0name)) / self.NUM_IMSET
        
        # Set EXPTIME for exposure
        self.exptime = self.getExptime(filename)
        
        # Build up list of chips in observation
        self.addMembers(filename)
                
        #print "*** Check"
        #pdb.set_trace()
        # Only need to worry about IDC tables for ACS...
        self.computeOffsets()

        # Set up the input members and create the product meta-chip
        self.buildProduct(filename, output, ref) #WZ


class STISObservation(Pattern):
    """This class defines an observation with information specific
       to STIS exposures.
    """

    # Default coefficients table to use for this instrument
    IDCKEY = 'cubic'
    # This parity is the multiplication of PC1 rotation matrix with
    # a flip in X for output image.
    __theta = 135.0
    __parity = drutil.buildRotMatrix(__theta) * N.array([[-1.,1.],[-1.,1.]])
    PARITY = {'CCD':__parity,'NUV_MAMA':__parity,'FUV_MAMA':__parity}

    # The dictionaries 'REFDATA' and 'REFPIX' are required for use with
    # cubic and Trauger coefficients tables in 'computeCubicCoeffs'.
    #
    # This information provides the absolute relationship between the chips
    # Latest plate scales: 0.05071, 0.0246
    REFDATA = {'CCD':{'psize':0.05,'xoff':0.0,'yoff':0.0,'theta':__theta},
              'NUV-MAMA':{'psize':0.024,'xoff':0.0,'yoff':0.0,'theta':__theta},
              'FUV-MAMA':{'psize':0.024,'xoff':0.0,'yoff':0.0,'theta':__theta}}
    REFPIX = {'x':512.,'y':512.}
                  
    def __init__(self, filename, output, pars=None):
       
        # Now initialize Pattern with all member exposures...
        Pattern.__init__(self,pars=pars)

        # build output rootnames here...
        self.setNames(filename,output)

        # Set EXPTIME for exposure
        self.exptime = self.getExptime(filename)

        self.nmembers = int(drutil.findNumExt(self.ext0name)) / self.NUM_IMSET
        
        # Now, build list of members and initialize them
        self.addMembers(filename)
        if self.members[0].chip == '1': self.members[0].chip = 'CCD'
        
        if self.members[0].geometry.ikey != 'idctab':
            # Correct distortion coefficients to match output pixel scale
            self.computeCubicCoeffs()
        else:
            self.computeOffsets()
            
        # Set up the input members and create the product meta-chip
        self.buildProduct(filename, output)

class WFPCObservation(Pattern):
    """This class defines an observation with information specific
       to WFPC2 exposures, including knowledge of how to mosaic the
       chips."""

    # Default coefficients table to use for this instrument
    IDCKEY = 'trauger'
    # This parity is the multiplication of PC1 rotation matrix with
    # a flip in X for output image.
    __theta = 224.908
    __parity = drutil.buildRotMatrix(__theta) * N.array([[-1.,1.],[-1.,1.]])
    PARITY = {'1':__parity,'2':__parity,'3':__parity,'4':__parity}
    
    # The dictionaries 'REFDATA' and 'REFPIX' are required for use with
    # cubic and Trauger coefficients tables in 'computeCubicCoeffs'.
    #
    # This information provides the absolute relationship between the chips
    #REFDATA ={'1':{'psize':0.04554,'xoff':354.356,'yoff':343.646,'theta':__theta},
    #          '2':{'psize':0.0996,'xoff':-371.27125,'yoff':350.50803,'theta':314.388},
    #          '3':{'psize':0.0996,'xoff':-369.014827,'yoff':-352.74708,'theta':44.698},
    #          '4':{'psize':0.0996,'xoff':353.39406,'yoff':-354.18689,'theta':135.258}}
    REFDATA ={'1':{'psize':0.04554,'xoff':354.356,'yoff':343.646,'theta':__theta},
              '2':{'psize':0.0996,'xoff':345.7481,'yoff':375.28818,'theta':__theta},
              '3':{'psize':0.0996,'xoff':366.56876,'yoff':354.79435,'theta':__theta},
              '4':{'psize':0.0996,'xoff':355.85016,'yoff':351.29183,'theta':__theta}}
    REFPIX = {'x':400.,'y':400.}
               
    def __init__(self, filename, output, pars=None):
       
        # Now initialize Pattern with all member exposures...
        Pattern.__init__(self,pars=pars)

        # Set this to None since WFPC2 does not have DQ arrays and
        # we will rely on default, uniform weighting for now...
        pars['bitvalue'] = None
        gcount = None

        # Determine syntax for accessing this observation based
        # on whether the input WFPC2 image is in GEIS or
        # multi-extension format
        # For GEIS image, we need to set values for 
        # 'ext0' and self.SCI_EXTN.
        if string.find(filename,'.fits') > 0:
            if drutil.getPrimaryKeyword(filename+self.EXT0,'EXTEND') == 'yes':
                # We have a 'waivered' FITS file...
                self.EXT0 = '[0]'
                self.SCI_EXTN = '[0][*,*,1]'
                self.SEPERATOR = '[0][*,*,'
                #gcount = drutil.getPrimaryKeyword(filename+self.EXT0,'SDASMGNU')                
                #print 'Trying to work with waivered WFPC2 FITS image...'
                
                _Error_msg = 'Can NOT use '+filename+'!'
                _Error_msg = _Error_msg + ' Not a GEIS or multi-extension FITS file!'
                print _Error_msg           
                raise TypeError            
        else:
            # We are working with a GEIS WFPC2 image
            gcount = drutil.getPrimaryKeyword(filename,'GCOUNT')
            self.EXT0 = '[1]'
            self.SCI_EXTN = '[1/'+gcount+']'
            self.SEPERATOR = '['

        # build output rootnames here...
        self.setNames(filename,output)

        # Determine how many 'chips' make up the observation	
        if gcount == None:
            self.nmembers = int(drutil.findNumExt(self.ext0name)) / self.NUM_IMSET
        else:
            self.nmembers = int(gcount)
        
        # Set EXPTIME for exposure
        self.exptime = self.getExptime(filename)

        # Now, build list of members and initialize them
        self.addMembers(filename)
        
        if self.members[0].geometry.ikey != 'idctab':
            # Correct distortion coefficients to match output pixel scale
            self.computeCubicCoeffs()
        else:
            self.computeOffsets()
            
        # Determine desired orientation of product
        self.setOrient()

        # Set up the input members and create the product meta-chip
        self.buildProduct(filename, output)
            
    def addMembers(self,filename):
        # Build rootname here for each SCI extension...
        #pdb.set_trace()
        for i in range(self.nmembers):
            extname = filename+self._setSciExtn(extn=i+1)
            _detnum = drutil.getPrimaryKeyword(extname,'DETECTOR')
            # Also, build mask file for this member chip
            maskname = buildmask.buildShadowMaskImage(filename,_detnum)
            self.members.append(Exposure(extname, idckey=self.idckey, mask=maskname,
                                nchips=self.nmembers))

    def setOrient(self):
        
        # Determine desired orientation of product
        meta_orient = None
        for exp in self.members:
            if int(exp.chip) == 1:
                meta_orient = exp.geometry.wcs.orient
        
        if meta_orient == None:
            meta_orient = self.members[0].geometry.wcs.orient
        
        # Set orient for all groups
        # Dither coefficients rotate all chips to chip 1 orientation
        # if not using Trauger coefficients
        for exp in self.members:
            exp.geometry.wcs.orient = meta_orient 
        
               
class DitherProduct(Pattern):
    """
    Builds an object for a set of dithered inputs, each of which
    will be one of the Observation objects.  
    """
    def __init__(self, prodlist,pars=None,ref=None): #WZ

        # Setup a default exposure to contain the results
        Pattern.__init__(self,pars=pars)

        print "Add product members"
        # Build temporary output drizzle product name
        #pdb.set_trace()
        output = drutil.buildRootname(string.lower(prodlist['output']),ext=['_drz.fits'])
        self.pars = prodlist['members']
        self.nmembers = len(prodlist['members'])
        self.ref = ref #WZ
                        
        self.addMembers(prodlist,pars,output)
        
        # Update member WCS info with any offsets given in ASN table
        self.applyAsnShifts()
        
        self.exptime = self.getExptime()
        
        _psize = pars['psize']
        output_wcs = self.buildMetachip(psize=_psize)

        self.size = (output_wcs.naxis1,output_wcs.naxis2)
        
        self.product = Exposure(output,wcs=output_wcs,new=yes)
        self.product.exptime = self.exptime
        
        # Compute the offsets relative to the product WCS
        self.computeOffsets()
    
    def getExptime(self):
        """
        Add up the exposure time for all the members in
        the pattern, since 'drizzle' doesn't have the necessary
        information to correctly set this itself.
        """
        _exptime = 0.
        _start = []
        _end = []
        for member in self.members:
            _exptime = _exptime + member.exptime[0]
            _start.append(member.exptime[1])
            _end.append(member.exptime[2])
        _expstart = min(_start)
        _expend = max(_end)
        
        return (_exptime,_expstart,_expend)
    
    def addMembers(self,prodlist,pars,output):
        """
        For each entry in prodlist, append the appropriate type
        of Observation to the members list.
        """
        print "prod image. Update ref"
        # pdb.set_trace()
        #pd = self.ref.prod
        #md = self.ref.model
        im = prodlist['members'].keys()[0]+'_flt.fits[1]'
        tem_wcs = wcsutil.WCSObject(im)  
        tem_wcs.prod = self.ref.prod
        tem_wcs.model = self.ref.model
        self.ref = tem_wcs.copy()
        del tem_wcs
        for memname in prodlist['members'].keys():
            pardict = self.pars[memname]
            pardict.update(pars)
            filename = drutil.buildRootname(memname)
            ref = self.ref #WZ
            self.members.append(selectInstrument(filename,output,ref,pars=pardict))
        
    def buildPars(self,ref=None, tweak=None): #WZ: This is an overall call

        if ref == None:
            #_field = None
            _field = SkyField()
        else:
            _field = ref

        # Remember the default situation for comparison with the
        # user-specified field given by 'ref'.
        print "Overall buildpars" # " Choose Default reference"
        #def_wcs = self.product.geometry.wcs
        def_wcs = self.ref #WZ Another key
        
        _field.mergeWCS(self.product.geometry.wcs)
        _field.exptime = self.exptime
                   
        parlist = []
        for member in self.members:
            parlist = parlist + member.buildPars(ref=_field,def_wcs=def_wcs, tweak=tweak)
        return parlist

    def buildMetachip(self,psize=None):
        """
        This method combines the results of the member's buildMetachip()
        methods into a composite WCS.
        (DitherProduct method)
        """
        prodlist = []

        _ref = self.members[0].product.geometry
        
        # Get pscale from model
        if psize == None:
            pscale = _ref.model.pscale
        else: 
            # User specified a different pixel scale to use...
            pscale = psize

        print "Making prodlist"
        #pdb.set_trace()
        for member in self.members:
            # merge the corrected shapes into a corrected meta-chip here
            # Start by computing the corected positions for reference points
            prodlist.append(member.product)
           
        # Use first input image WCS as initial reference WCS
        meta_wcs = _ref.wcs.copy()
        ref_crval = (meta_wcs.crval1,meta_wcs.crval2)
        ref_crpix = (meta_wcs.crpix1,meta_wcs.crpix2)
        
        # Specify the output orientation angle
        ref_orient = meta_wcs.orient
        
        _xr,_yr = [],[]
        # Compute offsets for each input relative to this reference WCS
        for member in prodlist:
            _wcs = member.geometry.wcs
            
            _crval = (_wcs.crval1,_wcs.crval2)
            # Determine where member CRVAL position falls in ref frame            
            npos = meta_wcs.rd2xy(_crval)
            
            # Shift of reference pixel position from member frame
            # into reference frame
            _dpos = (npos[0] - _wcs.crpix1, npos[1] - _wcs.crpix2)
            
            # Compute the range of pixels each input spans in the 
            # output meta-frame coordinates
            # apply
            #pdb.set_trace()
            _scale = meta_wcs.pscale / _wcs.pscale
            _range = (_wcs.naxis1 * _scale, _wcs.naxis2 * _scale)
                        
            _xr.append(_dpos[0])
            _xr.append(_dpos[0] + _range[0])
            _yr.append(_dpos[1])
            _yr.append(_dpos[1] + _range[1])        
                
        # Determine the full size of the metachip
        _xmin = N.minimum.reduce(_xr)
        _ymin = N.minimum.reduce(_yr)	
        _xmax = N.maximum.reduce(_xr)
        _ymax = N.maximum.reduce(_yr)

        if _xmin < 0.: _dx = _xmin
        else: _dx = 0.
        if _ymin < 0.: _dy = _ymin
        else: _dy = 0.

        xsize = int(_xmax+1.0) - int(_xmin-1.0)
        ysize = int(_ymax+1.0) - int(_ymin-1.0)

        # Determine offset of centers from center of image.
        _nref = (xsize/2. - (_xmax + _xmin )/2., ysize/2. - (_ymax + _ymin )/2.)
        
        # Update reference WCS for new size
        meta_wcs.naxis1 =  xsize
        meta_wcs.naxis2 =  ysize
        
        # Have to adjust the CRPIX by how much the center needs to shift
        # to fit into the reference frame.
        meta_wcs.crpix1 = meta_wcs.crpix1 + _nref[0]
        meta_wcs.crpix2 = meta_wcs.crpix2 + _nref[1]
        
        meta_wcs.crval1,meta_wcs.crval2 = meta_wcs.xy2rd((xsize/2.,ysize/2.))
        meta_wcs.crpix1 = xsize/2.
        meta_wcs.crpix2 = ysize/2.
        
        return meta_wcs

    def applyAsnShifts(self):
        """ This method updates each member's WCS to reflect any
            offsets/corrections specified in the ASN table.
        """
        print "Apply shifts"
        #pdb.set_trace()
        for member in self.members:
            in_wcs = member.product.geometry.wcs 
            in_wcslin = member.product.geometry.wcslin           
            # Match member name with name entry in dictionary
            for name in self.pars.keys():
                if string.find(member.ext0name,string.lower(name)) > -1:
                    memname = name
                    continue
            pars = self.pars[memname]
            
            # Check to see if there are any offsets given for this member...
            if pars['delta_x'] != 0. or pars['delta_y'] != 0.:
                _crval = (in_wcs.crval1,in_wcs.crval2)
                _crpix = (in_wcs.crpix1,in_wcs.crpix2)
                # convert it to units of degrees
                delta_ra = pars['delta_x'] * N.cos(_crval[1]) / 3600.
                delta_dec = pars['delta_y'] / 3600.
                # compute new CRVALs
                new_rd = (_crval[0]+delta_ra, _crval[1]+delta_dec)
                if pars['delta_rot'] != 0.:
                    new_orient = in_wcs.orient + pars['delta_rot']
                else:
                    new_orient = None
                # update WCS with new values
                in_wcs.updateWCS(orient=new_orient,refval=new_rd,refpos=_crpix)
                in_wcslin.updateWCS(orient=new_orient,refval=new_rd,refpos=_crpix)

    def computeOffsets(self):
        """
        This method will rely on final product's 'rd2xy' method
        to compute offsets between the different input chips.   key
        """
        ref_wcs = self.product.geometry.wcs.copy()
        print "Calculate offsets"
        #pdb.set_trace()
        ratio = self.ref.model.pscale/ref_wcs.pscale #WZ in native pixel scale
        if (ratio < 0.99 or ratio > 1.01):
            ref_wcs.cd11 = ref_wcs.cd11 * ratio
            ref_wcs.cd12 = ref_wcs.cd12 * ratio
            ref_wcs.cd21 = ref_wcs.cd21 * ratio
            ref_wcs.cd22 = ref_wcs.cd22 * ratio

        #ref_wcs = self.ref.prod #WZ Feb 11:  output frame
        ref_crv= (ref_wcs.crval1,ref_wcs.crval2)

        for member in self.members:
            in_wcs = member.product.geometry.wcs            
            # Match member name with name entry in dictionary
            for name in self.pars.keys():
                if string.find(member.ext0name,string.lower(name)) > -1:
                    memname = name
                    continue
                 
            #xc = in_wcs.naxis1/2.
            #yc = in_wcs.naxis2/2.
            #rac,decc=in_wcs.xy2rd([xc,yc])
            #x,y = ref_wcs.rd2xy([rac,decc])
            x,y = ref_wcs.rd2xy((in_wcs.crval1,in_wcs.crval2))
            
            xr,yr = in_wcs.rd2xy(ref_crv) #WZ Useful? Correct?
            xroff = ref_wcs.crpix1 - xr
            yroff = ref_wcs.crpix2 - yr
            
            xoff = x - ref_wcs.crpix1  # WZ in units of reference pixel
            yoff = y - ref_wcs.crpix2
            
            x1,y1 =self.xyrot([xoff,yoff],-ref_wcs.orient) # after rotation. Note the sign
            self.pars[memname]['xoff'] = x1 # xoff
            self.pars[memname]['yoff'] = y1 # yoff
            self.pars[memname]['rot'] = in_wcs.orient # - ref_wcs.orient #WZ
            #pdb.set_trace() #WZ
            print "Offsets: ",self.pars[memname]['xoff'],self.pars[memname]['yoff']
#
#
def selectInstrument(filename,output,ref=None,pars=None): #WZ
    """
    Method which encapsulates the logic for determining
    which class to instantiate for each file.

    """
    # Determine the instrument...
    instrument = drutil.getPrimaryKeyword(filename+'[0]','INSTRUME')
    print "Select instrument and set parameters"
    #pdb.set_trace()
    try:
        # ... then create an appropriate object.
        if instrument == INSTRUMENT[0]:
            member = ACSObservation(filename,output,ref,pars=pars) #WZ
        elif instrument == INSTRUMENT[1]:
            member = WFPCObservation(filename,output,pars=pars)
        elif instrument == INSTRUMENT[2]:
            member = STISObservation(filename,output,pars=pars)
        elif instrument == INSTRUMENT[3]:
            member = WFC3Observation(filename,output,ref,pars=pars)
        else:
            raise AttributeError, "Instrument '%s' not supported now."%instrument 
    except TypeError:
        return None
        
    return member

class SkyField:
    """
    This class will create a WCS object for a desired image, along 
    with all the necessary parameters for modifying that WCS.
    The user may optionally modify:
        ra,dec      - position of the center of the field (as a string)
        orient      - value of ORIENTAT for the field
        shape       - tuple containing the sizes of the field's x/y axes
        psize       - size of image's pixels in arc-seconds
        
    RA/Dec may be specified in either decimal degrees (eg.,124.5678) 
        or sexagesimal notation (hh:mm:ss.ss).

    """
    def __init__(self,shape=None,psize=None):

        self.shape = shape
        self.ra = None
        self.dec = None
        self.orient = None
        self.psize = psize

        self.crpix = (None,None)

        # Set up proper shape tuple for WCSObject
        if shape != None and psize != None:
            wshape = (shape[0],shape[1],psize)
        else:
            wshape = None
    
        self.wcs = wcsutil.WCSObject("SkyField",shape=wshape)
        
        # Set this to keep track of total exposure time for
        # output frame... Not user settable.
        self.exptime = None
        
    def mergeWCS(self,wcs):
        """ Sets up the WCS for this object based on another WCS.
            This method will NOT update object attributes other
            than WCS, as all other attributes reflect user-settings.           
        """
        if self.ra == None:
            ra = wcs.crval1
            dec = wcs.crval2
        else:
            ra = self.ra
            dec = self.dec
        _crval = (ra,dec)
        
        if self.shape == None:
            shape = (wcs.naxis1,wcs.naxis2)
            _crpix = (wcs.crpix1,wcs.crpix2)          
        else:
            shape = self.shape
            if self.crpix == None:
                _crpix = (self.shape[0]/2.,self.shape[1]/2.)
            else:
                _crpix = self.crpix
         
        if self.psize == None:
            _psize = wcs.pscale
        else:
            _psize = self.psize
        
        if self.orient == None:
            _orient = wcs.orient
        else:
            _orient = self.orient
              
        # Set up the new WCS based on values from old one.
        self.wcs.updateWCS(pixel_scale=_psize,orient=_orient,refpos=_crpix,refval=_crval)
        self.wcs.naxis1 =  shape[0]
        self.wcs.naxis2 =  shape[1]

    def set(self,psize=None,orient=None,ra=None,dec=None,
            shape=None,crpix=None):
        """ 
        Modifies the attributes of the SkyField and 
            updates it's WCS when appropriate.
        """
        # Converts(if necessary), then updates the RA and Dec.
        _ra,_dec = None,None
        if ra != None:
            if string.find(repr(ra),':') > 0: 
                _hms = string.split(repr(ra)[1:-1],':')
                if _hms[0][0] == '-': _sign = -1 
                else: _sign = 1
                
                for i in range(len(_hms)): _hms[i] = float(_hms[i])
                _ra = _sign * (_hms[0] + ((_hms[1] + _hms[2]/60.) / 60.)) * 15.
            else:
                _ra = float(ra)
            self.ra = _ra

        if dec != None:
            if string.find(repr(dec),':') > 0:
                _dms = string.split(repr(dec)[1:-1],':')
                if _dms[0][0] == '-': _sign = -1 
                else: _sign = 1
                
                for i in range(len(_dms)): _dms[i] = float(_dms[i])
                _dec = _sign * (_dms[0] + ((_dms[1] + _dms[2]/60.) / 60.))
            else:
                _dec = float(dec)
            self.dec = _dec

        if self.ra != None and self.dec != None:
            _crval = (self.ra,self.dec)
        else:
            _crval = None
            
        # Updates the shape, and reference position,
        # only if a new value is specified.
        _crpix = None
        if crpix == None:
            if shape != None:
                self.shape = shape
                _crpix = (self.shape[0]/2.,self.shape[1]/2.)
        else:
            _crpix = crpix

        self.crpix=_crpix
        
        if psize != None:
            self.psize = psize
        
        if orient != None:
            self.orient = orient
                    
        # Updates the WCS with all the changes, if there is enough info.
        self.wcs.updateWCS(pixel_scale=psize,orient=orient,refpos=_crpix,
                                refval=_crval,size=self.shape)

    def __str__(self):
        """ Prints the WCS information set for this object.
        """
        if self.psize != None and self.orient != None:
            block =  self.wcs.__str__()
        else:
            block = 'User parameters for SkyField object: \n'
            block = block + '    psize = '+repr(self.psize)+' \n'
            block = block + '   orient = '+repr(self.orient)+' \n'
            block = block + '    shape = '+repr(self.shape)+' \n'
            block = block + '       ra = '+repr(self.ra)+' \n'
            block = block + '      dec = '+repr(self.dec)+' \n'
            block = block + '   No WCS.\n'

        return block

    def help(self):
        """ Creates and prints usage information for this class.
        """
        _str = """
    An class for specifying the parameters and building a WCS object 
        for a user-specified drizzle product.
    The user may optionally modify the values for:
        psize       - size of image's pixels in arc-seconds
        orient      - value of ORIENTAT for the field
        shape       - tuple containing the sizes of the field's x/y axes
        ra,dec      - position of the center of the field
                      decimal (124.5678) or 
                      sexagesimal string _in quotes_ ('hh:mm:ss.ss')
        crpix       - tuple for pixel position of reference point in 
                        output image  
    Usage: 
        To specify a new field with a fixed output size of 1024x1024:
        --> field = xydrizzle.SkyField(shape=(1024,1024))
        
        The 'set()' method modifies one of the parameters listed above 
        without affecting the remainder of the parameters.
        --> field.set(psize=0.1,orient=0.0)        
        --> field.set(ra=123.45678,dec=0.1000,crpix=(521,576))
        
        View the WCS or user-specified values for this object:
        --> print field        
        """
        print _str

class PyDrizzle:		
    """
Program to process and/or dither-combine image(s) using (t)drizzle.
To create an object named 'test' that corresponds to a drizzle product:
    --> test = xydrizzle.PyDrizzle(input)
where input is the FULL filename of an ACS observation or ASN table.
This computes all the parameters necessary for running drizzle on all
the input images.  Once this object is created, you can run drizzle using:
    --> test.run()
Optional parameters:
    output      User-specified name for output products
    field       User-specified parameters for output image
                includes: psize, orient, ra, dec, shape
    units       Units for final product: 'counts' or 'cps'(DEFAULT)
    kernel      Specify which kernel to use in TDRIZZLE
                'square'(default),'point','gaussian','turbo','tophat'
    pixfrac     drizzle pixfrac value (Default: 1.0) 
    bits        Bit values considered good (Default: 0)
                'None' specifies that no input mask will be used
    idckey      User-specified keyword for determining IDCTAB filename
                'IDCTAB'(ACS default),'TRAUGER'(WFPC2),'CUBIC'(WFPC2)
    clean       Remove temp files? (Default: yes)              
Usage of optional parameters:    
    --> test = xydrizzle.PyDrizzle('test_asn.fits',bits=1894,units='counts')
To keep the individual 'drizzle' output products:
    --> test.run(save=yes)
    """
    def __init__(self, input, output=None, field=None, units=None,
                 kernel=None,pixfrac=None,bits=0,idckey=None,clean=yes,
                 filter1=None,filter2=None, ref=None, tweak=0): #WZ
                 # filter=None,ref=None, tweak=0): #WZ
        print 'Starting PyDrizzle Version ',version,' at ', _ptime()
        # GLOBALfilts[0] = filter
        GLOBALfilts[0] = filter1
        GLOBALfilts[1] = filter2
        #pdb.set_trace()
        
        # remember whether to clean out/delete temp files or not
        self.clean = clean
                
        # Do a quick sanity check on the input filename
        # Can it be found in the current directory?
        if fileutil.findFile(input) == no:
            raise IOError,'Can not find input file in current directory!'

        # Does it have a recognizable extension
        if string.find(input,'.') < 0:
            raise ValueError,"Please specify extension (i.e., .fits) for input '%s' "%input
        # These parameters are needed for buildPars()
        self.input = input
        
        # Extract user-specified parameters, if any have been set...
        # 'field' will be a SkyField object...
        if field != None:
            psize = field.psize
            orient = field.orient
        else:
            psize = None
            orient = None
        
        # These can also be set by the user.
        self.pars = {'psize':psize,'units':units,'kernel':kernel,'rot':orient,
            'pixfrac':pixfrac,'bitvalue':bits,'idckey':idckey}

        self.ref = ref #WZ
        #print "psize"
        #pdb.set_trace()
        #if field != None:
        #    self.ref.pscale=self.pars.psize
        
        # Check to see if user-supplied output name is complete
        # Append .FITS suffix to output name if necessary
        self.output = output
        if output != None:
            _indx = string.find(output,'.')
            if _indx < 0: output = output + '.fits'
        
        # Decipher input to determine whether we are working with
        # an ASN table or single image...
        #
        # Should this check for ASN table and 
        if string.find(input,'_asn') < 0:
            # Start with single image case...
            # Check to see if an output name was provided.
            if output == None:
                # We need to build a default output name...
                output = drutil.buildNewRootname(input,extn='_drz.fits')
                print 'Setting up default output name: ',output
                if fileutil.findFile(output): 
                    print 'Removing previous output product...'
                    os.remove(output)
                
            self.observation = selectInstrument(input,output,ref,pars=self.pars) #WZ
        else:
            # We are dealing with an ASN table...
            # 'input' - full filename for ASN table
            #pdb.set_trace()
            asndict = fileutil.readAsnTable(input,output)
            
            # Build output filename
            if output == None:
                output = drutil.buildRootname(string.lower(asndict['output']),ext=['_drz.fits'])
                print 'Setting up output name: ',output
            # Check for existance of output file.
            if fileutil.findFile(output): 
                print 'Removing previous output product...'
                os.remove(output)

            if len(asndict['members'].keys()) > 1:
                # pdb.set_trace()
                self.observation = DitherProduct(asndict,pars=self.pars,ref=self.ref) #WZ
            else:
                inroot = asndict['members'].keys()[0]
                pardict = asndict['members'][inroot]
                infile = drutil.buildRootname(inroot)
                # Append user specified parameters to shifts dictionary
                pardict.update(self.pars)
                self.observation = selectInstrument(infile,output,ref,pars=pardict) #WZ
        self.output = output

        # This call puts together the parameters for the input image
        # with those for the output to create a final parameter list
        # for running 'drizzle'.
        # It relies on the buildPars() methods for each exposure to
        # generate a complete set of parameters for all inputs
        #
        print "Main step to set pars"
        #pdb.set_trace() # WZ main step to set pars
        twk = {} #WZ >>
        twk['val']= tweak
        twk['name']=self.ref.tweak
        self.parlist = self.observation.buildPars(ref=field,tweak=twk)
        # Let the user know parameters have been successfully calculated
        print 'Drizzle parameters have been calculated. Ready to .run()...'
        print 'Finished calculating parameters at ',_ptime()
        
    # Run 'drizzle' here...
    #
    def run(self,save=no,build=yes):
        """Perform drizzle operation on input to create output.
         This method would rely on the buildPars() method for
         the output product to produce correct parameters
         based on the inputs. The output for buildPars() is a list 
         of dictionaries, one for each input, that matches the 
         primary parameters for an IRAF drizzle task.

         This method would then loop over all the entries in the 
         list and run 'drizzle' for each entry. """
        print "run Drizzle"
        #pdb.set_trace()
        #for i in range(len(self.parlist)): #WZ normalize WCS 
        #    substr = str.split(self.parlist[i]['data'],'[')
        #    im=substr[0]
        #    tem=str.split(substr[1],',')[1]
        #    fitsfile = pyfits.open(im,'update')
        #    j = 3*int(str.split(tem,']')[0])-2 # sci array
        #    fitsfile[j].header.update('CRPIX1',self.ref.crpix1)
        #    fitsfile[j].header.update('CRPIX2',self.ref.crpix2)
        #    fitsfile[j].header.update('CRVAL1',self.ref.crval1)
        #    fitsfile[j].header.update('CRVAL2',self.ref.crval2)
        #    fitsfile.close()
        #    del fitsfile
        runDrizzle(self.parlist)

        # Remove the drizzle executable from the process cache
        # for memory management.
        iraf.flprc('drizzle')
        
        # Update the EXPTIME keyword in the output SCI image
        _sci=self.parlist[0]['outdata']+'[0]'
        drutil.updateKeyword(_sci,'EXPTIME',self.parlist[0]['texptime'])
        
        # Now, package outputs into 'self.output' and clean-up temp files
        # Start by copying one of the inputs into the output file
        # Then, copy image data from temp files into new output
        # Next, update keywords CD*, CRPIX*, CRVAL*, ORIENTAT, 
        #	CTYPE*, and EXPTIME from self.parlist[0]['outdata']
        _wcs = self.observation.product.geometry.wcs
        _pardict=self.parlist[0]
        if build:
            fileutil.buildDthProduct(_pardict,self.output,wcs=_wcs)
        else:
            # If we do not create _drz product, then we need to keep 
            # the drizzle products no matter what...
            save = yes
            
        # Finally, delete seperate drizzle output files, if desired. 
        if not save:
            # Clean up directory by deleting dither products once they are copied...
            os.remove(_pardict['outdata'])
            os.remove(_pardict['outweight'])
            if fileutil.findFile(_pardict['outcontext']):
                os.remove(_pardict['outcontext'])
        
        # Remove temp files, if desired
        # Files to be removed are: 
        #   parlist['coeffs']
        #   parlist['in_mask']
        if self.clean:
            for img in self.parlist:
                os.remove(img['coeffs'])
                if img['in_mask'] != '':
                    os.remove(img['in_mask'])

    def help(self):
        """
	        Run the module level help function to provide syntax information.
        """
        help()


def runDrizzle(parlist,tdriz=None):
    """ This function sets up the parameters for the version of drizzle which
    has been tested with this wrapper, then starts the drizzle task.
    WZ modify input file WCS headers to that of ref, so that the output file will follow
    """
    # Check to see if drizzle is Version 2.6 or greater
    _tdriz = iraf.drizzle.getParDict().has_key('kernel')
    
    for plist in parlist:
        iraf.drizzle.data=plist['data']
        iraf.drizzle.outdata=plist['outdata']
        iraf.drizzle.outweig=plist['outweight']
        #if plist['in_mask'] != None:
        # pdb.set_trace() #WZ  >>
        if 'in_mask' in plist:
            # pdb.set_trace() #WZ  >>
            #if string.find(plist['in_mask'],'wt') > -1:
            #    plist['in_mask']= string.replace(plist['in_mask'],'wt','inmask')
            whtname=string.replace(plist['in_mask'],'inmask','wt'); #WZ <<
            if (os.path.exists(whtname)):
                #pdb.set_trace()
                mskname = plist['in_mask']
                iraf.drizzle.in_mask = whtname; #WZ <<
                iraf.drizzle.in_mask = plist['in_mask']
            else:
                iraf.drizzle.in_mask = plist['in_mask']
        else:
            iraf.drizzle.in_mask = ''
            
        # iraf.drizzle.wt_scl=1.0 # WZ

        if _tdriz:
            # New Version 2.6 DRIZZLE parameters
            iraf.drizzle.outcont=plist['outcontext']
            iraf.drizzle.kernel=plist['kernel']

        # Parameters generated by PyDrizzle
        print "** Drizzle **"
        iraf.drizzle.outnx=plist['outnx']
        iraf.drizzle.outny=plist['outny']
        iraf.drizzle.xsh=plist['xsh']
        iraf.drizzle.ysh=plist['ysh']
        iraf.drizzle.coeffs=plist['coeffs']
        iraf.drizzle.scale=plist['scale']
        # User selectable parameters
        iraf.drizzle.rot=plist['rot'] 
        print plist['xsh'],plist['ysh'],plist['rot'],plist['coeffs']
        iraf.drizzle.out_un=plist['units']
        iraf.drizzle.pixfrac=plist['pixfrac']
        # Fixed parameters
        iraf.drizzle.align='center'
        iraf.drizzle.shft_fr='output' # 'input' # WZ 
        iraf.drizzle.shft_un='output' # 'input' # WZ 
        # iraf.drizzle.expkey = DEXPTIME
        iraf.drizzle.expkey = "EXPTIME" #WZ
        iraf.drizzle.wt_scl = "exptime"
        iraf.drizzle.fillval=0.0
        iraf.drizzle.mode='h'
        # pdb.set_trace() #WZ
        iraf.drizzle()
        iraf.unlearn(iraf.drizzle)

        """
        # if _tdriz:
        if (string.find(plist['outdata'],'_dz') < 0):        
            iraf.drizzle.data=plist['data']
            iraf.drizzle.outdata=string.replace(plist['outdata'],'drz','drz2')
            iraf.drizzle.outweig=string.replace(plist['outweight'],'weight','exp')
            if string.find(plist['in_mask'],'wt') > -1:
                plist['in_mask']= string.replace(plist['in_mask'],'wt','inmask')
                iraf.drizzle.in_mask = plist['in_mask']
            iraf.drizzle.outcont=string.replace(plist['outcontext'],'context','context2')
            iraf.drizzle.kernel=plist['kernel']
            iraf.drizzle.outnx=plist['outnx']
            iraf.drizzle.outny=plist['outny']
            iraf.drizzle.xsh=plist['xsh']
            iraf.drizzle.ysh=plist['ysh']
            iraf.drizzle.coeffs=plist['coeffs']
            iraf.drizzle.scale=plist['scale']
            iraf.drizzle.rot=plist['rot'] 
            iraf.drizzle.out_un=plist['units']
            iraf.drizzle.pixfrac=plist['pixfrac']
            iraf.drizzle.align='center'
            iraf.drizzle.shft_fr='output' # 'input' # WZ 
            iraf.drizzle.shft_un='output' # 'input' # WZ 
            iraf.drizzle.expkey = "EXPTIME" #WZ
            iraf.drizzle.wt_scl = "exptime"
            iraf.drizzle.fillval=0.0
            iraf.drizzle.mode='h'
            iraf.drizzle()
            iraf.unlearn(iraf.drizzle)
        """
        # Once drizzle has completed its run, the WCS information
        # of the output image will need to be updated with the WCS
        # object from this wrapper.
        		
# End of 'runDrizzle'	
def help():
    help_str = """
Program to process and/or dither-combine image(s) using (t)drizzle.
To create an object named 'test' that corresponds to a drizzle product:
    --> test = xydrizzle.PyDrizzle(input)
where input is the FULL filename of an ACS observation or ASN table.
This computes all the parameters necessary for running drizzle on all
the input images.  Once this object is created, you can run drizzle using:
    --> test.run()
PyDrizzle object parameters:
    output      User-specified name for output products
    field       User-specified parameters for output image --
                psize, orient, ra, dec, shape
    units       Units for final product -- 'counts' or 'cps'(DEFAULT)
    kernel      Specify which kernel to use in TDRIZZLE--
                'square'(default),'point','gaussian','turbo','tophat'
    pixfrac     drizzle pixfrac value (Default: 1.0) 
    bits        Bit values considered good (Default: 0)
                'None' specifies that no input mask will be used
    idckey      User-specified keyword for determining IDCTAB filename--
                'IDCTAB'(ACS default),'TRAUGER'(WFPC2),'CUBIC'(WFPC2)
    clean       Remove temp files? -- yes(DEFAULT) or no              
.run() parameters:
    save        keep simple drizzle products? -- yes or no(DEFAULT)
    build       create multi-extension product? -- yes(DEFAULT) or no
Usage of optional parameters:    
    --> test = xydrizzle.PyDrizzle('test_asn.fits',bits=8194,units='counts')
	"""
    print help_str	

