#!/usr/local/bin/python -u

# $Id: wcsclass.py,v 1.3 2002/09/10 23:23:56 jpb Exp $
# ---------------------------------------------------------------------
# $Log: wcsclass.py,v $
# Revision 1.3  2002/09/10 23:23:56  jpb
# Added handling of 4th order IDCTABs.
#
# Revision 1.2  2002/05/13 16:44:21  anderson
# removed some header weirdness.  there appeared to be two places
# where __version__ and __version_date__ where being defined.
#
# K Anderson
#
# Adapted from the pipeline/lib/python module wcsclass to use numarray
# K Anderson. 22-Apr-2002
#
# wcsclass.py: Image FITS header WCS class and routines, mainly for ACS.
# jpb, 08-Jan-2001
#
# Two different classes are defined: BasicWCS and idcWCS.  Some of the
# BasicWCS methods are adapted or copied from functions in Warren Hack's
# pydrizzle code (but it's taxpayer's money, right?).  But here we don't
# use or import iraf and there's no duality between WCs attributes and a
# wcsinfo dictionary, which was confusing because you could change
# values in the dictionary without actually changing the wcs attributes.
# Here all the keywords only appear in the wcs dictionary (no other
# critical attributes).  Also, the idcWCS class implements the PVi_j
# FITS keywords coefficients, which it derives from orientation information
# in the header and the distortion coefficients in a supplied IDC table.
# This is a bit unfortunate, in that the PV coefficients become dependent
# on orientation, but I think it has to be this way to follow the
# Calabretta & Greisen definition (the problem is that the IDCTAB solution
# transform to the telescope V2/V3 coordinate system, not to RA/Dec).
#
# K Anderson 22-apr-2002

__version__      = '$Revision: 1.3 $ '[11:-3]
__version_date__ = '$Date: 2002/09/10 23:23:56 $ '[7:-3]

import math,sys
import matutil
from numpy import *

# define a few general conveniences
yes,no = 1,0
def deg2rad(deg):
    return (deg * pi / 180.) 
def rad2deg(rad):
    return (rad * 180. / pi)

def sgn(x):
    if x == 0: val = +1
    else:      val = x/abs(x)
    return val

def cubeval(x,y,coeff):
    """Evaluate a cubic geometric distortion of the form used by drizzle, etc."""
    eval = coeff[0] + coeff[1]*x + coeff[2]*y + \
           coeff[3]*x*x + coeff[4]*x*y + coeff[5]*y*y + \
           coeff[6]*x*x*x + coeff[7]*x*x*y + coeff[8]*x*y*y + coeff[9]*y*y*y 
    return eval

# List of the most basic WCS keys
basicWCSkeys = ['NAXIS1', 'NAXIS2', 'CTYPE1', 'CTYPE2',
		'CRPIX1', 'CRPIX2', 'CRVAL1', 'CRVAL2',
		'CD1_1',  'CD1_2',  'CD2_1',  'CD2_2', 'ORIENTAT']

class BasicWCS:
    """    A class to encapsulate the basic WCS information from an image header.
    The class contains the following methods:
         update(pscale,orient,refpos,refval,size)  # all pars default to None
         xy2rd(xypos)
         rd2xy(skypos,hour=0)
         rotateWCS(orient)
         wcsprint()
         copy2header(header)
    These allow the user to change wcs information, convert from x,y pos
    to sky pos and back again; rotate the wcs matrix; copy the wcs to another
    image header.
    """
    def __init__(self, header,shape=None): 
        # Initialize wcsinfo info from header.  If header uses old
        # CDELTi/CROTj convention (unlike HST instruments), convert
        # it to CDi_j matrix representation.

        self.wcs = {}
        self.wcskeys = []      # this just keeps names ordered

        # check to see if header uses old CDELTi/CROTj convention instead of CDi_j
        _old_system_ = 0
        if(not header.has_key('CD1_1')):
           if header.has_key('CDELT1') and header.has_key('CDELT2'):
               _old_system_ = 1
               CDELT1 = header['CDELT1']
               CDELT2 = header['CDELT2']
               if header.has_key('CROTA2'):
                   CROTA2=header['CROTA2']
               else:
                   CROTA2 = 0.0
               cosROT = math.cos(deg2rad(CROTA2))
               sinROT = math.sin(deg2rad(CROTA2))
           else:
               raise KeyError,"Header lacks both CDij matrix and CDELTi/CROTAj keys."

        for key in basicWCSkeys:
            self.wcskeys.append(key)
	    if not header.has_key(key):
                if _old_system_:
                    # following Calabretta & Greisen:
                    if   key == 'CD1_1':
                        self.wcs[key] =  CDELT1 * cosROT
                    elif key == 'CD1_2':
                        self.wcs[key] = -CDELT2 * sinROT
                    elif key == 'CD2_1':
                        self.wcs[key] =  CDELT1 * sinROT
                    elif key == 'CD2_2':
                        self.wcs[key] =  CDELT2 * cosROT
                    elif key == 'ORIENTAT':
                        # from the above, it seems CROTA must be -PA
                        self.wcs[key] = -CROTA2
                elif key == 'ORIENTAT':
                    if self.wcs['CD2_2'] == 0:
                        if self.wcs['CD2_1'] > 0:
                            self.wcs[key] =  90.0
                        else:
                            self.wcs[key] = -90.0
                    else:
                        self.wcs[key] = math.atan(self.wcs['CD2_1']/self.wcs['CD2_2']) 
                else:
                    raise KeyError,"Keyword "+key+" not found!"
            else:
                self.wcs[key] = header[key]

        if self.wcs['CTYPE1'][0:8] != 'RA---TAN' or self.wcs['CTYPE2'][0:8] != 'DEC--TAN':
            raise Exception,'Only handle RA/Dec tangent plane projections.'

        # then add one more, not officially a WCS keyword, but nice to have
        #self.pscale = pow((pow(self.cd11,2)+pow(self.cd12,2)),0.5) * 3600.
        self.wcskeys.append('PSCALE')
        self.wcs['PSCALE'] = 3600. * math.sqrt((self.wcs['CD1_1']**2 + self.wcs['CD1_2']**2 + \
                                        self.wcs['CD2_1']**2 + self.wcs['CD2_2']**2)/2.0)

    def update(self, pscale=None, orient=None, refpos=None, refval=None, size=None):
        """
        Update the WCS parameters: orient, refpos, refval, size
        Create a new CD Matrix from the absolute pixel scale 
        and reference image orientation.  The resulting CD matrix
        rotation is only right if there's no distortion, so the
        same pixel scale in both directions, etc.
        """
        # Set up parameters necessary for updating WCS
        # Check to see if new value is provided,
        # If not, fall back on old value as the default
        updateCD = no
        
        if orient != None:
	    pa = deg2rad(orient)
            self.wcs['ORIENTAT'] = orient
            updateCD = yes
        else:
            pa = deg2rad(self.wcs['ORIENTAT'])   # in case pscale was specified
        
        if pscale != None:
            self.wcs['PSCALE'] = pscale
            updateCD = yes
        else:
            pscale = self.wcs['PSCALE']         # in case orient was specified
            
        if size != None:
	        self.wcs['NAXIS1'] = size[0]
	        self.wcs['NAXIS2'] = size[1]

        if refpos != None:
	        self.wcs['CRPIX1'] = refpos[0]
	        self.wcs['CRPIX2'] = refpos[1]

        if refval != None:
	        self.wcs['CRVAL1'] = refval[0]
	        self.wcs['CRVAL2'] = refval[1]

        # Update CD matrix if needed
        # note, this is only right if image not flipped w.r.t. sky
        if updateCD:
            self.wcs['CD1_1'] = -pscale * cos(pa) / 3600.
            self.wcs['CD1_2'] =  pscale * sin(pa) / 3600.
            self.wcs['CD2_1'] =  pscale * sin(pa) / 3600.
            self.wcs['CD2_2'] =  pscale * cos(pa) / 3600. 

    def xy2rd(self,xypos):
        """ Apply the CDi_j keywords to a position to generate a new sky position.
        This is from Warren, who says it comes directly from 'imgtools.xy2rd'
        translate (x,y) to (ra, dec)
        """
        xi  = self.wcs['CD1_1'] * (xypos[0] - self.wcs['CRPIX1']) + \
              self.wcs['CD1_2'] * (xypos[1] - self.wcs['CRPIX2'])
        eta = self.wcs['CD2_1'] * (xypos[0] - self.wcs['CRPIX1']) + \
              self.wcs['CD2_2'] * (xypos[1] - self.wcs['CRPIX2'])

        xi = deg2rad(xi)
        eta = deg2rad(eta)
        ra0 = deg2rad(self.wcs['CRVAL1'])
        dec0 = deg2rad(self.wcs['CRVAL2'])

        ra = arctan((xi / (cos(dec0)-eta*sin(dec0)))) + ra0
        dec = arctan( ((eta*cos(dec0)+sin(dec0)) / 
                (sqrt((cos(dec0)-eta*sin(dec0))**2 + xi**2))) )

        ra = rad2deg(ra)
        dec = rad2deg(dec)
        ra = divmod (ra, 360.)[1]

        return ra,dec


    def rd2xy(self,skypos,hour=0):
        """
        This method would use the WCS keywords to compute the XY position
        from a given RA/Dec tuple (in deg).  Again, from Warren.
        """
        det = self.wcs['CD1_1']*self.wcs['CD2_2'] - self.wcs['CD1_2']*self.wcs['CD2_1']
        if det == 0.0: 
            raise ArithmeticError,"singular CD matrix!"

        cdinv11 =  self.wcs['CD2_2'] / det
        cdinv12 = -self.wcs['CD1_2'] / det
        cdinv21 = -self.wcs['CD2_1'] / det
        cdinv22 =  self.wcs['CD1_1'] / det

        # translate (ra, dec) to (x, y)

        ra0  = deg2rad(self.wcs['CRVAL1'])
        dec0 = deg2rad(self.wcs['CRVAL2'])
        if hour:
            skypos[0] = skypos[0] * 15.
        ra  = deg2rad(skypos[0])
        dec = deg2rad(skypos[1])

        bottom = sin(dec)*sin(dec0) + cos(dec)*cos(dec0)*cos(ra-ra0)
        if bottom == 0.0: 
            raise ArithmeticError,"Unreasonable RA/Dec range!"

        xi = cos(dec) * sin(ra-ra0) / bottom
        eta = (sin(dec)*cos(dec0) - cos(dec)*sin(dec0)*cos(ra-ra0)) / bottom 
        xi = rad2deg(xi)
        eta = rad2deg(eta)

        x = cdinv11 * xi + cdinv12 * eta + self.wcs['CRPIX1']
        y = cdinv21 * xi + cdinv22 * eta + self.wcs['CRPIX2']

        return x,y

    def copy2header(self,header):
        """Copies the wcs keywords from a wcs object to a given header""" 
        for key in self.wcs.keys():
            if header.has_key(key):
                header[key] = self.wcs[key]
            else:
                header.update(key, self.wcs[key])

    def wcsprint(self):
        for key in self.wcskeys:
            print "%12s:  %s" %(key,self.wcs[key])

    def help(self):
        print self.__doc__

class idcWCS:
    """    A WCS class containing the PVi_j keywords (among other stuff).
    This is not a subclass of BasicWCS because that had unwanted side effects.
    The constructor gets a BasicWCS object, rather than a header object.
    The MaxPV optional parameter allows the user to limit the number of PV,
    with the default being 11 for cubic transformations, as in the case of
    the IDC transformation.  The methods of this class are *not* general --
    they expect the terms found in the IDCTAB, thus the name of the class.

    The constructor does the following steps:
       1. Copies over the WCS keywords from the input baseWCS.
       2. Reads the distortion solution from the IDCTAB and checks
            that things are as expected.
       3. Determines the change in ORIENTAT produced in applying the
            de-distortion transformation; updates ORIENTAT key.
       4. Translates the reference pixel: set the new CRPIX to be
            the IDC REF pix, and determine it's RA,Dec from that of
            the old CRPIX; update CRVAL to be these RA,Dec coords.
       5. Constructs the PVi_k distortion keywords.

    The copy2header method then allows you to copy the revised WCS
    information to the image header.
    """
    def __init__(self, baseWCS, idcfile, chipID, MaxPV=None,identCD=0):
        # Initialize wcsinfo info from BasicWCS instance sent as baseWCS.
        # Note, this assumes the input WCS uses the modern CD-maxtrix convention
        # (as in ACS), which it will if it's a BasicWCS object (see above).

        self.wcs = {}
        self.wcskeys = []
        self.identCD = identCD
        
        # get distortion information
        CxMat,CyMat,refdict,self.Norder = matutil.getIDCinfo(idcfile,chipID)

        if self.Norder == 3:
            self.MaxPV = 10
        elif self.Norder == 4:
            self.MaxPV = 16
        else:
            print 'Error: Can only deal with cubic or quartic IDC transformations.'
            sys.exit()

        xco,yco = matutil.convertIDCtrans(CxMat,CyMat)

        for key in baseWCS.wcskeys:
            self.wcskeys.append(key)
            self.wcs[key] = baseWCS.wcs[key]

        # check that we have projections and distortion maps we can handle
        if self.wcs['CTYPE1'][0:8] != 'RA---TAN' or self.wcs['CTYPE2'][0:8] != 'DEC--TAN':
            raise Exception,'Only handle RA/Dec tangent plane projections.'
        if (xco[0] != yco[0] or xco[0] != 0.):
            raise Exception, "Zeroth order IDCtab terms should be 0."
        if (len(xco) != 10 and len(xco) != 15)  or (len(yco) != 10 and len(yco) != 15):
            raise Exception, 'Distortion coefficient array not of length 10 or 15.'

        # now create the distortion keywords
        for mu in range(1,3):
            for i in range(self.MaxPV+1):
                _newkey = 'PV'+str(mu)+'_'+str(i)
                self.wcskeys.append(_newkey)
                self.wcs[_newkey] = 0.0
	
        # by what angle does the IDC transformation change ORIENTAT?
        # Dtheta_y = math.atan(xco[2]/yco[2])
        # Dtheta_y = rad2deg(Dtheta_y)
        # that's not really appropropriate, a more accurate way is:
        Dtheta_y = self._getYrot(CxMat,CyMat)

        # make a new keyword for delta_ y-rotation 
        self.wcskeys.append('DTHETA_Y')
        self.wcs['DTHETA_Y'] = Dtheta_y

        # update WCS ORIENTAT for rotation resulting from dedistorting
        self.wcs['ORIENTAT'] += Dtheta_y

        # Ok, we need to find Ra,Dec of the refdict ref pixel, which will
        # be the new CRPIX. The idcWCS._translate function updates:
        # wcs['CRPIX1'], wcs['CRPIX2'], wcs['CRVAL1'], wcs['CRVAL2'] 
    
        self._translateCRPIXs(xco, yco, refdict)

        # rotate the 'rectified' x,y axes celestial axes and populate keywords
        self._construct(xco,yco, refdict)


    def _construct(self, Cx, Cy, refdict):
        " Populate the PVi_j keywords from the distortion coefficients."

        cosT = math.cos(self.wcs['ORIENTAT'] * math.pi/180.)
        sinT = math.sin(self.wcs['ORIENTAT'] * math.pi/180.)

        Xscale = math.sqrt(Cx[1] * Cx[1] + Cx[2] * Cx[2])
        Yscale = math.sqrt(Cy[2] * Cy[2] + Cy[1] * Cy[1])
        
        # make a new keywords for the scales
        self.wcskeys.append('PSCALE') 
        self.wcskeys.append('XSCALE')
        self.wcskeys.append('YSCALE') 
        self.wcs['XSCALE'] = Xscale
        self.wcs['YSCALE'] = Yscale
        self.wcs['PSCALE'] = sqrt(Xscale * Yscale)

        if self.identCD:
            Linv = 1.
        else:
            Linv  = 3600. / self.wcs['PSCALE']
        Linv2 = Linv * Linv
        Linv3 = Linv**3
        Linv4 = Linv**4
        
        # note, all CDi_j elements are positive because distortion
        # coefficients take care of orientation
        self.wcs['CD1_1'] =  1.0/Linv
        self.wcs['CD1_2'] =  0.0
        self.wcs['CD2_1'] =  0.0
        self.wcs['CD2_2'] =  1.0/Linv

        if (Cx[0] != Cy[0] or Cx[0] != 0.):
            raise Exception, "Zeroth order IDCtab terms should be 0."
        else:
            self.wcs['PV1_0'] = 0.
            self.wcs['PV2_0'] = 0.
	self.wcs['PV1_3']  = 0.
	self.wcs['PV2_3']  = 0.
        
	self.wcs['PV1_1']  = (-cosT * Cx[1] + sinT * Cy[1]) * Linv/3600.
	self.wcs['PV1_2']  = (-cosT * Cx[2] + sinT * Cy[2]) * Linv/3600.
	self.wcs['PV1_4']  = (-cosT * Cx[3] + sinT * Cy[3]) * Linv2/3600.
	self.wcs['PV1_5']  = (-cosT * Cx[4] + sinT * Cy[4]) * Linv2/3600.
	self.wcs['PV1_6']  = (-cosT * Cx[5] + sinT * Cy[5]) * Linv2/3600.
	self.wcs['PV1_7']  = (-cosT * Cx[6] + sinT * Cy[6]) * Linv3/3600.
	self.wcs['PV1_8']  = (-cosT * Cx[7] + sinT * Cy[7]) * Linv3/3600.
	self.wcs['PV1_9']  = (-cosT * Cx[8] + sinT * Cy[8]) * Linv3/3600.
	self.wcs['PV1_10'] = (-cosT * Cx[9] + sinT * Cy[9]) * Linv3/3600.
        if self.Norder == 4:
            self.wcs['PV1_11'] = 0.
            self.wcs['PV1_12'] = (-cosT * Cx[10] + sinT * Cy[10]) * Linv4/3600.
            self.wcs['PV1_13'] = (-cosT * Cx[11] + sinT * Cy[11]) * Linv4/3600.
            self.wcs['PV1_14'] = (-cosT * Cx[12] + sinT * Cy[12]) * Linv4/3600.
            self.wcs['PV1_15'] = (-cosT * Cx[13] + sinT * Cy[13]) * Linv4/3600.
            self.wcs['PV1_16'] = (-cosT * Cx[14] + sinT * Cy[14]) * Linv4/3600.
	
	self.wcs['PV2_1']  = ( sinT * Cx[2] + cosT * Cy[2]) * Linv/3600.
	self.wcs['PV2_2']  = ( sinT * Cx[1] + cosT * Cy[1]) * Linv/3600.
	self.wcs['PV2_4']  = ( sinT * Cx[5] + cosT * Cy[5]) * Linv2/3600.
	self.wcs['PV2_5']  = ( sinT * Cx[4] + cosT * Cy[4]) * Linv2/3600.
	self.wcs['PV2_6']  = ( sinT * Cx[3] + cosT * Cy[3]) * Linv2/3600.
	self.wcs['PV2_7']  = ( sinT * Cx[9] + cosT * Cy[9]) * Linv3/3600.
	self.wcs['PV2_8']  = ( sinT * Cx[8] + cosT * Cy[8]) * Linv3/3600.
	self.wcs['PV2_9']  = ( sinT * Cx[7] + cosT * Cy[7]) * Linv3/3600.
	self.wcs['PV2_10'] = ( sinT * Cx[6] + cosT * Cy[6]) * Linv3/3600.
        if self.Norder == 4:
            self.wcs['PV2_11'] = 0.
            self.wcs['PV2_12'] = (sinT * Cx[14] + cosT * Cy[14]) * Linv4/3600.
            self.wcs['PV2_13'] = (sinT * Cx[13] + cosT * Cy[13]) * Linv4/3600.
            self.wcs['PV2_14'] = (sinT * Cx[12] + cosT * Cy[12]) * Linv4/3600.
            self.wcs['PV2_15'] = (sinT * Cx[11] + cosT * Cy[11]) * Linv4/3600.
            self.wcs['PV2_16'] = (sinT * Cx[10] + cosT * Cy[10]) * Linv4/3600.
            

	#this would be better (reverse way around),
	#but not the dumb new convention 
        # self.wcs['CD1_1'] = -cosT * Xscale
        # self.wcs['CD1_2'] =  sinT * Yscale
        # self.wcs['CD2_1'] =  sinT * Xscale
        # self.wcs['CD2_2'] =  cosT * Yscale
	#
        # self.wcs['PV1_0']  = Cx[0] / Xscale
        # self.wcs['PV1_1']  = Cx[1] / Xscale
        # self.wcs['PV1_2']  = Cx[2] / Xscale
        # self.wcs['PV1_3']  = 0.
        # self.wcs['PV1_4']  = Cx[3] / Xscale
        # self.wcs['PV1_5']  = Cx[4] / Xscale
        # self.wcs['PV1_6']  = Cx[5] / Xscale
        # self.wcs['PV1_7']  = Cx[6] / Xscale
        # self.wcs['PV1_8']  = Cx[7] / Xscale
        # self.wcs['PV1_9']  = Cx[8] / Xscale
        # self.wcs['PV1_10'] = Cx[9] / Xscale
	# 
        # self.wcs['PV2_0']  = Cy[0] / Yscale
        # self.wcs['PV2_1']  = Cy[2] / Yscale
        # self.wcs['PV2_2']  = Cy[1] / Yscale
        # self.wcs['PV2_3']  = 0.
        # self.wcs['PV2_4']  = Cy[5] / Yscale
        # self.wcs['PV2_5']  = Cy[4] / Yscale
        # self.wcs['PV2_6']  = Cy[3] / Yscale
        # self.wcs['PV2_7']  = Cy[9] / Yscale
        # self.wcs['PV2_8']  = Cy[8] / Yscale
        # self.wcs['PV2_9']  = Cy[7] / Yscale
        # self.wcs['PV2_10'] = Cy[6] / Yscale
    
    def update(self, refpos=None, refval=None, size=None):
        """   Update the WCS parameters:  refpos, refval, size
        This is a much more limited update routine, basically allowing
        only changes in reference position and values, than that
        of the BasicWCS superclass.  The reason is that any change
        of scale or rotation changes all N header distortion keywords.
        """
        # Set up parameters necessary for updating WCS
        # Check to see if new value is provided,
        # If not, fall back on old value as the default

        if size != None:
	        self.wcs['NAXIS1'] = size[0]
	        self.wcs['NAXIS2'] = size[1]

        if refpos != None:
	        self.wcs['CRPIX1'] = refpos[0]
	        self.wcs['CRPIX2'] = refpos[1]

        if refval != None:
	        self.wcs['CRVAL1'] = refval[0]
	        self.wcs['CRVAL2'] = refval[1]

    def _getYrot(self,CxMat,CyMat):
        "Find Y-axis rotation caused by IDC matrix transformation."

        # want the angle of rotation of y-axis, "top" to "bottom" so:
        x_raw = 0.0
        y_raw = self.wcs['NAXIS2']/2.0
        
        # x_rect,y_rect = matutil.xyIDCtrans(x_raw, y_raw, CxMat, CyMat)

        x_top,y_top = matutil.xyIDCtrans(x_raw, y_raw, CxMat, CyMat)
        x_bot,y_bot = matutil.xyIDCtrans(x_raw, -y_raw, CxMat, CyMat)
        
        x_rect = x_top - x_bot
        y_rect = y_top - y_bot

        theta = math.atan(x_rect/y_rect) * 180./math.pi
        
        return theta
    
    def _translateCRPIXs(self, Cx, Cy, refdict):
        """Change CRPIX1,2 to be reference pixel from IDC map, and then
        find RA,Dec coordinates for that pixel and thus the new CRVAL's.
        """
        # first save coords of old CRPIX
        crpix1_old =  self.wcs['CRPIX1']
        crpix2_old =  self.wcs['CRPIX2']
        crval1_old =  self.wcs['CRVAL1']
        crval2_old =  self.wcs['CRVAL2']

        # Now, make the new rotated-frame distortion matrices,
        # which should take us to RA,DEC axes instead of |V2|,|V3|
        cosT = math.cos(self.wcs['ORIENTAT'] * math.pi/180.)
        sinT = math.sin(self.wcs['ORIENTAT'] * math.pi/180.)
        ABx = (-cosT * Cx + sinT * Cy) / 3600.
        ABy = ( sinT * Cx + cosT * Cy) / 3600.

        # evaluate the polynomial at the old crpix to find the offset in
        # tangential arcsec from the IDC reference pixel
        x = crpix1_old - refdict['XREF']
        y = crpix2_old - refdict['YREF']
        xi  = cubeval(x,y,ABx)             # these are now 
        eta = cubeval(x,y,ABy)             # in degrees
              
        # flip sign to make the direction from old_CRPIX1,2 to XREF,YREF
        # and change from deg to radians
        xi  = deg2rad((-xi)) 
        eta = deg2rad((-eta))

        # also convert units of old crval1,2 from degrees to radians
        ra0  = deg2rad(crval1_old)
        dec0 = deg2rad(crval2_old)

        # invert the gnomic projection to get true RA,Dec (in rad)

        ra = ra0 + arctan( xi / (cos(dec0) - eta*sin(dec0)) )
        
        dec = arctan( (eta*cos(dec0) + sin(dec0)) / 
                      sqrt( (cos(dec0) - eta*sin(dec0))**2 + xi**2 ) )

        # convert RA/Dec back to degrees
        ra = rad2deg(ra)
        ra = divmod (ra, 360.)[1]
        dec = rad2deg(dec)

        # finally, update the key WCS keys
        self.wcs['CRPIX1'] = refdict['XREF']
        self.wcs['CRPIX2'] = refdict['YREF']
        self.wcs['CRVAL1'] = ra
        self.wcs['CRVAL2'] = dec
        return
        
    def xy2rd(self, xypos):
        """  Apply the WCS keywords, INCLUDING PV DISTORTION COEFFICIENTS
        to an x,y position to generate a new sky position.
        """

        # first apply the CD matrix
        x  = self.wcs['CD1_1'] * (xypos[0] - self.wcs['CRPIX1']) + \
             self.wcs['CD1_2'] * (xypos[1] - self.wcs['CRPIX2'])
        y  = self.wcs['CD2_1'] * (xypos[0] - self.wcs['CRPIX1']) + \
             self.wcs['CD2_2'] * (xypos[1] - self.wcs['CRPIX2'])

        # can make it more general for other solutions, but for now...
        
        if (self.wcs['PV1_0'] != 0. or self.wcs['PV2_0'] != 0.):
            raise ValueError, 'Zero-order terms of projection should be zero!'

        if (self.wcs['PV1_3']  != 0. or	self.wcs['PV2_3'] != 0.):
            raise ValueError, 'radial terms of IDCTAB projection should be zero!!'

        # then the distortion coefficients (note: no PVi_3 terms!):
        xi  = self.wcs['PV1_1']  * x     + \
              self.wcs['PV1_2']  * y     + \
              self.wcs['PV1_4']  * x*x   + \
              self.wcs['PV1_5']  * x*y   + \
              self.wcs['PV1_6']  * y*y   + \
              self.wcs['PV1_7']  * x*x*x + \
              self.wcs['PV1_8']  * x*x*y + \
              self.wcs['PV1_9']  * x*y*y + \
              self.wcs['PV1_10'] * y*y*y    
                                 
        eta = self.wcs['PV2_1']  * y     + \
              self.wcs['PV2_2']  * x     + \
              self.wcs['PV2_4']  * y*y   + \
              self.wcs['PV2_5']  * y*x   + \
              self.wcs['PV2_6']  * x*x   + \
              self.wcs['PV2_7']  * y*y*y + \
              self.wcs['PV2_8']  * y*y*x + \
              self.wcs['PV2_9']  * y*x*x + \
              self.wcs['PV2_10'] * x*x*x

        if self.Norder == 4:
            xi  += (self.wcs['PV1_12'] * x*x*x*x + \
                    self.wcs['PV1_13'] * x*x*x*y + \
                    self.wcs['PV1_14'] * x*x*y*y + \
                    self.wcs['PV1_15'] * x*y*y*y + \
                    self.wcs['PV1_16'] * y*y*y*y)
            
            eta += (self.wcs['PV2_12'] * y*y*y*y + \
                    self.wcs['PV2_13'] * x*y*y*y + \
                    self.wcs['PV2_14'] * x*x*y*y + \
                    self.wcs['PV2_15'] * x*x*x*y + \
                    self.wcs['PV2_16'] * x*x*x*x)


        ## now proceed with the gnomic projection as usual
        # xi,eta and crvals are all in degrees now; make radians
        xi = deg2rad(xi)
        eta = deg2rad(eta)
        ra0 = deg2rad(self.wcs['CRVAL1'])
        dec0 = deg2rad(self.wcs['CRVAL2'])

        ra = ra0 + arctan( xi / (cos(dec0) - eta*sin(dec0)) )
        dec = arctan( (eta*cos(dec0) + sin(dec0)) / 
                      sqrt((cos(dec0) - eta*sin(dec0))**2 + xi**2) )

        ra = rad2deg(ra)
        ra = divmod (ra, 360.)[1]
        dec = rad2deg(dec)
        return ra,dec

    def wcsprint(self):
        for key in self.wcskeys:
            print "%12s:  %s" %(key,self.wcs[key])
    def help(self):
        print self.__doc__

    def copy2header(self,header,cnpix=0):
        "Copy the wcs keywords from a idcWCS object to a given header."
        # one more paranoic check
        if self.wcs['CTYPE1'][0:8] != 'RA---TAN' or self.wcs['CTYPE2'][0:8] != 'DEC--TAN':
            raise Exception,'Only handle RA/Dec tangent plane projections.'

        # copy all the standard WCS keys (should already be in header)
        for key in basicWCSkeys:
            header.update(key,self.wcs[key])

        if header.get('EQUINOX') != 2000.0:
            raise ValueError, 'EQUINOX needs to be 2000 for this to be valid.'
        
        if not header.has_key('RADESYS'):
            header.update('RADESYS','FK5', after='EQUINOX')

        # ok now for the new, distorted case keywords
        self.pvList = []
        for i in range(2,0,-1):
            for j in range(self.MaxPV,-1,-1):
                self.pvList.append('PV'+str(i)+'_'+str(j))

        for PVkey in self.pvList:
            header.update(PVkey,self.wcs[PVkey],after='CD2_2')

        # optionally add cnpix's so as to avoid possible confusion 
        if cnpix:
            header.update('CNPIX1',0, after=self.pvList[0])
            header.update('CNPIX2',0, after='CNPIX1')

        # and a few more Keywords that I'm making up
        header.update('DTHETA_Y',self.wcs['DTHETA_Y'],after=self.pvList[0])
        header.update('PSCALE_Y',self.wcs['YSCALE'], after=self.pvList[0])
        header.update('PSCALE_X',self.wcs['XSCALE'], after=self.pvList[0])

        # and I think 'LONGPOLE' should be set like this
        header.update('LONGPOLE',180.0, after=self.pvList[0])
        return
